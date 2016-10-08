#!/usr/bin/env python

import sys, os, subprocess, argparse, shutil, warnings, csv, multiprocessing
from Bio import SeqIO
from natsort import natsorted
import pandas as pd
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    from Bio import SearchIO

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='findHMM.py',
    description='''Find HMM model in GenBank flatfile genome.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', nargs='+', required=True, help='Genome (GBK) format')
parser.add_argument('-o','--out', default='findHMM', help='Output Basename')
parser.add_argument('-m','--hmm', required=True, help='HMM file')
parser.add_argument('--evalue', default='1e-50', help='HMM evalue')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs')
parser.add_argument('--debug', action='store_true', help='Debug')
parser.add_argument('--maxIntron', default='3000', help='Maximum intron length')
parser.add_argument('--split_hmms', action='store_true', help='Output multi-FASTA for each HMM model')
parser.add_argument('--reuse', help='Reuse data in tmp directory')
args=parser.parse_args()
FNULL = open(os.devnull, 'w')

def getSize(filename):
    st = os.stat(filename)
    return st.st_size

def countfasta(input):
    count = 0
    with open(input, 'rU') as f:
        for line in f:
            if line.startswith (">"):
                count += 1
    return count

def group_by_heading( some_source ):
    buffer= []
    for line in some_source:
        if line.startswith( ">" ):
            if buffer: yield buffer
            buffer= [ line ]
        else:
            buffer.append( line )
    yield buffer

def gb2output(input, output1, output3):
    with open(output1, 'w') as proteins:
        with open(output3, 'w') as scaffolds:
            with open(input, 'rU') as gbk:
                SeqRecords = SeqIO.parse(gbk, 'genbank')
                for record in SeqRecords:
                    scaffolds.write(">%s\n%s\n" % (record.id, record.seq))
                    for f in record.features:
                        if f.type == "CDS":
                            try:
                                locusID = f.qualifiers['locus_tag'][0]
                            except KeyError: #if no locus_id it isn't a real locus or is partial?
                                continue
                            try:
                                protID = f.qualifiers['protein_id'][0]
                            except KeyError:
                                protID = '???'
                            try:
                                protSeq = f.qualifiers['translation'][0]
                            except KeyError:
                                continue
                            proteins.write(">%s\n%s\n" % (locusID+'_'+protID, protSeq))

def tblastnFilter(input, query, cpus, output):
    global HitList, Scaffolds, tBlastN
    HitList = []
    Scaffolds = []
    tBlastN = {}
    #start by formatting blast db/dustmasker filtered format
    subprocess.call(['dustmasker', '-in', input, '-infmt', 'fasta', '-parse_seqids', '-outfmt', 'maskinfo_asn1_bin', '-out', os.path.join(output,'genome_dust.asnb')], stdout = FNULL, stderr = FNULL)
    subprocess.call(['makeblastdb', '-in', input, '-dbtype', 'nucl', '-parse_seqids', '-mask_data', os.path.join(output, 'genome_dust.asnb'), '-out', os.path.join(output, 'genome')], stdout = FNULL, stderr = FNULL)
    #okay, now run tblastn using uniprot proteins
    subprocess.call(['tblastn', '-num_threads', str(args.cpus), '-db', os.path.join(output, 'genome'), '-query', query, '-max_target_seqs', '1', '-db_soft_mask', '11', '-threshold', '999', '-max_intron_length', args.maxIntron, '-evalue', '1e-5', '-outfmt', '6', '-out', os.path.join(output,'filter.tblastn.tab')], stdout = FNULL, stderr = FNULL)
    #now parse through results, generating a list for exonerate function
    with open(os.path.join(output, 'filter.tblastn.tab')) as input:
        reader = csv.reader(input, delimiter='\t')
        for cols in reader:
            hit = cols[0] + '::' + cols[1]
            if hit not in HitList:
                HitList.append(hit)
            if cols[1] not in Scaffolds:
                Scaffolds.append(cols[1])
            if cols[0] not in tBlastN:
                tBlastN[cols[0]] = (cols[1]+":"+cols[8]+"-"+cols[9], cols[11], cols[10], 'NA', 'tBlastn')

def runExonerate(input):
    global Missing
    s = input.split('::')
    if s[0].startswith('sp|'):
        name = s[0].split("|")[1] + '_' + s[1]
    else:
        name = s[0].split()[0] + '_' + s[1]
    query = os.path.join(tmpdir, name+'.fa')
    with open(query, 'w') as output:
        rec = record_dict[s[0]]
        output.write(">%s\n%s\n" % (rec.id, rec.seq))
    scaffold = s[1] + '.fasta'
    scaffold = os.path.join(tmpdir, scaffold)
    exonerate_out = 'exonerate_' + name + '.out'
    exonerate_out = os.path.join(tmpdir, exonerate_out)
    ryo = ">%qi|pident=%pi|%ti:%tcb-%tce|Exonerate-Partial\n%tcs\n"
    with open(exonerate_out, 'w') as output4:
        subprocess.call(['exonerate', '--model', 'p2g', '--showvulgar', 'no', '--showalignment', 'no', '--showquerygff', 'no', '--showtargetgff', 'no', '--maxintron', args.maxIntron, '--percent', '25', '--ryo', ryo , query, scaffold], stdout = output4)
    os.remove(query)
    
def runHMMsearch(input, basename, tmpdir, cpus, evalue, hmm):
    Results = {}
    #load proteins into dictionary
    protein_dict = SeqIO.to_dict(SeqIO.parse(input, 'fasta'))
    #do hmmer search of proteins
    HMM = os.path.join(tmpdir, basename+'.hmmsearch.txt')
    subprocess.call(['hmmsearch', '-o', HMM, '--cpu', str(cpus), '-E', evalue, hmm, input], stdout = FNULL, stderr = FNULL)
    with open(HMM, 'rU') as results:
        for qresult in SearchIO.parse(results, "hmmer3-text"):
            query_length = qresult.seq_len #length of HMM model
            hits = qresult.hits
            num_hits = len(hits)
            if num_hits > 0:
                query = hits[0].id
                hit = hits[0].query_id
                score = hits[0].bitscore
                evalue = hits[0].evalue
                num_hsps = len(hits[0].hsps)
                aln_length = 0
                for x in range(0,num_hsps):
                    aln_length += hits[0].hsps[x].aln_span
                if hit not in Results:
                    Results[hit] = [query, score, evalue, aln_length, 'Hmmer3']
    for k,v in Results.items():
        description = base+'|'+k+"|"+v[0]+"|evalue="+str(v[2])+"|HMMer3-Complete"
        Results[k].append(description)
        Seq = str(protein_dict[v[0]].seq)
        Results[k].append(Seq)
    return Results


#set up tmpdir
tmpdir = 'tmp_'+str(os.getpid())
os.makedirs(tmpdir)

#setup dictionary to hold summary results
AllResults = []
labels = []
species = []
#setup output files    
FinalOut = args.out+'.proteins.fasta'
TextOut = args.out+'.hits.tsv'
SummaryOut = args.out+'.summary.csv'
with open(TextOut, 'w') as output:
    output.write('Genome\tHMM-Model\tHit\tBitScore\tEvalue\tAlign-Len\tMethod\n')

#make sure final output isn't there
if os.path.isfile(FinalOut):
    os.remove(FinalOut)
for file in args.input:
    #Split GBK into parts
    base = file.rsplit('.', 1)[0]
    if '/' in base:
        base = base.split('/') [-1]
    species.append(base)
    labels.append(base)
    Proteins = os.path.join(tmpdir, base+'.proteins.fa')
    Genome = os.path.join(tmpdir, base+'.genome.fa')
    gb2output(file, Proteins, Genome)
    
    #print status
    print '----------------------------------------------'
    print 'Working on %s' % base
    
    #check number of HMMer models
    Results = {}
    HMMstat = os.path.join(tmpdir, 'hmmstat.txt')
    if not os.path.isfile(HMMstat):
        with open(HMMstat, 'w') as output:
            subprocess.call(['hmmstat', args.hmm], stdout = output, stderr = FNULL)
        HMMmodels = []
        with open(HMMstat, 'rU') as input:
            for line in input:
                if line.startswith('\n'):
                    continue
                if not line.startswith('#'):
                    cols = line.split(' ')
                    cols = filter(None, cols)
                    if not cols[1] in HMMmodels:
                        HMMmodels.append(cols[1])
    print "Looking for %i protein HMM model(s)" % len(HMMmodels)
    
    #check for annotated genome
    Protsize = getSize(Proteins)
    if Protsize > 300:
        print "Scanning proteome using HMMsearch"       
        Results = runHMMsearch(Proteins, base, tmpdir, args.cpus, args.evalue, args.hmm)
    else:
        print "No annotation found in genome, will search DNA"
    
    #check for missing models
    notfound = []
    for i in HMMmodels:
        if not i in Results:
            notfound.append(i)

    if len(notfound) > 0: #have to do some more work here for these to be sure they really don't exist
        #get consensus from hmm model
        if args.debug:
            print "%i missing models [%s]" % (len(notfound), ', '.join(notfound))
        else:
            print "%i missing models" % (len(notfound))
        Consensus = os.path.join(tmpdir, 'missing.consensi.tmp')
        Consensi = os.path.join(tmpdir, 'missing.consensi.fa')
        with open(Consensus, 'w') as output1:
            subprocess.call(['hmmemit', '-c', args.hmm], stdout = output1, stderr = FNULL)
        with open(Consensi, 'w') as output2:
            with open(Consensus, 'rU') as input:
                for rec in SeqIO.parse(input, 'fasta'):
                    rec.id = rec.id.replace('-consensus', '')
                    rec.name = ''
                    rec.description = ''
                    if rec.id in notfound:
                        SeqIO.write(rec, output2, 'fasta')

        #now run tblastn against genome with those notfound
        Blast = os.path.join(tmpdir, 'tblastn.blast.tab')
        print "Try to recover models using tBlastn/Exonerate"
        tblastnFilter(Genome, Consensi, args.cpus, tmpdir)
        print "found %i preliminary tBlastn alignments" % (len(HitList))
        if len(HitList) != 0: 
            #split genome fasta into individual scaffolds
            with open(Genome, 'rU') as input:
                for record in SeqIO.parse(input, "fasta"):
                    if record.id in Scaffolds:
                        SeqIO.write(record, os.path.join(tmpdir, record.id + ".fasta"), "fasta")

            #Now run exonerate on hits
            print "Polishing alignments with Exonerate"
            record_dict = SeqIO.to_dict(SeqIO.parse(Consensi, 'fasta'))
            p = multiprocessing.Pool(args.cpus)
            rs = p.map_async(runExonerate, HitList)
            p.close()
            while (True):
                if (rs.ready()): break

            #now collect all exonerate results into one
            Exonerate = os.path.join(tmpdir, 'exonerate.output.txt')
            skip = ['Command line', '%', ' ','tmp_', '\n', '--', 'Hostname']
            with open(Exonerate, 'w') as output5:
                for root, dirs, files in os.walk(tmpdir):
                    for file in files:
                        if file.endswith('.out'):
                            filename = os.path.join(root, file)
                            with open(filename, 'rU') as readfile:
                                for line in group_by_heading(readfile):
                                    for i in line:
                                        if not any(i.startswith(x) for x in skip):
                                            i = i.replace('^>', '>'+base+'|')
                                            output5.write(i)
            if getSize(Exonerate) > 100:
                exonerate_tmp = os.path.join(tmpdir, base+'.exonerate.output.fa')
                with open(exonerate_tmp, 'w') as output6:
                    with open(Exonerate, 'rU') as input:
                        for rec in SeqIO.parse(input, 'fasta'):
                            ID = rec.id.split('|')
                            ExoID = ID[2]+'|'+ID[3]
                            rec.id = ExoID
                            rec.description = ''
                            rec.name = ''
                            output6.write('>%s\n%s\n' % (rec.id, rec.seq.translate()))
            
                #now check these proteins for matches with HMM models, getting best match
                print "Validating models using HMMsearch"
                ExoResults = runHMMsearch(exonerate_tmp, base+'.exonerate', tmpdir, args.cpus, args.evalue, args.hmm)
                for k,v in natsorted(ExoResults.items()):
                    if not k in Results:
                        Results[k] = v
                #clean up tmp folder
                for root, dirs, files in os.walk(tmpdir):
                    for file in files:
                        if file.endswith('.out') or file.endswith('.fasta'):
                            os.remove(os.path.join(root, file))
            
        else:
            print "No potential hits found"
            
        for i in HMMmodels:
            if not i in Results:
                if i in tBlastN:
                    Results[i] = tBlastN.get(i)
                else:
                    print 'HMM-model '+i+' not found'
                    Results[i] = ['None found', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']

    #all hits found via proteome, so write to output
    print "Saving %i hits to tsv file" % len(Results)
    goodCount = 0
    with open(FinalOut, 'ab') as output:
        for k,v in natsorted(Results.items()):
            if v[-2] != 'NA':
                goodCount += 1
                output.write(">%s\n%s\n" % (v[-2], v[-1]))
    print "Saving %i hits to FASTA file" % goodCount
    
    with open(TextOut, 'ab') as output:
        for k,v in natsorted(Results.items()):
            output.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (base, k, v[0], v[1], v[2], v[3], v[4]))
    SumResults = {}
    for k,v in Results.items():
        SumResults[k] = v[0]            
    AllResults.append(SumResults)

print '----------------------------------------------'
df = pd.DataFrame(AllResults, index=labels)
df.to_csv(SummaryOut)
print 'Summary table saved to %s' % SummaryOut
print '----------------------------------------------'
if args.split_hmms:
    complete = 0
    FNULL = open(os.devnull, 'w')
    concatout = args.out+'.concat.fa'
    concatseq = {}
    print 'Splitting into FASTA file for each HMM model'
    with open(FinalOut, 'ru') as inputfasta:
        for rec in SeqIO.parse(inputfasta, 'fasta'):
            model = rec.id.split('|')[1]
            with open(os.path.join(tmpdir, model+'.proteins.fa'), 'ab') as output:
                SeqIO.write(rec, output, 'fasta')
    for i in natsorted(HMMmodels):
        print "Working on "+i
        seen = []
        try:
            num_models = countfasta(os.path.join(tmpdir,i+'.proteins.fa'))
        except IOError:
            print "Skipping %s, protein fasta file not found, perhaps no models were found" % i
            continue
        if num_models == len(args.input):
            complete += 1
            mafftout = os.path.join(tmpdir, i+'.mafft.fa')
            trimalout = os.path.join(tmpdir, i+'.trimal.fa')
            with open(mafftout, 'w') as output:
                subprocess.call(['mafft', os.path.join(tmpdir,i+'.proteins.fa')], stdout = output, stderr = FNULL)
            trimalError = subprocess.Popen(['trimal', '-in', mafftout, '-out', trimalout, '-automated1', '-keepheader'], stderr=subprocess.STDOUT, stdout=subprocess.PIPE).communicate()[0].rstrip()
            if trimalError != '':  #if you get any warning or error from trimAl, then skip record.
                print trimalError
                print "Skipping %s, at least 1 seq composed only of gaps after trimAl"
                continue
            with open(trimalout, 'rU') as trimal:
                for rec in SeqIO.parse(trimal, 'fasta'):
                    ID = rec.id.split('|')[0]
                    if not ID in seen:
                        if not ID in concatseq:
                            concatseq[ID] = [str(rec.seq)]
                        else:
                            concatseq[ID].append(str(rec.seq))
                        seen.append(ID)
                    else:
                        print "already seen %s" % ID
        else:
            print "Skipping %s, as a single gene model was not found for all inputs (%i/%i found)" % (i, num_models, len(args.input))
    print '----------------------------------------------'
    print "Concatentating %i/%i models per genome into a single file" % (complete, len(HMMmodels))
    with open(concatout, 'w') as finalout:
        for k,v in natsorted(concatseq.items()):
            finalout.write('>%s\n%s\n' % (k, ''.join(v)))
    print '----------------------------------------------'    
if not args.debug:
    shutil.rmtree(tmpdir)
