#!/usr/bin/env python

import sys, argparse, os, subprocess, shutil, random
from natsort import natsorted
from Bio import SeqIO
import textwrap
import datetime

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

parser=argparse.ArgumentParser(prog='busco4phylogeny.py',
    description='''Script runs BUSCO2 on several GBK files, then parses results, finds all in common, and draws tree using RAxML''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)

parser.add_argument('-i','--input', nargs='+', required=True, help='Input GBK files')
parser.add_argument('-b','--busco', required=True, help='Path to BUSCO2')
parser.add_argument('-d','--busco_db', required=True, help='BUSCO2 DB to use')
parser.add_argument('-o','--out', required=True, help='Output name')
parser.add_argument('-n','--num', type=int, help='Number of BUSCO models to use for tree')
parser.add_argument('-c','--cpus', type=int, default=1, help='Number of CPUs')
parser.add_argument('--dir', help='Previously run directory')
args=parser.parse_args()

def printCMD(cmd):
    stringcmd = '{:}'.format(' '.join(cmd))
    prefix = '\033[96mCMD:\033[00m '
    wrapper = textwrap.TextWrapper(initial_indent=prefix, width=80, subsequent_indent=' '*8, break_long_words=False)
    print(wrapper.fill(stringcmd))

def status(string):
    print('\033[92m[{:}]\033[00m {:}'.format(datetime.datetime.now().strftime('%b %d %I:%M %p'), string))


def softwrap(string, every=80):
    lines = []
    for i in xrange(0, len(string), every):
        lines.append(string[i:i+every])
    return '\n'.join(lines)

def getID(input, type):
    #function to get ID from genbank record.features
    locusTag = None
    ID = None
    Parent = None
    if type == 'gene':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['gene'][0]
            except KeyError:
                pass
        else:
            try:
                ID = input.qualifiers['gene'][0]
            except KeyError:
                pass
        return locusTag, ID, locusTag
        
    elif type == 'mRNA' or type == 'tRNA' or type == 'ncRNA' or type == 'rRNA':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['transcript_id'][0]
                ID = locusTag
            except KeyError:
                pass
            try:
                Parent = input.qualifiers['gene'][0]
            except KeyError:
                pass
        else:
            try:
                ID = input.qualifiers['transcript_id'][0]
            except KeyError:
                pass
        if ID:
            if ':' in ID:
                ID = ID.split(':')[-1]
        return locusTag, ID, Parent
                   
    elif type == 'CDS':
        try:
            locusTag = input.qualifiers['locus_tag'][0]
            Parent = locusTag
        except KeyError:
            pass
        if not locusTag:
            try:
                locusTag = input.qualifiers['protein_id'][0]
            except KeyError:
                pass
            try:
                Parent = input.qualifiers['gene'][0]
            except KeyError:
                pass       
        else:
            try:
                ID = input.qualifiers['protein_id'][0]
            except KeyError:
                pass
        if ID:
            if ':' in ID:
                ID = ID.split(':')[-1]
        return locusTag, ID, Parent

def gb2name(input):
    with open(input, 'rU') as gbk:
        for record in SeqIO.parse(gbk, 'genbank'):
            for f in record.features:
                if f.type == "source":
                    organism = f.qualifiers.get("organism", ["???"])[0]
                    isolate = f.qualifiers.get("isolate", ["???"])[0]
                    if isolate == "???":
                        isolate = f.qualifiers.get("strain", ["???"])[0]
    if organism == '???': #default to file name
        organism = os.path.basename(input).split('.',-1)[0]
    else:
        organism = organism.replace(' ', '_').lower()
    if isolate == '???':
        return organism
    else:  
        return organism+'_'+isolate

def gb2prots(input, tmpdir):
    check = []
    basename = gb2name(input)
    outname = basename+'.prots.fa'
    if not os.path.isdir(os.path.join(tmpdir, 'transcripts')):
    	os.makedirs(os.path.join(tmpdir, 'transcripts'))
    if not os.path.isfile(os.path.join(tmpdir, 'transcripts', basename+'.cds.fa')):
		with open(os.path.join(tmpdir, outname), 'w') as proteins:
			with open(os.path.join(tmpdir, 'transcripts', basename+'.cds.fa'), 'w') as transcripts:
				with open(input, 'rU') as gbk:
					for record in SeqIO.parse(gbk, 'genbank'):
						for f in record.features:
							if f.type == "CDS":
								locusTag, ID, Parent = getID(f, f.type)
								if not ID:
									Name = locusTag
								else:
									Name = ID
								if not Name in check:
									check.append(Name)
								else: #duplicate locus tag which is so stupid, but apparently happens in GenBank files
									Name = Name+'_1'
									if Name in check:
										splitname = Name.split('_')
										num = int(splitname[1])+1
										Name = splitname[0]+'_'+str(num)
								proteins.write(">%s\n%s\n" % (Name, softwrap(f.qualifiers['translation'][0])))
								transcripts.write(">%s\n%s\n" % (Name, softwrap(str(f.extract(record.seq)))))

def parseBUSCO(input):
    #now parse output
    results = {}
    with open(input, 'rU') as busco:
        for line in busco:
            if line.startswith('#'):
                continue
            col = line.split('\t')
            if col[1] == 'Complete':
                ProtID = col[2]
                BuscoHit = col[0]
                score = col[3]
                if not BuscoHit in results:
                    results[BuscoHit] = (ProtID, score)
                else: #don't want any multiple hits, not sure that can ever happen but prevent here
                    del results[BuscoHit]
            elif col[1] == 'Duplicated': #genome might be diploid, so rescue the hit with better alignment score?
                ProtID = col[2]
                BuscoHit = col[0]
                score = col[3]
                if not BuscoHit in results:
                    results[BuscoHit] = (ProtID, score)
                else:
                    old_score = float(results.get(BuscoHit)[1])
                    if float(score) > old_score:
                        results[BuscoHit] = (ProtID, score)         
    return results

'''
first parse the input, if given proteins, i.e. ends with .fa, .fasta, .fna - then just copy to tempfolder
if passed GBK files, then convert to protein fasta files and put in tempfolder
then grab all protein fasta files in tempfolder and run BUSCO protein search
'''
print "----------------------------------"
#setup logfile
logfile = args.out+'.busco4phylogeny.log'
if os.path.isfile(logfile):
    os.remove(logfile)
    
#first parse the arguments
if not args.dir:
    tmpdir = 'phylomaBUSCO_'+str(os.getpid())
    os.makedirs(tmpdir)
else:
    tmpdir = args.dir

status('Loading {:,} species for analysis'.format(len(args.input)))
for i in args.input:
    if i.endswith('.gbk') or i.endswith('.gbf') or i.endswith('.gbff'):
        gb2prots(i, tmpdir)
    else:
        shutil.copy(i, tmpdir)

#now only files in tmpdir are protein fasta files, grab them all
file_list = []
for file in os.listdir(tmpdir):
    if os.path.isfile(os.path.join(tmpdir, file)):
        if not file.startswith('.'):
            file_list.append(file)

#now loop through each and run BUSCO, collect complete results in list of dictionaries
AllResults = []
status('Running iterative BUSCO analysis each proteome')
for x in file_list:
    name = os.path.basename(x).split('.',-1)[0]
    bs_results = os.path.join(tmpdir, 'run_'+name, 'full_table_'+name+'.tsv')
    if not os.path.isfile(bs_results):
        cmd = [os.path.abspath(args.busco), '-i', x, '-m', 'proteins', '-f', '-l', os.path.abspath(args.busco_db), '-o', name, '-c', str(args.cpus)]
        printCMD(cmd)
        with open(logfile, 'a') as log:
            subprocess.call(cmd, cwd=tmpdir, stdout = log, stderr = log)
    BuscoResults = parseBUSCO(bs_results)
    AllResults.append(BuscoResults)

#now get a list of all, just parsing keys here
status('Parsing BUSCO results, determining shared orthologs across all genomes')
AllBuscos = []
for x in AllResults:
    for k,v in x.items():
        if not k in AllBuscos:
            AllBuscos.append(k)

#loop through all buscos and determine if present in every dictionary
BadBuscos = []
for x in AllBuscos:
    for y in AllResults:
        if not x in y:
            BadBuscos.append(x)

#now find those that are found in all results
BadBuscos = set(natsorted(BadBuscos))
AllBuscos = natsorted(AllBuscos)
Keepers = [x for x in AllBuscos if x not in BadBuscos]

if args.num:
    if len(Keepers) <= args.num:
        Keepersfilt = Keepers
    else:
        Keepersfilt = random.sample(Keepers, args.num)
else:
    Keepersfilt = Keepers
print "BUSCO2 Results:"
print "----------------------------------"
#now loop through data and pull out the proteins and cds-transcripts for each
if '.fa' in args.out:
	outputbasename = args.out.rsplit('.fa', 1)[0]
else:
	outputbasename = args.out

protOut = outputbasename+'.proteins.fasta'
cdsOut = outputbasename+'.transcripts.fasta'
with open(protOut, 'w') as output:
	with open(cdsOut, 'w') as output2:
		for i in range(0,len(args.input)):
			SpeciesName = os.path.basename(file_list[i]).rsplit('.', 2)[0]
			Proteins = SeqIO.index(os.path.join(tmpdir, SpeciesName+'.prots.fa'), 'fasta')
			Transcripts = SeqIO.index(os.path.join(tmpdir, 'transcripts', SpeciesName+'.cds.fa'), 'fasta')
			Seq = []
			Seq2 = []
			print "%i BUSCOs found in %s" % (len(AllResults[i]), SpeciesName)
			with open(os.path.join(tmpdir, 'run_'+SpeciesName, SpeciesName+'.buscos.prots.fa'), 'w') as speciesout:
				with open(os.path.join(tmpdir, 'run_'+SpeciesName, SpeciesName+'.buscos.cds.fa'), 'w') as speciesout2:
					for y in AllBuscos:
						if y in AllResults[i]:
							ID = AllResults[i].get(y)[0]
							rec = Proteins[ID]
							rec.id = rec.id+'|'+y
							rec.name = ''
							rec.description = ''
							SeqIO.write(rec, speciesout, 'fasta')
							if y in Keepersfilt:
								Seq.append(str(rec.seq))
							rec = Transcripts[ID]
							rec.id = rec.id+'|'+y
							rec.name = ''
							rec.description = ''
							SeqIO.write(rec, speciesout2, 'fasta')
							if y in Keepersfilt:
								Seq2.append(str(rec.seq))
			output.write('>%s\n%s\n' % (SpeciesName, softwrap(''.join(Seq))))
			output2.write('>%s\n%s\n' % (SpeciesName, softwrap(''.join(Seq2))))

#finalize
print "\nFound %i BUSCOs conserved in all genomes, randomly chose %i" % (len(Keepers), len(Keepersfilt))
print("%s\n" %  ', '.join(Keepersfilt))
print "Concatenated protein sequences for all %i genomes located in: %s" % (len(args.input), protOut)
print "Concatenated transcript sequences for all %i genomes located in: %s" % (len(args.input), cdsOut)
print "----------------------------------"
