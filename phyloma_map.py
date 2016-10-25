#!/usr/bin/env python

#script to take BAM, GFF, and gene list to output consensus sequences
#dependencies:
#bwa, freebayes, bedtools, bcftools, samtools, tabix, bgzip

import sys, os, re, argparse, subprocess, shutil
from Bio import SeqIO

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=40)

parser = argparse.ArgumentParser(prog='phyloma_map.py', formatter_class = MyFormatter)
parser.add_argument('-r','--reads', nargs='+', help='Input FASTQ reads (single or paired)')
parser.add_argument('-g','--gff', required=True, help='GFF')
parser.add_argument('-s','--sequence', required=True, help='genome fasta file')
parser.add_argument('-b','--busco', help='BUSCO output (pulling SC-orthologs)')
parser.add_argument('-l','--list', help='Gene IDs (ID\tNewName\n)')
parser.add_argument('-o','--out', required=True, help='Output base folder')
parser.add_argument('-c','--cpus', default=6, type=int, help='Number of CPUS')
parser.add_argument('--force', action='store_true', help='use existing folder')
parser.add_argument('--bam', help='Precomputed BAM alignment (sorted)')
parser.add_argument('--variants', help='Precomputed FreeBayes variants (bgzipped)')
args = parser.parse_args()

FNULL = open(os.devnull, 'w')

#check inputs
if not args.variants:
    if not args.bam:
        if not args.reads:
            print "No valid input detected: -r, --bam, or --variants required."
            os._exit(1)

#setup tempdir
tmpdir = args.out
if not os.path.isdir(tmpdir):
    os.makedirs(tmpdir)
else:
    if not args.force:
        print "%s folder already exists, exiting" % tmpdir
        os._exit(1)

#check input, need busco or list but not both
geneDict = {}
if (args.busco and args.list) and not (args.busco or args.list):
    print "Error: must provide either --busco or --list, not both."
    os._exit(1)
else:
    print "------------------------------------"
    if args.busco:
        print "Parsing BUSCO report"
        with open(args.busco, 'rU') as input:
            for line in input:
                line = line.replace('\n', '')
                cols = line.split('\t')
                if cols[1] == 'Complete':
                    if not cols[2] in geneDict:
                        geneDict[cols[2]] = cols[0]
                    else:
                        print "Error: BUSCO mapping multiple models to %s, removing" % cols[2]
                        del geneDict[cols[2]]
    elif args.list:
        print "Parsing gene list file"
        with open(args.list, 'rU') as input:
            for line in input:
                if line.startswith('\n'):
                    continue
                line = line.replace('\n', '')
                if '\t' in line:
                    cols = line.split('\t')
                    geneDict[cols[0]] = cols[1]
                else:
                    geneDict[line] = ''
    else:
        print "Something is very wrong, should never get here....." 

print "Extracting coordinates for %i genes" % len(geneDict)
#get list of genes, check for rename if applicable

#get pattern match from the dictkeys
pattern = re.compile(r'\b(' + '|'.join(geneDict.keys()) + r')\b')

#now list of genes is located in genesDict
#slice GFF file to pull out coordinates
coordinates = os.path.join(tmpdir, 'coordinates.bed')
coordDict = {}
with open(coordinates, 'w') as output3:
    with open(args.gff, 'rU') as gff3:
        for line in gff3:
            line = line.replace('\n', '')
            if '\tgene\t' in line:
                if pattern.search(line):
                    match = pattern.search(line).group(0) #get regex match
                    cols = line.split('\t')
                    scaffold = cols[0]
                    if match in geneDict:
                        ID = geneDict.get(match)
                    else:
                        ID = match     
                    output3.write("%s\t%s\t%s\t%s\n" % (scaffold, cols[3], cols[4], ID))
                    if not ID in coordDict:
                        coordDict[ID] = scaffold+':'+cols[3]+'-'+cols[4]
if not args.variants:
    if not args.bam:
        #map reads to genome
        if len(args.reads) > 1:
            r1 = args.reads[0]
            r2 = args.reads[1]
        elif len(args.reads) > 2:
            print "Error: only one set of paired reads supported"
            os._exit(1)
        else:
            if args.reads[0].endswith('.bam'):
                basename = args.reads[0].split('.bam')[0]
                fastq = os.path.join(tmpdir, basename+'.fastq')
                print "Converting BAM to FASTQ"
                subprocess.call(['bedtools', 'bamtofastq', '-i', args.reads[0], '-fq', fastq])
            else:
                fastq = args.reads[0]
        
        samout = os.path.join(tmpdir, 'mapping.sam')
        bamout = os.path.join(tmpdir, 'mapping.bam')
        bamsort = os.path.join(tmpdir, 'mapping.sort.bam')

        if not os.path.isfile(bamsort):
            print "Building BWA index"
            if not os.path.isfile(os.path.join(tmpdir, 'genome.bwt')):
                subprocess.call(['bwa', 'index', '-p', 'genome', os.path.abspath(args.sequence)], cwd = tmpdir, stderr=FNULL, stdout=FNULL)
            if not os.path.isfile(samout):
                print "Mapping reads to genome using BWA"
                with open(samout, 'w') as output4:
                    if len(args.reads) > 1:
                        subprocess.call(['bwa', 'mem', '-t', str(args.cpus), os.path.join(tmpdir, 'genome'), r1, r2], stdout = output4, stderr=FNULL)
                    else:
                        subprocess.call(['bwa', 'mem', '-t', str(args.cpus), os.path.join(tmpdir, 'genome'), fastq], stdout = output4, stderr=FNULL)
            with open(bamout, 'w') as output5:
                subprocess.call(['samtools', 'view', '-bS', samout], stdout=output5)
            subprocess.call(['samtools', 'sort', '-o', bamsort, bamout], stderr=FNULL)
            subprocess.call(['samtools', 'index', bamsort])
            #cleanup
            os.remove(samout)
            os.remove(bamout)
        else:
            print "BWA mapping file found, skipping mapping"


    else:
        bamsort = os.path.abspath(args.bam)
        subprocess.call(['samtools', 'index', bamsort])
    
    #call variants using freebayes
    print "Calling variants using FreeBayes"
    variants = os.path.join(tmpdir, 'variants.vcf')
    if not os.path.isfile(variants+'.gz'):
        with open(variants, 'w') as output:
            subprocess.call(['freebayes', '-p', '1', '-q', '15', '-f', args.sequence, bamsort], stdout = output)
        #compress and index vcf file
        subprocess.call(['bgzip', '-f', variants])
    else:
        print "Variant calling output detected, using existing data"
    subprocess.call(['bcftools', 'index', '-f', variants+'.gz'])

else:
    if args.bam:
        bamsort = os.path.abspath(args.bam)
        subprocess.call(['samtools', 'index', bamsort])
    else:
        print "Sorted BAM file required (--bam) to map zero coverage regions"
        sys.exit(1)
    print "Pre-computed variants passed, re-indexing"
    variants = os.path.join(tmpdir, 'variants.vcf')
    shutil.copyfile(args.variants, variants+'.gz')
    subprocess.call(['bcftools', 'index', '-f', variants+'.gz'])
    
#get zero coverage information for each region
zeroCoverage = os.path.join(tmpdir, 'nocoverage.bed')
with open(zeroCoverage+'.tmp', 'w') as coverage:
    subprocess.call(['samtools', 'depth', '-aa', '-b', coordinates, bamsort], stdout=coverage)
with open(zeroCoverage, 'w') as coverage_out:
    with open(zeroCoverage+'.tmp', 'rU') as input:
        for line in input:
            line = line.replace('\n', '')
            cols = line.split('\t')
            if int(cols[2]) == 0:
                coverage_out.write('%s\t%i\t%s\n' % (cols[0], int(cols[1])-1, cols[1]))
os.remove(zeroCoverage+'.tmp')

#now get consensus sequence for each gene
print "Extracting consensus sequences for each gene"
for i in coordDict:
    location = coordDict.get(i)
    tmpseq = os.path.join(tmpdir, 'tmp.fa')
    with open(tmpseq, 'w') as temp:
        subprocess.call(['samtools', 'faidx', args.sequence, location], stdout = temp)
    outfasta = os.path.join(tmpdir, i+'.tmp.fa')
    with open(outfasta, 'w') as output5:
        subprocess.call(['bcftools', 'consensus', '-i', '-m', zeroCoverage, '-f', tmpseq, variants+'.gz'], stdout=output5, stderr=FNULL)
    finalfasta = os.path.join(tmpdir, i+'.fa')
    with open(finalfasta, 'w') as output6:
        with open(outfasta, 'rU') as input:
            for line in input:
                if line.startswith('>'):
                    line = '>'+i+'\n'
                output6.write(line)
    os.remove(outfasta)
    os.remove(tmpseq)

#loop through the results, printing out some stats for each locus
print "------------------------------------"
total = 0
good = 0
perfect = 0
missing = 0
missingList = []
coverage_stats = os.path.join(tmpdir, 'coverage_stats.txt')
with open(coverage_stats, 'w') as output:
    for file in os.listdir(tmpdir):
        if file.endswith('.fa'):
            name = file.replace('.fa', '')
            fasta = os.path.join(tmpdir, file)
            with open(fasta, 'rU') as input:
                for rec in SeqIO.parse(input, 'fasta'):
                    length = len(rec.seq)
                    ns = rec.seq.count('N')
                    try:
                        n_coverage = ns / float(length)
                    except ZeroDivisionError:
                        missing += 1
                        missingList.append(name)
                        continue
                    total +=1
                    coverage = 1 - n_coverage
                    if coverage >= 0.75:
                        good += 1
                    if coverage == 1:
                        perfect += 1
                    output.write('%s\t%f\n' % (name, coverage))
                    #print name+'\t'+'{0:.0f}% coverage'.format(coverage*100.0)       
#print "------------------------------------"
print "Found %i of %i total marker-genes" % (total, len(geneDict))
print "%i genes are atleast 75%% covered" % good
print "%i genes are 100%% covered" % perfect
print "%i genes had 0%% coverage: %s" % (missing, ', '.join(missingList))
print "------------------------------------"
