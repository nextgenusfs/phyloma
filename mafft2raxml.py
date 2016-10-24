#!/usr/bin/env python

import sys, argparse, subprocess, re, multiprocessing, os, inspect
from Bio import AlignIO
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

parser = argparse.ArgumentParser()
parser.add_argument('-f','--fasta', required=True, help='multi-fasta file')
parser.add_argument('-o','--out', default='out', help='Base name for output files')
parser.add_argument('-m','--raxml_method', default='GTRGAMMA', help='RAxML method')
parser.add_argument('--outgroup', help='Outgroup for RAxML')
parser.add_argument('--bootstrap', default='100', help='Num of Rapid Bootstraps for RAxML')
parser.add_argument('--cpus', default=6, help='Num of threads to use')
parser.add_argument('--skip_align', action='store_true', help='Skip alignment (seqs already aligned in FASTA format)')
args = parser.parse_args()

convertRAxML = os.path.join(currentdir, 'raxml_convert.pl')
newick2nexus = os.path.join(currentdir, 'newick2nexus.py')
FNULL = open(os.devnull, 'w')

def CleanAmbiguous(input, output):
    #make temporary file to replace nonGATC with gaps
    with open('mafft2raxml.tmp', 'w') as tmpout:
        with open(input, 'rU') as sequences:
            for rec in SeqIO.parse(sequences, 'fasta'):
                #replace non-GATC characters with dashes
                Seq = re.sub('[^GATC]', '-', str(rec.seq).upper())
                tmpout.write('>%s\n%s\n' % (rec.id, Seq))
    #now run trimal to remove all gaps
    subprocess.call(['trimal', '-in', 'mafft2raxml.tmp', '-nogaps', '-out', output, '-automated1'])
    #cleanup
    os.remove('mafft2raxml.tmp')

def NstoGaps(input, output):
    #convert N's to gaps as they are regions of low coverage
    with open('mafft2raxml.tmp', 'w') as tmpout:
        with open(input, 'rU') as filein:
            for rec in SeqIO.parse(filein, 'fasta'):
                Seq = str(rec.seq).upper()
                Seq = Seq.replace('N', '-')
                tmpout.write('>%s\n%s\n' % (rec.id, Seq))
    #now run trimal
    subprocess.call(['trimal', '-in', 'mafft2raxml.tmp', '-out', output, '-automated1'])
    #cleanup
    os.remove('mafft2raxml.tmp')

def RunRAxML(input):
    raxml_out = args.out + '.raxml.nwk'
    if args.outgroup:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', cores, '-f', 'a', '-m', args.raxml_method, '-p', '12345', '-x', '12345', '-o', args.outgroup, '-#', args.bootstrap,'-s', input, '-n', raxml_out], stdout = FNULL, stderr=FNULL)
    else:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', cores, '-f', 'a', '-m', args.raxml_method, '-p', '12345', '-x', '12345', '-#', args.bootstrap,'-s', input, '-n', raxml_out], stdout = FNULL, stderr=FNULL)
    #convert to proper newick format
    with open(args.out+'.nwk', 'w') as output:
        subprocess.call([convertRAxML, 'RAxML_bipartitionsBranchLabels.'+raxml_out], stdout = output)
    #convert to nexus format
    with open(args.out+'.nexus', 'w') as output2:
        subprocess.call([newick2nexus, args.out+'.nwk'], stdout = output2)
    for file in os.listdir("."):
        if file.startswith("RAxML_info"):
            os.rename(file, args.out+'.raxml.log')
        if file.startswith('RAxML'):
            try:
                os.remove(file)
            except:
                pass

cores = str(args.cpus)

if not args.skip_align:
    #first step run MAFFT alignment
    align = args.out + '.mafft.fa'
    with open(align, 'w') as align_handle:
        subprocess.call(['mafft','--quiet', '--thread', cores, args.fasta], stdout = align_handle)
else:
    align = args.fasta

#remove Ns as they are low coverage
clean = args.out + '.trimal.fa'
#CleanAmbiguous(align, clean)
NstoGaps(align, clean)

#run Raxml
RunRAxML(clean)