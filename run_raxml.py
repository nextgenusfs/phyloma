#!/usr/bin/env python

import sys, argparse, subprocess, re, multiprocessing, os, inspect
from Bio import AlignIO
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', required=True, help='multi-fasta file')
parser.add_argument('-o','--out', default='out', help='Base name for output files')
parser.add_argument('-m','--raxml_method', default='GTRGAMMA', help='RAxML method')
parser.add_argument('--outgroup', help='Outgroup for RAxML')
parser.add_argument('--bootstrap', default='100', help='Num of Rapid Bootstraps for RAxML')
parser.add_argument('--cpus', default=6, help='Num of threads to use')
parser.add_argument('--quiet', action='store_true', help='do not print any messages')
args = parser.parse_args()

convertRAxML = os.path.join(currentdir, 'raxml_convert.pl')
newick2nexus = os.path.join(currentdir, 'newick2nexus.py')
FNULL = open(os.devnull, 'w')

def getAlignLength(input):
    lengths = []
    with open(input, 'rU') as input:
        for rec in SeqIO.parse(input, 'fasta'):
            lengths.append(len(rec.seq))
    return lengths

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
if not args.quiet:
    print '----------------------------------------------'
#get cpus
cores = str(args.cpus)

#check alignment lengths
AlignLen = set(getAlignLength(args.input))
if len(AlignLen) > 1:
    if not args.quiet:
        print 'Different sequence lengths identified in alignment, re-running MAFFT/trimAl'
    #re-align with mafft
    align = args.out + '.mafft.fa'
    with open(align, 'w') as align_handle:
        subprocess.call(['mafft', '--quiet', '--thread', cores, input], stdout = align_handle)
else:
    align = args.input

#run trimal with automated ML settings
trimal = args.out + '.trimal.fa'
subprocess.call(['trimal', '-in', align, '-automated1', '-out', trimal])
if not args.quiet:
    print 'Running RAxML: %s method, %s bootstraps, %s CPUs' % (args.raxml_method, args.bootstrap, str(args.cpus))
RunRAxML(trimal)

#parse logfile and get some summary data
empiricalmatrix = ''
with open(args.out+'.raxml.log', 'rU') as logfile:
    for line in logfile:
        line = line.replace('\n', '')
        if line.startswith('This is RAxML version'):
            version = line.split(' released')[0]
            version = version.split('version ')[1]
        if line.startswith('Proportion of gaps'):
            gaps = line.split(': ')[1]
        if line.startswith('Alignment Patterns:'):
            patterns = line.split(': ')[1]
        if line.startswith('Substitution Matrix:'):
            submatrix = line.split(': ')[1]
        if 'best-scoring AA model:' in line:
            empiricalmatrix = line.split(' likelihood')[0]
            empiricalmatrix = empiricalmatrix.split(': ')[-1]

if not args.quiet:
    print 'RAxML version: %s' % version
    print 'alignment length: %i' % getAlignLength(trimal)[0]
    print 'characters with gaps: %s' % gaps
    print 'alignment patterns: %s' % patterns
    if empiricalmatrix == '':
        print 'substitution model: %s' % submatrix
    else:
        print 'substitution model: %s' % empiricalmatrix
    print '----------------------------------------------'
    print 'Newick tree saved as %s' % args.out+'.nwk'
    print 'Nexus tree saved as %s' % args.out+'.nexus'
    print '----------------------------------------------'
