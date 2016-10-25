#!/usr/bin/env python

#script to take phyloma mapped data and making individual gene trees and concatenated

import sys, os, re, argparse, subprocess, operator, inspect
from natsort import natsorted
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=40)

parser = argparse.ArgumentParser(prog='phyloma_concat.py', formatter_class = MyFormatter)
parser.add_argument('-i','--input', nargs='+', required=True, help='Input Folders from phyloma_map)')
parser.add_argument('-n','--num_concat', default=50, type=int, help='Number of genes to use')
parser.add_argument('-m','--min_cov', default=0.25, type=float, help='Minimum coverage')
parser.add_argument('-o','--out', required=True, help='Output base folder')
parser.add_argument('--single_trees', action='store_true', help='Run alignment/RAxML on individual genes')
parser.add_argument('--outgroup', help='Name of species to use as outgroup for concat tree')
parser.add_argument('--bootstrap', default=100, type=int, help='Number of bootstrap replicates')
parser.add_argument('--raxml', default='GTRGAMMA', help='RAxML method')
parser.add_argument('-c','--cpus', default=6, type=int, help='Number of CPUS')
args = parser.parse_args()

def runMAFFT(input, output):
    with open(output, 'w') as outfile:
        subprocess.call(['mafft','--quiet', '--thread', str(args.cpus), input], stdout = outfile)

def NstoGaps(input, output):
    #convert N's to gaps as they are regions of low coverage
    with open(output, 'w') as tmpout:
        with open(input, 'rU') as filein:
            for rec in SeqIO.parse(filein, 'fasta'):
                Seq = str(rec.seq).upper()
                Seq = Seq.replace('N', '-')
                tmpout.write('>%s\n%s\n' % (rec.id, Seq))

#check input
if len(args.input) < 4:
    print "Error: need more than %i species to run ML phylogeny" % len(args.input)
    os._exit(1)
coverage = {}
for folder in args.input:
    folder = folder.replace('/', '')
    if not os.path.isdir(folder):
        print "Error: %s is not a valid folder" % folder
        os._exit(1)
    #open coverage and parse results
    with open(os.path.join(folder, 'coverage_stats.txt'), 'rU') as input:
        for line in input:
            line = line.replace('\n', '')
            cols = line.split('\t')
            if not cols[0] in coverage:
                coverage[cols[0]] = [float(cols[1])]
            else:
                coverage[cols[0]].append(float(cols[1]))
print "------------------------------------"   
print "%i species loaded" % len(args.input)

single = args.out+'_singles' #maybe do PID here?
if not os.path.isdir(single):
    os.makedirs(single)

#find the "best" genes to use, i.e. highest average coverage?
results = {}
skipped = []
for k,v in coverage.items():
    if len(v) == len(args.input): #make sure that each gene is covered to be included
        Avg = sum(v)/float(len(v))
        results[k] = (Avg, min(v))
    else:
        skipping = k+" ("+str(len(v))+"/"+str(len(args.input))+")"
        skipped.append(skipping)

print "Skipping %i genes: %s" % (len(skipped), ', '.join(skipped))

sorted_results = sorted(results.items(), reverse=True, key=operator.itemgetter(1))

if len(results) < args.num_concat:
    cutoffnum = len(results)
else:
    cutoffnum = args.num_concat
coverage_output = args.out+'.gene_stats.csv'
final_list = []
with open(coverage_output, 'w') as cov_out:
    cov_out.write('GeneId,AvgCov,LowestCov\n')
    for i in sorted_results:
        cov_out.write('%s,%s,%s\n' % (i[0], i[1][0], i[1][1]))
        if i[1][1] < args.min_cov: #less than 25% coverage for minimum, move on?
            print "%s has at least 1 species less than %f coverage, skipping" % (i[0], args.min_cov)
            continue
        final_list.append(i[0])
        if len(final_list) == cutoffnum:
            print "Found best %i genes" % cutoffnum
            break

#sort alphabetically gene names
final_list = natsorted(final_list)

#now get each of these files
concat = os.path.join(args.out+'_concat.fa')
with open(concat, 'w') as concat_out:
    for x in final_list:
        single = os.path.join(args.out+'_singles', x+'.fa')
        with open(single, 'w') as output:
            for y in args.input:
                y = y.replace('/', '') #make sure no slashes
                for rec in SeqIO.parse(open(os.path.join(y,x+'.fa')), 'fasta'):
                    rec.id = y
                    rec.name = ''
                    rec.description = ''
                    SeqIO.write(rec, output, 'fasta')
    for z in args.input:
        z = z.replace('/', '') #make sure no slashes
        concat_out.write('>%s\n' % z)
        for w in final_list:
            for rec in SeqIO.parse(open(os.path.join(z, w+'.fa')), 'fasta'):
                concat_out.write(str(rec.seq))
        concat_out.write('\n')
        
align = os.path.join(currentdir, 'run_raxml.py')
if args.single_trees:
    if not os.path.isdir(args.out+'_align'):
        os.makedirs(args.out+'_align')
    if not os.path.isdir(args.out+'_trees'):
        os.makedirs(args.out+'_trees')
    for i in final_list:
        print "------------------------------------"
        print "Working on %s" % i
        seq = os.path.join(args.out+'_singles', i+'.fa')
        mafft_align = os.path.join(args.out+'_align', i+'.mafft.fa')
        clean_align = os.path.join(args.out+'_align', i+'.cleaned.fa')
        runMAFFT(os.path.abspath(seq), mafft_align)
        NstoGaps(mafft_align, clean_align)
        #then run trimal/raxml
        subprocess.call([sys.executable, align, '-i', os.path.abspath(clean_align), '-o', i, '--bootstrap', str(args.bootstrap), '-m', args.raxml, '--cpus', str(args.cpus), '--quiet'], cwd = args.out+'_align')
        os.rename(os.path.join(args.out+'_align', i+'.nwk'), os.path.join(args.out+'_trees', i+'.nwk'))
print "------------------------------------"
print "Running MAFFT alignment using %i cores" % args.cpus
concat_align = args.out + '.mafft.fa'
concat_clean = args.out + '.cleaned.fa'
runMAFFT(concat, concat_align)
NstoGaps(concat_align, concat_clean)
      
print "Running trimAl/RAxML on concatentated fasta using %i cores" % args.cpus
if not args.outgroup:
    subprocess.call([sys.executable, align, '-i', os.path.abspath(concat_clean), '-o', args.out+'_concat', '--bootstrap', str(args.bootstrap), '-m', args.raxml, '--cpus', str(args.cpus)])
else:
    subprocess.call([sys.executable, align, '-i', os.path.abspath(concat_clean), '-o', args.out+'_concat', '--bootstrap', str(args.bootstrap), '--outgroup', args.outgroup, '-m', args.raxml, '--cpus', str(args.cpus)])

