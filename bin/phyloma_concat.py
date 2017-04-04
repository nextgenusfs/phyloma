#!/usr/bin/env python

#script to take phyloma mapped data and making individual gene trees and concatenated

import sys, os, re, argparse, subprocess, operator, inspect, random
import pandas as pd
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
parser.add_argument('-f','--fasta_external', help='External Sequences, such as ITS, LSU, etc')
parser.add_argument('-g','--genes', nargs='+', help='Genes to use for concatenated tree, others will be dropped')
parser.add_argument('--single_trees', action='store_true', help='Run alignment/RAxML on individual genes')
parser.add_argument('--outgroup', help='Name of species to use as outgroup for concat tree')
parser.add_argument('--bootstrap', default=100, type=int, help='Number of bootstrap replicates')
parser.add_argument('--raxml', default='GTRGAMMA', help='RAxML method')
parser.add_argument('-c','--cpus', default=6, type=int, help='Number of CPUS')
parser.add_argument('--order', default='coverage', choices=['natural', 'coverage', 'random'], help='Order in which genes are selected')
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
    sys.exit(1)

#check if external fasta file is passed, i.e. Sanger data, parse and add to dataset
if args.fasta_external:
    args.input = [item.replace('/', '') for item in args.input]
    #load into dictionary
    with open(args.fasta_external, 'rU') as ExtFasta:
        for rec in SeqIO.parse(ExtFasta, 'fasta'):
            isolate = rec.id.rsplit('_',1)[0]
            geneID = rec.id.rsplit('_',1)[-1]
            match = next((x for x in args.input if isolate in x), None)
            if match:
                outfile = os.path.join(match, geneID+'.fa')
                with open(outfile, 'w') as newfile:
                    newfile.write('>%s\n%s\n' % (geneID, str(rec.seq)))
                #then also need to update coverage map, quickly load into list the genes already in coverage map
                coverage_data = os.path.join(match, 'coverage_stats.txt')
                present = []
                with open(coverage_data) as cov:
                    for line in cov:
                        gene = line.split('\t')[0]
                        if not gene in present:
                            present.append(gene)
                if not geneID in present:                
                    with open(coverage_data, 'ab') as cov:
                        cov.write('%s\t%s\n' % (geneID, '1.0000'))
            else:
                print '%s not found, skipping' % isolate
#create dictionary of the coverages for all isolates
coverage = {}
for folder in args.input:
    folder = folder.replace('/', '')
    if not os.path.isdir(folder):
        print "Error: %s is not a valid folder" % folder
        sys.exit(1)
    #open coverage and parse results
    with open(os.path.join(folder, 'coverage_stats.txt'), 'rU') as input:
        for line in input:
            line = line.replace('\n', '')
            cols = line.split('\t')
            if not cols[0] in coverage:
                coverage[cols[0]] = [float(cols[1])]
            else:
                coverage[cols[0]].append(float(cols[1]))
#convert to pandas dataframe
covdf = pd.DataFrame.from_dict(coverage, orient='index')
covdf.columns = args.input
covdf.fillna(value=0, inplace=True)
covdf = covdf.reindex(index=natsorted(covdf.index))
all_coverage = args.out+'.gene_coverage.csv'
covdf.to_csv(all_coverage)

print "------------------------------------"   
print "%i species loaded" % len(args.input)
print "%i genes in experiment" % len(covdf)

single = args.out+'_singles' #maybe do PID here?
if not os.path.isdir(single):
    os.makedirs(single)

#find the "best" genes to use, i.e. highest average coverage?
#first filter the pandas dataframe for lowest coverage filter, less than min_cov then drop
filtdf = covdf[~(covdf < float(args.min_cov)).any(1)]

#which were dropped and how many
original = set(list(covdf.index.values))
filtered = list(filtdf.index.values)
diff = [x for x in original if x not in filtered]
print "Minimum Coverage is %f, dropping %i genes: %s" % (args.min_cov, len(diff), ', '.join(diff))

#if pass a list of genes to use, grab those here
if args.genes:
    filtdf = filtdf.reindex(index=args.genes)

#now calculate average coverage for each gene
avgdf = filtdf.mean(axis=1)
if args.order == 'coverage':
    avgdf.sort_values(ascending=False, inplace=True)
elif args.order == 'natural':
    avgdf.index = natsorted(avgdf.index)
elif args.order == 'random':
    avgdf = avgdf.sample(frac=1)

#check number to use
if len(avgdf.index) < args.num_concat:
    cutoffnum = len(avgdf.index)
    print "Using all %i genes for concantenation" % len(avgdf.index)
else:
    cutoffnum = args.num_concat
    print "Using top %i genes for concantenation" % args.num_concat

#get list to keep
coverage_output = args.out+'.gene_concat_order.csv'
final_list = list(avgdf.index.values)[:cutoffnum]
finalcov = avgdf.reindex(index=final_list)
finalcov.to_csv(coverage_output)

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

