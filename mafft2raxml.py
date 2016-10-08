#!/usr/bin/env python

import sys, argparse, subprocess, re, multiprocessing, os, inspect
from Bio import AlignIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

parser = argparse.ArgumentParser()
parser.add_argument('-f','--fasta', required=True, help='multi-fasta file')
parser.add_argument('-o','--out', default='out', help='Base name for output files')
parser.add_argument('-m','--raxml_method', default='GTRGAMMA', help='RAxML method')
parser.add_argument('--outgroup', help='Outgroup for RAxML')
parser.add_argument('--bootstrap', default='100', help='Num of Rapid Bootstraps for RAxML')
parser.add_argument('--cpus', default=6, help='Num of threads to use')
args = parser.parse_args()

convertRAxML = os.path.join(currentdir, 'raxml_convert.pl')
newick2nexus = os.path.join(currentdir, 'newick2nexus.py')
FNULL = open(os.devnull, 'w')
def repl(m):
   return 'N' * len(m.group())

def flatten(l):
    flatList = []
    for elem in l:
        # if an element of a list is a list
        # iterate over this list and add elements to flatList
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList

def combinelists(x):
    y = []
    for i in range(len(x)):
        try:
            if x[i+1][0] < x[i][1]:
                if x[i+1][1] < x[i][1]:
                    del x[i+1]
                else:
                    del x[i][1]
                    x[i].append(x[i+1][1])
                    del x[i+1]
        except IndexError:
            continue
    for i in range(len(x)):
        try:
            if x[i+1][0] < x[i][1]:
                if x[i+1][1] < x[i][1]:
                    del x[i+1]
                else:
                    del x[i][1]
                    x[i].append(x[i+1][1])
                    del x[i+1]
        except IndexError:
            continue
    for i in range(len(x)):
        try:
            diff = x[i+1][0] - x[i][1]
            if diff < 3:
                y.append([x[i][0],x[i+1][1]])
            else:
                y.append(x[i])
        except IndexError:
            y.append(x[i])
    for i in range(len(y)):
        try:
            if y[i+1][1] == y[i][1]:
                del y[i+1]
        except IndexError:
            continue
    for i in range(len(y)):
        try:
            if y[i+1][0] < y[i][1]:
                del y[i+1]
        except IndexError:
            continue
    return y

def AlignClean(file, out):
    global countN, total_len
    #search pattern
    match = re.compile(r'(N)\1*')
    #create list to append start/stop to
    Ns = []
    handle = open(file, 'rU')
    outhandle = open(out, 'w')
    alignment = AlignIO.read(handle, 'fasta')
    for rec in alignment:
        total_len = len(rec.seq)
        string = str(rec.seq).upper()
        Seq = re.sub('N[-]*N', repl, string) #replace gaps between N's with N's for the next regex step
        for m in match.finditer(Seq):
            Ns.append( [m.start(),m.end()] )
    Ns.sort(key=lambda x: x[0])
    #now run the combinelist function as many times as necessary
    run1 = combinelists(Ns)
    flat = flatten(run1)
    flat.insert(0,0)
    flat.append(total_len)
    final = zip(*[iter(flat)] * 2)
    test = []
    for i in range(len(final)):
        cmd = "alignment[:, %i:%i]" % (final[i][0], final[i][1])
        test.append(cmd)

    edited = ' + '.join(test)
    AlignIO.write(eval(edited), outhandle, 'fasta')
    handle.close()
    countN = len(run1)

def RunRAxML(input):
    #re-align with mafft
    align2 = args.out + '.mafft2.fa'
    align2_handle = open(align2, 'w')
    subprocess.call(['mafft','--quiet', '--thread', cores, input], stdout = align2_handle)
    align2_handle.close()
    #run trimal
    trimal = args.out + '.trimal.phylip'
    subprocess.call(['trimal', '-in', align2, '-automated1', '-phylip', '-out', trimal])
    raxml_out = args.out + '.raxml.nwk'
    if args.outgroup:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', cores, '-f', 'a', '-m', args.raxml_method, '-p', '12345', '-x', '12345', '-o', args.outgroup, '-#', args.bootstrap,'-s', trimal, '-n', raxml_out], stdout = FNULL, stderr=FNULL)
    else:
        subprocess.call(['raxmlHPC-PTHREADS', '-T', cores, '-f', 'a', '-m', args.raxml_method, '-p', '12345', '-x', '12345', '-#', args.bootstrap,'-s', trimal, '-n', raxml_out], stdout = FNULL, stderr=FNULL)
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

#first step run MAFFT alignment
cores = str(args.cpus)
align = args.out + '.mafft.fa'
with open(align, 'w') as align_handle:
    subprocess.call(['mafft','--quiet', '--thread', cores, args.fasta], stdout = align_handle)

#Now run through Alignment clean-up to remove regions containing N's
clean_out = args.out + '.mafft.clean1.fa'
AlignClean(align, clean_out)
print("Length of alignment: %s\nN's in first pass: %s" % (str(total_len), str(countN)))
if countN > 0:
    clean2 = args.out + '.mafft.clean2.fa'
    AlignClean(clean_out, clean2)
    print("Length of alignment: %s\nN's in second pass: %s" % (str(total_len), str(countN)))
else:
    RunRAxML(clean_out)
    sys.exit(1)

if countN > 0:
    clean3 = args.out + '.mafft.clean3.fa'
    AlignClean(clean2, clean3)
    print("Length of alignment: %s\nN's in third pass: %s" % (str(total_len), str(countN)))
else:
    RunRAxML(clean2)
    sys.exit(1)

if countN > 0:
    clean4 = args.out + '.mafft.clean4.fa'
    AlignClean(clean3, clean4)
    print("Length of alignment: %s\nN's in fourth pass: %s" % (str(total_len), str(countN)))
else:
    RunRAxML(clean3)
    sys.exit(1)

if countN > 0:
    clean5 = args.out + '.mafft.clean5.fa'
    AlignClean(clean4, clean5)
    print("Length of alignment: %s\nN's in fifth pass: %s" % (str(total_len), str(countN)))
else:
    RunRAxML(clean4)
    sys.exit(1)

if countN > 0:
    clean6 = args.out + '.mafft.clean6.fa'
    AlignClean(clean5, clean6)
    print("Length of alignment: %s\nN's in sixth pass: %s" % (str(total_len), str(countN)))
    RunRAxML(clean6)
else:
    RunRAxML(clean5)
    sys.exit(1)

