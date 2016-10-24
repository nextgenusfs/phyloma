#!/usr/bin/env python

#Wrapper script for Phyloma package.

import sys, os, subprocess, inspect
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

version = '0.0.8'

default_help = """
Usage:       phyloma <command> <arguments>
version:     %s

Description: Phyloma: Phylogentic Marker-gene Analysis
    
Command:     map            Map reads to Reference Genome
             draw           Combine mapped data to draw trees
             hmm            Find HMM models in genomes (GBK format)
             raxml          Construct RAxML phylogeny.
             
Written by Jon Palmer (2016) nextgenusfs@gmail.com
        """ % version

if len(sys.argv) > 1:
    if sys.argv[1] == 'map':
        help = """
Usage:       phyloma %s <arguments>
version:     %s

Description: This scripts takes FASTQ sequencing reads (single or paired) and maps them
             to a reference genome using BWA.  Variants are called using FreeBayes and
             consensus sequence extracting using Samtools/BCFtools.
             Dependencies: BWA, Samtools, BCFtools, FreeBayes, Bedtools, bgzip 
    
Arguments:   -r, --reads    FASTQ sequence reads (Required)
             -g, --gff      Reference GFF3 annotation file (Required)
             -s, --sequence Reference genome FASTA. (Required)
             -l, --list     List of Consverved Gene models. LocusID\tNewName
             -b, --busco    BUSCO OGS results on proteome
             -o, --out      Base Output Name (i.e. isolate name)
             -c, --cpus     Number of CPUs. Default: 6
             --bam          Precomputed BAM alignment (sorted .bam)
             --variants     Precomputed FreeBayes output (vcf.gz)
             --force        Over-write folder
            
Written by Jon Palmer (2016) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'phyloma_map.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)
    elif sys.argv[1] == 'draw':
        help = """
Usage:       phyloma %s <arguments>
version:     %s

Description: This script takes phyloma folders as input, pulls out marker-genes, finds best
             marker-genes (coverage), creates alignments, trim alignments, and infers
             phylogeny using RAxML.
    
Arguments:   -i, --input      `phyloma map` output folders (Required)
             -o, --out         Base name for output (Required)
             -n, --num_concat  Number of marker-genes to concatenate for align/tree. Default: 50
             -m, --min_cov     Minimum coverage of marker-gene. Default: 0.25
             --bootstrap       Number of bootstrap replicates. Default: 100
             --outgroup        Outgroup species for RAxML.
             --single_trees    Create single gene trees
             --raxml           RAxML method. Default: GTRGAMMA
            
Written by Jon Palmer (2016) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'phyloma_concat.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)
    elif sys.argv[1] == 'hmm':
        help = """
Usage:       phyloma %s <arguments>
version:     %s

Description: This scripts takes a single or multiple genomes in GBK format and HMM file
             and searches for each protein model in the genome.  It will first search the 
             proteome, but if model is not found it will search DNA using tblastn/exonerate.
             The ouput is summary stats per genome as well can produce a concatenated multi-
             FASTA file for downstream phylogenies. 
    
Arguments:   -i, --input      Genome in GenBank format, list with spaces or *.gbk
             -m, --hmm        HMM file to search for (can be 1 or many)
             -o, --out        Basename for output files. Default: findHMM
             --evalue         Threshold for HMMScan. Default: 1e-50
             --cpus           Number of CPUs to use. Default: 2
             --split_hmms     Output multi-FASTA for each HMM and concatenated for each genome
             --maxIntron      Maxmimum intron length. Default: 3000
             --debug          Keep intermediate files, extra verbosity to output
           
Written by Jon Palmer (2016) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'findHMM.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)
    elif sys.argv[1] == 'raxml':
        help = """
Usage:       phyloma %s <arguments>
version:     %s

Description: This script is a wrapper for RAxML.
    
Arguments:   -i,--input         Input multi-FASTA file (Required)
             -o,--out           Base name for output. Default: out
             -m,--raxml_method  RAxML method. Default: GTRGAMMA
             --bootstrap        Number of bootstrap replicates. Default: 100
             --outgroup         Outgroup species for RAxML.
             --cpus             Number of CPUs to use. Default: 6

Written by Jon Palmer (2016) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'run_raxml.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)

    elif sys.argv[1] == 'version':
        print "phyloma v.%s" % version
    else:
        print "%s option not recognized" % sys.argv[1]
        print default_help
        sys.exit(1)

else:
    print default_help
