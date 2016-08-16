#!/usr/bin/env python

#Wrapper script for Phyloma package.

import sys, os, subprocess, inspect
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

version = '0.0.2'

default_help = """
Usage:       phyloma <command> <arguments>
version:     %s

Description: Phyloma: Phylogentic Marker-gene Analysis using Reference Genome
    
Command:     map            Map reads to Reference Genome
             draw           Combine mapped data to draw trees
             
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
             -l, --list     List of Consverved Gene models
             -b, --busco    BUSCO OGS results on proteome
             -o, --out      Base Output Name (i.e. isolate name)
             -c, --cpus     Number of CPUs. Default: 6
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
            os._exit(1)
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
             --single_trees    Create single gene trees
            
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
            os._exit(1)
    elif sys.argv[1] == 'version':
        print "phyloma v.%s" % version
    else:
        print "%s option not recognized" % sys.argv[1]
        print default_help
        os._exit(1)
    
else:
    print default_help
