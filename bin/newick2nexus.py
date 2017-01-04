#!/usr/bin/env python

import sys
from Bio import Phylo

Phylo.convert(sys.argv[1], 'newick', sys.stdout, 'nexus')