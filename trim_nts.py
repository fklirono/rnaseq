#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  cat reads.fastq | %prog [--min 0] n5p n3p > trimmed.fastq
"""

from Bio import SeqIO
from argparse import *
import sys

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('n5p', default=0, type=int, nargs='?', action='store', help="number of nucleotides to trim from the 5' end (default 0)")
parser.add_argument('n3p', default=0, type=int, nargs='?', action='store', help="number of nucleotides to trim from the 3' end (default 0)")
parser.add_argument('--min', default=0, type=int, nargs='?', action='store', help='filter our after trimming reads of length less than this')
parser.add_argument('--nospace', default=False, action='store_true', dest='nospace', help='rename reads by keeping the readID up to first space character?')
options = parser.parse_args()
OUT=sys.stdout
IN=sys.stdin


#import ipdb;  ipdb.set_trace()


for r in SeqIO.parse(IN, 'fastq'):
  l=len(r)
  if options.min!=0 and l-(options.n5p+options.n3p)<options.min:
    continue
  if options.nospace:
    r.description=r.id
  SeqIO.write(r[options.n5p:(l-options.n3p)], OUT, 'fastq')

