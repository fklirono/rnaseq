#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  this script retains reads from a FASTQ file with readIDs matching/containing the strings imported from a list

  keepIds.py --contains readids in.fastq out.fastq
"""
###################################################################################################

import sys
from argparse import *
from collections import defaultdict
from Bio import SeqIO

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('readids', default='', type=str, action='store', help='readIds to retain in single line per readId format')
parser.add_argument('infastq', default='', type=str, action='store', help='input FASTQ file of original reads')
parser.add_argument('outfastq', default='', type=str, action='store', help='output FASTQ file of the resulting reads')
parser.add_argument('--contains', default=False, action='store_true', help='look for string matches instead of expecting identical readIds')
options = parser.parse_args()


#import ipdb;  ipdb.set_trace()


#  import readIds to a dictionary
rids = defaultdict(int)  #  default values are zero
f=open(options.readids, 'r')
for r in f:
  rids[r.strip('\n')]+=1
f.close()

for read in SeqIO.parse(options.infastq, 'fastq'):
  write_read=False

  if not options.contains:  #  the readsIds list provided contains exact read identifiers
    if rids.get(read.id,-1)>0:
      write_read=True
  else:  #  the rids list provided contains substrings of the read identifiers
    for key in rids.keys():
      if key in read.id:
        write_read=True
        break
  
  if write_read:
    SeqIO.write(read, options.outfastq, 'fastq')


