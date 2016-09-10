#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  this program produces a pileup of aligned reads to STDOUT in the following format:
    
    READ_ID\\tREF_NAME\\tREF_POS\\tSTRAND\\tREAD_BASE_CALL\\tREAD_MULTIPLICITY
"""
###################################################################################################

import sys
import pysam
from argparse import *

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('-c', '--collapsed', default=False, action='store_true', help='do the aligned read IDs contain the multiplicity as a _x# suffix? If not, each read is reported with single multiplicity')
parser.add_argument('sbamfile', default='', type=str, action='store', help='SAM/BAM alignment file')
options = parser.parse_args()

#  try to open alignment file assumming it is SAM format, if it fails try BAM format or let it miserably die
try:
  f=pysam.Samfile(options.sbamfile, 'r')
except ValueError:
  f=pysam.Samfile(options.sbamfile, 'rb')

sys.stdout.write('READ_ID\tREF_NAME\tREF_POS\tSTRAND\tREAD_BASE_CALL\tREAD_MULTIPLICITY\n')
for pcolumn in f.pileup():
  for read in pcolumn.pileups:
    if (read.is_del !=0): 
      continue  #  deleted reference base should not be reported as covered
    
    if (options.collapsed):
      mult=read.alignment.qname.split('_x')[1]
    else:
      mult='1'
    
    sys.stdout.write('\t'.join(map(str, [read.alignment.qname, f.getrname(pcolumn.tid), pcolumn.pos+1, read.alignment.is_reverse and '-' or '+',
      read.alignment.seq[read.qpos], mult]))+'\n')

f.close()
