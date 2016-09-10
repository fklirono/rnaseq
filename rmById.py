#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  this script removes/retains reads from a BAM file with readIDs matching those imported from a list

  rmById.py [-v] readids in.bam out.bam
"""
###################################################################################################

import sys
import pysam
from argparse import *
from collections import defaultdict

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('readids', default='', type=str, action='store', help='readIds to remove/retain in single line per readId format')
parser.add_argument('inbam', default='', type=str, action='store', help='input BAM file of original reads')
parser.add_argument('outbam', default='', type=str, action='store', help='output BAM file of the resulting reads')
parser.add_argument('-v', '--invert', default=False, action='store_true', help='instead of filtering out, should we retain only those reads matching the readIds?')
options = parser.parse_args()


#import ipdb;  ipdb.set_trace()


#  import readIds to a dictionary
readId = defaultdict(int)  #  default values are zero
f=open(options.readids, 'r')
for r in f:
  readId[r.strip('\n')]+=1
f.close()


inbam=pysam.Samfile(options.inbam, 'rb')
outbam=pysam.Samfile(options.outbam, 'wb', header=inbam.header)


for read in inbam.fetch(until_eof=True):
  if not options.invert:
    if not read.qname in readId:  #  write the read if not found in the given list of readIds
      outbam.write(read) 
  else:
    if read.qname in readId:  #  write the read if found in the given list of readIds
      outbam.write(read) 


inbam.close()
outbam.close()
