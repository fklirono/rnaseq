#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  this script filters out the multiple alignments of reads, identified as the alignments with identical read-ids
  and reports ONLY THE FIRST ALIGNMENT ENCOUNTERED

  rmdup.py [-r] [--perfect-match] in.bam out.bam
"""
###################################################################################################

import sys
import pysam
from argparse import *
from collections import defaultdict

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('inbam', default='', type=str, action='store', help='input BAM file')
parser.add_argument('outbam', default='', type=str, action='store', help='output BAM file with filtered duplicates')
parser.add_argument('-r', '--report', default=False, action='store_true', help='should duplicate readIds with their multiplicities be reported to stderr?')
parser.add_argument('--perfect-match', default=False, action='store_true', help='should only perfectly matched reads (NM:i:0) be considered?')
options = parser.parse_args()

inbam=pysam.Samfile(options.inbam, 'rb')
outbam=pysam.Samfile(options.outbam, 'wb', header=inbam.header)

readId = defaultdict(int)  #  dictionary to hold readIDs with default values to zero


#import ipdb;  ipdb.set_trace()


for read in inbam.fetch(until_eof=True):
  #  skip mismatched alignments if only perfect desired
  if options.perfect_match and dict(read.tags).get('NM',0):
    continue
  
  if not read.qname in readId:  #  first time we encounter the read
    outbam.write(read)  #  write the read
  
  readId[read.qname]+=1  #  increase the read multiplicity (for later use)

inbam.close()
outbam.close()

if (options.report):
  for (r,m) in readId.iteritems():
    if (m!=1):
      sys.stderr.write('\t'.join(map(str, [r,m]))+'\n')


