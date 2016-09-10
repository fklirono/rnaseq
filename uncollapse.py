#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  this script reads collapsed BAM aligned reads and prints to STDOUT in SAM each alignment as many times
  as its multiplicity (_x###) found in the readID 

  uncollapse.py [-r] [--perfect-match] in.BAM
"""
###################################################################################################

import sys
import pysam
from argparse import *
from collections import defaultdict

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('inbam', default='', type=str, action='store', help='input BAM file')
parser.add_argument('outsam', default='', type=str, action='store', help='output SAM file')
parser.add_argument('-r', '--report', default=False, action='store_true', help='should the readIds with their multiplicities be reported to stderr?')
parser.add_argument('--perfect-match', default=False, action='store_true', help='should only perfectly matched reads (NM:i:0) be considered?')
options = parser.parse_args()

inbam=pysam.Samfile(options.inbam, 'rb')
outsam=pysam.Samfile(options.outsam, 'w', template=inbam)


#import ipdb;  ipdb.set_trace()


for read in inbam.fetch(until_eof=True):
  #  skip mismatched alignments if only perfect desired
  if options.perfect_match and dict(read.tags).get('NM',0):
    continue

  mult=int(read.qname.split('_x')[1])  #  read multiplicity

  for m in xrange(mult):
    outsam.write(read)
  
  if (options.report):
    sys.stderr.write('\t'.join([read.qname, str(mult)])+'\n')


inbam.close()

