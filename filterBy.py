#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################
usage="""
  this script isolates reads that fulfill a certain length criterion and have a certain 5' nucleotide (soft-clipped or not)
"""
###################################################################################################

import sys
import pysam
from argparse import *
from Rfunctions import reverse_complement

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('inbam', default='', type=str, action='store', help='input BAM file')
parser.add_argument('outbam', default='', type=str, action='store', help='output BAM file containing the filtered reads')
parser.add_argument('length', default=21, type=int, action='store', help='isolate reads of this length (including soft-clipped nucleotides)')
parser.add_argument('nt', default='', type=str, action='store', help="what the 5' nucleotide should be (no matter if it is soft-clipped or not)?")
options = parser.parse_args()


#  prepare input and output
inbam=pysam.Samfile(options.inbam, 'rb')
outbam=pysam.Samfile(options.outbam, 'wb', header=inbam.header)


#  make sure we stick to ACGT even if U is provided
NT=options.nt[0].upper()
NT=NT.replace('U', 'T')


#import ipdb;  ipdb.set_trace()


for read in inbam.fetch(until_eof=True):
    if read.rlen!=options.length:
        continue

    seq=read.seq               #  unclipped sequence
    strand=(read.flag & 0x16)  #  0 for plus, 16 for minus
    
    if strand==16:
        seq=reverse_complement(seq)
     
    if seq[0]==NT:
        outbam.write(read)


#  close files
inbam.close()
outbam.close()



