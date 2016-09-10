#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  this script filters out alignments from a BAM file according to demand 

  originally it was constructed for single-end reads, but its functionality has been extended to paired-end reads as well
"""
###################################################################################################

import sys
import pysam
from argparse import *
from collections import defaultdict

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('inbam', default='', type=str, action='store', help='input BAM file')
parser.add_argument('outbam', default='', type=str, action='store', help='output BAM file with filtered reads')
parser.add_argument('-r', '--report', default=False, action='store_true', help='should discarded reads be reported in STDERR?')
parser.add_argument('-m', '--min-mismatch', default=0, type=int, action='store', help='allow this minimum number of mismatches')
parser.add_argument('-M', '--max-mismatch', default=-1, type=int, action='store', help='allow this maximum number of mismatches')
parser.add_argument('-w', '--min-length', default=0, type=int, action='store', help='allow this minimum read length')
parser.add_argument('-W', '--max-length', default=0, type=int, action='store', help='allow this maximum read length. In order for this parameter to have an effect it needs to be non-zero and larger or equal than --min-length')
parser.add_argument('-d', '--with-deletions', default=False, action='store_true', help='allow for reference deletions')
parser.add_argument('-i', '--with-insertions', default=False, action='store_true', help='allow for reference insertions')
parser.add_argument('--no-chimeras', default=False, action='store_true', help='(only for paired-end) filter out chimeras where mates map to different chromosomes')
parser.add_argument('--no-splicing', default=False, action='store_true', help='filter out spliced alignments which report skipped reference bases (CIGAR: N operation)')
parser.add_argument('-f', '--flags', default=[], nargs='+', type=int, action='store', help='space-separated list of SAM FLAGs (followed by -- if it is that last flag passed) to accept, useful for paired-end alignments to pick (99, 147) or (83, 163) proper pairs and filter out everything else ')
options = parser.parse_args()

options.flags=map(long, options.flags)  #  convert to list of long integers (no error if empty)

inbam=pysam.Samfile(options.inbam, 'rb')
outbam=pysam.Samfile(options.outbam, 'wb', header=inbam.header)


if options.report:
  removedbam=pysam.Samfile('/dev/stderr', 'wb', header=inbam.header)


#import ipdb;  ipdb.set_trace()


for read in inbam.fetch(until_eof=True):
  tags = dict(read.tags)
  nm = tags.get('NM', 0)
  #  either no upper bound (max-mismatch<0) or there is an upper bound and we are below or up to it
  mismatch = ( nm>=options.min_mismatch and ( options.max_mismatch<0 or (options.max_mismatch>=0 and nm<=options.max_mismatch) ) )
  flags = len(options.flags)==0 or (read.flag in options.flags)
  not_chimera = True
  if options.no_chimeras:
    not_chimera = (read.rnext == -1) or (read.rnext == read.rname)  #  no need to convert to reference names: inbam.getrname(read.rnext)
  not_spliced = True
  if options.no_splicing:
    not_spliced = not ( 'XS' in tags )
  
  #  read.qlen = aligned nts excluding soft-clipped
  #  read.rlen = actual read length including soft-clipped nts
  #  read.alen = reference nts aligned with this read = (read.qlen + deletions, read.qlen - insertions)
  width=( read.qlen>=options.min_length )
  if options.max_length!=0: 
    width=( width and options.max_length>=options.min_length and read.qlen<=options.max_length )
  
  deletion=('^' in tags.get('MD', ''))
  insertion=('I' in read.cigarstring)
  
  if flags:
    if not_spliced:
      if not_chimera:
        if mismatch:
          if options.with_deletions or not deletion:
            if options.with_insertions or not insertion:
              if width:
                outbam.write(read)
              elif options.report:
                removedbam.write(read)  #  outside width range
            elif options.report:
              removedbam.write(read)  #  insertions not welcome and this is an insertion
          elif options.report:
            removedbam.write(read)  #  deletions not welcome and this is a deletion
        elif options.report:
          removedbam.write(read)   #  outside mismatch range
      elif options.report:
       removedbam.write(read)   #  chimeras not welcome and this is a chimera
    elif options.report:
      removedbam.write(read)   #  spliced read because the alignments skips reference bases
  elif options.report:
    removedbam.write(read)   #  FLAG is not part of our desired FLAG list


inbam.close()
outbam.close()
if options.report:
  removedbam.close()



