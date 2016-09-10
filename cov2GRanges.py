#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  this script expects a per-nucleotide coverage TSV file (usually the output of bedtools coverage -s -d ...)
  of the following BED-like format 
    
    seqnames feature_start feature_end feature_name score feature_strand nucleotide coverage

  and prints to STDOUT a GRanges-like TSV file of the following (1-based inclusive) format

    seqnames strand first_covered_nt_coordinate last_covered_nt_coordinate coverage feature_name

  Example:
    
    zcat TSS_coverage_sense.tsv.gz | cov2GRanges.py - > TSS_coverage_sense.GRanges.tsv.gz
"""
###################################################################################################

import sys
from argparse import *
from copy import deepcopy

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('infile', default='', type=str, action='store', help='input TSV file')
options = parser.parse_args()

if options.infile=='-':
  input_stream=sys.stdin
else:
  input_stream=open(options.infile, 'rU')  #  open as text file, expect any line-terminator


#import ipdb;  ipdb.set_trace()


first=['0',]*8
last=['0',]*8
for line in input_stream:
  line=line.strip('\n').split('\t')

  if len(line)!=8:
    sys.stderr.write('expected 8 fields in TSV input, found '+str(len(line))+'\n')
    sys.exit(123)
  
  if line[7]=='0':  #  ignore uncovered nucleotides, but check first if there is a GRanges() pending to close
    if last[7]!='0':
      #  1-based inclusive as in GRanges()
      sys.stdout.write('\t'.join([first[0], first[5], str(int(first[1])+int(first[6])), str(int(first[1])+int(last[6])), first[7], first[3]]) + '\n')
    first=deepcopy(line)
    last=deepcopy(line)
    continue  

  if first[7]=='0':  #  mark the start position for covered nucleotide of new feature
    first=deepcopy(line)
    last=deepcopy(line)
  else:
    #  quick feature comparison: join chromosome_start_end_strand and compare 
    if '_'.join([line[i] for i in [0,1,2,5]]) == '_'.join([last[i] for i in [0,1,2,5]]):  #  we are still walking over the same feature
      if line[7]==last[7]:  #  coverage has not changed, update end position
        last[6]=line[6]
        continue
      else:  #  coverage has changed on same feature, report the GRanges() with previous coverage and update start position
        #  1-based inclusive as in GRanges()
        sys.stdout.write('\t'.join([first[0], first[5], str(int(first[1])+int(first[6])), str(int(first[1])+int(last[6])), first[7], first[3]]) + '\n')
        first=deepcopy(line)
        last=deepcopy(line)
        continue
    else:  #  we are not walking over the same feature anymore, close GRanges() of previous feature
      #  1-based inclusive as in GRanges()
      sys.stdout.write('\t'.join([first[0], first[5], str(int(first[1])+int(first[6])), str(int(first[1])+int(last[6])), first[7], first[3]]) + '\n')
      first=deepcopy(line)
      last=deepcopy(line)
      continue


if last[7]!='0':  #  make sure to close any open GRanges()
  #  1-based inclusive as in GRanges()
  sys.stdout.write('\t'.join([first[0], first[5], str(int(first[1])+int(first[6])), str(int(first[1])+int(last[6])), first[7], first[3]]) + '\n')


input_stream.close()
