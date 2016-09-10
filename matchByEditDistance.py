#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  this script accepts a probe sequence [--seq] of given length 
  and looks for matches less than K-edit distance --distance away in a list of target sequences of identical length 
  provided by --targets

  when a target sequence has a number of mismatches <= --distance it is printed out in stdout

  sequences of different length than the probe are ignored

  if --strict is defined then targets with exactly --distance number of mismatches are printed out
  if --N_is_mismatch is defined then each N is considered a mismatch, rather than a match
"""
###################################################################################################

import sys
from argparse import *

parser = ArgumentParser(usage=usage)
parser.add_argument('--seq', default='', type=str, action='store', nargs='?', dest='seq', help='sequence to look for in targets')
parser.add_argument('--targets', default='', type=str, action='store', nargs='?', dest='targets',  help='target sequences to process')
parser.add_argument('--mm', default=0, action='store', type=int, dest='mm', help='edit distance')
parser.add_argument('--strict', default=False, action='store_true', dest='strict', help='look at exactly (not less than) number of mismatches == --mm?')
parser.add_argument('--N_is_mismatch', default=False, action='store_true', dest='Nmis', help='should N be considered a mismatch?')
options = parser.parse_args()

assert(options.mm>=0)

L=len(options.seq)  #  we are not going to look at sequences of different length

#import ipdb;  ipdb.set_trace()

with open(options.targets, 'r') as f:
  for line in f:
    l=line.strip()   

    if (len(l)!=L):
      continue

    matches=0
    for i in range(L):
      if ( (options.Nmis and l[i]=='N') or l[i]==options.seq[i] ):
        matches=matches+1

    if ( (options.strict and (L-matches)==options.mm) or (not options.strict and (L-matches)<=options.mm) ):
      sys.stdout.write(line)

