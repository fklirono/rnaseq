#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  this script:

    * uses the REVERSE-COMPLEMENT of the SS-SE sequences of miRs in order to group them as families
    * for each family it finds seed-matches with exactly --mm mismatches in ALL THE targets

  N.B. position counting is one-based
  N.B. overhangs are NOT considered mismatches UNLIKE vmatchPattern() in R
  N.B. multiple matches for the same pair are reported in multiple lines
  N.B. only matches are reported, unlike the `matchPattern.py` script which considers miRNA:target pairs individually

  Output is TSV:
  
    (target seed with no mismatch, i.e. the reverse-complemented SS-SE miRNA sequence)\t
    (target name)\t
    (target seed with --mm mismatches exactly)\t
    (start position relative to target 5')\t
    (end position relative to target 5')\t
"""
###################################################################################################

import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Rfunctions import *
import seqMotifs as sm
from argparse import *

parser = ArgumentParser(usage=usage)
parser.add_argument('--targets', default='targets.fa', action='store', dest='targets', help='FASTA file of the target sequences')
parser.add_argument('--mirs', default='mirs.fa', action='store', dest='mirs', help='FASTA file of the miRs')
parser.add_argument('--output', default='', action='store', dest='output', help='output file to store matches')
parser.add_argument('--mirS', default=2, action='store', type=int, dest='SS', help='starting nucleotide of miRs to use for the family seed')
parser.add_argument('--mirE', default=7, action='store', type=int, dest='SE', help='ending nucleotide of miRs to use for the family seed')
parser.add_argument('--mm', default=0, action='store', type=int, dest='MM', help='exact number of mismatches to look for')
options = parser.parse_args()


assert(options.MM>=0 and options.MM<=(options.SE-options.SS+1))


#  read targets
f=open(options.targets, 'r')
targets=[]
_=[targets.extend([s.id, str(s.seq)]) for s in SeqIO.parse(f, "fasta")]
f.close()


#  read miR sequences, convert to DNA just in case, reverse-complement 
f=open(options.mirs, 'r')
mirs=[str(s.seq[(options.SS-1):options.SE].back_transcribe().reverse_complement()) for s in SeqIO.parse(f, 'fasta')]
f.close()


#  classify miRNA REVERSE-COMPLEMENTED seeds as families
fam={}
for m in mirs:
  if m not in fam:
    fam[m]=1
  else:
    fam[m]+=1


if (options.output==''):
  f=sys.stdout
else:
  f=open(options.output, 'w')


#import ipdb;  ipdb.set_trace()


for m in fam.keys():
  for (j,t) in enumerate(targets[1::2]):
    #  compile all mismatch strings as regular expressions
    regs=[re.compile(mm, re.IGNORECASE) for mm in sm.seqs_with_mms(m, options.MM)]

    #  list of non-empty matches
    matches = unlist([sm.myfinditer(r,t) for r in regs])

    if len(matches)!=0:
      for l in matches:
        start, end =  (str(l[0]+1) , str(l[1]))  #  convert to one-based 
        f.writelines('\t'.join([m, targets[2*j], t[l[0]:l[1]], start, end])+'\n')


f.close()



