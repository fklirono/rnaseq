#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  this script generates all possible pentamers HHHHH (H={A,C,T}) and counts the amount of times they are 
  found in the target sequences imported.

  The pentamers are printed in the first line of STDOUT and their corresponding frequencies in the second line.
"""
###################################################################################################

import sys
import re
import random
from Bio import SeqIO
from Rfunctions import *
import seqMotifs as sm
import itertools
from argparse import *

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('targets', default='targets.fa', type=str, action='store', help='FASTA sequences of targets')
options = parser.parse_args()


#  read target sequences
f=open(options.targets, 'r')
targets=[]
_=[targets.extend([s.id, str(s.seq)]) for s in SeqIO.parse(f, "fasta")]
f.close()


#  generate all HHHHH
H=['A','C','T']
p=[''.join(r) for r in itertools.product(H,H,H,H,H)]
regs=[ re.compile(r, re.IGNORECASE) for r in p ]  #  compile motifs to regular expressions

sys.stdout.write('\t'.join(p)+'\n')

matches=[0]*len(p)  #  matches per motif


#import ipdb;  ipdb.set_trace()


for (m,motif) in enumerate(p):
  matches[m]=0
  for t in targets[1::2]:
    matches[m]+=len(sm.myfindall(regs[m],t))  #  count matches including overlaps


sys.stdout.write('\t'.join(map(str, matches))+'\n')



