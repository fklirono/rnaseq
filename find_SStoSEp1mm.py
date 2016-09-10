#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  %prog [options]

  this script uses the SS-SE sequences of miRs and tries to find a SS-SE+1mm seed in the corresponding targets

  the output file has the format:
  
    (miR name) (target name) (miR sequence) (target sequence) (SS-SE reverse complemented) (start match 1) (end match 1) 
  
  N.B. position counting starts from 1, not 0 and it is relative to 5' end of target

  N.B. the SS-SE+1mm motif is explored by introducing the mismatch from the last miR nucleotids to first 
       but ONLY the first successful SS-SE+1mm is reported, all others are ignored
"""
###################################################################################################

import sys
#sys.path.append('/Users/fkliron/bio/lib')
import re, copy
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Rfunctions import *
import seqMotifs as sm
from argparse import *

parser = ArgumentParser(usage=usage)
parser.add_argument('--targets', default='targets.fa', action='store', dest='targets', help='FASTA file of the target sequences')
parser.add_argument('--mirs', default='mirs.fa', action='store', dest='mirs', help='FASTA file of the miRs')
parser.add_argument('--output', default='matches.bam', action='store', dest='output', help='output file to store alignment details')
parser.add_argument('--mirS', default=2, action='store', type=int, dest='SS', help='starting nucleotide of miRs to use for the sequence matching')
parser.add_argument('--mirE', default=7, action='store', type=int, dest='SE', help='ending nucleotide of miRs to use for the sequence matching')
options = parser.parse_args()

#  read taget sequences
print 'loading the target sequences...'
f=open(options.targets, 'r')
targets=[]
_=[targets.extend([s.id, str(s.seq)]) for s in SeqIO.parse(f, "fasta")]
f.close()
print 'done!\n'

#  read miR sequences
print 'loading the miR sequences...'
f=open(options.mirs, 'r')
mirs_original=[]
_=[mirs_original.extend([s.id, str(s.seq)]) for s in SeqIO.parse(f, "fasta")]
f.close()
print 'done!\n'

#  copy the originals, replace U with T, keep only SS-SE and reverse complement
mirs=copy.deepcopy(mirs_original)
mirs[1::2] = [ str( Seq(s[(options.SS-1):options.SE].replace('U','T')).reverse_complement() ) for s in mirs[1::2] ]


print 'starting the matching...'
f=open(options.output, "w")
for (n,(t,m)) in enumerate(zip(targets[1::2], mirs[1::2])):
  regs=sm.s_plus_1mm(m)  #  generate all regular expressions of one nt mismatch
  match = [ first for first in [ re.search(r, t) for r in regs ] if not first==None ]  #  return the first match if any
  if len(match)!=0:
    start, end = match[0].span()
    start, end = (str(start+1), str(end))  #  end is defined inclusive here, no need to add one
  else:
    start, end = ('NA', 'NA')
  f.writelines('\t'.join([mirs[2*n], targets[2*n], mirs_original[2*n+1], targets[2*n+1], m, start, end])+'\n')

f.close()
print 'done!\n'


