#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  this script uses the SS-SE sequences of miRs and tries to find REVERSE-COMPLEMENTED matches [--norc to turn off] 
  in ALL THE targets with mismatches according to the --mm flag 

  N.B. position counting is 1-based
  N.B. overhangs are NOT considered mismatches UNLIKE vmatchPattern() in R
  N.B. multiple matches for the same pair are reported in multiple lines
  N.B. only matches are reported, unlike the `matchPattern.py` script which considers miRNA:target pairs individually
  N.B. the REVERSE-COMPLEMENT of the SS-SE sequence of miRs is matched UNLIKE vmatchPattern() in R 
       pass the flag --norc to match the exact sequence

  Output is TSV:
  
    (miR name)\t
    (target name)\t
    (miR original sequence)\t  #  only if --with-seqs
    (target original sequence)\t  #  only if --with-seqs 
    (target seed that matched the --mm criterion)\t
    (start position relative to target 5')\t
    (end position relative to target 5')\t

  Example:

    vmatchPattern.py --targets targets.fa --mirs mature.miRBase20.fa --output 2_7+0mm.tsv --mirS 2 --mirE 7 --mm 0 --norc
"""
###################################################################################################

import sys
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
parser.add_argument('--output', default='matches.tsv', action='store', dest='output', help='output file to store matches')
parser.add_argument('--mirS', default=2, action='store', type=int, dest='SS', help='starting nucleotide of miRs to use for the sequence matching')
parser.add_argument('--mirE', default=7, action='store', type=int, dest='SE', help='ending nucleotide of miRs to use for the sequence matching')
parser.add_argument('--mm', default=0, action='store', type=int, dest='MM', help='exact number of mismatches to look for')
parser.add_argument('--with-seqs', default=False, action='store_true', dest='withseqs', help='should the miRNA and TARGET sequences be included?')
parser.add_argument('--norc', default=False, action='store_true', dest='norc', help='do not look for reverse-complement (seed) matches but identical sequence matches')
options = parser.parse_args()


assert(options.MM>=0 and options.MM<=(options.SE-options.SS+1))


#  read target sequences
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

#  copy the originals, replace U with T, keep only SS-SE and reverse complement if desired
mirs=copy.deepcopy(mirs_original)
if (options.norc):
  mirs[1::2] = [ str( Seq(s[(options.SS-1):options.SE].replace('U','T')) ) for s in mirs[1::2] ]
else:
  mirs[1::2] = [ str( Seq(s[(options.SS-1):options.SE].replace('U','T')).reverse_complement() ) for s in mirs[1::2] ]

print 'starting the matching...'
f=open(options.output, "w")

#import ipdb;  ipdb.set_trace()

for (i,m) in enumerate(mirs[1::2]):
  for (j,t) in enumerate(targets[1::2]):
    #  compile all mismatch strings as regular expressions
    regs=[re.compile(mm, re.IGNORECASE) for mm in sm.seqs_with_mms(m, options.MM)]

    #  list of non-empty matches
    matches = unlist([sm.myfinditer(r,t) for r in regs])

    if len(matches)!=0:
      for l in matches:
        start, end =  (str(l[0]+1) , str(l[1]))  #  convert to one-based 
        if (options.withseqs):
          f.writelines('\t'.join([mirs[2*i], targets[2*j], mirs_original[2*i+1], targets[2*j+1], t[l[0]:l[1]], start, end])+'\n')
        else:
          f.writelines('\t'.join([mirs[2*i], targets[2*j], t[l[0]:l[1]], start, end])+'\n')

f.close()
print 'done!\n'



