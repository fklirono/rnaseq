#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  this script loads FASTA sequences,
  shuffles them at the dinucleotide level using a HMM, thus preserving the dinucleotide frequencies,
  saves the shuffled sequences as FASTA.
  
  There is the posibility to add a prefix to the sequence names by passing the --prefix flag
  
  e.g.

  shuffleDinucleotides.py --prefix random_ infasta outfasta
"""
###################################################################################################

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import seqMotifs as sm
import ghmm
from argparse import *

parser = ArgumentParser(usage=usage)
parser.add_argument('infasta', default='', action='store', help='input FASTA file of sequences we want to shuffle their dinucleotides')
parser.add_argument('outfasta', default='-', action='store', help='output FASTA file of shuffled sequences')
parser.add_argument('--prefix', default='', action='store', help='add this prefix to the sequence ids')
options = parser.parse_args()

fin=open(options.infasta, 'r')

if (options.outfasta!='-'):
  fout=open(options.outfasta, 'w')
else:
  fout=sys.stdout


sigma=ghmm.Alphabet(['A','C','G','T'])


#import ipdb;  ipdb.set_trace()


#  HMM with transition probabilities learned from the dinucleotide frequencies
#
#    A->A, A->C, A->G, A->T
#    C->A, C->C, C->G, C->T
#    G->A, G->C, G->G, G->T
#    T->A, T->C, T->G, T->T
#
#  emission probabilities
#
#    A->A (1.0), A->C (0.0), A->G (0.0), A->T (0.0)
#    C->A (0.0), C->C (1.0), C->G (0.0), C->T (0.0)
#    G->A (0.0), G->C (0.0), G->G (1.0), G->T (0.0)
#    T->A (0.0), T->C (0.0), T->G (0.0), T->T (1.0)
# 
#  and probabilities for the initial states learned from the nucleotide probabilities
#
#    A , C, G, T
#
emission_probs = [ [1.0, 0.0, 0.0, 0.0], 
      [0.0, 1.0, 0.0, 0.0], 
      [0.0, 0.0, 1.0, 0.0], 
      [0.0, 0.0, 0.0, 1.0] ]  #  emission probabilities (A emits A, C emits C, etc...)

for s in SeqIO.parse(fin, 'fasta'):

  #  print empty line if sequence is empty
  if len(str(s.seq))==0:
    fout.write('>'+s.description+'\n\n')  #  SeqIO.write() fails to print the newline, we must do it ourselves
    continue

  #  frequencies of dinucleotides with add-one smoothing
  freqs=sm.dinucleotideFrequency(str(s.seq), add_one=True, as_prob=False)  
  
  #  transition frequencies
  transition_probs = [[freqs['AA'],freqs['AC'],freqs['AG'],freqs['AT']],
       [freqs['CA'],freqs['CC'],freqs['CG'],freqs['CT']],
       [freqs['GA'],freqs['GC'],freqs['GG'],freqs['GT']],
       [freqs['TA'],freqs['TC'],freqs['TG'],freqs['TT']]]  
  
  #  normalize marginals to get transition probabilities
  transition_probs = map( lambda x: [float(y)/float(sum(x)) for y in x], transition_probs )  
  
  #  nucleotide frequencies with add-one smoothing
  freqs_initial=sm.nucleotideFrequency(str(s.seq), add_one=True, as_prob=False)  

  #  normalize initial probabilities
  initial_probs=map( float , [ 1+freqs_initial['A'], 1+freqs_initial['C'], 1+freqs_initial['G'], 1+freqs_initial['T'] ])
  initial_probs=[ x/sum(initial_probs) for x in initial_probs]
  
  #  initialize the HMM
  hmm=ghmm.HMMFromMatrices(sigma, ghmm.DiscreteDistribution(sigma), transition_probs, emission_probs, initial_probs)
  
  #  generate the random sequence
  s.seq=Seq(''.join(map(sigma.external, hmm.sampleSingle(len(s.seq)))))  
  
  if (options.prefix!=''):
    s.id=options.prefix+s.id
    s.name=options.prefix+s.name
    s.description=options.prefix+s.description
  
  SeqIO.write(s, fout, 'fasta')

fin.close()
fout.close()


