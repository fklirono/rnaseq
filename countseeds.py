#!/usr/bin/env python
# -*- coding: utf-8 -*-
###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################
usage = """
    this script takes FASTA sequences of miRNA seeds, groups them by miRNA name, and counts the number of perfect
    matches on the corresponding target FASTA sequences.
    All matches of seeds of the same miRNA group that overlap by at least one nucleotide are counted as single matches.
    
    Output is in TSV format:

        (miRNA name)\t(target name)\t(non-zero counts)
"""
###################################################################################################

import sys,re
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from itertools import izip
from Rfunctions import unlist
import seqMotifs as sm
from argparse import *

parser = ArgumentParser(usage=usage)
parser.add_argument('--targets', default='targets.fa', action='store', dest='targets', help='FASTA file of the named target sequences')
parser.add_argument('--seeds', default='seeds.fa', action='store', dest='seeds', help='FASTA file of the named seed sequences')
parser.add_argument('output', default='counts.tsv', action='store', help='output file where non-zero counts are stored')
options = parser.parse_args()


#  read target sequences
print 'loading the target sequences...'
f=open(options.targets, 'r')
targets=[]
_=[targets.extend([s.id, str(s.seq)]) for s in SeqIO.parse(f, 'fasta')]
f.close()
print 'done!\n'


#  read seed sequences
print 'loading the seed sequences...'
f=open(options.seeds, 'r')
seeds=[]
_=[seeds.extend([s.id, str(s.seq)]) for s in SeqIO.parse(f, 'fasta')]
f.close()
print 'done!\n'


print 'starting the matching...'
f=open(options.output, 'w')


#  capitalize sequences and replace U with T if applicable
seeds[1::2] = [ s.upper().replace('U','T') for s in seeds[1::2] ]
targets[1::2] = [ s.upper().replace('U','T') for s in targets[1::2] ]


#import ipdb;  ipdb.set_trace()
#from IPython.core.debugger import Pdb; ipdb=Pdb(); ipdb.set_trace()


#  make a dictionary with keys the miRNA names and values the corresponding seeds
names=defaultdict(list)
for n,s in izip( seeds[0::2], seeds[1::2] ):
    names[n].append(s)


#  go over all seeds of each miRNA and match them against all targets
for (seeds_name,all_seeds) in names.items():
    
    #  compile seeds as regular expressions (case-sensitive since we have capitalized everything)
    regs=[re.compile(s) for s in all_seeds]  
    

    #  go over all targets and find matches for the given seeds
    for (target_name, target_seq) in izip( targets[0::2], targets[1::2] ):
        
        #  find possible matches
        matches = unlist([sm.myfinditer(r,target_seq) for r in regs]) 
    
        
        #  merge overlapping ranges if we have matches
        if len(matches)!=0:
            matches = sm.reduce_ranges(matches)


            #  write the entry
            f.writelines('\t'.join([seeds_name, target_name, str(len(matches))])+'\n')


f.close()
print 'done!\n'



