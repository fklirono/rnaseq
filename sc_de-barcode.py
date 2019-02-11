#!/usr/bin/env python2
# -*- coding: utf-8 -*-
###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################
usage = """
    this script takes from STDIN the uncompressed FASTQ sequences of a single cell sequencing experiment for 
    the mates that contain the UMIs, the barcodes and the fillers and identifies each one with either exact 
    matches or up to a single mismatch.

    It is also possible to define barcode pairs, so if a certain barcode mate is found in particular place it 
    should be considered as if the paired mate barcode was found. This can happen if barcode pairs are used in 
    wells during a barcoding run of SplitSeq.

    The sequence order of the different elements is:

        UMI + barcode1 + filler1 + barcode2 + filler2 + barcode3

    The output to STDOUT is in TSV format:

        (read_id)\t(actual UMI)\t(actual bc1)\t(actual bc2)\t(actual bc3)\t(assigned bc1)\t(assigned bc2)\t(assigned bc3)

    Typical run:

      zcat mate2.fastq.gz | %(prog)s --barcodes bc1.fa,bc2.fa,bc3.fa --pairs ,,bc3_equivalent.fa --lengths 10,8,30,8,30,8 --with-mismatch --discard-uncalled > read_ids.tsv 2> log.txt
"""
###################################################################################################

import sys,re
import seqMotifs as sm
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from argparse import *
from numpy import cumsum

parser = ArgumentParser(usage=usage)
parser.add_argument('--barcodes', default='', action='store', dest='barcodes', help='comma-separated FASTA files that contain the barcode sequences')
parser.add_argument('--pairs', default='', action='store', dest='pairs', help='comma-separated FASTA files that define the corresponding paired barcode sequences. Order of the files matters.')
parser.add_argument('--lengths', default='', action='store', dest='lengths', help='comma-separated list of 6 lenghts for each of the elements on the sequence including the UMIs and the fillers: (UMI length) (barcode1 length) (filler1 length) (barcode2 length) (filler2 length) (barcode3 length)')
parser.add_argument('--with-mismatch', default=False, action='store_true', help='allow in barcode sequences a single mismatch?')
parser.add_argument('--discard-uncalled', default=False, action='store_true', help='discard reads where not all barcodes have been identified?')
options = parser.parse_args()


#parser = ArgumentParser(prog='sc_de-barcode.py')
#options = parser.parse_args('--barcodes bc1.fa,bc2.fa,bc3.fa --pairs ,,bc3_equivalent.fa --lengths 10,8,30,8,30,8'.split())


assert( options.lengths != '' and options.barcodes != '' )


#  import barcode sequences into a dictionary and make sure they are capitalized
barcodes = defaultdict(list)
for n,barcode_file in enumerate(options.barcodes.split(',')):
    f = open(barcode_file, 'r')
    bc = []
    _ = [bc.append(str(s.seq).upper()) for s in SeqIO.parse(f, 'fasta')]
    f.close()
    barcodes[n] = bc


#  import barcode pairs if defined
if options.pairs != '':
    pairs = defaultdict(list)
    for n,pair_file in enumerate(options.pairs.split(',')):
        #  import particular paired barcode only if defined
        if pair_file =='':
            continue
        f = open(pair_file, 'r')
        p = []
        _ = [p.append(str(s.seq).upper()) for s in SeqIO.parse(f, 'fasta')]
        f.close()
        pairs[n] = p


#import ipdb;  ipdb.set_trace()
#from IPython.core.debugger import Pdb; ipdb = Pdb(); ipdb.set_trace()


#  convert lengths to integers and compute the cumulative sum that will indicate the breakpoint positions
POS = cumsum(map(int, options.lengths.split(','))).tolist()


ntotal = 0
ncalled = 0
for read in SeqIO.parse(sys.stdin, 'fastq'):
    ntotal += 1
    bc_guess = ['']*3
    bc = ['']*3


    #  isolate the read sequence
    s = str(read.seq)


    #  identify the different elements by their declared lengths
    umi = s[0:POS[0]]
    bc[0] = s[POS[0]:POS[1]]
    #fill1 = s[POS[1]:POS[2]]
    bc[1] = s[POS[2]:POS[3]]
    #fill2 = s[POS[3]:POS[4]]
    bc[2] = s[POS[4]:POS[5]]

    
    #  check each barcode for a perfect match
    for p in range(len(bc)):
        if bc[p] in barcodes[p]:
            bc_guess[p] = barcodes[p][ barcodes[p].index(bc[p]) ] 

    
    #  check on unassigned barcodes if single mismatches are wished 
    #  join with commas all unique matches one Hamming distance away into a single string
    if options.with_mismatch:
        for p in range(len(bc)):
            if bc_guess[p] == '':
                bc_guess[p] = ','.join(set(sm.find_all_1mm(bc[p], barcodes[p])))

    
    #  check for equivalent barcode pairs if they are defined
    if options.pairs != '':
        for p in pairs.keys():
            #  identical matches
            if bc_guess[p] == '':
                if bc[p] in pairs[p]:  
                    bc_guess[p] = barcodes[p][ pairs[p].index(bc[p]) ]  #  replace
                elif options.with_mismatch:
                    #  replace with comma-separated matches 
                    bc_guess[p] = ','.join(set(sm.find_all_1mm(bc[p], pairs[p])))


    #  if at this point we have any barcode uncalled and we want to keep only fully recognized barcodes this it the time to check
    if options.discard_uncalled and '' in bc_guess:
        continue


    #  count this read into the pool of those with all barcodes called if applicable
    if not '' in bc_guess:
        ncalled +=1


    #  write everything to the output
    sys.stdout.write('\t'.join([read.id, umi, bc[0], bc[1], bc[2], bc_guess[0], bc_guess[1], bc_guess[2]])+'\n')


#  report to STDERR total counts and total number of kept reads
sys.stderr.write(''.join(['Number of reads = ', str(ntotal), '\nnumber of reads with all barcodes identified = ', str(ncalled), ' (', str(round(float(ncalled)/float(ntotal)*100.0, 1)), '%)\n']))




