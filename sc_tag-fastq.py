#!/usr/bin/env python2
# -*- coding: utf-8 -*-
###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################
usage = """
    this script takes from STDIN the uncompressed FASTQ sequences of a single cell sequencing experiment for the 
    mates that contain the transcripts and the called barcodes from a run of 'sc_de-barcode.py' and then:

        1) discards all input reads without all three barcodes called
        2) tags with complete cell barcodes (bc1+bc2+bc3) and UMIs the read_ids outputing them to STDOUT

    N.B. it expects and checks that the reads and the barcode results are in identical order!

    Typical usage:
    
        zcat R1.fastq.gz | %(prog)s R2.fastq.gz.barcodes 
"""
###################################################################################################

import sys,csv
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from argparse import *
from itertools import izip

parser = ArgumentParser(usage=usage)
parser.add_argument('barcodes', default='', type=str, action='store', help='sc_de-barcode.py result file with identified UMIs and barcodes per read_id')
parser.add_argument('--output', default='', type=str, action='store', dest='output', help='output file of tagged sequences')
options = parser.parse_args()


#import ipdb;  ipdb.set_trace()
#from IPython.core.debugger import Pdb; ipdb = Pdb(); ipdb.set_trace()


ntotal = 0
ntagged = 0
with open(options.barcodes, 'r') as BARCODES:
    for read,barcodes in izip(SeqIO.parse(sys.stdin, 'fastq'), csv.reader(BARCODES, delimiter='\t')):

        #  count the read for the total
        ntotal += 1

        #  make sure read order for transcripts and barcodes is identical!
        assert( read.name==barcodes[0] )

        #  check if barcodes are incomplete and skip in this case
        if '' in barcodes[5:8]:  
            continue
        
        #  tag the read_id with cell barcode and UMI
        read_id = '_'.join([ read.id, ''.join(barcodes[5:8]), barcodes[1] ])
        read.id = read_id
        read.description = ' '.join([ read_id, read.description.split(' ')[1] ])  #  both id and description need to be updated
        SeqIO.write(read, sys.stdout, 'fastq')
        ntagged += 1


#  report to STDERR total counts and total number of tagged reads
sys.stderr.write(''.join(['number of reads = ', str(ntotal), '\nnumber of tagged reads = ', str(ntagged), ' (', str(round(float(ntagged)/float(ntotal)*100.0, 1)), '%)\n']))




