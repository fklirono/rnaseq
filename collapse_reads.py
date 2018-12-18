#!/usr/bin/env python2
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
#
#  this code collapses identical reads to a single one and appends the multiplicity count at the
#  read_id converting essentially also the FASTQ to FASTA 
#
#  N.B. it must be fed GZIPPED FASTA/FASTQ files and it will output GZIPPED FASTA
#
###################################################################################################
usage = """
  %(prog)s reads1.fastq.gz reads2.fastq.gz ...

  Output:

      each input file is processed separately and unique reads are saved as prefix.collapsed.fasta.gz
"""

from argparse import *
from collections import OrderedDict
from Bio import SeqIO, bgzf
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import sys, gzip


parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('reads', metavar='input', type=str, nargs='+', help='space-separated FASTQ files to process')
options = parser.parse_args()


#import ipdb;  ipdb.set_trace()

#sys.stderr.write('\nprocessing the reads...')


for reads_file in options.reads:

    #  output filename
    output = reads_file.split('.')[0]+'.collapsed.fasta.gz'


    #  dictionary to count multiplicities of unique read sequences as key values, preserving the order of appearance
    unique_seqs = OrderedDict()  

    with gzip.open(reads_file, 'rt') as fd:
        for r in SeqIO.parse(fd, 'fastq'):
            seq = str(r.seq)
            _=unique_seqs.setdefault(seq, int(0))  #  it the key is there the count is returned, otherwise the count is set to zero
            unique_seqs[seq]+=1
    

    #  save FASTA records one at a time, let the OS deal with the IO...
    with bgzf.BgzfWriter(output, 'wb') as fd:
        for n,(seq,count) in enumerate(unique_seqs.items()):
            SeqIO.write(SeqRecord(Seq(seq, IUPAC.ambiguous_dna), id='read_'+str(n+1)+'_x'+str(count), description=''), fd, 'fasta')


#sys.stderr.write('done!\n')

