#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  This script counts the nucleotide occurrence in alignments for a given genomic range of nucleotides of interest.
  Nucleotide counts are split in forward (+) and reverse (-) strand based on the alignment orientation.
  Both forward and reverse stand use the same nucleotide index, i.e. coordinates are always on the forward strand
  just like in GRanges() in R.
  
  The genomic region defines the nucleotides but alignments are treated as edited-reads generating from the region and 
  soft-clipped nucleotides that overlap the region are included in the counts so that the edits can be tabulated.
  
  N.B. we are counting nucleotide occurrences on the alignments, not the reference, so that soft-clipped
       parts of alignments are included in the tabs if they fall within the genomic range
  
  Return value is a tab-separated array with strand-split nucleotide occurences for the nucleotides of the genomic region:

       N1  N2  N3  N4 ...
   A+  0   0   0   0
   C+  0   0   0   0
   G+  0   0   0   0
   T+  0   0   0   0
   N+  0   0   0   0
   *+  0   0   0   0
   A-  0   0   0   0
   C-  0   0   0   0
   G-  0   0   0   0
   T-  0   0   0   0
   N-  0   0   0   0
   *-  0   0   0   0

"""
###################################################################################################

import sys
import re
import pysam
from argparse import *
import pandas
import numpy
from Rfunctions import reverse_complement
from expand_cigar import *


parser = ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('region', default='', type=str, action='store', help='genomic region 1-based inclusive, e.g. chr4:48045153-48086447')
parser.add_argument('inbam', default='', type=str, action='store', help='input BAM file')
parser.add_argument('--l', default=1, type=int, action='store', help='controls the minimum nucleotide overlap SOFT-CLIPPED NUCLEOTIDES INCLUDED that an alignment should have with the genomic region in order to be counted')
parser.add_argument('--within', default=False, action='store_true', help='demands that the entire alignment SOFT-CLIPPED NUCLEOTIDES INCLUDED is contained within the genomic region')
parser.add_argument('--no-deletions', default=False, action='store_true', help='skips alignments with deletions in the reference')
parser.add_argument('--with-insertions', default=False, action='store_true', help='allows alignments with insertions to the reference to count, this can be problematic because the insertion introduces a frame shift to the remaining nucleotides, however neither the insertion nor the frame shift are taken into account, instead the remaining nucleotides are counted as if there was no insertion')
parser.add_argument('--debug', default=False, action='store_true', help='show on stderr the accepted read ids')


options = parser.parse_args()
#options = parser.parse_args(['chr13:58392841-58392863', '/home/fklirono/Downloads/Aligned.sortedByCoord.out.bam'])
#options = parser.parse_args(['chr13:58392841-58392863', '/home/fklirono/bio/lib/test/reads.bam'])


#import ipdb;  ipdb.set_trace()


#  deparse genomic position using <chromosome>:<start>-<end> assumption 1-based inclusive
chr, start, end = re.match(r'(^[^:]*):([0-9]*)-([0-9]*)$', options.region).group(1,2,3)
start = int(start) - 1  #  0-based left inclusive
end = int(end)  #  0-based right exclusive but last nucleotide of region should be included


#  basic checks
assert options.l > 0  #  overlap needs to be at least 1 nucleotide
assert options.l <= end - start  #  we cannot ask for more overlap than available nucleotides


#  define the array that will keep the tabs as a nice DataFrame
counts = pandas.DataFrame( numpy.zeros(shape=(12, end - start), dtype=numpy.int64), 
         index=['A+', 'C+', 'G+', 'T+', 'N+', '*+', 'A-', 'C-', 'G-', 'T-', 'N-', '*-'])


#  open BAM file for reading
fin = pysam.Samfile(options.inbam, 'rb')


#  go over the alignments one by one
for read in fin.fetch(reference=chr, start=start, end=end, until_eof=True):


    #  skip alignments with insertions unless instructed not to
    if not options.with_insertions and 'I' in read.cigarstring: 
        continue


    #  skip alignments with deletions if instructed to
    if options.no_deletions and 'D' in read.cigarstring:  
        continue


    #  expand the CIGAR string
    cigar = expand_cigar(read.cigartuples)


    #  genomic position of first nucleotide of read with soft-clipped nucleotides included
    r_start = read.reference_start - ( read.cigartuples[0][1] if read.cigartuples[0][0]==4 else 0 ) 


    #  genomic position of last nucleotide of read with soft-clipped nucleotides included
    r_end = r_start + len(read.query_sequence) - 1


    #  if instructed to do so, skip alignment not entirely within region with soft-clipped nucleotides included in the counting
    if options.within and ( r_start < start or r_end >= end ):  #  end is 0-based exclusive
        continue


    #  skip read if the region of overlap is smaller than desired  
    if r_start < start and r_end - start + 1 < options.l:  #  +1 because we are counting nucleotides not distances
        continue  
    elif r_start >= start and end - r_start < options.l:  #  end is 0-based exclusive but there is a +1 because we are counting nucleotides
        continue 


    if options.debug:
        sys.stderr.write('[DEBUG] accepted read = %s\n' % (read.query_name))


    #  iterate over the nucleotides of the read in absolute genomic coordinates
    for pos in xrange(r_start, r_end + 1):

        #  keep moving if we are to the left of the region
        if pos < start:
            continue

        #  stop moving if we are already out of the region
        if pos >= end:  
            break

        #  ignore insertions to the reference
        if cigar[pos - r_start] == 'I':
            continue

        if cigar[pos - r_start] == 'D':  #  deletion to the reference
            if read.flag & 0x16:
                row = '*-'
            else:
                row = '*+'
        else:  #  match or mismatch
            if read.flag & 0x16:
                row = ''.join( [ reverse_complement(read.query_sequence[pos - r_start]), '-' ] )
            else:
                row = ''.join( [ read.query_sequence[pos - r_start], '+' ] )
         
        #  identify the column
        col = pos - start

        #  increment the corresponding count
        counts.loc[ row, col ] += 1


fin.close()


#  spit out the counts to stdout
counts.to_csv(sys.stdout, sep='\t')



