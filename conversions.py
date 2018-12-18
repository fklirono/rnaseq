#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  version 2

  this script reports the positions on the genome (0-based) in STDOUT of the conversions encountered in single-end sequenced RNA alignments
  and their summaries at the end in STDERR

  the script can optionally filter to [--tc output.bam] the reads with at least one RNA T:C conversion, for multi-mappers it adds the 
  corresponding SEQ and QUAL SAM fields as well

  reads with insertions are ignored 

  secondary alignments (e.g. produced by 'bwa mem -a' with FLAG=={256,272}) are included in the analysis ONLY IF THEY FOLLOW 
  DIRECTLY UNDER THE PRIMARY ALIGNMENT, that means the alignments need to be sorted by read-ID and FLAG (0,16,256,272) before
  fed to this script with a command like:

    samtools view INPUT.bam | sort -S 1G -k1,1 -k2,2n | samtools view -S -bt REFERENCE.fa.fai -

  ALL CONVERSION COUNTS REFER TO THE CHANGE OF THE RNA MOLECULE, for example:

    an alignment on the reverse-strand with a G:A conversion means that the actual RNA molecule had a C:T conversion, consequently
    the C:T conversion count is incremented in this case, not the G:A count

  Example how to run:

  conversions.py <(samtools view input.bam | sort -S 1G -k1,1 -k2,2n | samtools view -S -bt hg19.fa.fai -) > positions.tsv 2> stats.tsv

  results:

    stats.tsv:

      conversion  number_from_primary_alignments  number_from_secondary_alignments 

    positions.tsv:

      chromosome  position  conversion(A:G,A:C,...)  supported_by_primary_alignment(T(rue)/F(alse))
  
  KNOWN BUGS:
    
    NONE
"""
###################################################################################################

import sys
import re
import pysam
from argparse import *
from Rfunctions import reverse_complement

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('inbam', default='', type=str, action='store', help='input BAM file')
parser.add_argument('--tc', default='', type=str, action='store', help='save reads with RNA T:C conversions separately under the given BAM filename')


options = parser.parse_args()


#  define the primary and secondary conversions counts
conv_pri={ 'A:G': 0, 'A:C': 0, 'A:T': 0, 'A:*':0, 
           'G:A': 0, 'G:C': 0, 'G:T': 0, 'G:*':0,
           'T:A': 0, 'T:C': 0, 'T:G': 0, 'T:*':0, 
           'C:A': 0, 'C:T': 0, 'C:G': 0, 'C:*':0}
conv_sec={ 'A:G': 0, 'A:C': 0, 'A:T': 0, 'A:*':0, 
           'G:A': 0, 'G:C': 0, 'G:T': 0, 'G:*':0,
           'T:A': 0, 'T:C': 0, 'T:G': 0, 'T:*':0, 
           'C:A': 0, 'C:T': 0, 'C:G': 0, 'C:*':0}


fin=pysam.Samfile(options.inbam, 'rb')
if options.tc!='':
  ftc=pysam.Samfile(options.tc, 'wb', header=fin.header)



#import ipdb;  ipdb.set_trace()


#  CIGAR (all operations are on the reference):
#    0 -> M
#    1 -> I
#    2 -> D
#    3 -> N
#    4 -> S
#    5 -> H
#    6 -> P
#    7 -> =
#    8 -> X
#
#  read.rlen = actual read length including soft-clipped nts
#  read.qlen = aligned nts excluding soft-clipped
#  read.alen = reference nts aligned with this read = (read.qlen + deletions - insertions)
#    e.g. 
#      101M -> [(0, 101)] -> read.rlen = read.qlen = read.alen = 101
#      28S43M1D30M -> [(4, 28), (0, 43), (2, 1), (0, 30)] -> read.rlen = 101 (D is on the reference), read.qlen = 73, read.alen = 74
#
#  read.query = soft/hard clipped read sequence 
#  read.seq = unclipped read sequence 
#  read.qstart = number of nts clipped from 5' 
rid=''
for read in fin.fetch(until_eof=True):
  #  check if the read-id has changed and update the primary read sequence and strand 
  if read.qname!=rid:
    rid=read.qname
    assert read.flag in (0, 16)  #  primary aligment should always be first!
    strand=(read.flag & 0x16)  #  0 for plus, 16 for minus
    seq=read.seq  #  use full read sequence (no clippling) because supplementary aligments might not clip
    qual=read.qual  #  nucleotide quality scores to report in filtered RNA T:C reads
    seq_rc=reverse_complement(seq)
   
  #  after having updated the new rid, check if the primary or secondary alignment has insertion or is of a weird FLAG and ignore
  if (read.flag not in (0,16,256,272)) or ('I' in read.cigarstring):
    continue
  
  #  walk through the MD string to identify mismatches/deletions
  pos=read.qstart
  gpos=0  #  position on the read is different than position on the genome when soft/hard-clipping is present
  one_tc=False  #  activate flag if at least one T:C RNA conversion is encountered for this read in order to print this read out if required
  for md in re.finditer(r"(\d+)|([ACGT]+)|\^([ACGT]+)",read.get_tag('MD')):
    m,mm,d = md.groups()  #  match, mismatch, deletion
     
    if m:
      pos+=int(m)  #  move on the aligned read sequence by the number of matches found
      gpos+=int(m)  #  move on the genome by the number of matches found
    
    elif mm:
      mm=list(mm)
      for i in mm:  #  process consecutive nucleotide mismatches one at a time
        if read.flag in (16, 272):  #  current alignment on - strand
          i=reverse_complement(i)  #  RNA edit is the reverse-complement of the reported
          if strand==0:  #  primary alignment on + strand
            s=reverse_complement(seq_rc[pos])  #  RNA nucleotide is the reverse-comlement of the primary reverse-complement SEQ nucleotide
          else:  #  primary alignment on - strand
            s=reverse_complement(seq[pos])  #  RNA nucleotide is the reverse-comlement of the primary SEQ nucleotide
        else:  #  current alignment on + strand
          if strand==0:  #  primary alignment on + strand
            s=seq[pos]  #  RNA nucleotide is the primary SEQ nucleotide
          else:  #  primary aligment on - strand
            s=seq_rc[pos]  #  RNA nucleotide is the reverse-comlement of the primary SEQ nucleotide
        
        k=':'.join([i,s])  
        if options.tc!='' and k=='T:C':
          one_tc=True
        
        if read.flag in (0,16):
          conv_pri[k]+=1
          sys.stdout.write('\t'.join(map(str, [fin.getrname(read.tid), read.pos+gpos, k, 'T']))+'\n')
        else:
          conv_sec[k]+=1
          sys.stdout.write('\t'.join(map(str, [fin.getrname(read.tid), read.pos+gpos, k, 'F']))+'\n')
        
        pos+=1  #  move past the mismatch on the read
        gpos+=1  #  move past the mismatch on the genome as well
     
    elif d:
      d=list(d)
      for i in d:  #  process consecutive nucleotide deletions one at a time
        gpos+=1  #  advance one the genome position over the deleted nucleotide
        if read.flag in (16,272):
          i=reverse_complement(i)  #  RNA deletion is the reverse-complement of the reported deletion
        k=':'.join([i,'*'])
        
        if read.flag in (0,16):
          conv_pri[k]+=1
          sys.stdout.write('\t'.join(map(str, [fin.getrname(read.tid), read.pos+gpos, k, 'T']))+'\n')
        else:
          conv_sec[k]+=1
          sys.stdout.write('\t'.join(map(str, [fin.getrname(read.tid), read.pos+gpos, k, 'F']))+'\n')
        
        #  for a deletion there is no need to move to the next nucleotide 
   
  if one_tc and options.tc!='':  #  we need to define SEQ and QUAL in case this is not a primary alignment
    #  current alignment on opposite strand from primary alignment, we need to reverse-complement SEQ and QUAL
    if (read.flag in (16, 272) and strand==0) or (read.flag in (0, 256) and strand==16):  
      read.seq=seq_rc
      read.qual=qual[::-1]
    #  current alignment on identical strand as primary alignment
    else:
      read.seq=seq
      read.qual=qual
    ftc.write(read)  #  write read that contains at least one T:C conversion


fin.close()
if options.tc!='':
  ftc.close()


sys.stderr.write('\t'.join(['conversion', 'primary', 'secondary'])+'\n')
for c in conv_pri:  #  iterate over keys and print out the counts
  sys.stderr.write('\t'.join(map(str, [c, conv_pri[c], conv_sec[c]]))+'\n')


