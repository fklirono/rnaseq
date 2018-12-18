#!/usr/bin/env python
#
#  version modified by FDK of the original /home/fklirono/pcp/sequence_data/bin/bam2fastq.py
#
#  how to run:
#
#    cat aligned_reads.bam | BAM2fastq.py > reads.fastq
#
#  N.B. samtools bam2fq ought to be faster...
import pysam,re,sys,os
from Rfunctions import reverse_complement

bam = pysam.Samfile('-',"rb")
for a in bam:
  seq = a.seq
  qual = a.qual
  if a.is_reverse:
    seq = reverse_complement(seq)
    qual = qual[::-1]
  
  sys.stdout.write('@'+a.qname+'\n'+seq+'\n+\n'+qual+'\n')

