#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
#
#  maps old miRBase annotations in the chimeras to new -5p, -3p and replaces all mapped
#  sequences of miRs with the corresponding new ones
#
#  !!!INPUT BED/TSV/whatever file needs to be tab-separated!!!
#
###################################################################################################
usage = """
  %prog --bed file.bed --cols 7,8 --mirbase miRBase.fa > output.log

  in the end a file.bed.map is generated with the mappings

  if XXXX are present in the miR sequence they are removed

  columns need to be tab-separated
"""

import re
from optparse import *
import seqMotifs as sm
from Rfunctions import *

parser = OptionParser(usage=usage)
parser.add_option('--bed', default='interactions.bed', action='store', dest='bedfile', help='chimeras in BED format')
parser.add_option('--cols', default='7,8', action='store', dest='cols', help='columns in BED file where miR name and miR sequence (in that order) are to be found [first column = column number 1]')
parser.add_option('--mirbase', default='miRBase.rel20.fa', action='store', dest='mirbase', help='reference miRBase FASTA file to map miRs to')
options,args = parser.parse_args()


#  turn string to integers for the columns where miR names and sequences are
#  convert to first column = column number 0
names_col,seq_col = tuple([ int(x)-1 for x in re.split(',', options.cols)])

#  load the chimera miRs
f=open(options.bedfile, 'r')
chimeras=f.readlines()
f.close()
print 'Number of chimeras = %d\n' % (len(chimeras))

#  get the miRs and the corresponding sequences 
mirs=unlist([re.split(r'\t', w)[names_col:seq_col+1] for w in chimeras])

#  strip possible *, \n from names
mirs[0::2]=[m.strip('*\n') for m in mirs[0::2]]

#  strip possible X,\n from sequences
mirs[1::2]=[m.strip('X\n') for m in mirs[1::2]]

#  load the reference
f=open(options.mirbase, "r")
mirbase=f.readlines()
f.close()
print "Number of reference miRs = %d\n" % (len(mirbase)/2)
# strip newlines, > and *
mirbase=[s.strip('\n>*') for s in mirbase]

#  find miRNA by name and compare sequences of -5p, -3p
for k,m in enumerate(mirs[0::2]):
    seq=mirs[2*k+1]
    if mirbase.count(m)!=0:  #  does it exist already without -5p, -3p annotation?
        i=mirbase.index(m)
        mirs[2*k+1]=mirbase[i+1]  #  replace sequence with miRBase sequence just in case
        continue
     
    if mirbase.count(m+'-3p')!=0:  #  try to look for the -3p
        i=mirbase.index(m+'-3p')
        if seq in mirbase[i+1]:  #  found a -3p, accept it if the sequence fragment matches some part of the -3p sequence
            mirs[2*k]=mirbase[i]
            mirs[2*k+1]=mirbase[i+1]  #  replace sequence with miRBase sequence
            continue
     
    if mirbase.count(m+'-5p')!=0:  #  -3p is not found, try to look for the -5p
        i=mirbase.index(m+'-5p')
        if seq in mirbase[i+1]:  #  found a -5p, accept it if the sequence fragment matches some part of the -5p sequence
            mirs[2*k]=mirbase[i]
            mirs[2*k+1]=mirbase[i+1]  #  replace sequence with miRBase sequence
            continue
     
    #  no -3p, or -5p was found, look up for first match of sequence fragment in whole miRBase (lower miRBase numbers are better)
    first=[ seq in s for s in mirbase[1::2] ]  
    if first.count(True)!=0:
      first=first.index(True)*2  #  index of mirbase[1::2] times 2 is the index of the matching miR
      print "%s (%s) cannot be mapped to -5p, -3p but to %s" % (mirs[2*k], seq, mirbase[first])
      mirs[2*k]=mirbase[first]
      mirs[2*k+1]=mirbase[first+1]
    else:  #  this sequence fragment does not match anywhere in miRBase 20
      print '********** %s (%s) cannot be mapped by sequence to anything **********' % (mirs[2*k], seq)

#  add tabs and newlines back
mirs[0::2]=[ s+'\t' for s in mirs[0::2]]
mirs[1::2]=[ s+'\n' for s in mirs[1::2]]

#  write the mappings to file
f=open(options.bedfile+'.map', 'w')
f.writelines(mirs)
f.close()
