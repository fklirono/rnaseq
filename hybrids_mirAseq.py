#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  script isolates the alignment pattern of miRNAs in the hybrids by removing the gaps ('-') 
  and saves the new alignments so that they can be loaded in R  
"""

import os
import re 
from optparse import *

parser = OptionParser(usage=usage)
parser.add_option("-o","--output", default="mirAseq", action="store", dest="output", help="output file to save MIRZA miRNA alignments")
options,args = parser.parse_args()


print 'reading the hybrids...'
#  read the hybrids
f=open(os.path.dirname(options.output)+'/hybrids.processed', "r")
hybrids=f.readlines()
f.close()
print 'done! Number of hybrids found = %d' % (len(hybrids))

print 'saving mirAseq conversions concurrently at:\n %s' % (options.output)
f=open(options.output, "w")
#  get the corresponding mirAseq and tA columns
seqs=[ [re.split(r' ', w)[3], re.split(r' ', w)[5]] for w in hybrids]

#  save memory
del hybrids

#  iterate over all pairs and isolate the sequence chunks in the 
#  mirAseq (anything not '-') and their coordinates

for s in seqs:
  #  find indices in s[0] of all non '-' characters and return the corresponding s[1] character 
  f.write(''.join([ s[1][n] for n in [m for m,x in enumerate(s[0]) if x!='-'] ])+'\n')

f.close()
print 'done!\n\n'
