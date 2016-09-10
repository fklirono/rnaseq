#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  %prog [options]

  this script finds the index in the HYBRIDS of all pairs found in the CHIMERAS 

"""
###################################################################################################

import re
from Rfunctions import *
from argparse import *


parser=ArgumentParser(usage=usage)
parser.add_argument("--hybrids", default="hybrids.txt", action="store", dest="hybrids", help="hybrids file with (target, mir) columns")
parser.add_argument("--chimeras", default="chimeras.txt", action="store", dest="chimeras", help="chimeras file with (target, mir) columns")
options=parser.parse_args()


#  read chimeras
print 'loading the hybrids and the chimeras...'
f=open(options.hybrids, 'r')
hybrids=f.readlines()
f.close()
f=open(options.chimeras, 'r')
chimeras=f.readlines()
f.close()
print 'done! hybrids found = %d, chimeras found = %d\n' % (len(hybrids), len(chimeras))

#  remove \n from everywhere
print 'processing them...'
hybrids=[p.strip('\n') for p in hybrids]
chimeras=[p.strip('\n') for p in chimeras]
print 'done!\n'

print 'starting pulling down the pairs and saving them concurrently...'
f=open(options.chimeras+'_with_indices', "w")
for pair in chimeras:
  c=hybrids.count(pair)  #  get the count
  if c==1:
    ind=hybrids.index(pair)  #  the first one returned by .index() is the last one also
    f.writelines(' '.join([pair,str(ind+1)])+'\n')  #  +1 for R objects, don't use TAB but space, since target-miR have space between them
  elif c==0:  #  it is possible that a given miR (e.g. < 21nts) has dropped out from MIRZA
    print 'interaction pair (%s) cannot be found in hybrids' % (','.join(re.split(' ', pair)))
    f.writelines(' '.join([pair,'NA'])+'\n')
  else:
    print '\n***** more than one hybrids found for pair (%s), this should not happen *****\n' % (','.join(re.split(' ', pair)))

f.close()
print 'finished!\n'


