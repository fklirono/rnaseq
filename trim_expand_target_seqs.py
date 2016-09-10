#!/usr/bin/env python
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage = """
  %prog [options]

  script that trims or expands sequences mapped to a reference transcriptome 
  with a choice of which end to add to and remove from

  the default is to add at the 3' end and remove at the 5' end
"""
###################################################################################################

import re
from Rfunctions import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from argparse import *


parser=ArgumentParser(usage=usage)
parser.add_argument('--input', default='targets.bed', action='store', dest='input', help='file with miR and target pairs')
parser.add_argument('--cols', default='4,9', action='store', dest='cols', help='target name, target sequence columns in --input file in that order [first column = column number 1]')
parser.add_argument('--reference', default='reference.fa', action='store', dest='reference', help='reference transcriptome in FASTA format')
parser.add_argument('--output', default='sequences.procesed.', action='store', dest='output', help='root name for output file to concatenate --readsize --add5p --remove3p values as strings to it')
parser.add_argument("--add5p", default=True, action="store_false", dest="add3p", help="if a sequence needs to be augmented, always augment the 3' end of it if possible")  #  --add5p : add3p == False
parser.add_argument("--remove3p", default=True, action="store_false", dest="remove5p", help="if a sequence needs to be shrunk, always shrink the 5' end of it if possible")  #  --remove3p : remove5p == False
parser.add_argument('--readsize', type=int, default=51, required=True, action='store', dest='readsize', help='root name for output file to concatenate options.readsiez _add3p _remove5p etc...')
options = parser.parse_args()


#  target names, target sequences columns, convert to first column = column number 0
c_target, c_targetseq = tuple([ int(x)-1 for x in re.split(',', options.cols)])


#  output FASTA file
OUTFILE=options.output+str(options.readsize)+'nts'
if options.add3p:
  OUTFILE+='_add3p'
else:
  OUTFILE+='_add5p'
if options.remove5p:
  OUTFILE+='_remove5p'
else:
  OUTFILE+='_remove3p'

OUTFILE+='.fa'


#  read chimeras
print 'loading chimeras...'
f=open(options.input, 'r')
chimeras=f.readlines()
f.close()
print 'done! Found = %d\n' % (len(chimeras))

#  get target id and sequence
targets=unlist([ [m[c_target], m[c_targetseq]] for m in [re.split('\t', c) for c in chimeras] ])

#  strip possible spaces, \n from everywhere
targets=[r.strip(' \n') for r in targets]


#  read the transcriptome
print 'loading reference...'
f=open(options.reference, 'r')
txs=[]
_=[txs.extend([s.id, str(s.seq)]) for s in SeqIO.parse(f, 'fasta')]
f.close()
print 'done! Number of transcripts in reference = %d\n' % (len(txs))


if options.add3p:
  print "sequences smaller than %d will be augmented at the 3' end if possible" % (options.readsize)
else:
  print "sequences smaller than %d will be augmented by the 5' end if possible" % (options.readsize)

if options.remove5p:
  print "sequences larger than %d will be trimmed from the 5' end if possible" % (options.readsize)
else:
  print "sequences larger than %d will be trimmed from the 3' end if possible" % (options.readsize)

print '\noutput file will be = %s\n' % (OUTFILE)


print 'Starting to map the target sequences...'
#  iterate over the ligated targets
for k,r in enumerate(targets[1::2]):
  r_len=len(r)
  if r_len < options.readsize:
    matches=[ r in s for s in txs[1::2]]  #  find all transcripts that match
    if matches.count(True)==0:  #  ligated read does not map to hg19 transcriptome!
      print 'target = %s (%s) is not found in the reference' % (targets[2*k], targets[2*k+1])
      continue
    else:
      m=matches.index(True)*2+1  #  pick the first transcript match
      coords=re.search(r, txs[m])  #  find the coordinates 
      width=len(coords.string)  #  width of the transcript 
      start,end=coords.span()  #  in python: end-start = width of ligated read
      possible_to_add3p=(width-end > options.readsize - r_len)
      possible_to_add5p= (start - options.readsize + r_len > 0)
      if possible_to_add3p and options.add3p:  #  plenty of space to extend in the 3' direction provided we wish for it
        targets[2*k+1]=txs[m][start:(end+options.readsize-r_len)]
      elif possible_to_add5p and not options.add3p :  #  plenty of space to extend in the 5' direction provided we wish for it
        targets[2*k+1]=txs[m][(start-options.readsize+r_len):end]
      elif possible_to_add3p:  #  forced to extend in the direction that is available
        targets[2*k+1]=txs[m][start:(end+options.readsize-r_len)]
      elif possible_to_add5p:  #  forced to extend in the direction that is available
        targets[2*k+1]=txs[m][(start-options.readsize+r_len):end]
      else:  #  not enough exons to extend in either direction for the first picked transcript
        print "cannot extend ligated read = %s (%s)" % (targets[2*k], targets[2*k+1])
        continue
   
  elif r_len > options.readsize:  
    if options.remove5p:  #  remove from 5' end
      targets[2*k+1]=r[(r_len-options.readsize):]
    else:  #  remove from 3' end
      targets[2*k+1]=r[0:options.readsize]

print 'finished!\n'

#  write trimmed ligated targets to FASTA file
print 'saving sequences...'
reads_fa=[]
for (r,i) in (zip(targets[1::2], targets[0::2])):
  reads_fa.append(SeqRecord(Seq(r), id=i, description=''))

SeqIO.write(reads_fa,OUTFILE, 'fasta')
print 'done! Make sure to remove targets that did not map to the reference (if any)!\n'

