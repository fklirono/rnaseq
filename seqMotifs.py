# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
#
#  set of functions used in sequence data analysis
#
###################################################################################################


def occurrences(string, pattern):
  '''
  overlapping count of occurrences of pattern in string

  stolen from: http://stackoverflow.com/a/2970542
  '''
  count = start = 0
  while True:
    start = string.find(pattern, start) + 1
    if start > 0:
      count+=1
    else:
      return count


def dinucleotideFrequency(seq, add_one=False, as_prob=False):
  '''
  compute the single-step sliding window dinucleotide frequency and return a dictionary of them

  as_prob : shall the frequencies be returned as probabilities?
  add_one : add one to the counts?

  THERE IS NO CHECK IF SEQUENCE CONSISTS OF A,C,G,T NUCLEOTIDES ONLY
  '''
  seq=seq.upper()

  dinucleotides = set(['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'CT', 'GA', 'GC', 'GG', 'GT', 'TA', 'TC', 'TG', 'TT'])

  frequencies = dict()
  for d in dinucleotides:
    frequencies[d]=int(add_one)
  
  total=0
  for d in dinucleotides:
    frequencies[d]+=occurrences(seq,d)
    total+=frequencies[d]

  if as_prob:
    for k in frequencies.keys():
      frequencies[k]=float(frequencies[k])/float(total)
  
  return(frequencies)


def nucleotideFrequency(seq, add_one=False, as_prob=False):
  '''
  compute nucleotide frequencies and return a dictionary of them

  as_prob : shall the frequencies be returned as probabilities?
  add_one : add one to the counts?

  THERE IS NO CHECK IF SEQUENCE CONSISTS OF A,C,G,T NUCLEOTIDES ONLY
  '''
  seq=seq.upper()

  nucleotides = set(['A','C','G','T'])

  frequencies = dict()
  for n in nucleotides:
    frequencies[n]=int(add_one)
  
  total=0
  for n in nucleotides:
    frequencies[n]+=occurrences(seq,n)
    total+=frequencies[n]

  if as_prob:
    for k in frequencies.keys():
      frequencies[k]=float(frequencies[k])/float(total)
  
  return(frequencies)


def diShuffle(seq):
  '''
  shuffle sequence "seq" at the dinucleotide level WITHOUT PRESERVING THE DINUCLEOTIDE FREQUENCY!
  '''
  from random import shuffle
  
  y=[seq[i:i+2] for i in range(0, len(seq), 2)]

  shuffle(y)

  return ''.join(y)

  
def seqs_with_mms(seq, mms=1):
  '''
  [Andrei's code]

  returns an iterator of all possible sequences with exactly "mms" mismatches to the original "seq"
  '''
  from itertools import combinations, product

  nucleotides = set("ACTG")

  #import ipdb;  ipdb.set_trace()

  for positions in combinations(range(len(seq)), mms):  #  generate all possible mismatch positions
    variants = []
    for p in positions:
      nt = set(seq[p])
      variants.append(nucleotides-nt)  #  for each nucleotide generate list of mismatches
    mms_nts = product(*variants)  #  for each set of mismatch positions generate all possible nucleotide mismatches
    for nts in mms_nts:
      l = list(seq)
      for i, m in enumerate(nts):
         l[positions[i]] = m  #  replace original sequence at mismatch position with given mismatch nucleotide
      yield "".join(l)


def insertions(sequence, insert):
  """  
  returns an iterator of all possible sequences generated by inserting "insert" at all possible places in "sequence"

  Multiplicities might arise, e.g all possible insertions for 

    sequence=ACG, insert=AA 
    
  are 
  
    AAACG , AAACG , ACAAG , ACGAA 

  """
  for s in range(len(sequence)):
    yield sequence[:s] + insert + sequence[s:] 


def deletions(sequence, n=1):
  """  
  returns an iterator of all possible sequences generated by deleting "n" characters from "sequence"
  """
  from itertools import combinations

  if (len(sequence)<=n):
    yield ''
  else:
    for r in combinations(sequence, len(sequence)-n):
      yield ''.join(r) 


def s_plus_1d(s):  
  """  
  find motif of sequence + 1 deletion by sequentially inserting ? next to each character of the sequence 
  and finally returning a list of the regex compiles 
  """

  import re
  return([ re.compile(s[0:(n+1)]+'?'+s[(n+1):], re.IGNORECASE) for n in range(len(s)) ])


def s_plus_1mm(s):  
  """  
  generate regular expressions for +1 mismatch in sequence 

  each sequence character s[i] is sequentially replaced by [^s[i]] regular expression

  attention needs to be paid at beginning and end of reference sequence where removing 
  first s[0] or last s[len(s)-1] character of pattern sequence also counts as +1 mismatch
  """

  import re
  regs=[ re.compile(s[:n]+'[^'+s[n]+']'+s[n+1:], re.IGNORECASE) for n in range(len(s)) ]

  #  include mismatches at beginning or end in reference sequence
  regs.extend([re.compile('^'+s[1:], re.IGNORECASE), re.compile(s[0:len(s)-1]+'$', re.IGNORECASE)])  
  return(regs)


def myfindall(rgx, seq):
  """
  re.findall() replacement reporting overlapping matches as well

  template taken from https://mail.python.org/pipermail/tutor/2005-September/041126.html
  """

  import re

  results=[]
  pos=0
  while True:
    r = rgx.search(seq, pos)
    if r is None:
      break
    results.append(seq[r.start():r.end()])  #  .extend() will break the sequence to letters
    pos = r.start()+1
  return results


def myfinditer(rgx, seq):
  """
  This function is not a direct replacement of re.finditer() since it does not return an iterator.
  Instead, it returns as a list all (start, end) tuples of matches.
  """

  import re

  results=[]
  pos=0
  while True:
    r = rgx.search(seq, pos)
    if r is None:
      break
    results.append(r.span())  #  .extend() will break the tuples up
    pos = r.start()+1
  return results

