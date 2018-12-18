# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
#
#  set of functions used in sequence data analysis
#
###################################################################################################



def occurrences(pattern, string, no_overlaps=False):
  '''
  1-step sliding-window count of occurrences of pattern in string 
  if no overlaps are allowed then the sliding-window becomes as large as the length of the string

  original code stolen from: http://stackoverflow.com/a/2970542
  '''

  count=start=0
  step=1
  if no_overlaps:
      step=len(pattern)


  while True:
    start=string.find(pattern, start) + step
    if start>0:
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

  returns an iterator of all possible sequences with exactly mms mismatches to the original 'seq'.

  N.B. N's are replaced with ACGT but only when their turn comes, i.e. for mms=1 the sequence TNC will generate all possible
       mutations for T=[ACG] and pass them down as ANC, CNC, GNC and only when N is considered will it be passed down as N=[ACGT]
       however most of the times this is not what we want...

  '''
  from itertools import combinations, product

  nucleotides = set('ACGT')

  #import ipdb;  ipdb.set_trace()

  for positions in combinations(range(len(seq)), mms):  #  generate all possible mismatch positions
    variants = []

    for p in positions:
      nt = set(seq[p])
      variants.append(nucleotides-nt)  #  for each nucleotide generate list of mismatches, if nucleotide is N it will generate ACGT mismatches

    mms_nts = product(*variants)  #  for each set of mismatch positions generate all possible nucleotide mismatches
    for nts in mms_nts:
      l = list(seq)
      for i, m in enumerate(nts):
         l[positions[i]] = m  #  replace original sequence at mismatch position with given mismatch nucleotide
      yield "".join(l)



def seqs_with_1mm_and_N(seq):
    '''
    [Andrei's code, modified]
    
    wrapper function of seqs_with_mms(mms=1) which treats N as the single mismatch and replaces it with ACGT. 
    If more that one N's are present it returns the original sequence.
    If no N's are present it calls seqs_with_mms(mms=1).
    '''
    from itertools import combinations, product
    

    #  if we are looking for one mismatch but we have already more than one N's return empty string instead of error
    #  if we have single N then return the 4 different sequences
    if seq.count('N')>1:
        yield seq
    elif seq.count('N')==1:
        i = seq.find('N')
        variants = []
        for s in 'ACGT':
            x = list(seq)
            x[i] = s
            yield ''.join(x)

    
    #  if we made it to here there are no N's present so just call seqs_with_mms()
    for s in seqs_with_mms(seq, mms=1):
        yield s



def find_all_1mm(bc, bclist):
    '''
    return all matches in 'bclist' of all sequences generated from 'bc' with a single nucleotide mismatch
    '''

    bc_1mm = []
    for mm in seqs_with_1mm_and_N(bc):
        if mm in bclist:
            bc_1mm.append(mm)

    return bc_1mm



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



def reduce_ranges(tuples):
    """
    ranges provided as a list of tuples of (start, end) are checked for overlaps and are reduced similarly to the reduce() R function of IRanges
    """
    from Rfunctions import unlist
    from copy import deepcopy


    #  sanity check
    if (len(tuples)<2):
        return tuples
    
    
    #  sort the ranges by start position and for ties by end position
    tuples.sort()


    #import ipdb;  ipdb.set_trace()
    #from IPython.core.debugger import Pdb; ipdb=Pdb(); ipdb.set_trace()

    
    i = 0
    j = i + 1
    while True:
        #  convert tuples to inclusive ranges that are in turn converted to sets for easy overlap comparisons
        ranges = [ set(range(start, end+1)) for start, end in tuples ]


        #  go over all downstream ranges and identify the overlaps with the current range
        if len( ranges[i] & ranges[j] )!=0:
            u = unlist([ ranges[i], ranges[j] ])
            start = min(u)
            end = max(u)


            #  replace the upstream tuple with the new (start,end) values 
            tuples[i] = (start, end)


            #  delete the merged downstream tuple
            del tuples[j]


            #  update the ranges to reflect the merging
            ranges = [ set(range(start, end+1)) for start, end in tuples ]  
        else:
            j+=1


        #  last tuple to be compared with current tuple has been reached, time to move one step down
        if j>=len(ranges):
            i+=1
            j = i + 1

        #  all comparisons have been made, we are sitting at the last tuple with no other downstream of it
        if i>=(len(ranges)-1):
            break


    return tuples

    


