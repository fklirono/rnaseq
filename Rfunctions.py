# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
#
#  set of R functions to have in python
#
###################################################################################################

def head(x, n=10):
  return x[0:n+1]
    

def tail(x, n=10):
  return x[len(x)-n:]
    

def unlist(l):
  tmp=list()
  _=[tmp.extend(s) for s in l]
  return tmp

def which(lst, what):
  """look into list lst for elements in the list what and return indices"""
  #  'lst' has to be iterable (list)
  #  'what' has to be iterable (list)
  return [ m for m,x in enumerate(lst) if x in what]

