# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################



def expand_cigar(cigartuples):
    """
    expand the CIGAR string by converting tuples to per-nucleotide operation string
    """
    
    cigar = ''
    for t,n in cigartuples:
        if t == 0:  #  M
            cigar = cigar + 'M' * n
        elif t == 1:  #  I
            cigar = cigar + 'I' * n
        elif t == 2:  #  D
            cigar = cigar + 'D' * n
        elif t == 3:  #  N
            cigar = cigar + 'N' * n
        elif t == 4:  #  S
            cigar = cigar + 'S' *n
        elif t == 5:  #  H
            cigar = cigar + 'H' *n
        elif t == 6:  #  P
            cigar = cigar + 'P' *n
        elif t == 7:  #  =
            cigar = cigar + '=' *n
        elif t == 8:  #  X
            cigar = cigar + 'X' *n
        else:
            raise ValueError('unknown CIGAR operation')

    return cigar



