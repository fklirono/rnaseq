#!/usr/bin/env python3
# -*- coding: utf-8 -*-
###################################################################################################
#  Filippos Klironomos, Berlin Institute for Medical Systems Biology, Max-Delbrück Center, Berlin
###################################################################################################
usage="""
  This script filters out alignments from a BAM file according to demand.

  N.B. Certain filters are aligner-specific! For example, spliced alignments are detected 
       by the presence of the XS:A tag which STAR uses to report the strand for spliced reads.

  Typical use: %(prog)s in.bam out.sam --hard-clip --sam --report 2> discarded.bam
"""
###################################################################################################

import sys
import pysam
from argparse import *
from collections import defaultdict

parser=ArgumentParser(description=usage, formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('inbam', default='', type=str, action='store', help='input BAM file')
parser.add_argument('outbam', default='', type=str, action='store', help='output BAM file with filtered reads')
parser.add_argument('--report', default=False, action='store_true', help='should discarded reads be reported in STDERR?')
parser.add_argument('--min-mismatch', default=0, type=int, action='store', help='allow this minimum number of mismatches')
parser.add_argument('--max-mismatch', default=-1, type=int, action='store', help='allow this maximum number of mismatches')
parser.add_argument('--min-length', default=0, type=int, action='store', help='allow this minimum read alignment length, i.e. the length of the part of the read that aligned')
parser.add_argument('--max-length', default=0, type=int, action='store', help='allow this maximum read alignment length, i.e. the length of the part of the read that aligned. In order for this parameter to have an effect it needs to be non-zero and larger or equal than --min-length')
parser.add_argument('--with-deletions', default=False, action='store_true', help='allow for reference deletions')
parser.add_argument('--with-insertions', default=False, action='store_true', help='allow for reference insertions')
parser.add_argument('--no-chimeras', default=False, action='store_true', help='(only for paired-end) filter out chimeras where mates map to different chromosomes')
parser.add_argument('--no-splicing', default=False, action='store_true', help='filter out spliced alignments that have the XS:A strand tag defined (STAR aligner-specific)')
parser.add_argument('--hard-clip', default=False, action='store_true', help='should we hard-clip any soft-clipped reads?')
parser.add_argument('--sam', default=False, action='store_true', help='shall we output in SAM format?')
parser.add_argument('--flags', default=[], nargs='+', type=int, action='store', help='space-separated list of SAM FLAGs (followed by -- if it is that last flag passed) to accept, useful for paired-end alignments to pick (99, 147) or (83, 163) proper pairs and filter out everything else ')
options = parser.parse_args()


#  convert to list of integers (no error if empty)
options.flags=list(map(int, options.flags))  


inbam=pysam.Samfile(options.inbam, 'rb')
if options.sam:
    outbam=pysam.Samfile(options.outbam, 'w', template=inbam, add_sam_header=False)
else:
    outbam=pysam.Samfile(options.outbam, 'wb', template=inbam)


if options.report:
  removedbam=pysam.Samfile('/dev/stderr', 'wb', template=inbam)


#import ipdb;  ipdb.set_trace()


for read in inbam.fetch(until_eof=True):

    #  convert list of tuples to dictionary for easy lookup with no idiotic KeyError exceptions raised if tags are not present
    tags = dict(read.tags)
    nm = tags.get('NM', 0)  #  read.get_tag('NM') would raise KeyError if 'NM' tag is not present...

    #  either no upper bound (max_mismatch<0) or there is an upper bound and we are below or up to it
    mismatches_ok = ( nm>=options.min_mismatch and ( options.max_mismatch<0 or (options.max_mismatch>=0 and nm<=options.max_mismatch) ) )

    #  either no flags were provided or we have the right flags
    flags_ok = len(options.flags)==0 or (read.flag in options.flags)

    chimeras_ok = True
    if options.no_chimeras:
        chimeras_ok = (read.rnext == -1) or (read.rnext == read.rname)  #  no need to convert to reference names: inbam.getrname(read.rnext)

    #  STAR-aligner specific
    splicing_ok = True
    if options.no_splicing:
        splicing_ok = not ( 'XS' in tags )
    
    #  read.qlen = aligned nts excluding soft-clipped (or full read length for unmapped reads)
    #  read.rlen = actual read length including soft-clipped nts
    #  read.alen = reference nts aligned with this read = (read.qlen + deletions, read.qlen - insertions)
    width_ok=( read.qlen>=options.min_length )  #  the alignment part needs to fulfill the length requirement
    if options.max_length!=0: 
        width_ok=( width_ok and options.max_length>=options.min_length and read.qlen<=options.max_length )

    if read.qlen != read.rlen and options.hard_clip:
        op,n = read.cigartuples[0]  
        if op == 4:  #  is there any left soft-clipping?
            read.seq = read.seq[n:read.rlen]
            read.cigarstring = read.cigarstring.replace('S', 'H', 1)  #  replace only first occurrence 
        op,n = read.cigartuples[len(read.cigartuples)-1]  
        if op == 4:  #  is there any right soft-clipping?
            read.seq = read.seq[0:(read.rlen-n)]
            read.cigarstring = read.cigarstring.replace('S', 'H', 1)  #  there should be only one occurrence anyway...
    
    #sys.stderr.write(read.qname + '\n')  #  DEBUG

    #  make sure you check for indels only if the read is aligned
    deletion = False
    insertion = False
    if read.flag != 4:
        deletion=('D' in read.cigarstring)
        insertion=('I' in read.cigarstring)
    
    if flags_ok:
        if splicing_ok:
            if chimeras_ok:
                if mismatches_ok:
                    if options.with_deletions or not deletion:
                        if options.with_insertions or not insertion:
                            if width_ok:
                                outbam.write(read)
                            elif options.report:
                                removedbam.write(read)  #  outside width range
                        elif options.report:
                          removedbam.write(read)  #  insertions not welcome and this is an insertion
                    elif options.report:
                      removedbam.write(read)  #  deletions not welcome and this is a deletion
                elif options.report:
                  removedbam.write(read)   #  outside mismatch range
            elif options.report:
             removedbam.write(read)   #  chimeras not welcome and this is a chimera
        elif options.report:
          removedbam.write(read)   #  spliced read because the alignments skips reference bases
    elif options.report:
        removedbam.write(read)   #  FLAG is not part of our desired FLAG list


inbam.close()
outbam.close()
if options.report:
    removedbam.close()



