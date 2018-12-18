#!/usr/bin/env python2
# -*- coding: utf-8 -*-
###############################################################################################################
#
#  Filippos Klironomos, Department of Pediatric Hematology, Oncology and SCT Charité University Hospital Berlin
#
###############################################################################################################
usage = """
    this script takes BAM alignments with read_ids that contain the cell barcode and the UMI according to:

        read_id_(cell barcode)_(UMI)

    and introduces the XC:Z:(cell barcode) and XM:Z:(UMI) tags instead restoring the original read_id.

    Typical usage:
    
        %(prog)s input.bam output.tagged.bam
"""
###################################################################################################

import pysam
from argparse import *

parser = ArgumentParser(usage=usage)
parser.add_argument('input', default='', type=str, action='store', help='BAM file with cell barcode and UMI embeded in the read_id')
parser.add_argument('output', default='', type=str, action='store', help='output BAM file with XC:Z:(cell barcode) and XM:Z:(UMI) tags introduced and with restored read_ids.')
options = parser.parse_args()


#import ipdb;  ipdb.set_trace()
#from IPython.core.debugger import Pdb; ipdb = Pdb(); ipdb.set_trace()


#  open BAM files for reading and writing
IN = pysam.AlignmentFile(options.input, 'rb')
OUT = pysam.AlignmentFile(options.output, 'wb', template=IN)


for read in IN:

    #  identify cell barcode and UMI from the read_id from the last two elements (in case the original read_id contains _ as well)
    #  reconstruct if necessary the original read_id
    read_id = read.query_name.split('_')
    barcode, umi = read_id[-2:]
    read_id = '_'.join(read_id[:-2])
    
    
    #  add the tags and reinstate the read_id
    read.tags += [('XC', barcode)]
    read.tags += [('XM', umi)]
    read.query_name = read_id


    #  write the read
    OUT.write(read)



