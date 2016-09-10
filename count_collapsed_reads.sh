#!/bin/bash
#  $1 : BAM file of collapsed reads where multiplicities are appended at the end of the readID as: _xNNN
#  $2 : the genomic coordinate in samtools format to subset the collapsed reads, e.g. chr2:127951774-127988283

if [[ $# -lt 1 ]]
then
  echo -ne "Please provide the BAM file to count collapsed reads from\nand optionally the genomic range to focus on,  e.g.\n\n"
  echo -ne "  $(basename $0) some_alignment.bam chr2:127951774-127988283\n\n"
  exit 255
fi

samtools view $1 $2 | cut -f1 | gawk 'BEGIN{sum=0} {match($1, "_x([[:digit:]]+)$", a); sum+=a[1]}END{print sum}'
