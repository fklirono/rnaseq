#!/bin/bash
#
#  example:
#
#    t2c.sh Ableson_Bim_WT_WT /local_data/bim_data/par-clip/indices/WT.fa

READS=$1
REF=$2

if [[ ! -e "$2".fai ]]; then
  samtools faidx ${REF}
fi

#  filter in only reads with T:C conversion and mapped sense to our features and print out the position of the T:C conversion at STDERR
samtools view -F16 ${READS}.bam | t2c_sam.awk 2> ${READS}_T\:C.positions > ${READS}_T\:C.sam

#  SAM to BAM conversion
samtools view -bt ${REF}.fai ${READS}_T\:C.sam > ${READS}_T\:C.bam
samtools index ${READS}_T\:C.bam
rm ${READS}_T\:C.sam
