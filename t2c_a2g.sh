#!/bin/bash
#
#  the difference with t2c.sh is that this script considers both T:C for plus-strand and A:G for minus-strand conversions
#
#  example:
#
#    t2c_a2g.sh Abelson_Bim_mm9 /local_data/fklirono/systems/mm9/mm9.fa

READS=$1
REF=$2

if [[ ! -e "$2".fai ]]; then
  samtools faidx ${REF}
fi

#  filter in reads with T:C and A:G conversions and print out their positions
samtools view  ${READS}.bam | t2c_sam.awk 2> ${READS}_T\:C_A\:G.positions > ${READS}_T\:C_A\:G.sam

#  SAM to BAM conversion
samtools view -bt ${REF}.fai ${READS}_T\:C_A\:G.sam > ${READS}_T\:C_A\:G.bam
samtools index ${READS}_T\:C_A\:G.bam
rm ${READS}_T\:C_A\:G.sam
