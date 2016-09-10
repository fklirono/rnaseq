#!/usr/bin/awk -f
#
#  this script prints out to STDOUT all reads from STDIN that
#
#    are plus-strand mapped uniquely with a single T:C mismatch 
#    are negative-strand mapped uniquely with a single plus-strand A:G mismatch
#
#  additionally, it prints to STDERR the reference, strand, and position of the corresponding T:C or A:G conversion
#
#  Running example: 
#
#    samtools view <READS.bam> | t2c_sam.awk
#    samtools -F16 view <READS.bam> | t2c_sam.awk  #  only plus-strand T:C conversions

!($12=="X1:i:0" && $13=="XM:i:1" && $19~/[^^][TA]/) {next}
{
  N="T"  #  preset nucleotide choice to T
  M="C"  #  only mismatch to look for
  S="+"  #  preset strand to +
  x=substr($19, 6)  #  remove MD:Z: part
}
and($2, 16) {
  N="A"  #  negative strand mapped-read needs to reset nucleotide to A
  M="G"  #  only mismatch to look for in this case
  S="-"  #  strand is -
}
split(x, y, N)!=2{next}  #  no plus-strand mapped reads with A mismatch, or minus-strand reads with T mismatch, or deletion of C or G
{
  nucl=substr($10,y[1]+1,1)  #  even if N is first letter this will work and set nucl to the read nucleotide that N converted to
}
nucl==M {
  print $0 #  print the whole line if T:C of A:G conversion
  print $3"\t"S"\t"$4+y[1] > "/dev/stderr"  #  print to STDERR: reference name, strand, position of T:C or A:G conversion
}
