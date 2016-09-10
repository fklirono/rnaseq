#!/usr/bin/perl
#  reports in comma-separated entries the positions of occurrences of [ACGT] letters in a string, if no occurrence NA is produced

use strict;
use warnings;

while (<>){
  my $string=$_;
  my @a;
  while ($string =~ /[ACGT]/g){
    push @a, pos($string);  #  record all positions of matches
  }
  if ($#a==-1){
    print "NA\n";
  } else {
    ($string = "@a") =~ s/ /,/g;  #  collapse positions into space separaed string and replace the spaces with commas
    print $string."\n";
  }
}
