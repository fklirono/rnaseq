#!/usr/bin/perl

use Getopt::Long;

my $usage="
tabulate values extracted from given column of file1, ...

usage: $0 --simple -c/--column= [default:1] -s/--separator= [default:'\\t'] file1, ...

Example: $0 -c 5 -s ',' file.csv
";

my $col=1;  #  column to read from [one-based]
my $sep="\t";  #  separator to use to determine the column
my $simple="";  #  report only 'value TAB multiplicity' without header or : separator
GetOptions('c|column=i' => \$col, 's|separator=s' => \$sep, 'simple' => \$simple);

if ( -t STDIN && !@ARGV ){die $usage}  #  no input and no pipe is provided, STDIN will wait for input from the terminal

my %v;  #  hash to hold the values extracted as keys

while(<>){  #  open piped STDIN or @ARGV 
  chomp;
  $v{(split(/$sep/, $_))[$col-1]}++;  #  split line by $sep and extract $col-1 (zero-based) as key, increasing its value by one in the hash
}

my @v_s = sort { $a <=> $b } keys %v;  

printf "%6s  : %6s\n", "value", "multiplicity" unless ($simple);
print "--------:-------------\n" unless ($simple);
foreach (@v_s){ 
  if ($simple){
    printf "%s\t%s\n", $_, $v{$_}
  } else {
    printf "%1s\t:\t%1s\n", $_, $v{$_}
  }
}

