#!/usr/bin/perl
# USAGE:
# unix command line:
use strict;
use warnings;

my $inFile = shift @ARGV;    # read file names from terminal input
open IN, '<', $inFile or die "Cannot open '$inFile' because: $!";

my $outFile = shift @ARGV;    # read file names from terminal input
open OUT, '>', $outFile or die "Cannot open '$inFile' because: $!";


while (<IN>){
chomp;

my @tmpo=split("\t", $_);
print OUT "$tmpo[0]\t$tmpo[1]\t";
for (my $i=2; $i<=$#tmpo; $i++){
my $geno;

if ($tmpo[$i] =~ /^0\/0/){
$geno=0;
}#inner if
elsif ($tmpo[$i] =~ /^1\/1/){
$geno=2;
}#
elsif ($tmpo[$i] =~ /^0\/1/){
$geno=1;
}#
elsif ($tmpo[$i] =~ /^\.\/\./){
$geno=9;
}#

else{
print "dummy, you are making mistakes!";
}
print OUT "$geno\t";
} #for

print OUT "\n";


} #while
   
close IN;
close OUT;





