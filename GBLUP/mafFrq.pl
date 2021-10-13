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
my $allele1 = $tmpo[2];
my $allele1_frq = $tmpo[3];
my $allele2 = $tmpo[4];
my $allele2_frq = $tmpo[5];




my $frq;

if($allele1_frq > $allele2_frq){
$frq = $allele2_frq;
}#out if

elsif($allele1_frq < $allele2_frq){
$frq =  $allele1_frq;
}#elsif

elsif($allele1_frq = $allele2_frq){
$frq =  $allele1_frq;
}#elsif

else{
$frq = "NA";
}
print OUT "$frq\t";

for (my $i=6; $i<=$#tmpo; $i++){
print OUT "$tmpo[$i]\t";
}#for
print OUT "\n";
} #while
   
close IN;
close OUT;





