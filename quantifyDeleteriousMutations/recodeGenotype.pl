#!/usr/bin/perl
# USAGE:
# unix command line:
use strict;
use warnings;

my $inFile = shift @ARGV;    # read file names from terminal input
open IN, '<', $inFile or die "Cannot open '$inFile' because: $!";

my $outFile = shift @ARGV;    # read file names from termoutal output
open OUT, '>', $outFile or die "Cannot open '$outFile' because: $!";


while (<IN>){
chomp;
next if ($_ =~ m/^CHROM/);
my $delAlleleCount;

my @tmpo=split("\t", $_);
next if ($tmpo[1] <= 0 );
print OUT "$tmpo[0]\t$tmpo[1]\t$tmpo[2]\t";

for (my $i=3; $i<=$#tmpo; $i++){

if ($tmpo[$i] eq "x") {
$delAlleleCount="NA";
}#if

elsif ($tmpo[$i] eq "X") {
$delAlleleCount=1;
}#elsif

elsif ($tmpo[$i] eq $tmpo[5]) {
$delAlleleCount=0;
}#elsif

else {
$delAlleleCount=2;
}#elsif

#else{
#print "line $. abnormal\n";
#}
print OUT "$delAlleleCount\t";
}#for

print OUT "\n";


}#while

close IN;
close OUT;




