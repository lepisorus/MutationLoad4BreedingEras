




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
if ($_ =~ m/^POS/){
print OUT "$_\t\n";
}

else{
my @tmpo=split("\t", $_);
next if ($tmpo[1] ne $tmpo[3] && $tmpo[1] ne $tmpo[4]);
my $delAlleleCount;

print OUT "$tmpo[0]\t$tmpo[1]\t$tmpo[2]\t$tmpo[3]\t$tmpo[4]\t";

for (my $i=5; $i<=$#tmpo; $i++){

if ($tmpo[$i] eq "x") {
$delAlleleCount="NA";
}#if

elsif ($tmpo[$i] eq "X") {
$delAlleleCount=1;
}#elsif

elsif ($tmpo[$i] eq $tmpo[1]) {
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


}#else



}#while

close IN;
close OUT;



