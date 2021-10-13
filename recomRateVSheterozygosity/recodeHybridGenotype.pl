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
if ($_ =~ m/^chr/){
print OUT $_"\n";
}

else{
my @tmpo=split("\t", $_);
next if ($tmpo[3] ne $tmpo[5] && $tmpo[3] ne $tmpo[6]);

my $delAlleleCount;

print OUT "$tmpo[0]\t$tmpo[1]\t$tmpo[2]\t$tmpo[3]\t$tmpo[4]\t$tmpo[5]\t$tmpo[6]\t";

for (my $i=7; $i<=$#tmpo; $i++){

if ($tmpo[$i] eq "x") {
$delAlleleCount="NA";
}#if

elsif ($tmpo[$i] eq $tmpo[6].$tmpo[6]) {
$delAlleleCount=0;
}#elsif

elsif ($tmpo[$i] eq $tmpo[5].$tmpo[5]) {
$delAlleleCount=0;
}#elsif

elsif ($tmpo[$i] eq $tmpo[5].$tmpo[6]) {
$delAlleleCount=1;
}#elsif

elsif ($tmpo[$i] eq $tmpo[6].$tmpo[5]) {
$delAlleleCount=1;
}#elsif


else {
print "something is wrong\n";
}#elsif


print OUT "$delAlleleCount\t";
}#for

print OUT "\n";


}#else



}#while

close IN;
close OUT;



