#!/usr/bin/perl
# USAGE:
# unix command line:
use strict;
use warnings;

my $inFile1 = shift @ARGV;    # read file names from terminal input
open IN1, '<', $inFile1 or die "Cannot open '$inFile1' because: $!";

my $inFile2 = shift @ARGV;    # read file names from terminal input
open IN2, '<', $inFile2 or die "Cannot open '$inFile2' because: $!";

my $outFile = shift @ARGV;    # read file names from termoutal output
open OUT, '>', $outFile or die "Cannot open '$outFile' because: $!";

my @maternalID;
my @paternalID;
my @hybridID;
my @parentID;

while (<IN2>){
chomp;
next if ($_ =~ m/^Variety/);
my @tmp=split("\t", $_);

push @hybridID, $tmp[0];
push @maternalID, $tmp[2];
push @paternalID, $tmp[4];
}# while

while (<IN1>){
chomp;
next if ($_ =~ /^#/);
my @tmpo=split("\t", $_);
my $i;
my $j;

#print OUT "$tmpo[0]\t$tmpo[1]\t$tmpo[2]\t";

if ($_ =~ m/^CHROM/) {
for ($i=3; $i<=$#tmpo; $i++){
push @parentID, $tmpo[$i];
}#for
print OUT "$tmpo[0]\t$tmpo[1]\t$tmpo[2]\t";
for ($j=0; $j<=$#hybridID; $j++){
print OUT "$hybridID[$j]\t";
}#for
print OUT "\n";
}#if

elsif ($_ !~ m/^CHROM/) {

print OUT "$tmpo[0]\t$tmpo[1]\t$tmpo[2]\t";
for ($j=0; $j<=$#hybridID; $j++){
my $maternalGenotype;
my $paternalGenotype;
my $hybridGenotype;

for ($i=3; $i<=$#tmpo; $i++){
if ($maternalID[$j] eq $parentID[$i-3]) {
$maternalGenotype = $tmpo[$i];
}#if
elsif ($paternalID[$j] eq $parentID[$i-3]) {
$paternalGenotype = $tmpo[$i];
}#elsif
}#for

if ($maternalGenotype eq "x" || $maternalGenotype eq "X" || $paternalGenotype eq "x" || $paternalGenotype eq "X"){
$hybridGenotype = "X";
}#if
else{
$hybridGenotype=$maternalGenotype.$paternalGenotype;
}#else

print OUT "$hybridGenotype\t";
}#for
print OUT "\n";

}#elsif
}#while

close IN1;
close IN2;
close OUT;



