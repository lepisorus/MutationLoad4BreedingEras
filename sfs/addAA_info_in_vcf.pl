
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



my @tmpo=split("\t", $_);

my $INFO;

print OUT "$tmpo[2]\t$tmpo[3]\t$tmpo[4]\t$tmpo[0]\t$tmpo[5]\t$tmpo[6]\t$tmpo[7]\t$tmpo[8]\t";

$INFO="AA=".$tmpo[1].";".$tmpo[9];
print OUT "$INFO\t";

for (my $i=10; $i<=$#tmpo; $i++){

print OUT "$tmpo[$i]\t";
}#for

print OUT "\n";


}#while

close IN;
close OUT;



