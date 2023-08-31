#!/bin/env perl
# need a file, one marker per line, tab seperated, e.g.:
#primer_ID	Left_primer_seq	Right_primer_seq

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<__EOUSAGE__;

############################################################
#
# Usage:  $0 --input primer.txt --output output.txt
#
# Required:
#
#	--input <string>			input filename, one pair primer per line, tab seperated, e.g.:
#								primerID	Left_primer_seq	Right_primer_seq
#
#	--output <string>			output filename.
#
############################################################


__EOUSAGE__

;

my $help_flag;
my $input;
my $output;

&GetOptions('help|h' => \$help_flag,
            'input|i=s' => \$input,
            'output|o=s' => \$output,
            );

unless ($input && $output) {
	die $usage;
}
open IN, "$input";
open OUT, ">$output";

print OUT "Pimer\thit_num\tChrom\tstrand\tfrom\tto\tmism\tgaps\tact_len/exp_len\n";
while(<IN>){
	chomp;
	my @a = split /\t/, $_;
	my @epcr = readpipe("re-PCR -s genome.hash -n 1 -g 1 $a[1] $a[2] 50-1000");
	my $hit = @epcr;
	$hit -= 2;
	print OUT "$a[0] has no hit\n" if $hit == 0;
	if($hit > 20){
		print OUT "$a[0] has to $hit hits!!! It's not specific\n";
		next;
	}
	my $i = 0;
	while($i < $hit){
		chomp $epcr[$i+1];
		my @b = split /\t/, $epcr[$i+1];
		print OUT "$a[0]\t$hit\t$b[1]\t$b[2]\t$b[3]\t$b[4]\t$b[5]\t$b[6]\t$b[7]\n";
		$i++;
	}
}

close IN;
