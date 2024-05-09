#!/bin/env perl
#####################################################################
# zcat pop.recode.vcf.gz | grep -v "##" | cut -f1,2,10- > geno.txt  #
# perl geno2Numeric.pl geno.txt GD.tmp GM.txt                       #
# csvtk transpose -t -T GD.tmp -o GD.txt                            #
#####################################################################
use strict;
use warnings;

open IN, "$ARGV[0]";
open GD, ">$ARGV[1]";
open GM, ">$ARGV[2]";
my $header = <IN>;
$header =~ s/^#CHROM\tPOS/taxa/;
print GD "$header";
print GM "SNP\tChromosome\tPosition\n";

while(<IN>){
	chomp;
	my @a = split /\t/, $_;
	my $snp = join "_", $a[0], $a[1];
	print GM "$snp\t$a[0]\t$a[1]\n";
	shift @a for (1, 2);
	print GD "$snp";
	foreach my $geno (@a){
		my @b = split /\|/, $geno;
		my $code;
		if(($b[0] == 0 || $b[0] == 1) && ($b[1] == 0 || $b[1] == 1)) {
			$code = ($b[0] + $b[1]);
			print GD "\t$code";
		} else {
			print GD "\t";
		}
	}
	print GD "\n";
}

close IN;
close GD;
close GM;
