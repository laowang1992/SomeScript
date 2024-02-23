#!/usr/bin/env perl
use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long qw(:config no_ignore_case bundling);

my $usage = <<__EOUSAGE__;

############################################################
#
# Usage:  $0 --fasta <cds_or_pep.fa> --gff <genes.gff> --out <outprefix>
#
# Required:
#
#	--fasta <string>			CDS or pep fasta file.
#
#	--gff   <string>			gff file.
#
#	--out   <string>			output prefix.
#
############################################################


__EOUSAGE__

    ;

my $help_flag;
my $fasta;
my $gff;
my $out;

&GetOptions('help|h' => \$help_flag,
            'fasta|f=s' => \$fasta,
            'gff|g=s' => \$gff,
            'out|o=s' => \$out,
            );

unless ($fasta && $gff && $out) {
	die $usage;
}

open GFF, $gff;

# 读取genome文件
my %mrna;
my $fa = Bio::SeqIO->new (-file =>$fasta, -f =>'fasta');
while (my $seq_obj = $fa->next_seq) {
	my $id = $seq_obj->id;
	my $seq = $seq_obj->seq;
	$mrna{$id}{'sequence'} = $seq;
	$mrna{$id}{'length'} = length $seq;
}

my %genes;
while(<GFF>){
	chomp;
	next unless $_;
	my @a = split /\t/, $_;
	next unless $a[2] eq "mRNA";
	my ($mrnaID, $geneID);
	if(/ID=(.+?);/ || /ID=(.+?)$/){
		$mrnaID = $1;
		#print "$mrnaID\t";
	}else{
		die "Error: There is no mRNA ID in line: $_\n";
	}
	if(/Parent=(.+?)\;/ || /Parent=(.+?)$/){
		$geneID = $1;
		#print "$geneID\n";
	}else{
		die "Error: There is no gene ID in line: $_\n";
	}
	if(!exists($genes{$geneID}) || $genes{$geneID}{'length'} < $mrna{$mrnaID}{'length'}){
		$genes{$geneID} = {
			"mrnaID" => $mrnaID,
			"length" => $mrna{$mrnaID}{'length'},
		};
	}
}
close GFF;

my $outfasta = join ".", $out, "longest.fa";
my $outgene = join ".", $out, "longest.list";
open OUTFA, ">$outfasta";
open OUTGENE, ">$outgene";

foreach my $gene (keys %genes) {  
	print OUTFA ">$gene\n$mrna{$genes{$gene}{'mrnaID'}}{'sequence'}\n";
	print OUTGENE "$gene\t$genes{$gene}{'mrnaID'}\t$genes{$gene}{'length'}\n";
}

close OUTFA;
close OUTGENE;
