import sys
import os
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def find_gaps(fasta_file, output_file):
    with open(fasta_file, "r") as input_handle, open(output_file, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            sequence = str(record.seq)
            chrom = record.id
            gap_start = None

            for i, base in enumerate(sequence):
                if base == 'N':
                    if gap_start is None:
                        gap_start = i
                else:
                    if gap_start is not None:
                        output_handle.write(f"{chrom}\t{gap_start+1}\t{i}\n")
                        gap_start = None

            if gap_start is not None:
                output_handle.write(f"{chrom}\t{gap_start+1}\t{len(sequence)}\n")

def split_fasta(fasta_file, contig_file):
    with open(fasta_file, "r") as input_handle, open(contig_file, "w") as contig_handle:
        for record in SeqIO.parse(input_handle, "fasta"):
            sequence = str(record.seq)
            chrom = record.id
            gap_start = None
            contig_start = 0

            for i, base in enumerate(sequence):
                if base == 'N':
                    if gap_start is None:
                        gap_start = i
                else:
                    if gap_start is not None:
                        if contig_start < gap_start:
                            new_seq = sequence[contig_start:gap_start]
                            new_header = f"{chrom}_{contig_start+1}_{gap_start}"
                            new_record = SeqRecord(Seq(new_seq), id=new_header, description="")
                            SeqIO.write(new_record, contig_handle, "fasta")
                        contig_start = i + 1
                        gap_start = None

            if contig_start < len(sequence):
                new_seq = sequence[contig_start:]
                new_header = f"{chrom}_{contig_start+1}_{len(sequence)}"
                new_record = SeqRecord(Seq(new_seq), id=new_header, description="")
                SeqIO.write(new_record, contig_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Find gaps (N regions) in a genome FASTA file and optionally split at gaps.")
    parser.add_argument("-i", "--input", required=True, help="Input genome FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output file for gap positions")
    parser.add_argument("-s", "--split", action="store_true", help="Split the FASTA file at gaps")
    parser.add_argument("-c", "--contig", help="Output file for split contigs (default: split.fa)")

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)

    find_gaps(args.input, args.output)
    print(f"Gap positions have been written to {args.output}")

    if args.split:
        contig_file = args.contig if args.contig else "split.fa"
        split_fasta(args.input, contig_file)
        print(f"Split contigs have been written to {contig_file}")

if __name__ == "__main__":
    main()
