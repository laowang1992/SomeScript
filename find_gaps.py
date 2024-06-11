import sys
import os
import argparse
from Bio import SeqIO

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

def main():
    parser = argparse.ArgumentParser(description="Find gaps (N regions) in a genome FASTA file.")
    parser.add_argument("-i", "--input", required=True, help="Input genome FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output file for gap positions")

    args = parser.parse_args()

    if not os.path.isfile(args.input):
        print(f"Error: Input file '{args.input}' not found.")
        sys.exit(1)

    find_gaps(args.input, args.output)
    print(f"Gap positions have been written to {args.output}")

if __name__ == "__main__":
    main()
