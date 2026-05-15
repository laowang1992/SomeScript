#!/usr/bin/env python3

import re
import argparse
import math


def revcomp(seq):
    comp = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
    return "".join(comp.get(b, "N") for b in seq[::-1])


def rotations(seq):
    return {seq[i:] + seq[:i] for i in range(len(seq))}


def parse_attr(attr, key):
    m = re.search(rf"{key}=([^;]+)", attr)
    return m.group(1) if m else None


def main():
    parser = argparse.ArgumentParser(
        description="Filter telomere tandem repeats from TRF gff3"
    )

    parser.add_argument("-i", "--input", required=True, help="Input GFF3")
    parser.add_argument("-o", "--output", required=True, help="Output file")

    parser.add_argument("-m", "--motif", default="TTTAGGG",
                        help="Telomere motif (default: TTTAGGG)")

    parser.add_argument("-c", "--copy", type=float, default=50,
                        help="Minimum CopyNumber (default: 50)")

    parser.add_argument("--match", type=float, default=0,
                        help="Minimum PercentMatches (default: 0)")

    parser.add_argument("--indel", type=float, default=math.inf,
                        help="Maximum PercentIndels (default: inf)")

    args = parser.parse_args()

    motif = args.motif.upper()

    # motif + reverse complement + all rotations
    motifs = set()
    for seq in [motif, revcomp(motif)]:
        motifs |= rotations(seq)

    with open(args.input) as fin, open(args.output, "w") as fout:

        for line in fin:
            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")
            if len(cols) < 9:
                continue

            attr = cols[8]

            consensus = parse_attr(attr, "Consensus")
            copy = parse_attr(attr, "CopyNumber")
            match = parse_attr(attr, "PercentMatches")
            indel = parse_attr(attr, "PercentIndels")

            if not (consensus and copy and match and indel):
                continue

            consensus = consensus.upper()

            try:
                copy = float(copy)
                match = float(match)
                indel = float(indel)
            except:
                continue

            if (consensus in motifs and
                copy >= args.copy and
                match >= args.match and
                indel <= args.indel):

                fout.write(line)


if __name__ == "__main__":
    main()
