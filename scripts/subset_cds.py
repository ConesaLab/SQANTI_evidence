#!/usr/bin/env python3

import sys

def subset_cds(input_gtf, output_gtf):
    with open(input_gtf, "r") as infile, open(output_gtf, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                continue  # skip headers/comments
            fields = line.strip().split("\t")
            if len(fields) > 2 and fields[2] == "CDS":
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input.gtf> <output.gtf>")
        sys.exit(1)

    input_gtf = sys.argv[1]
    output_gtf = sys.argv[2]
    subset_cds(input_gtf, output_gtf)

