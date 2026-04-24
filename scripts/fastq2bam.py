#!/usr/bin/env python3
"""
convert_to_pacbio_bam.py

Optimized script to convert a FASTA or FASTQ file with reads
to a flnc.bam generated with IsoSeq3
"""

import os
import sys
import multiprocessing
from Bio import SeqIO
import pysam

print(r"""
    ____                          _     _          ____            ____  _         ____    _    __  __
   / ___|___  _ ____  _____ _  __| |_  | |_ ___   |  _ \ __ _  ___| __ )(_) ___   | __ )  / \  |  \/  |
  | |   / _ \| '_ \ \ / / _ \ '__| __| | __/ _ \  | |_) / _` |/ __|  _ \| |/ _ \  |  _ \ / _ \ | |\/| |
  | |__| (_) | | | \ V /  __/ |  | |_  | || (_) | |  __/ (_| | (__| |_) | | (_) | | |_) / ___ \| |  | |
   \____\___/|_| |_|\_/ \___|_|   \__|  \__\___/  |_|   \__,_|\___|____/|_|\___/  |____/_/   \_\_|  |_|
""")

def detect_format(file_path):
    """Detects if a file is FASTA or FASTQ based on the first character."""
    with open(file_path, "r") as f:
        first_char = f.readline().strip()[0]
        if first_char == ">":
            return "fasta"
        elif first_char == "@":
            return "fastq"
        else:
            raise ValueError(f"Unknown format: First character '{first_char}' is not FASTA (>) or FASTQ (@)")

def template_bam_generator(bam_path):
    """Yields alignment records from the template BAM, looping indefinitely."""
    while True:
        with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as in_sam:
            for aln in in_sam:
                yield aln

reads_file = snakemake.input[0]
bam_file = snakemake.params.bam
out_bam = snakemake.output[0]
threads = snakemake.threads
keep_files = True

if not os.path.exists(reads_file):
    print("ERROR: reads file does not exist. Provide a valid path", file=sys.stderr)
    sys.exit(1)

if not os.path.exists(bam_file):
    print("ERROR: bam file does not exist. Provide a valid path", file=sys.stderr)
    sys.exit(1)

if os.path.exists(out_bam):
    print("WARNING: output file already exists. Overwriting!", file=sys.stderr)

fmt = detect_format(reads_file)
print(f"Format detected: {fmt}")

# Initialize the infinite template generator
template_gen = template_bam_generator(bam_file)

# Open output SAM with threads for parallel compression
print(f"Converting reads using {threads} threads for BAM compression...")

# We grab the header from the template file once
with pysam.AlignmentFile(bam_file, "rb", check_sq=False) as header_source:
    header = header_source.header

with pysam.AlignmentFile(out_bam, "wb", header=header, threads=threads) as out_sam:
    # Stream the sequences one by one (No memory overload)
    for record in SeqIO.parse(reads_file, fmt):
        
        # Fetch the next template alignment
        aln = next(template_gen)

        seq_str = str(record.seq)
        seq_len = len(seq_str)

        # Replace ID and sequence
        aln.query_name = record.id
        aln.query_sequence = seq_str

        # Handle qualities
        if fmt == "fastq":
            aln.query_qualities = record.letter_annotations["phred_quality"]
        else:
            aln.query_qualities = pysam.qualitystring_to_array("~" * seq_len)

        # Update tags
        tags = dict(aln.tags)
        if "qe" in tags:
            tags["qe"] = seq_len
        aln.tags = list(tags.items())

        # Write out parallelized by pysam
        out_sam.write(aln)

print("SAM file created")
print("Completed")
