#!/usr/bin/env python3
"""Wrap FASTA/FASTQ reads in an unmapped PacBio BAM, using a template BAM for the record structure.

`isoseq cluster2` only reads PacBio BAM, so a FASTA/FASTQ input to the IsoSeq reconstruction route has
to be converted first. Records are taken from a template (``envs/pacbio_mock.bam``) in a loop, and for
each read the identifier, sequence and base qualities are replaced.

**What is real and what is not.** Read names, sequences and (for FASTQ) base qualities come from the
input. Everything else - the read group, movie name, instrument, and the per-read PacBio tags such as
``zm`` (ZMW), ``np`` (passes) and ``rq`` (predicted accuracy) - is copied from the template and is
therefore *not* provenance. Two consequences worth knowing:

* Any tool that filters on ``rq``/``np``, or that assumes one record per ZMW, operates on template
  values and will behave uniformly rather than per read.
* The header declares ``READTYPE=CCS`` whatever the reads actually are, so subreads converted this way
  are accepted silently downstream. Check the reads are really CCS/HiFi first (``samtools stats``:
  HiFi is below ~2 % error). `scripts/input_check.py` warns about this at configuration time.

Adapted from `scripts/fastq2bam.py` of the pre-refactor `dev-IsoSeq` branch.

Runs as a Snakemake ``script:`` (uses the ``snakemake`` object) or standalone::

    fastq2bam.py --reads reads.fastq --template envs/pacbio_mock.bam --out sample.bam [--threads 4]
"""
import argparse
import sys

FASTA_START, FASTQ_START = ">", "@"


def detect_format(path):
    """'fasta' or 'fastq', from the first non-empty character of the file."""
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == FASTA_START:
                return "fasta"
            if line[0] == FASTQ_START:
                return "fastq"
            raise ValueError(
                f"{path}: first character {line[0]!r} is neither '{FASTA_START}' (FASTA) nor "
                f"'{FASTQ_START}' (FASTQ)."
            )
    raise ValueError(f"{path}: file is empty.")


def template_records(bam_path):
    """Yield alignment records from the template BAM, looping for as long as reads keep coming."""
    import pysam

    while True:
        with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as template:
            empty = True
            for record in template:
                empty = False
                yield record
            if empty:
                raise ValueError(f"{bam_path}: template BAM holds no records.")


def convert(reads, template, out_bam, threads=1, log=sys.stderr):
    """Write the reads as an unmapped PacBio BAM; returns the number of records written."""
    import pysam
    from Bio import SeqIO

    fmt = detect_format(reads)
    print(f"Input format: {fmt}", file=log)
    print(f"Template records: {template}", file=log)

    records = template_records(template)
    written = 0
    with pysam.AlignmentFile(template, "rb", check_sq=False) as header_source:
        header = header_source.header
    with pysam.AlignmentFile(out_bam, "wb", header=header, threads=threads) as out:
        for read in SeqIO.parse(reads, fmt):
            aln = next(records)
            sequence = str(read.seq)
            aln.query_name = read.id
            aln.query_sequence = sequence
            if fmt == "fastq":
                aln.query_qualities = read.letter_annotations["phred_quality"]
            else:
                aln.query_qualities = pysam.qualitystring_to_array("~" * len(sequence))
            tags = dict(aln.tags)
            if "qe" in tags:          # query end: the only template tag that must track the sequence
                tags["qe"] = len(sequence)
            aln.tags = list(tags.items())
            out.write(aln)
            written += 1

    if not written:
        raise ValueError(f"{reads}: no reads found, nothing was written to {out_bam}.")
    print(f"Records written: {written}", file=log)
    print("NOTE: read names, sequences and qualities are from the input; read group, movie, instrument "
          "and the per-read PacBio tags (zm/np/rq) come from the template and are not provenance.", file=log)
    return written


def main_snakemake():
    with open(snakemake.log[0], "w") as log:  # noqa: F821 (injected by Snakemake)
        convert(reads=snakemake.input[0],  # noqa: F821
                template=snakemake.params.template,  # noqa: F821
                out_bam=snakemake.output[0],  # noqa: F821
                threads=snakemake.threads,  # noqa: F821
                log=log)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--reads", required=True, help="FASTA or FASTQ file of full-length reads")
    p.add_argument("--template", required=True, help="PacBio BAM whose records/header are reused")
    p.add_argument("--out", required=True, help="output unmapped PacBio BAM")
    p.add_argument("--threads", type=int, default=1, help="BAM compression threads")
    a = p.parse_args(argv)
    convert(a.reads, a.template, a.out, a.threads)


if __name__ == "__main__":
    if "snakemake" in globals():
        main_snakemake()
    else:
        main()
