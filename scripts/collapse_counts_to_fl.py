#!/usr/bin/env python3
"""Convert the abundance table of `isoseq collapse` into the FL-count matrix SQANTI3 expects.

SQANTI3's ``-fl`` option takes a two-column table, ``feature_id<TAB>count``, which is the shape of
IsoQuant's ``*.discovered_transcript_counts.tsv``. ``isoseq collapse`` instead writes
``<prefix>.abundance.txt``: three ``#`` comment lines, then a header
``pbid, count_fl, fl_assoc, cell_barcodes``.

The count taken is **fl_assoc**, not ``count_fl``. With clustering enabled (always, in this pipeline)
``count_fl`` counts the clusters supporting a collapsed model, while ``fl_assoc`` counts the
underlying full-length reads, which is what IsoQuant's column means and therefore what keeps the two
reconstruction routes comparable.

Runs as a Snakemake ``script:`` (uses the ``snakemake`` object) or standalone::

    collapse_counts_to_fl.py --abundance sample.collapsed.abundance.txt --out sample.fl_counts.tsv
"""
import argparse
import sys

COUNT_COLUMN = "fl_assoc"
ID_COLUMN = "pbid"


def read_abundance(path):
    """Yield (transcript_id, count) from an `isoseq collapse` abundance table."""
    with open(path) as fh:
        header = None
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                for column in (ID_COLUMN, COUNT_COLUMN):
                    if column not in header:
                        raise ValueError(
                            f"{path}: missing column '{column}' in the abundance header {header}. "
                            "Is this the output of `isoseq collapse`?"
                        )
                id_i, count_i = header.index(ID_COLUMN), header.index(COUNT_COLUMN)
                continue
            if len(fields) <= max(id_i, count_i):
                raise ValueError(f"{path}: truncated row {fields!r}")
            yield fields[id_i], fields[count_i]
    if header is None:
        raise ValueError(f"{path}: no header found; the file holds only comments or is empty.")


def convert(abundance, out_path, log=sys.stderr):
    """Write the two-column FL matrix; returns (transcripts, total counted reads)."""
    rows = list(read_abundance(abundance))
    if not rows:
        raise ValueError(f"{abundance}: no transcripts listed; `isoseq collapse` produced nothing to curate.")
    total = 0
    with open(out_path, "w") as out:
        out.write("feature_id\tcount\n")
        for transcript_id, count in rows:
            out.write(f"{transcript_id}\t{count}\n")
            try:
                total += float(count)
            except ValueError:
                raise ValueError(f"{abundance}: non-numeric {COUNT_COLUMN} {count!r} for {transcript_id}")
    print(f"Collapsed transcripts: {len(rows)}", file=log)
    print(f"Full-length reads assigned ({COUNT_COLUMN}): {total:.0f}", file=log)
    print("Compare that total with the number of input reads: a large shortfall means the reads were "
          "deduplicated upstream (for example by a shared ZMW in converted BAMs).", file=log)
    return len(rows), total


def main_snakemake():
    with open(snakemake.log[0], "w") as log:  # noqa: F821 (injected by Snakemake)
        convert(abundance=snakemake.input.abundance,  # noqa: F821
                out_path=snakemake.output[0],  # noqa: F821
                log=log)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--abundance", required=True, help="<prefix>.collapsed.abundance.txt from `isoseq collapse`")
    p.add_argument("--out", required=True, help="two-column FL matrix for SQANTI3 -fl")
    a = p.parse_args(argv)
    convert(a.abundance, a.out)


if __name__ == "__main__":
    if "snakemake" in globals():
        main_snakemake()
    else:
        main()
