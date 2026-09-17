#!/usr/bin/env python3
"""Reference-free noise diagnostic of the reconstructed transcriptome.

Computed from the SQANTI3 classification alone, before the rules filter runs, so it is
available for any genome (no reference annotation needed). The four fractions are the ones
that tracked the reference-based junk fraction (class codes ``u`` + ``i`` against a curated
reference) on the eight benchmark species; the non-coding fraction is the best single
predictor (Spearman rho 0.88), and per isoform 82-91 % of the junk models are non-coding.

Metrics (share of transcript models):
  noncoding_frac   ``coding`` column != ``coding``
  fl2_frac         FL <= 2 (minimal read support)
  monoexon_frac    exons == 1
  multifail_frac   failing >= 3 of the six strict rules (intra-priming, canonical junctions,
                   coding, CDS >= 100 aa, not NMD, psauron >= 0.7); NA counts as pass for the
                   string columns and as fail for the numeric ones, as SQANTI3 does.

Benchmark scale for ``noncoding_frac`` (default-filter runs, September 2026):
  Arabidopsis 2.6 %, fly 12.0 %, yeast 14.1 %, human 16.7 %, worm 24.2 %, chicken 25.3 %,
  mouse 35.3 %, zebrafish 60.8 % (zebrafish: 56 % of its models were junk against the reference).

Why it matters: the default rules filter rescues isoforms that fail one of the intra-priming,
psauron or canonical-junction checks when they carry >= 10 full-length reads (``FL`` in
``envs/filter_rules.json``). In a noisy transcriptome that gate lets more junk through
(zebrafish: +973 ``u``/``i`` isoforms for +891 reference-exact ones); raising the ``FL`` value of
the three rescue paths is the knob.

Runs as a Snakemake ``script:`` (uses the ``snakemake`` object) or standalone::

    noise_metrics.py --classification sp_classification.txt --out noise_metrics.tsv [--sample sp]
"""
import argparse
import csv
import sys

SCALE = [("arabidopsis", 2.6), ("drosophila", 12.0), ("yeast", 14.1), ("human", 16.7),
         ("cElegans", 24.2), ("gallus", 25.3), ("mouse", 35.3), ("danio", 60.8)]
CLEAN_MAX, NOISY_MIN = 20.0, 40.0   # % non-coding: below = clean, above = noisy


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _in_range(v, lo, hi):
    x = _num(v)
    return x is not None and lo <= x <= hi


# The six strict rules with SQANTI3's NA semantics (blank string column passes, numeric NA fails).
STRICT_RULES = {
    "intrapriming": lambda r: _in_range(r.get("perc_A_downstream_TTS"), 0, 59),
    "canonical":    lambda r: r.get("all_canonical") in ("canonical", "", "NA", None),
    "coding":       lambda r: r.get("coding") == "coding",
    "CDS_length":   lambda r: (_num(r.get("CDS_length")) or -1) >= 100,
    "not_NMD":      lambda r: r.get("predicted_NMD") in ("FALSE", "False", "", "NA", None),
    "psauron":      lambda r: (_num(r.get("psauron_score")) or -1) >= 0.7,
}


def compute_metrics(rows):
    """Return the metrics dict for a list of classification rows (dicts)."""
    n = len(rows)
    if n == 0:
        return {"models": 0, "noncoding_frac": 0.0, "fl2_frac": 0.0, "monoexon_frac": 0.0, "multifail_frac": 0.0}
    noncoding = sum(r.get("coding") != "coding" for r in rows)
    fl2 = sum((_num(r.get("FL")) or 0) <= 2 for r in rows)
    mono = sum((_num(r.get("exons")) or 0) == 1 for r in rows)
    multifail = sum(sum(not rule(r) for rule in STRICT_RULES.values()) >= 3 for r in rows)
    return {"models": n, "noncoding_frac": 100 * noncoding / n, "fl2_frac": 100 * fl2 / n,
            "monoexon_frac": 100 * mono / n, "multifail_frac": 100 * multifail / n}


def verdict(noncoding_pct):
    if noncoding_pct < CLEAN_MAX:
        return "clean"
    if noncoding_pct < NOISY_MIN:
        return "moderate"
    return "noisy"


def message(sample, m):
    v = verdict(m["noncoding_frac"])
    nearest = min(SCALE, key=lambda s: abs(s[1] - m["noncoding_frac"]))
    lines = [
        f"[noise diagnostic] {sample}: {m['models']:,} transcript models; non-coding {m['noncoding_frac']:.1f} %, "
        f"FL<=2 {m['fl2_frac']:.1f} %, mono-exon {m['monoexon_frac']:.1f} %, failing >=3 strict rules {m['multifail_frac']:.1f} %.",
        f"[noise diagnostic] verdict: {v} transcriptome (benchmark scale: "
        + ", ".join(f"{s} {p:.1f} %" for s, p in SCALE) + f"; closest: {nearest[0]}).",
    ]
    if v == "noisy":
        lines.append("[noise diagnostic] The FL>=10 rescue paths of the default filter will admit more unsupported "
                     "isoforms here; consider raising the 'FL' value of the three rescue paths in the rules JSON "
                     "(curation.filter_rules), or use envs/filter_rules.strict.json.")
    elif v == "moderate":
        lines.append("[noise diagnostic] Default filter appropriate; inspect the rescued isoforms if precision matters more than sensitivity.")
    else:
        lines.append("[noise diagnostic] Default filter appropriate.")
    return "\n".join(lines)


def write_table(path, sample, m):
    with open(path, "w") as fh:
        fh.write("sample\tmodels\tnoncoding_frac\tfl2_frac\tmonoexon_frac\tmultifail_frac\tverdict\n")
        fh.write("\t".join([sample, str(m["models"])] + [f"{m[k]:.1f}" for k in
                 ("noncoding_frac", "fl2_frac", "monoexon_frac", "multifail_frac")] + [verdict(m["noncoding_frac"])]) + "\n")


def run(classification, out, sample, log=None):
    with open(classification) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    m = compute_metrics(rows)
    write_table(out, sample, m)
    msg = message(sample, m)
    if log:
        with open(log, "w") as fh:
            fh.write(msg + "\n")
    print(msg, file=sys.stderr)
    return m


def main_snakemake():
    smk = globals()["snakemake"]
    run(smk.input.classification, smk.output[0], smk.params.sample, smk.log[0] if smk.log else None)


def main_cli(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--classification", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", default="sample")
    ap.add_argument("--log")
    a = ap.parse_args(argv)
    run(a.classification, a.out, a.sample, a.log)


if __name__ == "__main__":
    if "snakemake" in globals():
        main_snakemake()
    else:
        main_cli()
