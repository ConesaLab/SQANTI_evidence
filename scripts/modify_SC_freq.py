"""
Set the stop-codon frequencies of an Augustus species from the etraining output.

`etraining` reports the observed frequency of the three stop codons at the end of its output;
`rule extract_stop_codon_freq` keeps those three lines (SC_freq.txt). This script rewrites the
corresponding `/Constant/{amber,ochre,opal}prob` parameters of
`$AUGUSTUS_CONFIG_PATH/species/<species>/<species>_parameters.cfg`.

The parameter lines are matched by name and their value column is replaced whatever its
current content, so the script gives the same result on a first run and on a rerun (the old
version substituted the literal template defaults "0.33"/"0.34" and silently did nothing once
they had been changed).

Runs as a Snakemake `script:` (uses the `snakemake` object) or standalone:
    python modify_SC_freq.py <SC_freq.txt> <parameters.cfg> [--copy-to <file>]
"""
import argparse
import os
import re
import shutil
import sys

# Augustus parameter -> codon name as printed by etraining
PARAM_TO_CODON = {
    "/Constant/amberprob": "tag",
    "/Constant/ochreprob": "taa",
    "/Constant/opalprob": "tga",
}

# e.g. "tag:   313 (0.207)"  ->  codon "tag", freq "0.207"
FREQ_LINE_RE = re.compile(r"^\s*(?P<codon>tag|taa|tga)\s*:\s*\d+\s*\(\s*(?P<freq>[0-9.]+)\s*\)", re.IGNORECASE)

# e.g. "/Constant/amberprob                   0.33   # Prob(stop codon = tag) ..."
PARAM_LINE_RE = re.compile(
    r"^(?P<key>/Constant/(?:amberprob|ochreprob|opalprob))(?P<ws>\s+)(?P<value>\S+)(?P<rest>.*)$"
)


def parse_stop_codon_freqs(path):
    """Returns {"tag": "0.207", "taa": ..., "tga": ...} from the etraining excerpt."""
    freqs = {}
    with open(path) as f_in:
        for line in f_in:
            m = FREQ_LINE_RE.match(line)
            if m:
                freqs[m.group("codon").lower()] = m.group("freq")
    missing = sorted(set(PARAM_TO_CODON.values()) - set(freqs))
    if missing:
        raise ValueError(
            f"Could not find stop codon frequencies for {missing} in {path}. "
            "Expected three lines like 'tag:   313 (0.207)' from the end of the etraining output."
        )
    return freqs


def rewrite_line(line, freqs):
    """Returns (new_line, changed). Replaces only the value column of the three stop-codon parameters."""
    m = PARAM_LINE_RE.match(line.rstrip("\n"))
    if not m:
        return line, False
    codon = PARAM_TO_CODON[m.group("key")]
    return f"{m.group('key')}{m.group('ws')}{freqs[codon]}{m.group('rest')}\n", True


def rewrite_parameters(config_file, freqs):
    """Rewrites config_file in place (via a temp file). Raises if not exactly the three parameters were found."""
    tmp_file = config_file + ".tmp"
    changed = 0
    with open(config_file) as f_in, open(tmp_file, "w") as f_out:
        for line in f_in:
            new_line, was_changed = rewrite_line(line, freqs)
            changed += was_changed
            f_out.write(new_line)
    if changed != len(PARAM_TO_CODON):
        os.remove(tmp_file)
        raise ValueError(
            f"Expected to rewrite {len(PARAM_TO_CODON)} stop-codon parameters in {config_file}, "
            f"found {changed}. Is this a valid Augustus species parameters file?"
        )
    os.replace(tmp_file, config_file)
    return changed


def species_parameters_file(species):
    config_path = os.environ.get("AUGUSTUS_CONFIG_PATH")
    if not config_path:
        raise EnvironmentError("AUGUSTUS_CONFIG_PATH is not set; run inside the Augustus conda environment.")
    return os.path.join(config_path, "species", species, f"{species}_parameters.cfg")


def run(freq_file, config_file, copy_to=None, sentinel=None, log=sys.stderr):
    freqs = parse_stop_codon_freqs(freq_file)
    rewrite_parameters(config_file, freqs)
    log.write(f"Stop codon frequencies written to {config_file}: "
              + ", ".join(f"{c}={freqs[c]}" for c in ("tag", "taa", "tga")) + "\n")
    if copy_to:
        os.makedirs(os.path.dirname(os.path.abspath(copy_to)), exist_ok=True)
        shutil.copyfile(config_file, copy_to)
    if sentinel:
        with open(sentinel, "w") as f_out:
            f_out.write(f"{config_file} modified: " + " ".join(f"{c}={freqs[c]}" for c in freqs) + "\n")


def main_snakemake():
    species = snakemake.params.name  # noqa: F821 (injected by Snakemake)
    config_file = species_parameters_file(species)
    with open(snakemake.log[0], "w") as log:  # noqa: F821
        run(
            freq_file=snakemake.input.train,  # noqa: F821
            config_file=config_file,
            copy_to=snakemake.output.params_cfg,  # noqa: F821
            sentinel=snakemake.output.mod,  # noqa: F821
            log=log,
        )


def main_cli():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("freq_file", help="SC_freq.txt (three stop-codon lines from etraining)")
    parser.add_argument("config_file", help="<species>_parameters.cfg to modify in place")
    parser.add_argument("--copy-to", help="also write a copy of the modified cfg here")
    args = parser.parse_args()
    run(args.freq_file, args.config_file, copy_to=args.copy_to)


if __name__ == "__main__":
    if "snakemake" in globals():
        main_snakemake()
    else:
        main_cli()
