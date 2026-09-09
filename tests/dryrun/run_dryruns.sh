#!/usr/bin/env bash
# DAG regression test: dry-run the workflow for every training.mode x prediction.mode
# combination using the tiny fixtures in tests/dryrun/. Exits non-zero on the first failure.
#
# Usage (from anywhere):  bash tests/dryrun/run_dryruns.sh [-v]
#   -v  print Snakemake's job table for each combination
set -uo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$ROOT"

if command -v snakemake >/dev/null 2>&1; then
    SMK=(snakemake)
elif command -v conda >/dev/null 2>&1 && conda env list | grep -q '^snakemake '; then
    SMK=(conda run -n snakemake snakemake)
else
    echo "snakemake not found on PATH nor in a conda env named 'snakemake'" >&2
    exit 2
fi

verbose=0; [[ "${1:-}" == "-v" ]] && verbose=1
status=0
for cfg in tests/dryrun/configs/*.yaml; do
    name=$(basename "$cfg" .yaml)
    out=$("${SMK[@]}" -n -s snakefile --configfile "$cfg" 2>&1); rc=$?
    if [[ $rc -eq 0 ]]; then
        printf "%-20s OK\n" "$name"
        [[ $verbose -eq 1 ]] && echo "$out" | sed -n '/^Job stats/,/^total/p'
    else
        printf "%-20s FAIL\n" "$name"
        echo "$out" | tail -30
        status=1
    fi
done
exit $status
