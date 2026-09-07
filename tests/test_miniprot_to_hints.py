#!/usr/bin/env python3
"""
Unit tests for scripts/miniprot_to_hints.py

Run:  pytest tests/test_miniprot_to_hints.py
"""
import os
import sys

import pytest

SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
sys.path.insert(0, SCRIPTS_DIR)

import miniprot_to_hints as mth  # noqa: E402


def run(tmp_path, gff_text, priority=2):
    gff = tmp_path / "miniprot.gff"
    gff.write_text(gff_text)
    out = tmp_path / "protein.hints.gff"
    alignments = mth.parse_miniprot_gff(str(gff))
    mth.write_augustus_protein_hints(alignments, str(out), src="P", priority=priority)
    lines = [l.strip() for l in out.read_text().splitlines() if l.strip()]
    return alignments, lines


def feats(lines, ftype):
    return [l for l in lines if l.split("\t")[2] == ftype]


# Miniprot CDS features include the stop codon; stop_codon coincides with the last codon.
PLUS_COMPLETE = (
    "chr1\tminiprot\tmRNA\t1000\t3000\t150\t+\t.\tID=MP001;Rank=1;Identity=0.98;Target=prot_gene1 1 267\n"
    "chr1\tminiprot\tCDS\t1000\t1200\t.\t+\t0\tParent=MP001;Rank=1;Target=prot_gene1 1 67\n"
    "chr1\tminiprot\tCDS\t1800\t2200\t.\t+\t1\tParent=MP001;Rank=1;Target=prot_gene1 68 200\n"
    "chr1\tminiprot\tCDS\t2800\t3000\t.\t+\t2\tParent=MP001;Rank=1;Target=prot_gene1 201 267\n"
    "chr1\tminiprot\tstop_codon\t2998\t3000\t.\t+\t0\tParent=MP001;Rank=1\n"
)

MINUS_COMPLETE = (
    "chr2\tminiprot\tmRNA\t5000\t7000\t200\t-\t.\tID=MP002;Rank=1;Identity=0.95;Target=prot_gene2 1 367\n"
    "chr2\tminiprot\tCDS\t6500\t7000\t.\t-\t0\tParent=MP002;Rank=1;Target=prot_gene2 1 167\n"
    "chr2\tminiprot\tCDS\t5000\t5500\t.\t-\t0\tParent=MP002;Rank=1;Target=prot_gene2 168 367\n"
    "chr2\tminiprot\tstop_codon\t5000\t5002\t.\t-\t0\tParent=MP002;Rank=1\n"
)


def test_complete_alignment_plus_strand(tmp_path):
    alignments, lines = run(tmp_path, PLUS_COMPLETE)
    assert len(alignments) == 1
    a = alignments[0]
    assert (a["id"], a["target"], a["target_start"], a["target_end"], a["rank"]) == ("MP001", "prot_gene1", 1, 267, 1)
    assert a["stop_codon"] == (2998, 3000)
    # 3 CDSpart + 2 introns + start + stop
    assert len(lines) == 7
    attr = "grp=MP001;pri=2;src=P;target=prot_gene1"
    assert f"chr1\tminiprot\tCDSpart\t1000\t1200\t0\t+\t0\t{attr}" in lines
    assert f"chr1\tminiprot\tCDSpart\t1800\t2200\t0\t+\t1\t{attr}" in lines
    assert f"chr1\tminiprot\tCDSpart\t2800\t3000\t0\t+\t2\t{attr}" in lines
    assert f"chr1\tminiprot\tintron\t1201\t1799\t0\t+\t.\t{attr}" in lines
    assert f"chr1\tminiprot\tintron\t2201\t2799\t0\t+\t.\t{attr}" in lines
    term = "grp=MP001;pri=3;src=P;target=prot_gene1"
    assert f"chr1\tminiprot\tstart\t1000\t1002\t0\t+\t0\t{term}" in lines
    assert f"chr1\tminiprot\tstop\t2998\t3000\t0\t+\t0\t{term}" in lines


def test_complete_alignment_minus_strand(tmp_path):
    _, lines = run(tmp_path, MINUS_COMPLETE)
    assert len(lines) == 5  # 2 CDSpart + 1 intron + start + stop
    attr = "grp=MP002;pri=2;src=P;target=prot_gene2"
    assert f"chr2\tminiprot\tintron\t5501\t6499\t0\t-\t.\t{attr}" in lines
    term = "grp=MP002;pri=3;src=P;target=prot_gene2"
    # '-' strand: start at the highest coordinate, stop at the stop_codon feature (lowest)
    assert f"chr2\tminiprot\tstart\t6998\t7000\t0\t-\t0\t{term}" in lines
    assert f"chr2\tminiprot\tstop\t5000\t5002\t0\t-\t0\t{term}" in lines


def test_partial_secondary_alignment_gets_no_start_or_stop(tmp_path):
    """A paralog hit covering residues 80-300 of a 400 aa protein must not assert gene boundaries."""
    gff = (
        "chr3\tminiprot\tmRNA\t100\t900\t80\t+\t.\tID=MP010;Rank=2;Identity=0.61;Target=prot_gene9 80 300\n"
        "chr3\tminiprot\tCDS\t100\t400\t.\t+\t0\tParent=MP010;Rank=2;Target=prot_gene9 80 180\n"
        "chr3\tminiprot\tCDS\t600\t900\t.\t+\t0\tParent=MP010;Rank=2;Target=prot_gene9 181 300\n"
    )
    alignments, lines = run(tmp_path, gff)
    assert alignments[0]["rank"] == 2 and alignments[0]["target_start"] == 80
    assert len(feats(lines, "CDSpart")) == 2
    assert len(feats(lines, "intron")) == 1
    assert feats(lines, "start") == []
    assert feats(lines, "stop") == []


def test_n_terminal_complete_without_stop_codon_gets_start_only(tmp_path):
    gff = (
        "chr3\tminiprot\tmRNA\t100\t900\t80\t+\t.\tID=MP011;Rank=2;Identity=0.7;Target=prot_gene9 1 250\n"
        "chr3\tminiprot\tCDS\t100\t400\t.\t+\t0\tParent=MP011;Rank=2;Target=prot_gene9 1 100\n"
        "chr3\tminiprot\tCDS\t600\t900\t.\t+\t0\tParent=MP011;Rank=2;Target=prot_gene9 101 250\n"
    )
    _, lines = run(tmp_path, gff)
    assert len(feats(lines, "start")) == 1
    assert feats(lines, "stop") == []


def test_two_alignments_of_same_protein_get_distinct_groups(tmp_path):
    """Primary + paralog alignment of one protein: distinct grp=, shared target=."""
    gff = PLUS_COMPLETE + (
        "chr7\tminiprot\tmRNA\t100\t900\t80\t+\t.\tID=MP099;Rank=2;Identity=0.58;Target=prot_gene1 1 267\n"
        "chr7\tminiprot\tCDS\t100\t400\t.\t+\t0\tParent=MP099;Rank=2;Target=prot_gene1 1 100\n"
        "chr7\tminiprot\tCDS\t600\t900\t.\t+\t0\tParent=MP099;Rank=2;Target=prot_gene1 101 267\n"
        "chr7\tminiprot\tstop_codon\t898\t900\t.\t+\t0\tParent=MP099;Rank=2\n"
    )
    alignments, lines = run(tmp_path, gff)
    assert len(alignments) == 2
    groups = {l.split("\t")[8].split(";")[0] for l in lines}
    assert groups == {"grp=MP001", "grp=MP099"}
    assert all("target=prot_gene1" in l for l in lines)
    # both are complete, so both get start and stop
    assert len(feats(lines, "start")) == 2 and len(feats(lines, "stop")) == 2


def test_score_column_is_numeric_and_phase_kept(tmp_path):
    _, lines = run(tmp_path, PLUS_COMPLETE)
    for l in lines:
        p = l.split("\t")
        float(p[5])  # column 6 must be numeric, never '.'
        if p[2] == "CDSpart":
            assert p[7] in ("0", "1", "2")
        else:
            assert p[7] in (".", "0")


def test_missing_or_empty_input_returns_no_alignments(tmp_path):
    assert mth.parse_miniprot_gff(str(tmp_path / "nope.gff")) == []
    empty = tmp_path / "empty.gff"
    empty.write_text("")
    assert mth.parse_miniprot_gff(str(empty)) == []


def test_snakemake_entry_point(tmp_path, monkeypatch):
    gff = tmp_path / "mp.gff"
    gff.write_text(PLUS_COMPLETE)
    out = tmp_path / "out" / "protein.hints.gff"
    from types import SimpleNamespace
    monkeypatch.setattr(mth, "snakemake", SimpleNamespace(
        input=SimpleNamespace(miniprot_gff=str(gff)),
        output=SimpleNamespace(protein_hints=str(out)),
        params=SimpleNamespace(src="P", priority=2),
    ), raising=False)
    mth.main()
    assert out.exists() and len(out.read_text().splitlines()) == 7
