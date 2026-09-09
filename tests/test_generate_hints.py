import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import generate_hints as gh  # noqa: E402

DEFAULT_CFG = {"exon": False, "exonpart": False, "intron": True, "CDS": True, "start": True, "stop": True}
ATTR = "grp={t};pri=1;src=lrRNA"


def types(lines):
    return [l.split("\t")[2] for l in lines]


def test_complete_orf_plus_strand():
    # exons 100-200, 300-400, 500-600 ; CDS 150-200, 300-400, 500-550 (stop codon inside CDS at 548-550)
    lines = gh.transcript_hints("t1", "chr1", "+", [(100, 200), (300, 400), (500, 600)],
                                [(150, 200), (300, 400), (500, 550)], {"t1": "complete"}, DEFAULT_CFG)
    a = ATTR.format(t="t1")
    assert f"chr1\tHints\tintron\t201\t299\t0\t+\t.\t{a}\n" in lines
    assert f"chr1\tHints\tintron\t401\t499\t0\t+\t.\t{a}\n" in lines
    assert types(lines).count("CDS") == 3 and "CDSpart" not in types(lines)
    assert f"chr1\tHints\tstart\t150\t152\t0\t+\t0\t{a}\n" in lines
    assert f"chr1\tHints\tstop\t548\t550\t0\t+\t0\t{a}\n" in lines
    assert "exon" not in types(lines) and "exonpart" not in types(lines)


def test_complete_orf_minus_strand():
    lines = gh.transcript_hints("t2", "chr1", "-", [(100, 200), (300, 400)], [(150, 200), (300, 350)],
                                {"t2": "complete"}, DEFAULT_CFG)
    a = ATTR.format(t="t2")
    # '-' strand: start at the highest CDS coordinate, stop at the lowest
    assert f"chr1\tHints\tstart\t348\t350\t0\t-\t0\t{a}\n" in lines
    assert f"chr1\tHints\tstop\t150\t152\t0\t-\t0\t{a}\n" in lines


@pytest.mark.parametrize("cds_type,strand,expect_start,expect_stop,cdspart_index", [
    ("5prime_partial", "+", False, True, 0),    # missing start: first CDS block on + is partial
    ("3prime_partial", "+", True, False, 2),    # missing stop: last CDS block on + is partial
    ("5prime_partial", "-", False, True, 2),    # on '-', the start is the last block
    ("3prime_partial", "-", True, False, 0),    # on '-', the stop is the first block
    ("internal", "+", False, False, None),
])
def test_partial_orfs(cds_type, strand, expect_start, expect_stop, cdspart_index):
    cds = [(150, 200), (300, 400), (500, 550)]
    lines = gh.transcript_hints("t", "c", strand, [(100, 200), (300, 400), (500, 600)], cds, {"t": cds_type}, DEFAULT_CFG)
    tps = types(lines)
    assert ("start" in tps) == expect_start
    assert ("stop" in tps) == expect_stop
    cds_lines = [l for l in lines if l.split("\t")[2] in ("CDS", "CDSpart")]
    if cdspart_index is None:
        assert [l.split("\t")[2] for l in cds_lines] == ["CDSpart", "CDS", "CDSpart"]
    else:
        assert cds_lines[cdspart_index].split("\t")[2] == "CDSpart"
        assert sum(1 for l in cds_lines if l.split("\t")[2] == "CDSpart") == 1


def test_introns_restricted_to_cds_span_when_exon_hints_off():
    # UTR intron 201-299 lies outside the CDS span (300-550) and must be dropped; 401-499 kept
    lines = gh.transcript_hints("t", "c", "+", [(100, 200), (300, 400), (500, 600)], [(300, 400), (500, 550)],
                                {"t": "complete"}, DEFAULT_CFG)
    introns = [l for l in lines if "\tintron\t" in l]
    assert len(introns) == 1 and "\t401\t499\t" in introns[0]


def test_all_introns_when_exon_hints_on():
    cfg = dict(DEFAULT_CFG, exon=True)
    lines = gh.transcript_hints("t", "c", "+", [(100, 200), (300, 400), (500, 600)], [(300, 400), (500, 550)],
                                {"t": "complete"}, cfg)
    assert sum(1 for l in lines if "\tintron\t" in l) == 2
    assert types(lines).count("exon") == 3


def test_exonpart_marks_terminal_exons():
    cfg = dict(DEFAULT_CFG, exon=True, exonpart=True)
    lines = gh.transcript_hints("t", "c", "+", [(100, 200), (300, 400), (500, 600)], [], {}, cfg)
    assert [l.split("\t")[2] for l in lines if l.split("\t")[2] in ("exon", "exonpart")] == ["exonpart", "exon", "exonpart"]


def test_score_column_numeric_and_unknown_transcript_defaults_to_internal():
    lines = gh.transcript_hints("nope", "c", "+", [(1, 10), (20, 30)], [(1, 10), (20, 30)], {}, DEFAULT_CFG)
    for l in lines:
        float(l.split("\t")[5])
    assert "start" not in types(lines) and "stop" not in types(lines)  # 'internal' by default


def test_ungrouped_gtf_is_collected_per_transcript(tmp_path):
    """Lines of one transcript interleaved with another must still form a single hint group."""
    gtf = tmp_path / "x.gtf"
    gtf.write_text(
        'c\tS\texon\t100\t200\t.\t+\t.\ttranscript_id "A"; gene_id "gA";\n'
        'c\tS\texon\t1000\t1100\t.\t-\t.\ttranscript_id "B"; gene_id "gB";\n'
        'c\tS\texon\t300\t400\t.\t+\t.\ttranscript_id "A"; gene_id "gA";\n'
        'c\tS\tCDS\t150\t200\t.\t+\t0\ttranscript_id "A"; gene_id "gA";\n'
        'c\tS\tCDS\t300\t350\t.\t+\t2\ttranscript_id "A"; gene_id "gA";\n'
        'c\tS\texon\t1200\t1300\t.\t-\t.\ttranscript_id "B"; gene_id "gB";\n'
    )
    tx = gh.collect_transcripts(str(gtf))
    assert list(tx) == ["A", "B"]                       # first-appearance order preserved
    assert tx["A"]["exons"] == [(100, 200), (300, 400)] and tx["A"]["cds"] == [(150, 200), (300, 350)]
    assert tx["B"]["exons"] == [(1000, 1100), (1200, 1300)]


def test_hint_config_parsing(tmp_path):
    cfg = tmp_path / "hint_config.tsv"
    cfg.write_text("feature\tenabled\nexon\tfalse\nintron\tTrue\nCDS\ttrue\n\nstart\ttrue\nstop\tfalse\n")
    parsed = gh.load_hint_config(str(cfg))
    assert parsed == {"exon": False, "exonpart": False, "intron": True, "CDS": True, "start": True, "stop": False}


@pytest.mark.parametrize("bad,msg", [
    ("feature\tenabled\nexonparts\ttrue\n", "unknown hint feature"),
    ("feature\tenabled\nintron\tyes\n", "must be true or false"),
    ("feature\tenabled\nintron\n", "expected 'feature<TAB>true|false'"),
])
def test_hint_config_errors_name_the_line(tmp_path, bad, msg):
    cfg = tmp_path / "hint_config.tsv"
    cfg.write_text(bad)
    with pytest.raises(ValueError, match=msg):
        gh.load_hint_config(str(cfg))


def test_end_to_end_matches_repo_default_config(tmp_path):
    repo_cfg = os.path.join(os.path.dirname(__file__), "..", "envs", "hint_config.tsv")
    gtf = tmp_path / "x.gtf"
    gtf.write_text(
        'c\tS\ttranscript\t100\t600\t.\t+\t.\ttranscript_id "A"; gene_id "gA";\n'
        'c\tS\texon\t100\t200\t.\t+\t.\ttranscript_id "A"; gene_id "gA";\n'
        'c\tS\texon\t300\t600\t.\t+\t.\ttranscript_id "A"; gene_id "gA";\n'
        'c\tS\tCDS\t150\t200\t.\t+\t0\ttranscript_id "A"; gene_id "gA";\n'
        'c\tS\tCDS\t300\t450\t.\t+\t2\ttranscript_id "A"; gene_id "gA";\n'
    )
    cls = tmp_path / "cls.txt"
    cls.write_text("isoform\tlength\tCDS_type\nA\t500\tcomplete\n")
    out = tmp_path / "hints.gff"
    n = gh.generate_hints(str(gtf), str(cls), repo_cfg, str(out))
    lines = out.read_text().splitlines()
    assert n == len(lines) == 5  # intron, 2 CDS, start, stop with the default config
    assert {l.split("\t")[2] for l in lines} == {"intron", "CDS", "start", "stop"}
