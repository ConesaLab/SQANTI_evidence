import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import noise_metrics as nm  # noqa: E402

HDR = "isoform\texons\tFL\tcoding\tperc_A_downstream_TTS\tall_canonical\tCDS_length\tpredicted_NMD\tpsauron_score\n"
# 1 clean coding multi-exon, 1 non-coding mono-exon FL=2 failing 4 rules, 1 coding failing only intra-priming,
# 1 non-coding with NA everywhere (numeric NA fails, string NA passes)
ROWS = (
    "a\t3\t25\tcoding\t10\tcanonical\t300\tFALSE\t0.95\n"
    "b\t1\t2\tnon_coding\t80\tnon_canonical\tNA\tFALSE\tNA\n"
    "c\t4\t12\tcoding\t75\tcanonical\t250\tFALSE\t0.9\n"
    "d\t2\t1\tnon_coding\tNA\tNA\tNA\tNA\tNA\n"
)


def rows():
    import csv, io
    return list(csv.DictReader(io.StringIO(HDR + ROWS), delimiter="\t"))


def test_metrics_fractions():
    m = nm.compute_metrics(rows())
    assert m["models"] == 4
    assert m["noncoding_frac"] == 50.0
    assert m["fl2_frac"] == 50.0          # b (FL 2) and d (FL 1)
    assert m["monoexon_frac"] == 25.0     # b
    # b fails intrapriming, canonical, coding, CDS_length, psauron (5); d fails coding, CDS_length, psauron,
    # intrapriming (NA numeric) = 4; a 0; c 1  -> 2 of 4 fail >= 3
    assert m["multifail_frac"] == 50.0


def test_strict_rules_na_semantics():
    r = {"perc_A_downstream_TTS": "NA", "all_canonical": "NA", "coding": "coding", "CDS_length": "NA",
         "predicted_NMD": "NA", "psauron_score": "NA"}
    assert nm.STRICT_RULES["canonical"](r) and nm.STRICT_RULES["not_NMD"](r)
    assert not nm.STRICT_RULES["intrapriming"](r) and not nm.STRICT_RULES["CDS_length"](r) and not nm.STRICT_RULES["psauron"](r)


def test_verdict_thresholds():
    assert nm.verdict(2.6) == "clean" and nm.verdict(19.9) == "clean"
    assert nm.verdict(25.3) == "moderate"
    assert nm.verdict(60.8) == "noisy"


def test_run_writes_table_and_log(tmp_path):
    c = tmp_path / "sp_classification.txt"
    c.write_text(HDR + ROWS)
    out, log = tmp_path / "noise_metrics.tsv", tmp_path / "noise.log"
    m = nm.run(str(c), str(out), "sp", str(log))
    lines = out.read_text().splitlines()
    assert lines[0].split("\t") == ["sample", "models", "noncoding_frac", "fl2_frac", "monoexon_frac", "multifail_frac", "verdict"]
    assert lines[1].split("\t")[:3] == ["sp", "4", "50.0"] and lines[1].endswith("noisy")
    text = log.read_text()
    assert "verdict: noisy" in text and "raising the 'FL' value" in text
    assert m["models"] == 4


def test_empty_classification(tmp_path):
    c = tmp_path / "empty.txt"
    c.write_text(HDR)
    m = nm.run(str(c), str(tmp_path / "o.tsv"), "e")
    assert m["models"] == 0 and nm.verdict(m["noncoding_frac"]) == "clean"
