import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import modify_SC_freq as msc  # noqa: E402

TEMPLATE = """# comment line
/Constant/almost_identical_a          1.0    # unrelated
/Constant/amberprob                   0.33   # Prob(stop codon = tag), if 0 tag is assumed to code for amino acid
/Constant/ochreprob                   0.33   # Prob(stop codon = taa), if 0 taa is assumed to code for amino acid
/Constant/opalprob                    0.34   # Prob(stop codon = tga), if 0 tga is assumed to code for amino acid
/Constant/decomp_num_at               1
"""

ETRAIN_TAIL = """tag:   313 (0.207)
taa:   521 (0.345)
tga:   676 (0.448)
"""


@pytest.fixture
def files(tmp_path):
    cfg = tmp_path / "sp_parameters.cfg"
    cfg.write_text(TEMPLATE)
    freq = tmp_path / "SC_freq.txt"
    freq.write_text(ETRAIN_TAIL)
    return cfg, freq


def values(cfg_text):
    out = {}
    for line in cfg_text.splitlines():
        for key in msc.PARAM_TO_CODON:
            if line.startswith(key):
                out[key] = line.split()[1]
    return out


def test_parse_freqs(files):
    _, freq = files
    assert msc.parse_stop_codon_freqs(str(freq)) == {"tag": "0.207", "taa": "0.345", "tga": "0.448"}


def test_parse_freqs_missing_codon_raises(tmp_path):
    f = tmp_path / "bad.txt"
    f.write_text("tag: 1 (0.5)\ntaa: 1 (0.5)\n")
    with pytest.raises(ValueError, match="tga"):
        msc.parse_stop_codon_freqs(str(f))


def test_first_run_sets_values_and_keeps_comments(files):
    cfg, freq = files
    msc.run(str(freq), str(cfg), log=open(os.devnull, "w"))
    text = cfg.read_text()
    assert values(text) == {"/Constant/amberprob": "0.207", "/Constant/ochreprob": "0.345", "/Constant/opalprob": "0.448"}
    assert "# Prob(stop codon = tag)" in text
    assert "/Constant/almost_identical_a          1.0    # unrelated" in text  # untouched
    assert "/Constant/decomp_num_at               1" in text


def test_rerun_with_new_frequencies_overwrites_previous_values(files):
    """The old implementation replaced the literal '0.33'/'0.34' and was a silent no-op on rerun."""
    cfg, freq = files
    msc.run(str(freq), str(cfg), log=open(os.devnull, "w"))
    freq.write_text("tag:   10 (0.100)\ntaa:   20 (0.200)\ntga:   70 (0.700)\n")
    msc.run(str(freq), str(cfg), log=open(os.devnull, "w"))
    assert values(cfg.read_text()) == {"/Constant/amberprob": "0.100", "/Constant/ochreprob": "0.200", "/Constant/opalprob": "0.700"}


def test_non_default_template_values_are_replaced(files):
    cfg, freq = files
    cfg.write_text(TEMPLATE.replace("0.33", "0.3333").replace("0.34", "0.3334"))
    msc.run(str(freq), str(cfg), log=open(os.devnull, "w"))
    assert values(cfg.read_text())["/Constant/amberprob"] == "0.207"


def test_missing_parameter_line_raises_and_leaves_file_intact(files):
    cfg, freq = files
    cfg.write_text(TEMPLATE.replace("/Constant/opalprob", "/Constant/renamed"))
    before = cfg.read_text()
    with pytest.raises(ValueError, match="found 2"):
        msc.run(str(freq), str(cfg), log=open(os.devnull, "w"))
    assert cfg.read_text() == before
    assert not os.path.exists(str(cfg) + ".tmp")


def test_copy_and_sentinel_are_written(files, tmp_path):
    cfg, freq = files
    copy = tmp_path / "out" / "sp_parameters.cfg"
    sentinel = tmp_path / "out" / "SC_freq_mod.done"
    msc.run(str(freq), str(cfg), copy_to=str(copy), sentinel=str(sentinel), log=open(os.devnull, "w"))
    assert copy.read_text() == cfg.read_text()
    assert "tag=0.207" in sentinel.read_text()
