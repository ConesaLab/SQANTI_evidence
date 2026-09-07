import copy
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import input_check  # noqa: E402


@pytest.fixture
def base_config(tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGT\n")
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r\nACGT\n+\nIIII\n")
    tools = tmp_path / "tools"
    tools.mkdir()
    return {
        "project": {"genome": str(genome), "input": str(reads), "toolsdir": str(tools),
                    "outdir": str(tmp_path / "out")},
        "training": {"mode": "mixed", "lineage": "dummy_odb12", "miniprot_threshold": 0.95,
                     "flanking_region": 1000, "test_size": 5000},
        "prediction": {"species": "sp", "mode": "split", "filter_mode": "medium"},
        "curation": {"mode": "placebo"},
        "evaluation": {"reference_gtf": ""},
    }


def expect_die(config, fragment):
    with pytest.raises(SystemExit) as exc:
        input_check.check_inputs(config)
    assert exc.value.code == 1
    return fragment


def test_valid_config_passes(base_config, monkeypatch):
    monkeypatch.delenv("AUGUSTUS_CONFIG_PATH", raising=False)
    input_check.check_inputs(base_config)


def test_ab_initio_requires_busco_only(base_config, caplog):
    cfg = copy.deepcopy(base_config)
    cfg["curation"]["mode"] = "ab_initio"
    for mode in ["mixed", "sqanti_only"]:
        cfg["training"]["mode"] = mode
        caplog.clear()
        expect_die(cfg, "busco_only")
        assert "must be trained using only BUSCO genes" in caplog.text


def test_ab_initio_with_busco_only_passes(base_config, monkeypatch):
    monkeypatch.delenv("AUGUSTUS_CONFIG_PATH", raising=False)
    cfg = copy.deepcopy(base_config)
    cfg["curation"]["mode"] = "ab_initio"
    cfg["training"]["mode"] = "busco_only"
    input_check.check_inputs(cfg)


@pytest.mark.parametrize("section,key,bad", [
    ("training", "mode", "all"),
    ("prediction", "mode", "chromosome"),
    ("prediction", "filter_mode", "monoexon"),
    ("curation", "mode", "reference"),
    ("curation", "data_type", "illumina"),
])
def test_invalid_enum_dies(base_config, section, key, bad, caplog):
    cfg = copy.deepcopy(base_config)
    cfg[section][key] = bad
    expect_die(cfg, key)
    assert f"'{section}.{key}' must be one of" in caplog.text


def test_missing_lineage_for_busco_modes(base_config):
    cfg = copy.deepcopy(base_config)
    cfg["training"]["lineage"] = ""
    expect_die(cfg, "lineage")


def test_sqanti_only_does_not_need_lineage(base_config, monkeypatch):
    monkeypatch.delenv("AUGUSTUS_CONFIG_PATH", raising=False)
    cfg = copy.deepcopy(base_config)
    cfg["training"] = {"mode": "sqanti_only"}
    input_check.check_inputs(cfg)


def test_optional_file_set_but_missing_dies(base_config, caplog):
    cfg = copy.deepcopy(base_config)
    cfg["prediction"]["hint_weights"] = "envs/does_not_exist.cfg"
    expect_die(cfg, "hint_weights")
    assert "prediction.hint_weights" in caplog.text


def test_non_numeric_values_are_reported_not_traceback(base_config, caplog):
    cfg = copy.deepcopy(base_config)
    cfg["training"]["miniprot_threshold"] = None
    expect_die(cfg, "miniprot_threshold")
    assert "must be a number" in caplog.text

    cfg = copy.deepcopy(base_config)
    cfg["training"]["test_size"] = [5000]
    expect_die(cfg, "test_size")
    assert "must be an integer" in caplog.text


def test_user_mode_requires_existing_gtf(base_config, tmp_path, monkeypatch):
    monkeypatch.delenv("AUGUSTUS_CONFIG_PATH", raising=False)
    cfg = copy.deepcopy(base_config)
    cfg["curation"]["mode"] = "user"
    expect_die(cfg, "user_gtf")
    gtf = tmp_path / "ref.gtf"
    gtf.write_text("chr1\tx\texon\t1\t10\t.\t+\t.\ttranscript_id \"t\"; gene_id \"g\";\n")
    cfg["curation"]["user_gtf"] = str(gtf)
    input_check.check_inputs(cfg)


def test_unwritable_toolsdir_dies(base_config, tmp_path, caplog):
    if os.geteuid() == 0:
        pytest.skip("root ignores directory permissions")
    cfg = copy.deepcopy(base_config)
    ro = tmp_path / "ro_tools"
    ro.mkdir()
    ro.chmod(0o555)
    cfg["project"]["toolsdir"] = str(ro)
    try:
        expect_die(cfg, "toolsdir")
        assert "not writable" in caplog.text
    finally:
        ro.chmod(0o755)


def test_unwritable_augustus_config_path_dies(base_config, tmp_path, monkeypatch, caplog):
    if os.geteuid() == 0:
        pytest.skip("root ignores directory permissions")
    ro = tmp_path / "aug_config"
    (ro / "species").mkdir(parents=True)
    (ro / "species").chmod(0o555)
    monkeypatch.setenv("AUGUSTUS_CONFIG_PATH", str(ro))
    try:
        expect_die(base_config, "AUGUSTUS_CONFIG_PATH")
        assert "AUGUSTUS_CONFIG_PATH" in caplog.text
    finally:
        (ro / "species").chmod(0o755)
