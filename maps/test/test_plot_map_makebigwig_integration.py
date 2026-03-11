import importlib
import hashlib
import os
import shutil
import sys
import types

import pytest


def _install_density_stubs(monkeypatch):
    density_module = types.ModuleType("density")
    density_peak_module = types.ModuleType("density.Peak")
    density_read_density_module = types.ModuleType("density.ReadDensity")
    density_norm_module = types.ModuleType("density.normalization_functions")
    density_map_module = types.ModuleType("density.Map")

    class _Dummy(object):
        pass

    density_peak_module.Peak = _Dummy
    density_read_density_module.ReadDensity = _Dummy
    density_norm_module.normalize_and_per_region_subtract = lambda *a, **k: None
    density_norm_module.get_density = lambda *a, **k: None
    density_norm_module.per_region_subtract_and_normalize = lambda *a, **k: None
    density_norm_module.read_entropy = lambda *a, **k: None
    density_norm_module.get_input = lambda *a, **k: None

    density_module.Peak = density_peak_module
    density_module.ReadDensity = density_read_density_module
    density_module.normalization_functions = density_norm_module
    density_module.Map = density_map_module

    monkeypatch.setitem(sys.modules, "density", density_module)
    monkeypatch.setitem(sys.modules, "density.Peak", density_peak_module)
    monkeypatch.setitem(sys.modules, "density.ReadDensity", density_read_density_module)
    monkeypatch.setitem(
        sys.modules, "density.normalization_functions", density_norm_module
    )
    monkeypatch.setitem(sys.modules, "density.Map", density_map_module)


def _load_plot_map(monkeypatch):
    _install_density_stubs(monkeypatch)
    sys.modules.pop("maps.plot_map", None)
    return importlib.import_module("maps.plot_map")


def _parse_bg(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            chrom, start, end, value = line.split("\t")
            rows.append((chrom, int(start), int(end), float(value)))
    return rows


def _assert_bg_equal(generated, expected):
    gen_rows = _parse_bg(generated)
    exp_rows = _parse_bg(expected)
    assert len(gen_rows) == len(exp_rows)
    for g, e in zip(gen_rows, exp_rows):
        assert g[0] == e[0]
        assert g[1] == e[1]
        assert g[2] == e[2]
        assert abs(g[3] - e[3]) < 1e-9


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def _assert_binary_equal(generated, expected):
    assert _sha256(generated) == _sha256(expected)


def _require_bigwig_toolchain():
    required = ("samtools", "bedtools", "bedGraphToBigWig")
    missing = [cmd for cmd in required if shutil.which(cmd) is None]
    if missing:
        pytest.skip("Missing required external tools: {}".format(", ".join(missing)))


@pytest.mark.integration
def test_run_makebigwigfiles_matches_testbam_references(monkeypatch, tmp_path):
    _require_bigwig_toolchain()
    plot_map = _load_plot_map(monkeypatch)

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fixtures = os.path.join(repo_root, "tests")

    bam = os.path.join(fixtures, "test.bam")
    genome = os.path.join(fixtures, "GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes")

    pos_bw = str(tmp_path / "test.pos.bw")
    neg_bw = str(tmp_path / "test.neg.bw")
    workdir = str(tmp_path / "work")

    plot_map.run_makebigwigfiles(
        bam=bam,
        pos_bw=pos_bw,
        neg_bw=neg_bw,
        genome_file=genome,
        direction="f",
        workdir=workdir
    )

    assert os.path.exists(pos_bw)
    assert os.path.exists(neg_bw)
    _assert_binary_equal(pos_bw, os.path.join(fixtures, "test.pos.bw"))
    _assert_binary_equal(neg_bw, os.path.join(fixtures, "test.neg.bw"))

    _assert_bg_equal(
        os.path.join(workdir, "test.norm.pos.bg"),
        os.path.join(fixtures, "test.sorted.norm.pos.bg")
    )
    _assert_bg_equal(
        os.path.join(workdir, "test.norm.neg.bg"),
        os.path.join(fixtures, "test.sorted.norm.neg.bg")
    )


@pytest.mark.integration
def test_run_makebigwigfiles_matches_rbfox2_references(monkeypatch, tmp_path):
    _require_bigwig_toolchain()
    plot_map = _load_plot_map(monkeypatch)

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fixtures = os.path.join(repo_root, "tests")

    bam = os.path.join(fixtures, "204_01_RBFOX2.merged.r2.bam")
    genome = os.path.join(fixtures, "hg19.chrom.sizes")

    pos_bw = str(tmp_path / "rbfox2.norm.pos.bw")
    neg_bw = str(tmp_path / "rbfox2.norm.neg.bw")
    workdir = str(tmp_path / "work")

    plot_map.run_makebigwigfiles(
        bam=bam,
        pos_bw=pos_bw,
        neg_bw=neg_bw,
        genome_file=genome,
        direction="r",
        workdir=workdir
    )

    assert os.path.exists(pos_bw)
    assert os.path.exists(neg_bw)
    _assert_binary_equal(
        pos_bw, os.path.join(fixtures, "204_01_RBFOX2.merged.r2.norm.pos.bw")
    )
    _assert_binary_equal(
        neg_bw, os.path.join(fixtures, "204_01_RBFOX2.merged.r2.norm.neg.bw")
    )

    _assert_bg_equal(
        os.path.join(workdir, "204_01_RBFOX2.merged.r2.norm.pos.bg"),
        os.path.join(fixtures, "204_01_RBFOX2.merged.r2.sorted.norm.pos.bg")
    )
    _assert_bg_equal(
        os.path.join(workdir, "204_01_RBFOX2.merged.r2.norm.neg.bg"),
        os.path.join(fixtures, "204_01_RBFOX2.merged.r2.sorted.norm.neg.bg")
    )


@pytest.mark.integration
def test_main_end_to_end_density_generation_testbam(monkeypatch, tmp_path):
    _require_bigwig_toolchain()
    plot_map = _load_plot_map(monkeypatch)

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fixtures = os.path.join(repo_root, "tests")
    test_feature = os.path.join(
        repo_root,
        "maps",
        "density",
        "test",
        "test_features",
        "RBFOX2-BGHLV26-HepG2-SE.MATS.JunctionCountOnly.txt"
    )

    ip_bam = os.path.join(fixtures, "test.bam")
    input_bam = os.path.join(fixtures, "test.bam")
    genome = os.path.join(fixtures, "GRCh38_no_alt_analysis_set_GCA_000001405.15.chrom.sizes")
    ip_pos_bw = str(tmp_path / "ip.pos.bw")
    ip_neg_bw = str(tmp_path / "ip.neg.bw")
    input_pos_bw = str(tmp_path / "input.pos.bw")
    input_neg_bw = str(tmp_path / "input.neg.bw")
    workdir = str(tmp_path / "work")
    out_plot = str(tmp_path / "out.svg")

    captured = {}
    monkeypatch.setattr(
        plot_map,
        "run_make_density",
        lambda *args: captured.setdefault("run_make_density_args", args)
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_map",
            "--ip", ip_bam,
            "--input", input_bam,
            "--ip_pos_bw", ip_pos_bw,
            "--ip_neg_bw", ip_neg_bw,
            "--input_pos_bw", input_pos_bw,
            "--input_neg_bw", input_neg_bw,
            "--genome", genome,
            "--make_bigwig_files_direction", "f",
            "--make_bigwig_files_workdir", workdir,
            "--annotations", test_feature,
            "--annotation_type", "rmats",
            "--output", out_plot,
            "--event", "se",
        ],
    )
    plot_map.main()

    assert "run_make_density_args" in captured
    _assert_binary_equal(ip_pos_bw, os.path.join(fixtures, "test.pos.bw"))
    _assert_binary_equal(ip_neg_bw, os.path.join(fixtures, "test.neg.bw"))
    _assert_binary_equal(input_pos_bw, os.path.join(fixtures, "test.pos.bw"))
    _assert_binary_equal(input_neg_bw, os.path.join(fixtures, "test.neg.bw"))


@pytest.mark.integration
def test_main_end_to_end_density_generation_rbfox2(monkeypatch, tmp_path):
    _require_bigwig_toolchain()
    plot_map = _load_plot_map(monkeypatch)

    repo_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    fixtures = os.path.join(repo_root, "tests")
    test_feature = os.path.join(
        repo_root,
        "maps",
        "density",
        "test",
        "test_features",
        "RBFOX2-BGHLV26-HepG2-SE.MATS.JunctionCountOnly.txt"
    )

    ip_bam = os.path.join(fixtures, "204_01_RBFOX2.merged.r2.bam")
    input_bam = os.path.join(fixtures, "204_01_RBFOX2.merged.r2.bam")
    genome = os.path.join(fixtures, "hg19.chrom.sizes")
    ip_pos_bw = str(tmp_path / "ip.pos.bw")
    ip_neg_bw = str(tmp_path / "ip.neg.bw")
    input_pos_bw = str(tmp_path / "input.pos.bw")
    input_neg_bw = str(tmp_path / "input.neg.bw")
    workdir = str(tmp_path / "work")
    out_plot = str(tmp_path / "out.svg")

    captured = {}
    monkeypatch.setattr(
        plot_map,
        "run_make_density",
        lambda *args: captured.setdefault("run_make_density_args", args)
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_map",
            "--ip", ip_bam,
            "--input", input_bam,
            "--ip_pos_bw", ip_pos_bw,
            "--ip_neg_bw", ip_neg_bw,
            "--input_pos_bw", input_pos_bw,
            "--input_neg_bw", input_neg_bw,
            "--genome", genome,
            "--make_bigwig_files_direction", "r",
            "--make_bigwig_files_workdir", workdir,
            "--annotations", test_feature,
            "--annotation_type", "rmats",
            "--output", out_plot,
            "--event", "se",
        ],
    )
    plot_map.main()

    assert "run_make_density_args" in captured
    _assert_binary_equal(
        ip_pos_bw, os.path.join(fixtures, "204_01_RBFOX2.merged.r2.norm.pos.bw")
    )
    _assert_binary_equal(
        ip_neg_bw, os.path.join(fixtures, "204_01_RBFOX2.merged.r2.norm.neg.bw")
    )
    _assert_binary_equal(
        input_pos_bw, os.path.join(fixtures, "204_01_RBFOX2.merged.r2.norm.pos.bw")
    )
    _assert_binary_equal(
        input_neg_bw, os.path.join(fixtures, "204_01_RBFOX2.merged.r2.norm.neg.bw")
    )
