import importlib
import sys
import types
import argparse
import subprocess


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


def test_ensure_density_bigwigs_calls_makebigwigfiles_when_missing(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)

    existing_files = {"ip.bam"}
    generation_calls = []

    def fake_isfile(path):
        return path in existing_files

    def fake_run_makebigwigfiles(
            bam, pos_bw, neg_bw, genome_file, direction, makebigwigfiles_cmd,
            extra_args, workdir=None
    ):
        generation_calls.append(
            {
                "bam": bam,
                "pos_bw": pos_bw,
                "neg_bw": neg_bw,
                "genome_file": genome_file,
                "direction": direction,
                "makebigwigfiles_cmd": makebigwigfiles_cmd,
                "extra_args": extra_args,
                "workdir": workdir,
            }
        )
        existing_files.update({pos_bw, neg_bw})

    monkeypatch.setattr(plot_map.os.path, "isfile", fake_isfile)
    monkeypatch.setattr(plot_map, "run_makebigwigfiles", fake_run_makebigwigfiles)

    plot_map.ensure_density_bigwigs(
        bam="ip.bam",
        pos_bw="ip.norm.pos.bw",
        neg_bw="ip.norm.neg.bw",
        genome_file="hg19.chrom.sizes",
        direction="r",
        makebigwigfiles_cmd="makebigwigfiles",
        makebigwigfiles_extra_args="--threads 8",
    )

    assert len(generation_calls) == 1
    assert generation_calls[0]["bam"] == "ip.bam"
    assert generation_calls[0]["pos_bw"] == "ip.norm.pos.bw"
    assert generation_calls[0]["neg_bw"] == "ip.norm.neg.bw"
    assert generation_calls[0]["genome_file"] == "hg19.chrom.sizes"
    assert generation_calls[0]["direction"] == "r"
    assert generation_calls[0]["makebigwigfiles_cmd"] == "makebigwigfiles"
    assert generation_calls[0]["extra_args"] == "--threads 8"


def test_ensure_density_bigwigs_forwards_makebigwigfiles_workdir(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    existing_files = {"ip.bam"}
    generation_calls = []

    monkeypatch.setattr(plot_map.os.path, "isfile", lambda path: path in existing_files)

    def fake_run_makebigwigfiles(
            bam, pos_bw, neg_bw, genome_file, direction, makebigwigfiles_cmd,
            extra_args, workdir=None
    ):
        generation_calls.append(workdir)
        existing_files.update({pos_bw, neg_bw})

    monkeypatch.setattr(plot_map, "run_makebigwigfiles", fake_run_makebigwigfiles)

    plot_map.ensure_density_bigwigs(
        bam="ip.bam",
        pos_bw="ip.norm.pos.bw",
        neg_bw="ip.norm.neg.bw",
        genome_file="hg19.chrom.sizes",
        direction="r",
        makebigwigfiles_cmd="makebigwigfiles",
        makebigwigfiles_extra_args="",
        makebigwigfiles_workdir="/tmp/signal-work",
    )

    assert generation_calls == ["/tmp/signal-work"]


def test_ensure_density_bigwigs_skips_generation_when_outputs_exist(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    existing_files = {"ip.norm.pos.bw", "ip.norm.neg.bw"}
    generation_calls = []

    monkeypatch.setattr(plot_map.os.path, "isfile", lambda path: path in existing_files)
    monkeypatch.setattr(
        plot_map,
        "run_makebigwigfiles",
        lambda *a, **k: generation_calls.append((a, k)),
    )

    plot_map.ensure_density_bigwigs(
        bam="ip.bam",
        pos_bw="ip.norm.pos.bw",
        neg_bw="ip.norm.neg.bw",
        genome_file="hg19.chrom.sizes",
        direction=None,
        makebigwigfiles_cmd=None,
        makebigwigfiles_extra_args="",
    )

    assert generation_calls == []


def test_resolve_density_input_paths_uses_generated_signal_dir(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)

    class Args(object):
        ipbam = "/readonly/ip.bam"
        inputbam = "/readonly/input.bam"
        ip_pos_bw = None
        ip_neg_bw = None
        input_pos_bw = None
        input_neg_bw = None
        generated_signal_dir = "/tmp/generated-signal"

    (
        ip_bam, input_bam, ip_pos_bw, ip_neg_bw, input_pos_bw, input_neg_bw
    ) = plot_map.resolve_density_input_paths(Args())

    assert ip_bam == "/readonly/ip.bam"
    assert input_bam == "/readonly/input.bam"
    assert ip_pos_bw == "/tmp/generated-signal/ip.norm.pos.bw"
    assert ip_neg_bw == "/tmp/generated-signal/ip.norm.neg.bw"
    assert input_pos_bw == "/tmp/generated-signal/input.norm.pos.bw"
    assert input_neg_bw == "/tmp/generated-signal/input.norm.neg.bw"


def test_maybe_subset_rmats_annotations_subsets_only_rmats(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    calls = []

    def fake_subset(annotation_path, event, output_dir=None, force=False, runner=None):
        calls.append((annotation_path, event, output_dir, force))
        return annotation_path + ".nr.txt"

    monkeypatch.setattr(plot_map, "subset_rmats_annotation_file", fake_subset)

    annotations = ["a.rmats.txt", "b.miso.txt", "c.rmats.txt"]
    annotation_types = ["rmats", "miso", "rmats"]
    updated = plot_map.maybe_subset_rmats_annotations(
        annotations=annotations,
        annotation_types=annotation_types,
        event="se",
        auto_subset_rmats=True,
        subset_rmats_dir="/tmp/subset",
        subset_rmats_force=True
    )

    assert updated == ["a.rmats.txt.nr.txt", "b.miso.txt", "c.rmats.txt.nr.txt"]
    assert calls == [
        ("a.rmats.txt", "se", "/tmp/subset", True),
        ("c.rmats.txt", "se", "/tmp/subset", True),
    ]


def test_maybe_subset_rmats_annotations_rejects_unsupported_event(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    try:
        plot_map.maybe_subset_rmats_annotations(
            annotations=["a.rmats.txt"],
            annotation_types=["rmats"],
            event="bed",
            auto_subset_rmats=True
        )
        assert False, "Expected ValueError for unsupported event"
    except ValueError as exc:
        assert "only supported for events" in str(exc)


def test_resolve_density_input_paths_defaults_to_bam_directory(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)

    class Args(object):
        ipbam = "/a/ip.bam"
        inputbam = "/b/input.bam"
        ip_pos_bw = None
        ip_neg_bw = None
        input_pos_bw = None
        input_neg_bw = None
        generated_signal_dir = None

    (_, _, ip_pos_bw, ip_neg_bw, input_pos_bw, input_neg_bw) = \
        plot_map.resolve_density_input_paths(Args())
    assert ip_pos_bw == "/a/ip.norm.pos.bw"
    assert ip_neg_bw == "/a/ip.norm.neg.bw"
    assert input_pos_bw == "/b/input.norm.pos.bw"
    assert input_neg_bw == "/b/input.norm.neg.bw"


def test_resolve_density_input_paths_prefers_explicit_bigwig_paths(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)

    class Args(object):
        ipbam = "/a/ip.bam"
        inputbam = "/b/input.bam"
        ip_pos_bw = "/x/ip.pos.bw"
        ip_neg_bw = "/x/ip.neg.bw"
        input_pos_bw = "/y/in.pos.bw"
        input_neg_bw = "/y/in.neg.bw"
        generated_signal_dir = "/tmp/ignored"

    (_, _, ip_pos_bw, ip_neg_bw, input_pos_bw, input_neg_bw) = \
        plot_map.resolve_density_input_paths(Args())
    assert ip_pos_bw == "/x/ip.pos.bw"
    assert ip_neg_bw == "/x/ip.neg.bw"
    assert input_pos_bw == "/y/in.pos.bw"
    assert input_neg_bw == "/y/in.neg.bw"


def test_resolve_density_input_paths_requires_density_bams(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)

    class Args(object):
        ipbam = None
        inputbam = "/b/input.bam"
        ip_pos_bw = None
        ip_neg_bw = None
        input_pos_bw = None
        input_neg_bw = None
        generated_signal_dir = None

    try:
        plot_map.resolve_density_input_paths(Args())
        assert False, "Expected ValueError for missing BAMs"
    except ValueError as exc:
        assert "--ip/--ipbam" in str(exc)


def test_ensure_density_bigwigs_raises_if_outputs_still_missing(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    monkeypatch.setattr(plot_map.os.path, "isfile", lambda path: False)
    monkeypatch.setattr(plot_map, "run_makebigwigfiles", lambda *a, **k: None)

    try:
        plot_map.ensure_density_bigwigs(
            bam="ip.bam",
            pos_bw="ip.norm.pos.bw",
            neg_bw="ip.norm.neg.bw",
            genome_file="hg19.chrom.sizes",
            direction=None,
            makebigwigfiles_cmd=None,
            makebigwigfiles_extra_args="",
        )
        assert False, "Expected RuntimeError when outputs remain missing"
    except RuntimeError as exc:
        assert "still missing" in str(exc)


def test_subset_rmats_annotation_file_reuses_existing_file(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    runner_calls = []

    monkeypatch.setattr(plot_map.os.path, "exists", lambda p: p.endswith(".nr.txt"))
    monkeypatch.setattr(plot_map.os, "makedirs", lambda *a, **k: None)

    output = plot_map.subset_rmats_annotation_file(
        annotation_path="/d/a.txt",
        event="se",
        output_dir="/out",
        force=False,
        runner=lambda *a, **k: runner_calls.append((a, k))
    )
    assert output == "/out/a.txt.nr.txt"
    assert runner_calls == []


def test_subset_rmats_annotation_file_calls_runner_when_needed(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    runner_calls = []

    monkeypatch.setattr(plot_map.os.path, "exists", lambda p: False)
    monkeypatch.setattr(plot_map.os, "makedirs", lambda *a, **k: None)

    output = plot_map.subset_rmats_annotation_file(
        annotation_path="/d/a.txt",
        event="a3ss",
        output_dir="/out",
        force=True,
        runner=lambda *a: runner_calls.append(a)
    )
    assert output == "/out/a.txt.nr.txt"
    assert runner_calls == [("/d/a.txt", "/out/a.txt.nr.txt", "a3ss", "rmats")]


def test_maybe_subset_rmats_annotations_noop_when_disabled(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    annotations = ["a", "b"]
    out = plot_map.maybe_subset_rmats_annotations(
        annotations=annotations,
        annotation_types=["rmats", "miso"],
        event="se",
        auto_subset_rmats=False
    )
    assert out == annotations


def test_maybe_subset_rmats_annotations_is_case_insensitive(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    monkeypatch.setattr(
        plot_map,
        "subset_rmats_annotation_file",
        lambda annotation_path, **kwargs: annotation_path + ".nr.txt",
    )
    out = plot_map.maybe_subset_rmats_annotations(
        annotations=["a", "b"],
        annotation_types=["RMATS", "MiSo"],
        event="se",
        auto_subset_rmats=True
    )
    assert out == ["a.nr.txt", "b"]


def test_run_makebigwigfiles_requires_genome(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    try:
        plot_map.run_makebigwigfiles("a.bam", "a.pos.bw", "a.neg.bw", None)
        assert False, "Expected ValueError for missing genome file"
    except ValueError as exc:
        assert "--genome" in str(exc)


def test_run_makebigwigfiles_requires_tools_on_path(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    monkeypatch.setattr(plot_map.shutil, "which", lambda cmd: None)
    try:
        plot_map.run_makebigwigfiles(
            bam="a.bam",
            pos_bw="a.pos.bw",
            neg_bw="a.neg.bw",
            genome_file="g.sizes"
        )
        assert False, "Expected RuntimeError for missing external tool"
    except RuntimeError as exc:
        assert "not found in PATH" in str(exc)


def test_run_makebigwigfiles_rejects_invalid_direction(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    monkeypatch.setattr(plot_map.shutil, "which", lambda cmd: "/usr/bin/" + cmd)
    monkeypatch.setattr(plot_map.subprocess, "check_output", lambda *a, **k: b"100")
    monkeypatch.setattr(plot_map.os, "makedirs", lambda *a, **k: None)
    try:
        plot_map.run_makebigwigfiles(
            bam="a.bam",
            pos_bw="/tmp/a.pos.bw",
            neg_bw="/tmp/a.neg.bw",
            genome_file="g.sizes",
            direction="bad"
        )
        assert False, "Expected ValueError for invalid direction"
    except ValueError as exc:
        assert "Invalid direction" in str(exc)


def test_run_makebigwigfiles_requires_mapped_reads(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    monkeypatch.setattr(plot_map.shutil, "which", lambda cmd: "/usr/bin/" + cmd)
    monkeypatch.setattr(plot_map.subprocess, "check_output", lambda *a, **k: b"0")
    monkeypatch.setattr(plot_map.os, "makedirs", lambda *a, **k: None)
    try:
        plot_map.run_makebigwigfiles(
            bam="a.bam",
            pos_bw="/tmp/a.pos.bw",
            neg_bw="/tmp/a.neg.bw",
            genome_file="g.sizes",
            direction="r"
        )
        assert False, "Expected RuntimeError for zero mapped reads"
    except RuntimeError as exc:
        assert "no mapped reads" in str(exc)


def test_run_makebigwigfiles_success_executes_pipeline(monkeypatch, tmp_path):
    plot_map = _load_plot_map(monkeypatch)
    monkeypatch.setattr(plot_map.shutil, "which", lambda cmd: "/usr/bin/" + cmd)

    popen_calls = []
    check_call_calls = []

    class _FakeStdout(object):
        def close(self):
            return None

    class _FakeProc(object):
        def __init__(self, has_stdout=False):
            self.stdout = _FakeStdout() if has_stdout else None

        def wait(self):
            return 0

    def fake_popen(cmd, stdout=None, stdin=None):
        popen_calls.append(cmd)
        if cmd[0] == "bedtools":
            return _FakeProc(has_stdout=True)
        return _FakeProc(has_stdout=False)

    monkeypatch.setattr(plot_map.subprocess, "Popen", fake_popen)
    monkeypatch.setattr(plot_map.subprocess, "check_output", lambda *a, **k: b"100")
    monkeypatch.setattr(
        plot_map.subprocess,
        "check_call",
        lambda cmd: check_call_calls.append(cmd)
    )

    pos_bw = str(tmp_path / "ip.norm.pos.bw")
    neg_bw = str(tmp_path / "ip.norm.neg.bw")
    workdir = str(tmp_path / "work")
    plot_map.run_makebigwigfiles(
        bam="/data/ip.bam",
        pos_bw=pos_bw,
        neg_bw=neg_bw,
        genome_file="/ref/genome.sizes",
        direction="f",
        workdir=workdir
    )

    assert any(cmd[0] == "bedtools" for cmd in popen_calls)
    assert any(cmd[0] == "sort" for cmd in popen_calls)
    assert check_call_calls == [
        ["bedGraphToBigWig", workdir + "/ip.norm.pos.bg", "/ref/genome.sizes", pos_bw],
        ["bedGraphToBigWig", workdir + "/ip.norm.neg.bg", "/ref/genome.sizes", neg_bw],
    ]


def test_run_make_density_invalid_event_raises_systemexit(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    monkeypatch.setattr(plot_map.density.ReadDensity, "ReadDensity", lambda **k: object())
    try:
        plot_map.run_make_density(
            outfile="o.svg",
            ip_pos_bw="a", ip_neg_bw="b", ip_bam="c",
            input_pos_bw="d", input_neg_bw="e", input_bam="f",
            norm_func=lambda *a, **k: None,
            event="not-an-event",
            exon_or_upstream_offset=50,
            intron_or_downstream_offset=300,
            confidence=0.95,
            annotation_dict={},
            condition_list=[],
            bg_filename=None,
            test_method="permutation",
            scale=False,
        )
        assert False, "Expected SystemExit for invalid event"
    except SystemExit:
        assert True


def test_run_make_peak_calls_significance_when_conditions_and_bg(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    method_calls = []

    class _FakeMap(object):
        def create_matrices(self):
            method_calls.append("create_matrices")

        def normalize_matrix(self):
            method_calls.append("normalize_matrix")

        def create_lines(self):
            method_calls.append("create_lines")

        def set_background_and_calculate_significance(self, c, b, t):
            method_calls.append(("sig", c, b, t))

        def write_intermediates_to_csv(self):
            method_calls.append("write_csv")

        def plot(self, conditions):
            method_calls.append(("plot", conditions))

    monkeypatch.setattr(plot_map.density.Peak, "Peak", lambda peaks: ("peak", peaks))
    monkeypatch.setattr(
        plot_map.Map, "SkippedExon", lambda *a, **k: _FakeMap(), raising=False
    )

    plot_map.run_make_peak(
        outfile="o.svg",
        peak_file="p.bed",
        norm_func=lambda *a, **k: None,
        event="se",
        exon_or_upstream_offset=50,
        intron_or_downstream_offset=300,
        confidence=0.95,
        annotation_dict={"a": "rmats"},
        condition_list=["a"],
        bg_filename="bg",
        test_method="fisher",
        scale=False,
    )
    assert ("sig", ["a"], "bg", "fisher") in method_calls
    assert ("plot", ["a"]) in method_calls


def test_main_requires_density_bams_when_peak_not_set(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    monkeypatch.setattr(
        sys, "argv",
        [
            "plot_map",
            "--output", "o.svg",
            "--annotations", "a",
            "--annotation_type", "rmats",
        ]
    )
    try:
        plot_map.main()
        assert False, "Expected SystemExit parser error for missing --ip/--input"
    except SystemExit as exc:
        assert exc.code == 2


def test_main_peak_mode_uses_subsetted_annotations_for_bg_and_tests(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)
    captured = {}

    monkeypatch.setattr(
        plot_map,
        "maybe_subset_rmats_annotations",
        lambda **k: ["a.nr.txt", "b.nr.txt", "c.nr.txt"]
    )

    def fake_run_make_peak(
            outfile, peak_file, norm_func, event, exon_offset, intron_offset,
            confidence, annotation_dict, files_to_test, background_file,
            test_method, scale
    ):
        captured["annotation_dict"] = annotation_dict
        captured["files_to_test"] = files_to_test
        captured["background_file"] = background_file
        captured["peak_file"] = peak_file

    monkeypatch.setattr(plot_map, "run_make_peak", fake_run_make_peak)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_map",
            "--peak", "p.bed",
            "--output", "o.svg",
            "--event", "se",
            "--annotations", "a.txt", "b.txt", "c.txt",
            "--annotation_type", "rmats", "rmats", "rmats",
            "--testnums", "0", "2",
            "--bgnum", "1",
            "--auto_subset_rmats",
        ],
    )
    plot_map.main()

    assert captured["peak_file"] == "p.bed"
    assert list(captured["annotation_dict"].keys()) == ["a.nr.txt", "b.nr.txt", "c.nr.txt"]
    assert captured["files_to_test"] == ["a.nr.txt", "c.nr.txt"]
    assert captured["background_file"] == "b.nr.txt"


def test_main_density_mode_requires_genome_when_bigwigs_missing(monkeypatch):
    plot_map = _load_plot_map(monkeypatch)

    monkeypatch.setattr(plot_map, "check_for_index", lambda bam: None)
    monkeypatch.setattr(plot_map.os.path, "isfile", lambda p: False)
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "plot_map",
            "--ip", "ip.bam",
            "--input", "in.bam",
            "--output", "o.svg",
            "--annotations", "a.txt",
            "--annotation_type", "rmats",
        ],
    )
    try:
        plot_map.main()
        assert False, "Expected parser error when genome is required"
    except SystemExit as exc:
        assert exc.code == 2
