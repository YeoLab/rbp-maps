import importlib
import sys
import types


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
