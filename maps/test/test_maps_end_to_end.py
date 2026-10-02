"""End-to-end maps on a 2 kb synthetic chromosome: matrices, lines, files and figure."""

import importlib
import os
import sys
from collections import OrderedDict

import pandas as pd
import pytest

from maps.conftest import rmats_line
from density import Map, Peak
from density import normalization_functions as norm

EXON, INTRON = 5, 10
REGION = EXON + INTRON
PEAKS = os.path.join(os.path.dirname(__file__), '..', 'density', 'test', 'test_Peak', 'a_pos_chr1_0_10.bed.bb')

# event -> (rMATS coordinate columns, number of plotted regions)
SPLICE_EVENTS = {
    'se': ([400, 500, 100, 200, 700, 800], 4),
    'ri': ([100, 800, 100, 200, 700, 800], 2),
    'a5ss': ([100, 250, 100, 200, 700, 800], 3),
    'a3ss': ([650, 800, 700, 800, 100, 200], 3),
    'mxe': ([400, 500, 700, 800, 100, 200, 1000, 1100], 6),
}


@pytest.fixture()
def plot_map():
    sys.modules.pop('maps.plot_map', None)
    return importlib.import_module('maps.plot_map')


@pytest.fixture()
def bigwigs(make_bigwig, make_bam):
    """Paths for an IP at 2 RPM and an input at 1 RPM across chr1."""
    reads = [('chr1', 10 * i, 20, False) for i in range(4)]
    return dict(
        ip_pos_bw=make_bigwig('ip.norm.pos.bw', [('chr1', 0, 2000, 2)]),
        ip_neg_bw=make_bigwig('ip.norm.neg.bw', [('chr1', 0, 2000, -2)]),
        ip_bam=make_bam('ip.bam', reads),
        input_pos_bw=make_bigwig('input.norm.pos.bw', [('chr1', 0, 2000, 1)]),
        input_neg_bw=make_bigwig('input.norm.neg.bw', [('chr1', 0, 2000, -1)]),
        input_bam=make_bam('input.bam', reads),
    )


def run_density(plot_map, bigwigs, outfile, event, annotations, conditions=(), background=None,
                test='mannwhitneyu', norm_func=norm.per_region_subtract_and_normalize,
                upstream=EXON, downstream=INTRON):
    plot_map.run_make_density(
        outfile=outfile, norm_func=norm_func, event=event,
        exon_or_upstream_offset=upstream, intron_or_downstream_offset=downstream,
        confidence=1, annotation_dict=annotations, condition_list=list(conditions),
        bg_filename=background, test_method=test, scale=False, **bigwigs
    )


def read_means(outfile, annotation):
    base = os.path.splitext(outfile)[0]
    return [float(line) for line in open('{}.{}.means.txt'.format(base, os.path.basename(annotation)))]


### bed maps ###

def test_bed_map_writes_figure_and_intermediates(plot_map, bigwigs, bed_annotation, tmp_path):
    bed = bed_annotation('sites.bed')
    outfile = str(tmp_path / 'map.png')
    run_density(plot_map, bigwigs, outfile, 'bed', OrderedDict([(bed, 'bed')]), upstream=3, downstream=2)

    assert os.path.getsize(outfile) > 0
    base = str(tmp_path / 'map.sites.bed')
    for clip, rpm in (('ip', 2), ('input', 1)):
        raw = pd.read_csv('{}.{}.raw_density.txt'.format(base, clip), index_col=0)
        assert raw.shape == (20, 6)
        assert (raw == rpm).all().all()
    normed = pd.read_csv(str(tmp_path / 'map.sites.bed.normed_matrix.txt'), index_col=0)
    means = read_means(outfile, bed)
    assert len(means) == 6
    assert means == pytest.approx(list(normed.mean()))
    assert means[0] > 0 and means == pytest.approx([means[0]] * 6)


def test_bed_map_normalization_levels(plot_map, bigwigs, bed_annotation, tmp_path):
    """Raw IP and raw input maps report the RPM of each track."""
    bed = bed_annotation('sites.bed')
    for func, rpm in ((norm.get_density, 2), (norm.get_input, 1)):
        outfile = str(tmp_path / '{}.png'.format(func.__name__))
        run_density(plot_map, bigwigs, outfile, 'bed', OrderedDict([(bed, 'bed')]), norm_func=func)
        assert read_means(outfile, bed) == [rpm] * (EXON + INTRON + 1)


def test_bed_map_negative_strand_sites_read_the_negative_bigwig(plot_map, bigwigs, bed_annotation, tmp_path):
    bed = bed_annotation('sites.bed', strand='-')
    outfile = str(tmp_path / 'map.png')
    run_density(plot_map, bigwigs, outfile, 'bed', OrderedDict([(bed, 'bed')]), norm_func=norm.get_density)
    assert read_means(outfile, bed) == [2] * (EXON + INTRON + 1)


@pytest.mark.parametrize('test', ['mannwhitneyu', 'ks', 'zscore'])
def test_bed_map_significance_writes_pvalues(plot_map, bigwigs, bed_annotation, tmp_path, test):
    condition, background = bed_annotation('condition.bed'), bed_annotation('background.bed', first=500)
    outfile = str(tmp_path / 'map.png')
    run_density(
        plot_map, bigwigs, outfile, 'bed',
        OrderedDict([(condition, 'bed'), (background, 'bed')]),
        conditions=[condition], background=background, test=test,
    )
    pvalues = pd.read_csv(str(tmp_path / 'map.condition.bed.pvalues.txt'), sep='\t', header=None)
    assert pvalues.shape == (EXON + INTRON + 1, 2)
    assert not os.path.exists(str(tmp_path / 'map.background.bed.pvalues.txt'))


def test_permutation_test_sets_the_background_error_band(bigwigs, bed_annotation, tmp_path):
    """The band around the background line is the envelope of resampled background means."""
    condition, background = bed_annotation('condition.bed'), bed_annotation('background.bed', first=500)
    outfile = str(tmp_path / 'map.png')
    density_map = Map.Bed(
        *density_objects(bigwigs), output_filename=outfile,
        norm_function=norm.get_density,
        annotation=OrderedDict([(condition, 'bed'), (background, 'bed')]),
        upstream_offset=3, downstream_offset=2, conf=1
    )
    density_map.create_matrices()
    density_map.normalize_matrix()
    density_map.create_lines()
    density_map.set_background_and_calculate_significance(
        [condition], background, test='permutation', num_permutations=200
    )
    condition_line, background_line = density_map.lines
    assert background_line.label.endswith('*')
    assert background_line.error_pos == [2] * 6
    assert background_line.error_neg == [2] * 6
    assert condition_line.error_pos == condition_line.means
    samples = pd.read_csv(str(tmp_path / 'map.png.background.bed.condition.bed.randsample.tsv'), sep='\t', index_col=0)
    assert samples.shape == (200, 6)


def density_objects(bigwigs):
    from density import ReadDensity
    return (
        ReadDensity.ReadDensity(pos=bigwigs['ip_pos_bw'], neg=bigwigs['ip_neg_bw'], bam=bigwigs['ip_bam']),
        ReadDensity.ReadDensity(pos=bigwigs['input_pos_bw'], neg=bigwigs['input_neg_bw'], bam=bigwigs['input_bam']),
    )


def test_multi_length_bed_map(plot_map, bigwigs, write_lines, tmp_path):
    bed = write_lines('features.bed', [
        'chr1\t{}\t{}\tf{}\t0\t+'.format(100 * i, 100 * i + 40, i) for i in range(1, 6)
    ])
    outfile = str(tmp_path / 'map.png')
    run_density(plot_map, bigwigs, outfile, 'multi-length-bed', OrderedDict([(bed, 'bed')]),
                norm_func=norm.get_density)
    assert os.path.getsize(outfile) > 0
    assert read_means(outfile, bed) == [2] * (2 * (EXON + INTRON))


### splice maps ###

@pytest.mark.parametrize('event', sorted(SPLICE_EVENTS))
def test_splice_map(plot_map, bigwigs, write_lines, tmp_path, event):
    coords, regions = SPLICE_EVENTS[event]
    annotation = write_lines(event + '.txt', [rmats_line('chr1', '+', coords, i) for i in range(3)])
    outfile = str(tmp_path / 'map.png')
    run_density(plot_map, bigwigs, outfile, event, OrderedDict([(annotation, 'rmats')]),
                norm_func=norm.get_density)
    assert os.path.getsize(outfile) > 0
    assert read_means(outfile, annotation) == [2] * (regions * REGION)


def test_atac_map_uses_the_retained_intron_layout(plot_map, bigwigs, write_lines, tmp_path):
    coords, regions = SPLICE_EVENTS['ri']
    annotation = write_lines('ri.txt', [rmats_line('chr1', '+', coords)])
    outfile = str(tmp_path / 'map.png')
    run_density(plot_map, bigwigs, outfile, 'atac', OrderedDict([(annotation, 'rmats')]),
                norm_func=norm.get_density)
    assert len(read_means(outfile, annotation)) == regions * REGION


### metagene and cds ###

def transcripts(write_lines, name, start):
    return write_lines(name, [
        'chr1\t{}\t{}\ttx{}\t0\t+'.format(start + 300 * i, start + 300 * i + 60, i) for i in range(4)
    ])


def test_metagene_map_joins_utr5_cds_and_utr3(plot_map, bigwigs, write_lines, tmp_path):
    """5'UTR, CDS and 3'UTR are scaled to 7, 52 and 41 positions."""
    annotations = OrderedDict([
        (transcripts(write_lines, 'utr5.bed', 100), 'utr5'),
        (transcripts(write_lines, 'cds.bed', 160), 'cds'),
        (transcripts(write_lines, 'utr3.bed', 220), 'utr3'),
    ])
    outfile = str(tmp_path / 'map.png')
    run_density(plot_map, bigwigs, outfile, 'metagene', annotations, norm_func=norm.get_density,
                upstream=0, downstream=0)
    assert os.path.getsize(outfile) > 0
    assert read_means(outfile, 'meta') == [2] * 100


### peak maps ###

def test_peak_map_reports_the_fraction_of_events_with_a_peak(plot_map, write_lines, tmp_path):
    """The peak covers chr1:0-10 (+); one of two sites is under it."""
    bed = write_lines('sites.bed', ['chr1\t5\t6\tunder\t0\t+', 'chr1\t500\t501\taway\t0\t+'])
    background = write_lines('background.bed', ['chr1\t700\t701\tbg{}\t0\t+'.format(i) for i in range(4)])
    outfile = str(tmp_path / 'peaks.png')
    plot_map.run_make_peak(
        outfile=outfile, peak_file=PEAKS, norm_func=norm.get_density, event='bed',
        exon_or_upstream_offset=2, intron_or_downstream_offset=2, confidence=0.95,
        annotation_dict=OrderedDict([(bed, 'bed'), (background, 'bed')]),
        condition_list=[bed], bg_filename=background, test_method='fisher', scale=False
    )
    assert os.path.getsize(outfile) > 0
    hist = [float(line) for line in open(str(tmp_path / 'peaks.sites.bed.hist.txt'))]
    assert hist == [0.5] * 5
    pvalues = pd.read_csv(str(tmp_path / 'peaks.sites.bed.pvalues.txt'), sep='\t', header=None)
    assert pvalues.shape == (5, 2)


def test_map_type_follows_the_signal_source(bigwigs):
    ip, _ = density_objects(bigwigs)
    assert Map.Map(ip, 'out.png', norm.get_density, {}).map_type == 'density'
    assert Map.Map(Peak.Peak(PEAKS), 'out.png', norm.get_density, {}).map_type == 'peak'


### command line ###

def test_main_density_bed_map(plot_map, bigwigs, bed_annotation, tmp_path, monkeypatch):
    condition, background = bed_annotation('condition.bed'), bed_annotation('background.bed', first=500)
    outfile = str(tmp_path / 'cli.svg')
    monkeypatch.setattr(sys, 'argv', [
        'plot_map', '--ip', bigwigs['ip_bam'], '--input', bigwigs['input_bam'],
        '--output', outfile, '--event', 'bed',
        '--annotations', condition, background, '--annotation_type', 'bed', 'bed',
        '--exon_offset', '3', '--intron_offset', '2', '--confidence', '1',
        '--testnums', '0', '--bgnum', '1', '--sigtest', 'ks',
    ])
    plot_map.main()
    assert os.path.getsize(outfile) > 0
    assert len(read_means(outfile, condition)) == 6
    assert os.path.exists(str(tmp_path / 'cli.condition.bed.pvalues.txt'))


def test_main_flip_swaps_the_strand_bigwigs(plot_map, bigwigs, bed_annotation, tmp_path, monkeypatch, make_bigwig):
    """With --flip a (+) site reads the *.neg.bw track."""
    make_bigwig('ip.norm.neg.bw', [('chr1', 0, 2000, 7)])
    bed = bed_annotation('sites.bed')
    outfile = str(tmp_path / 'cli.png')
    monkeypatch.setattr(sys, 'argv', [
        'plot_map', '--ip', bigwigs['ip_bam'], '--input', bigwigs['input_bam'],
        '--output', outfile, '--event', 'bed', '--annotations', bed, '--annotation_type', 'bed',
        '--normalization_level', '0', '--confidence', '1', '--flip',
    ])
    plot_map.main()
    assert set(read_means(outfile, bed)) == {7}


def test_main_peak_map(plot_map, write_lines, tmp_path, monkeypatch):
    bed = write_lines('sites.bed', ['chr1\t5\t6\tunder\t0\t+'])
    outfile = str(tmp_path / 'cli.png')
    monkeypatch.setattr(sys, 'argv', [
        'plot_map', '--peak', PEAKS, '--output', outfile, '--event', 'bed',
        '--annotations', bed, '--annotation_type', 'bed',
        '--exon_offset', '2', '--intron_offset', '2', '--normalization_level', '0',
    ])
    plot_map.main()
    assert os.path.getsize(outfile) > 0


def test_main_rejects_mismatched_annotation_types(plot_map, bigwigs, bed_annotation, tmp_path, monkeypatch):
    monkeypatch.setattr(sys, 'argv', [
        'plot_map', '--ip', bigwigs['ip_bam'], '--input', bigwigs['input_bam'],
        '--output', str(tmp_path / 'cli.png'), '--event', 'bed',
        '--annotations', bed_annotation('sites.bed'), '--annotation_type', 'bed', 'bed',
    ])
    with pytest.raises(SystemExit):
        plot_map.main()
