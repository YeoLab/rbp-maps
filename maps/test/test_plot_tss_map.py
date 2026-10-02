"""End-to-end stranded TSS maps on a 2 kb synthetic chromosome."""

import os
import sys

import pandas as pd
import pytest

from maps import plot_tss_map

SLOP = 50


@pytest.fixture()
def inputs(make_bigwig, make_bam, write_lines):
    """IP: 4 RPM on (+) at 200-210 and 2 RPM on (-) at 180-190. Input: empty apart from one far base.
    TSS: (+) at 200 and 1000, (-) at 240 (its window overlaps the window at 200)."""
    reads = [('chr1', 10 * i, 20, False) for i in range(4)]
    return dict(
        ip=make_bam('ip.bam', reads),
        input=make_bam('input.bam', reads),
        ip_pos_bw=make_bigwig('ip.norm.pos.bw', [('chr1', 200, 210, 4)]),
        ip_neg_bw=make_bigwig('ip.norm.neg.bw', [('chr1', 180, 190, 2)]),
        input_pos_bw=make_bigwig('input.norm.pos.bw', [('chr1', 1900, 1901, 1)]),
        input_neg_bw=make_bigwig('input.norm.neg.bw', [('chr1', 1900, 1901, 1)]),
        tss=write_lines('tss.bed', [
            'chr1\t200\t201\ta\t10\t+', 'chr1\t1000\t1001\tb\t2\t+', 'chr1\t240\t241\tc\t5\t-',
        ]),
    )


def run(monkeypatch, tmp_path, inputs, *extra):
    outfile = str(tmp_path / 'tss.png')
    monkeypatch.setattr(sys, 'argv', [
        'plot_tss_map', '--ip', inputs['ip'], '--input', inputs['input'],
        '--tss', inputs['tss'], '--output', outfile, '--slop', str(SLOP),
    ] + list(extra))
    plot_tss_map.main()
    base = os.path.splitext(outfile)[0]
    return (
        pd.read_csv(base + '.windows.tsv', sep='\t'),
        pd.read_csv(base + '.profiles.csv', index_col=0),
        pd.read_csv(base + '.summary.csv'),
    )


def test_tss_map_writes_figure_and_tables(monkeypatch, tmp_path, inputs):
    windows, profiles, summary = run(monkeypatch, tmp_path, inputs, '--shift', '0')
    assert os.path.getsize(str(tmp_path / 'tss.png')) > 0
    assert list(windows.x) == [200, 1000, 240]
    assert list(windows.overlaps_opposite) == [True, False, True]
    assert windows.kept.all()
    assert list(profiles.columns) == ['sense', 'antisense', 'total']
    assert list(profiles.index) == list(range(-SLOP, SLOP + 1))
    assert list(summary.regions.unique()) == ['TSS']


def test_tss_map_separates_sense_from_antisense(monkeypatch, tmp_path, inputs):
    """Only the (+) window at 200 is kept: (+) reads at 0..+9 are sense, (-) reads at -20..-11 antisense."""
    _, profiles, _ = run(monkeypatch, tmp_path, inputs, '--shift', '0', '--min_height', '10')
    assert (profiles.sense.loc[0:9] > 0).all() and profiles.sense.drop(range(0, 10)).eq(0).all()
    assert (profiles.antisense.loc[-20:-11] > 0).all() and profiles.antisense.drop(range(-20, -10)).eq(0).all()
    assert profiles.sense.loc[0] == pytest.approx(2 * profiles.antisense.loc[-20])
    assert (profiles.total == profiles.sense + profiles.antisense).all()


def test_tss_map_reads_a_window_on_the_minus_strand_5p_to_3p(monkeypatch, tmp_path, inputs):
    """For the (-) TSS at 240, the (+) reads at 200-210 are antisense, 31 to 40 nt downstream."""
    windows, profiles, _ = run(monkeypatch, tmp_path, inputs, '--shift', '0', '--min_height', '5')
    assert list(windows.kept) == [True, False, True]
    antisense_only = profiles.antisense.loc[31:40]
    assert (antisense_only > 0).all()


def test_tss_map_filters(monkeypatch, tmp_path, inputs):
    windows, _, _ = run(monkeypatch, tmp_path, inputs, '--shift', '0', '--exclude_opposite_overlaps')
    assert list(windows.kept) == [False, True, False]


def test_tss_map_shifted_control(monkeypatch, tmp_path, inputs):
    """200 nt downstream of the windows there is no signal."""
    _, profiles, summary = run(monkeypatch, tmp_path, inputs, '--shift', '200')
    assert not profiles[['shifted_sense', 'shifted_antisense', 'shifted_total']].any().any()
    assert list(summary.regions.unique()) == ['TSS', 'shifted 200 nt']


def test_tss_map_fails_when_no_window_is_left(monkeypatch, tmp_path, inputs):
    with pytest.raises(ValueError):
        run(monkeypatch, tmp_path, inputs, '--min_files', '2')
