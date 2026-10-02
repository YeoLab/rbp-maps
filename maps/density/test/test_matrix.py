import os

import numpy as np
import pandas as pd
import pytest

from maps.conftest import rmats_line
from density import matrix

curdir = os.path.dirname(__file__)

EXON, INTRON = 5, 10
REGION = EXON + INTRON


def rng(start, end, step=1):
    return [float(v) for v in range(start, end, step)]


def features(name):
    return os.path.join(
        curdir, 'test_features/RBFOX2-BGHLV26-HepG2-{}.MATS.JunctionCountOnly.txt'.format(name)
    )


### same_length_region (bed) ###

def test_same_length_region_positive_strand(fake_density, write_lines):
    bed = write_lines('sites.bed', ['chr1\t100\t101\tsite\t0\t+'])
    df = matrix.same_length_region(bed, fake_density, 'bed', 3, 2, False)
    assert list(df.index) == ['chr1:100-101:site:+']
    assert list(df.iloc[0]) == rng(97, 103)


def test_same_length_region_negative_strand_is_5p_to_3p(fake_density, write_lines):
    """Upstream of a (-) site is at higher coordinates."""
    bed = write_lines('sites.bed', ['chr1\t100\t101\tsite\t0\t-'])
    df = matrix.same_length_region(bed, fake_density, 'bed', 3, 2, False)
    assert list(df.iloc[0]) == rng(103, 97, -1)


def test_same_length_region_skips_header_lines(fake_density, write_lines):
    bed = write_lines('sites.bed', ['ID\tsomething', 'chr1\t100\t101\tsite\t0\t+'])
    assert matrix.same_length_region(bed, fake_density, 'bed', 1, 1, False).shape == (1, 3)


def test_same_length_region_scale_gives_100_columns(fake_density, write_lines):
    bed = write_lines('sites.bed', ['chr1\t100\t350\ta\t0\t+', 'chr1\t400\t450\tb\t0\t+'])
    df = matrix.same_length_region(bed, fake_density, 'bed', 0, 0, True)
    assert df.shape == (2, 100)
    assert df.loc['chr1:400-450:b:+', 0] == 400
    assert df.loc['chr1:400-450:b:+', 99] == 449


def test_same_length_region_scales_when_lengths_differ(fake_density, write_lines):
    bed = write_lines('sites.bed', ['chr1\t100\t350\ta\t0\t+', 'chr1\t400\t450\tb\t0\t+'])
    assert matrix.same_length_region(bed, fake_density, 'bed', 0, 0, False).shape == (2, 100)


### skipped_exon ###

def test_skipped_exon_positive_strand_regions(fake_density, write_lines):
    """[exon|intron] 3' site, [intron|exon] 5' site, for each junction 5'->3'."""
    event = write_lines('se.txt', [rmats_line('chr1', '+', [400, 500, 100, 200, 700, 800])])
    row = list(matrix.skipped_exon(event, fake_density, EXON, INTRON).iloc[0])
    assert row == rng(195, 210) + rng(390, 405) + rng(495, 510) + rng(690, 705)


def test_skipped_exon_negative_strand_regions(fake_density, write_lines):
    event = write_lines('se.txt', [rmats_line('chr1', '-', [400, 500, 100, 200, 700, 800])])
    row = list(matrix.skipped_exon(event, fake_density, EXON, INTRON).iloc[0])
    assert row == rng(704, 689, -1) + rng(509, 494, -1) + rng(404, 389, -1) + rng(209, 194, -1)


def test_skipped_exon_pads_with_minus_one_at_the_next_exon(fake_density, write_lines):
    """A 4-nt intron and a 3-nt upstream exon: positions beyond them are -1."""
    event = write_lines('se.txt', [rmats_line('chr1', '+', [204, 500, 197, 200, 700, 800])])
    row = list(matrix.skipped_exon(event, fake_density, EXON, INTRON).iloc[0])
    assert row[:REGION] == [-1, -1] + rng(197, 204) + [-1] * 6
    assert row[REGION:2 * REGION] == [-1] * 6 + rng(200, 209)


def test_skipped_exon_real_rmats_file_shape(fake_density):
    n_events = sum(1 for line in open(features('SE')) if not line.startswith('ID'))
    df = matrix.skipped_exon(features('SE'), fake_density, 50, 300)
    assert df.shape == (n_events, 4 * 350)


### other splice events ###

@pytest.mark.parametrize("name,function,regions", [
    ('A3SS', matrix.alt_3p_splice_site, 3),
    ('A5SS', matrix.alt_5p_splice_site, 3),
    ('RI', matrix.retained_intron, 2),
    ('MXE', matrix.mutually_exc_exon, 6),
])
def test_splice_event_real_rmats_file_shape(fake_density, name, function, regions):
    n_events = sum(1 for line in open(features(name)) if not line.startswith('ID'))
    df = function(features(name), fake_density, 50, 300)
    assert df.shape == (n_events, regions * 350)
    assert list(df.columns) == list(range(regions * 350))


def test_retained_intron_regions(fake_density, write_lines):
    event = write_lines('ri.txt', [rmats_line('chr1', '+', [100, 800, 100, 200, 700, 800])])
    row = list(matrix.retained_intron(event, fake_density, EXON, INTRON).iloc[0])
    assert row == rng(195, 210) + rng(690, 705)


def test_alt_5p_splice_site_regions(fake_density, write_lines):
    """short (alt2) 3' site, long (alt1) 3' site, downstream 5' site."""
    event = write_lines('a5ss.txt', [rmats_line('chr1', '+', [100, 250, 100, 200, 700, 800])])
    row = list(matrix.alt_5p_splice_site(event, fake_density, EXON, INTRON).iloc[0])
    assert row == rng(195, 210) + rng(245, 260) + rng(690, 705)


def test_alt_3p_splice_site_regions(fake_density, write_lines):
    """upstream 3' site, long (alt1) 5' site, short (alt2) 5' site."""
    event = write_lines('a3ss.txt', [rmats_line('chr1', '+', [650, 800, 700, 800, 100, 200])])
    row = list(matrix.alt_3p_splice_site(event, fake_density, EXON, INTRON).iloc[0])
    assert row == rng(195, 210) + rng(640, 655) + rng(690, 705)


def test_mutually_exc_exon_regions(fake_density, write_lines):
    event = write_lines('mxe.txt', [
        rmats_line('chr1', '+', [400, 500, 700, 800, 100, 200, 1000, 1100])
    ])
    row = list(matrix.mutually_exc_exon(event, fake_density, EXON, INTRON).iloc[0])
    assert row == (rng(195, 210) + rng(390, 405) + rng(495, 510)
                   + rng(690, 705) + rng(795, 810) + rng(990, 1005))


### multi_length_regions ###

def test_multi_length_regions_stops_at_the_feature_midpoint(fake_density, write_lines):
    """A 10-nt feature: only 5 nt on each side of the midpoint are reported, the rest is -1."""
    bed = write_lines('features.bed', ['chr1\t100\t110\tf\t0\t+'])
    row = list(matrix.multi_length_regions(bed, fake_density, 'bed', 8, 3).iloc[0])
    five_prime, three_prime = row[:11], row[11:]
    assert five_prime == rng(97, 105) + [-1] * 3
    assert three_prime == [-1] * 3 + rng(105, 113)


### meta ###

def test_meta_concatenates_exons_and_scales(fake_density, write_lines):
    bed = write_lines('cds.bed', [
        'chr1\t100\t150\ttxA\t0\t+',
        'chr1\t200\t250\ttxA\t0\t+',
        'chr1\t500\t600\ttxB\t0\t-',
    ])
    df = matrix.meta(bed, fake_density, 0, 0, scale_to=100)
    assert df.shape == (2, 100)
    assert list(df.loc['txA']) == rng(100, 150) + rng(200, 250)
    assert list(df.loc['txB']) == rng(599, 499, -1)


def test_meta_uses_every_exon_of_a_transcript(fake_density, write_lines):
    bed = write_lines('cds.bed', [
        'chr1\t{}\t{}\ttx\t0\t+'.format(100 * i, 100 * i + 10) for i in range(1, 11)
    ])
    df = matrix.meta(bed, fake_density, 0, 0, scale_to=100)
    assert df.loc['tx', 99] == 1009


def test_meta_negative_strand_exons_are_read_3p_exon_last(fake_density, write_lines):
    bed = write_lines('cds.bed', ['chr1\t100\t150\ttx\t0\t-', 'chr1\t200\t250\ttx\t0\t-'])
    df = matrix.meta(bed, fake_density, 0, 0, scale_to=100)
    assert list(df.loc['tx']) == rng(249, 199, -1) + rng(149, 99, -1)
