import numpy as np
import pandas as pd
import pytest

from density import ReadDensity, tss
from density import normalization_functions as norm

SIZES = {'chr1': 2000}


def calls(*rows):
    """rows of (x, strand, h) on chr1, each from its own file unless a file is given."""
    return pd.DataFrame([
        ('chr1', r[0], r[1], float(r[2]), float(r[2]), r[3] if len(r) > 3 else 'file{}'.format(i))
        for i, r in enumerate(rows)
    ], columns=['chrom', 'x', 'strand', 'h', 'c', 'file'])


### read_tss_calls ###

def test_read_tss_calls_takes_the_leftmost_highest_base_of_a_cage_peak(write_lines):
    peaks = write_lines('cage.bed', [
        'chr1\t100\t104\tpk1\t1000\t+\t13\tg\tg\tpk1\t1.0,5.0,5.0,2.0',
        'chr1\t300\t302\tpk2\t1000\t-\t4\tg\tg\tpk2\t1.0,3.0',
    ])
    df = tss.read_tss_calls([peaks])
    assert list(df.x) == [101, 301]
    assert list(df.strand) == ['+', '-']
    assert list(df.h) == [5, 3]
    assert list(df.c) == [3.25, 2]
    assert set(df.file) == {'cage.bed'}


def test_read_tss_calls_uses_the_5p_end_of_bed6_intervals(write_lines):
    bed = write_lines('tss.bed', ['chr1\t100\t150\ta\t7\t+', 'chr1\t100\t150\tb\t9\t-'])
    df = tss.read_tss_calls([bed])
    assert list(df.x) == [100, 149]
    assert list(df.h) == [7, 9]


def test_read_tss_calls_pools_files(write_lines):
    one = write_lines('one.bed', ['chr1\t100\t101\ta\t1\t+'])
    two = write_lines('two.bed', ['chr1\t500\t501\tb\t1\t+'])
    assert list(tss.read_tss_calls([one, two]).file) == ['one.bed', 'two.bed']


### build_windows ###

def test_build_windows_are_centered_on_the_tss():
    windows = tss.build_windows(calls((100, '+', 1)), 10, SIZES)
    assert (windows.start[0], windows.end[0]) == (90, 111)
    assert windows.region_id[0] == 'chr1:90-111:+'


def test_build_windows_keeps_the_highest_call_of_overlapping_same_strand_windows():
    windows = tss.build_windows(calls((100, '+', 5), (105, '+', 9), (500, '+', 1)), 10, SIZES)
    assert list(windows.x) == [105, 500]
    assert list(windows.n_files) == [2, 1]


def test_build_windows_breaks_height_ties_by_mean_signal():
    df = calls((100, '+', 5), (105, '+', 5))
    df.loc[0, 'c'] = 4
    df.loc[1, 'c'] = 2
    assert list(tss.build_windows(df, 10, SIZES).x) == [100]


def test_build_windows_clusters_windows_that_only_touch():
    """[90, 111) and [111, 132) are book-ended, which bedtools merge joins."""
    assert len(tss.build_windows(calls((100, '+', 5), (121, '+', 9)), 10, SIZES)) == 1
    assert len(tss.build_windows(calls((100, '+', 5), (122, '+', 9)), 10, SIZES)) == 2


def test_build_windows_counts_distinct_files_per_cluster():
    windows = tss.build_windows(calls((100, '+', 5, 'a.bed'), (105, '+', 9, 'a.bed')), 10, SIZES)
    assert list(windows.n_files) == [1]


def test_build_windows_does_not_merge_across_strands_but_flags_the_overlap():
    windows = tss.build_windows(calls((100, '+', 5), (110, '-', 9), (500, '-', 1)), 10, SIZES)
    assert list(windows.x) == [100, 110, 500]
    assert list(windows.overlaps_opposite) == [True, True, False]


@pytest.mark.parametrize('distance,overlap', [(20, True), (21, False)])
def test_build_windows_opposite_overlap_needs_one_shared_base(distance, overlap):
    windows = tss.build_windows(calls((100, '+', 1), (100 + distance, '-', 1)), 10, SIZES)
    assert list(windows.overlaps_opposite) == [overlap, overlap]


def test_build_windows_drops_windows_off_the_chromosome_and_unknown_chromosomes():
    df = calls((5, '+', 1), (1995, '+', 1), (100, '+', 1))
    df.loc[len(df)] = ('chrUn', 100, '+', 1.0, 1.0, 'file')
    assert list(tss.build_windows(df, 10, SIZES).x) == [100]


### shift_windows ###

def test_shift_windows_moves_downstream_in_the_direction_of_transcription():
    windows = tss.build_windows(calls((100, '+', 1), (500, '-', 1)), 10, SIZES)
    shifted = tss.shift_windows(windows, 50, SIZES)
    assert list(shifted.start) == [140, 440]
    assert list(shifted.end) == [161, 461]
    assert list(shifted.x) == [150, 450]
    assert list(shifted.region_id) == ['chr1:90-111:+_shift50', 'chr1:490-511:-_shift50']


def test_shift_windows_drops_windows_that_leave_the_chromosome():
    windows = tss.build_windows(calls((100, '-', 1), (1900, '+', 1), (1000, '+', 1)), 10, SIZES)
    assert list(tss.shift_windows(windows, 95, SIZES).x) == [1095]


### stranded_matrices ###

@pytest.fixture()
def density(make_bigwig):
    """(+) track: 3 at chr1:95-100. (-) track: -7 at chr1:102-104 (legacy negative values)."""
    return ReadDensity.ReadDensity(
        pos=make_bigwig('pos.bw', [('chr1', 95, 100, 3)]),
        neg=make_bigwig('neg.bw', [('chr1', 102, 104, -7)]),
    )


def test_stranded_matrices_positive_window(density):
    windows = tss.build_windows(calls((100, '+', 1)), 5, SIZES)
    sense, antisense = tss.stranded_matrices(density, windows)
    assert list(sense[0]) == [3] * 5 + [0] * 6
    assert list(antisense[0]) == [0] * 7 + [7, 7] + [0] * 2


def test_stranded_matrices_negative_window_is_reversed(density):
    """Position 0 of a (-) window is its highest coordinate."""
    windows = tss.build_windows(calls((100, '-', 1)), 5, SIZES)
    sense, antisense = tss.stranded_matrices(density, windows)
    assert list(sense[0]) == [0] * 2 + [7, 7] + [0] * 7
    assert list(antisense[0]) == [0] * 6 + [3] * 5


def test_stranded_matrices_sum_to_the_unstranded_coverage(density):
    windows = tss.build_windows(calls((100, '+', 1), (100, '-', 1)), 5, SIZES)
    sense, antisense = tss.stranded_matrices(density, windows)
    plus, minus = sense + antisense
    assert list(plus) == list(minus[::-1])


### normalize, profiles, summarize ###

def test_normalize_divides_both_orientations_by_one_denominator():
    """|1| + |-2| + 2 positions x 0.5 = 4."""
    ip = (np.array([[2., 0.]]), np.array([[0., 1.]]))
    inp = (np.array([[1., 0.]]), np.array([[0., 3.]]))
    normalized = tss.normalize(ip, inp, 0.5)
    assert list(normalized['sense'][0]) == [0.25, 0]
    assert list(normalized['antisense'][0]) == [0, -0.5]
    assert list(normalized['total'][0]) == [0.25, -0.5]


def test_normalize_total_is_the_unstranded_normalization():
    rng = np.random.RandomState(0)
    ip = (rng.rand(5, 11), rng.rand(5, 11))
    inp = (rng.rand(5, 11), rng.rand(5, 11))
    normalized = tss.normalize(ip, inp, 0.5)
    unstranded = norm.per_region_subtract_and_normalize(
        pd.DataFrame(ip[0] + ip[1]), pd.DataFrame(inp[0] + inp[1]), 0.5, 0.5
    )
    assert np.allclose(normalized['total'], unstranded.values)
    assert np.allclose(normalized['sense'] + normalized['antisense'], normalized['total'])


def test_normalize_window_without_signal_is_zero():
    zeros = (np.zeros((1, 3)), np.zeros((1, 3)))
    assert not tss.normalize(zeros, zeros, 0.5)['total'].any()


def test_profiles_are_indexed_by_position_relative_to_the_tss():
    normalized = {o: np.array([[1., 2., 3.], [3., 4., 5.]]) for o in tss.ORIENTATIONS}
    profile = tss.profiles(normalized)
    assert list(profile.index) == [-1, 0, 1]
    assert list(profile.sense) == [2, 3, 4]
    assert list(tss.profiles(normalized, keep=np.array([False, True])).sense) == [3, 4, 5]


def profile_of(sense, antisense, slop=300):
    index = pd.Index(np.arange(-slop, slop + 1), name='position')
    return pd.DataFrame({'sense': sense, 'antisense': antisense, 'total': sense + antisense}, index=index)


def test_summarize_sums_windows_centered_on_the_tss():
    summary = tss.summarize(profile_of(3., 1.))
    assert list(summary.window) == ['+/-50', '+/-100', '+/-250', '+/-300']
    assert list(summary.sense) == [303, 603, 1503, 1803]
    assert list(summary.antisense) == [101, 201, 501, 601]
    assert list(summary.sense_frac) == [0.75] * 4


def test_summarize_skips_windows_wider_than_the_profile():
    assert list(tss.summarize(profile_of(1., 1., slop=100)).window) == ['+/-50', '+/-100']


def test_summarize_gives_no_fraction_when_a_sum_is_not_positive():
    assert tss.summarize(profile_of(1., -1.)).sense_frac.isna().all()
