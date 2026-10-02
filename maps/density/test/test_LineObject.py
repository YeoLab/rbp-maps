import numpy as np
import pandas as pd
import pytest

from density import LineObject


def density_line(matrix, name='condition.rmats', conf=1, threshold=100):
    return LineObject.create_line(
        event_matrix=matrix, annotation_src_file=name, conf=conf, color='red',
        min_event_threshold=threshold, map_type='density',
        num_events=[matrix.shape[0]] * matrix.shape[1]
    )


def peak_line(matrix, name='condition.rmats'):
    return LineObject.create_line(
        event_matrix=matrix, annotation_src_file=name, conf=0.95, color='red',
        min_event_threshold=100, map_type='peak',
        num_events=[matrix.shape[0]] * matrix.shape[1]
    )


@pytest.fixture()
def matrix():
    """200 events, 3 positions: values are 0/2, 10/12, 20/22."""
    return pd.DataFrame({0: [0., 2.] * 100, 1: [10., 12.] * 100, 2: [20., 22.] * 100})


### density lines ###

def test_create_line_returns_density_line(matrix):
    line = density_line(matrix)
    assert isinstance(line, LineObject.DensityLine)
    assert line.has_mean() and not line.has_hist() and not line.has_pvalues()


def test_density_line_means_and_error_boundaries(matrix):
    line = density_line(matrix)
    sem = matrix[0].sem()
    assert line.means == [1, 11, 21]
    assert line.values == line.means
    assert line.error_pos == pytest.approx([1 + sem, 11 + sem, 21 + sem])
    assert line.error_neg == pytest.approx([1 - sem, 11 - sem, 21 - sem])
    assert line.max == pytest.approx(21 + sem)
    assert line.min == pytest.approx(1 - sem)


def test_density_line_with_one_event_has_zero_error():
    line = density_line(pd.DataFrame([[1., 2.]]))
    assert line.sems == [0, 0]
    assert line.std == [0, 0]
    assert line.error_pos == [1, 2]


def test_line_is_dimmed_below_the_event_threshold(matrix):
    assert not density_line(matrix, threshold=100).dim
    assert density_line(matrix, threshold=200).dim


@pytest.mark.parametrize("name,label", [
    ('/path/RBFOX2-HepG2-included-upon-knockdown', 'Included upon knockdown (200 events)'),
    ('/path/RBFOX2-HepG2-excluded-upon-knockdown', 'Excluded upon knockdown (200 events)'),
    ('/path/HepG2_native_cassette_exons_all', 'Native cassette exons all (200 events)'),
    ('/path/shorter-isoform-controls', 'Short isoform in >50% controls (200 events)'),
    ('/path/longer-isoform', 'Long isoform upon kd (200 events)'),
    ('/path/positive.se.txt', 'Positive.se (200 events)'),
])
def test_line_label(matrix, name, label):
    assert density_line(matrix, name=name).label == label


def test_line_file_label_is_the_file_basename(matrix):
    assert density_line(matrix, name='/path/a.b.rmats').file_label == 'a.b.rmats'


def test_set_std_error_boundaries_replaces_nan_with_the_mean(matrix):
    line = density_line(matrix)
    line._set_std_error_boundaries([0, np.nan, 20], [2, np.nan, 22])
    assert line.error_neg == [0, 11, 20]
    assert line.error_pos == [2, 11, 22]
    assert (line.min, line.max) == (0, 22)


### density significance ###

@pytest.fixture()
def high_and_background():
    rng = np.random.RandomState(0)
    high = pd.DataFrame(rng.normal(10, 1, (60, 2)))
    background = pd.DataFrame(rng.normal(0, 1, (60, 2)))
    return density_line(high, conf=0.95), background


def test_mannwhitneyu_reports_minus_log10_p(high_and_background):
    line, background = high_and_background
    line.calculate_and_set_significance(background, 'mannwhitneyu')
    assert line.has_pvalues()
    assert len(line.p_values) == 2
    assert all(p > 10 for p in line.p_values)


def test_mannwhitneyu_is_one_sided(high_and_background):
    """A line below its background is not significant."""
    line, background = high_and_background
    low = density_line(background, conf=0.95)
    low.calculate_and_set_significance(line.event_matrix, 'mannwhitneyu')
    assert all(p < 0.01 for p in low.p_values)


def test_ks_reports_minus_log10_p(high_and_background):
    line, background = high_and_background
    line.calculate_and_set_significance(background, 'ks')
    assert all(p > 10 for p in line.p_values)


def test_zscore_against_background(matrix):
    line = density_line(matrix)
    background = pd.DataFrame({0: [0., 2.] * 50, 1: [0., 2.] * 50, 2: [0., 2.] * 50})
    line.calculate_and_set_significance(background, 'zscore')
    std = background[0].std()
    assert line.p_values == pytest.approx([0, 10 / std, 20 / std])


def test_unknown_test_leaves_line_without_pvalues(matrix):
    line = density_line(matrix)
    line.calculate_and_set_significance(matrix, 'fisher')
    assert not line.has_pvalues()


### peak lines ###

@pytest.fixture()
def peaks():
    """4 events, 3 positions; peaks overlap 4, 2 and 0 events."""
    return pd.DataFrame([[1, 1, 0], [1, 1, 0], [1, 0, 0], [1, 0, 0]])


def test_create_line_returns_peak_line(peaks):
    line = peak_line(peaks)
    assert isinstance(line, LineObject.PeakLine)
    assert line.has_hist() and not line.has_mean()
    assert line.conf == 1


def test_peak_line_reports_fraction_of_events_with_a_peak(peaks):
    line = peak_line(peaks)
    assert line.hist == [4, 2, 0]
    assert line.values == [1, 0.5, 0]
    assert line.means == line.values


def test_peak_line_error_is_binomial(peaks):
    line = peak_line(peaks)
    assert line.error_pos == pytest.approx([1, 0.75, 0])
    assert line.error_neg == pytest.approx([1, 0.25, 0])


def test_peak_line_fisher_against_background(peaks):
    background = pd.DataFrame(np.zeros((40, 3), dtype=int))
    line = peak_line(peaks)
    line.calculate_and_set_significance(background, 'fisher')
    assert line.p_values[0] > line.p_values[1] > 1
    assert line.p_values[2] == 0
