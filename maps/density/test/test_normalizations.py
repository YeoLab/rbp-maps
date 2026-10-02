import numpy as np
import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from density import normalization_functions as norm


def frame(rows, index=None):
    return pd.DataFrame(rows, index=index, dtype=float)


### clean / add_missing_events ###

def test_clean_turns_nan_into_zero_and_minus_one_into_nan():
    cleaned = norm.clean(frame([[1, np.nan, -1]]))
    assert cleaned.iloc[0, 0] == 1
    assert cleaned.iloc[0, 1] == 0
    assert np.isnan(cleaned.iloc[0, 2])


def test_add_missing_events_returns_input_when_nothing_is_missing():
    ip = frame([[1, 2]], index=['a'])
    inp = frame([[3, 4]], index=['a'])
    assert norm.add_missing_events(ip, inp) is inp


def test_add_missing_events_appends_nan_rows():
    ip = frame([[1, 2], [5, 6]], index=['a', 'b'])
    inp = frame([[3, 4]], index=['a'])
    added = norm.add_missing_events(ip, inp)
    assert list(added.index) == ['a', 'b']
    assert added.loc['b'].isna().all()


### pdf ###

def test_calculate_pdf_rows_sum_to_one():
    pdf = norm.calculate_pdf(frame([[1, 2, 3, 4], [0, 0, 0, 0]]), pseudocount=1)
    assert list(pdf.iloc[0]) == pytest.approx([2 / 14., 3 / 14., 4 / 14., 5 / 14.])
    assert list(pdf.iloc[1]) == pytest.approx([0.25] * 4)


def test_calculate_pdf_ignores_nan_positions():
    pdf = norm.calculate_pdf(frame([[1, np.nan]]), pseudocount=1)
    assert pdf.iloc[0, 0] == 1
    assert np.isnan(pdf.iloc[0, 1])


def test_calculate_pdf_without_pseudocount_uses_smallest_positive_value():
    pdf = norm.calculate_pdf(frame([[0, 2], [4, 4]]))
    assert list(pdf.iloc[0]) == pytest.approx([2 / 6., 4 / 6.])


def test_calculate_pdf_refuses_negative_values():
    assert norm.calculate_pdf(frame([[1, -2]]), pseudocount=1) == 1


def test_calculate_abs_pdf_docstring_example():
    df = frame([[0, 1, 2, 3], [3, 4, 5, -1], [0, 0, 0, 0]])
    pdf = norm.calculate_abs_pdf(df, pseudocount=1)
    assert list(pdf.iloc[0]) == pytest.approx([0, 0.1, 0.2, 0.3])
    assert list(pdf.iloc[1]) == pytest.approx([3 / 17., 4 / 17., 5 / 17., -1 / 17.])
    assert list(pdf.iloc[2]) == [0, 0, 0, 0]


def test_get_abs_sum_counts_pseudocount_only_at_non_nan_positions():
    assert norm.get_abs_sum(pd.Series([1, -2, np.nan]), 1) == 5


def test_get_abs_sum_of_an_empty_row_is_one():
    assert norm.get_abs_sum(pd.Series([np.nan, np.nan]), 1) == 1


### ip vs input ###

def test_per_region_subtract_and_normalize():
    """(ip - input) / (sum|ip - input| + pseudocount * positions)."""
    normed = norm.per_region_subtract_and_normalize(
        frame([[1, 2, 3, 4]]), frame([[1, 1, 1, 1]]), 1, 1
    )
    assert list(normed.iloc[0]) == pytest.approx([0, 0.1, 0.2, 0.3])


def test_per_region_subtract_and_normalize_keeps_negative_signal():
    normed = norm.per_region_subtract_and_normalize(
        frame([[0, 0]]), frame([[1, 3]]), 1, 1
    )
    assert list(normed.iloc[0]) == pytest.approx([-1 / 6., -3 / 6.])


def test_per_region_subtract_and_normalize_masks_padding():
    """-1 marks positions beyond the feature: not counted, reported as NaN."""
    normed = norm.per_region_subtract_and_normalize(
        frame([[2, -1]]), frame([[1, -1]]), 1, 1
    )
    assert normed.iloc[0, 0] == pytest.approx(0.5)
    assert np.isnan(normed.iloc[0, 1])


def test_per_region_subtract_and_normalize_treats_missing_input_events_as_zero():
    normed = norm.per_region_subtract_and_normalize(
        frame([[1, 1]], index=['only_in_ip']), frame([[5, 5]], index=['other']), 1, 1
    )
    assert list(normed.loc['only_in_ip']) == pytest.approx([0.25, 0.25])


def test_normalize_and_per_region_subtract():
    normed = norm.normalize_and_per_region_subtract(
        frame([[1, 2, 3, 4]]), frame([[1, 1, 1, 1]]), 1, 1
    )
    assert list(normed.iloc[0]) == pytest.approx(
        [2 / 14. - 0.25, 3 / 14. - 0.25, 4 / 14. - 0.25, 5 / 14. - 0.25]
    )


def test_read_entropy_is_zero_when_ip_equals_input():
    df = frame([[1, 2, 3]])
    entropy = norm.read_entropy(df, df, 1, 1)
    assert list(entropy.iloc[0]) == [0, 0, 0]


def test_read_entropy_value():
    """p = rpm / 1e6 + one read; entropy = p_ip * log2(p_ip / p_input)."""
    entropy = norm.read_entropy(frame([[1, 0]]), frame([[0, 0]]), 1, 1)
    assert entropy.iloc[0, 0] == pytest.approx(2e-6)
    assert entropy.iloc[0, 1] == 0


def test_read_entropy_drops_events_at_or_below_the_density_threshold():
    entropy = norm.read_entropy(
        frame([[0, 0], [1, 1]], index=['empty', 'covered']),
        frame([[0, 0], [1, 1]], index=['empty', 'covered']), 1, 1,
        min_density_threshold=0
    )
    assert entropy.loc['empty'].isna().all()
    assert list(entropy.loc['covered']) == [0, 0]


def test_get_density_and_get_input_return_the_cleaned_matrix():
    ip, inp = frame([[1, -1]]), frame([[np.nan, 2]])
    assert np.isnan(norm.get_density(ip, inp, 1, 1).iloc[0, 1])
    assert list(norm.get_input(ip, inp, 1, 1).iloc[0]) == [0, 2]


### means, errors, outlier removal ###

def test_get_means_and_sems_keeps_every_value_at_conf_1():
    means, sems, stds, merged = norm.get_means_and_sems(frame([[1, 10], [3, 30]]), conf=1)
    assert means == [2, 20]
    assert sems == pytest.approx([1, 10])
    assert stds == pytest.approx([np.sqrt(2), np.sqrt(200)])
    assert merged is None


def test_get_means_and_sems_drops_the_tails_per_position():
    """40 events at conf 0.95: one value is dropped from each tail."""
    df = frame({0: list(range(40)), 1: [0] * 39 + [1000]})
    means, _, _, _ = norm.get_means_and_sems(df, conf=0.95)
    assert means == [19.5, 0]


def test_get_means_and_sems_ignores_nan():
    means, _, _, _ = norm.get_means_and_sems(frame([[1], [np.nan], [3]]), conf=1)
    assert means == [2]


def test_get_means_and_sems_with_merged_masks_outliers():
    df = frame({0: list(range(40))})
    means, _, _, merged = norm.get_means_and_sems_with_merged(df, conf=0.95)
    assert means == [19.5]
    assert merged.shape == (40, 1)
    assert list(merged.index[merged[0].isna()]) == [0, 39]


def test_bottom_top_values_from_dataframe():
    """0.5% of 1,000 values: the 5th smallest and 5th largest."""
    df = frame({0: list(range(1000))})
    assert norm.bottom_top_values_from_dataframe(df, 0.5, 0.5) == ([4], [995])


### peaks ###

def test_divide_by_num_events():
    assert norm.divide_by_num_events([10, 20, 30, 40, 50], [10, 10, 10, 10, 5]) == [1, 2, 3, 4, 10]


def test_dev_is_binomial_standard_error():
    assert norm.dev(0.5, 0.5, 100) == pytest.approx(0.05)


def test_std_error_per_position():
    assert norm.std_error([50, 0], [100, 100]) == pytest.approx([0.05, 0])
