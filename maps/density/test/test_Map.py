#!/usr/env python

import os
import pytest
import pybedtools
import pandas as pd
from density import Peak
from pandas.testing import assert_series_equal

### Fixtures ###
curdir = os.path.dirname(__file__)

### These are for testing the interval/boundary regions.

@pytest.fixture()
def get_peak_pos_chr1_0_10():
    return Peak.Peak(os.path.join(curdir, 'test_Peak/a_pos_chr1_0_10.bed.bb'))

@pytest.fixture()
def get_peak_neg_chr1_0_10():
    return Peak.Peak(os.path.join(curdir, 'test_Peak/b_neg_chr1_0_10.bed.bb'))

@pytest.fixture()
def get_peak_pos_chr1_5_15():
    return Peak.Peak(os.path.join(curdir, 'test_Peak/c_pos_chr1_5_15.bed.bb'))

@pytest.fixture()
def get_peak_neg_chr1_5_15():
    return Peak.Peak(os.path.join(curdir, 'test_Peak/d_neg_chr1_5_15.bed.bb'))

@pytest.fixture()
def get_a_and_c():
    return Peak.Peak(os.path.join(curdir, 'test_Peak/a_and_c.bed.bb'))

@pytest.fixture()
def get_e_and_f():
    return Peak.Peak(os.path.join(curdir, 'test_Peak/e_and_f.bed.bb'))

### Test single overlap functions ###

def test_peak_values_p1():
    print("Tests the basic functionality of the values() function. "
          "peak and region should totally overlap")
    expect_series = pd.Series([1,1,1,1,1,1,1,1,1,1])
    test_peak = get_peak_pos_chr1_0_10()
    test_series = test_peak.values(
        'chr1', 0, 10, '+'
    )
    assert_series_equal(test_series, expect_series)
def test_peak_values_n1():
    print("Tests the basic functionality of the values() function. "
          "peak and region should totally overlap")
    expect_series = pd.Series([1,1,1,1,1,1,1,1,1,1])
    test_peak = get_peak_neg_chr1_0_10()
    test_series = test_peak.values(
        'chr1', 0, 10, '-'
    )
    assert_series_equal(test_series, expect_series)

def test_peak_values_p2():
    print("Tests the lower boundary of the region, and upper boundary "
          "of the peak (region and peak should overlap by just one)")
    expect_series = pd.Series([1,0,0,0,0])
    test_peak = get_peak_pos_chr1_0_10()
    test_series = test_peak.values(
        'chr1', 9, 14, '+'
    )
    assert_series_equal(test_series, expect_series)
def test_peak_values_n2():
    print("Tests the lower boundary of the region, and upper boundary "
          "of the peak (region and peak should overlap by just one)")
    expect_series = pd.Series([0,0,0,0,1])
    test_peak = get_peak_neg_chr1_0_10()
    test_series = test_peak.values(
        'chr1', 9, 14, '-'
    )
    assert_series_equal(test_series, expect_series)

def test_peak_values_p3():
    print("Tests the lower boundary of the region, and upper boundary "
          "of the peak (region and peak should NOT overlap)")
    expect_series = pd.Series([0,0,0,0,0])
    test_peak = get_peak_pos_chr1_0_10()
    test_series = test_peak.values(
        'chr1', 10, 15, '+'
    )
    assert_series_equal(test_series, expect_series)
def test_peak_values_n3():
    print("Tests the lower boundary of the region, and upper boundary "
          "of the peak (region and peak should NOT overlap)")
    expect_series = pd.Series([0,0,0,0,0])
    test_peak = get_peak_neg_chr1_0_10()
    test_series = test_peak.values(
        'chr1', 10, 15, '-'
    )
    assert_series_equal(test_series, expect_series)

def test_peak_values_p4():
    print("Tests the upper boundary of the region, and lower boundary "
          "of the peak (region and peak should overlap by just one)")
    expect_series = pd.Series([0,0,0,0,0,1])
    test_peak = get_peak_pos_chr1_5_15()
    test_series = test_peak.values(
        'chr1', 0, 6, '+'
    )
    assert_series_equal(test_series, expect_series)
def test_peak_values_n4():
    print("Tests the upper boundary of the region, and lower boundary "
          "of the peak (region and peak should overlap by just one)")
    expect_series = pd.Series([1,0,0,0,0,0])
    test_peak = get_peak_neg_chr1_5_15()
    test_series = test_peak.values(
        'chr1', 0, 6, '-'
    )
    assert_series_equal(test_series, expect_series)

def test_peak_values_p5():
    print("Tests the upper boundary of the region, and lower boundary "
          "of the peak (region and peak should NOT overlap)")
    expect_series = pd.Series([0,0,0,0,0])
    test_peak = get_peak_pos_chr1_5_15()
    test_series = test_peak.values(
        'chr1', 0, 5, '+'
    )
    assert_series_equal(test_series, expect_series)
def test_peak_values_n5():
    print("Tests the upper boundary of the region, and lower boundary "
          "of the peak (region and peak should NOT overlap)")
    expect_series = pd.Series([0,0,0,0,0])
    test_peak = get_peak_neg_chr1_5_15()
    test_series = test_peak.values(
        'chr1', 0, 5, '-'
    )
    assert_series_equal(test_series, expect_series)

### test bigBed with multiple peaks ###

def test_peak_values_p6():
    print("tests the values function when multiple overlapping peaks "
          "are present. This shouldn't happen with compressed.bed files "
          "since those don't have any overlapping peaks. For now, "
          "the expeted behavior is to add regions with overlapping peaks.")
    expect_series = pd.Series([1,1,1,1,1,2,2,2,2,2,1,1,1,1,1,0,0,0,0,0])
    test_peak = get_a_and_c()
    test_series = test_peak.values(
        'chr1', 0, 20, '+'
    )
    assert_series_equal(test_series, expect_series)

def test_peak_values_p7():
    print("tests the values function when multiple non-overlapping peaks "
          "are present.")
    expect_series = pd.Series([0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,1,1,1,1,1])
    test_peak = get_e_and_f()
    test_series = test_peak.values(
        'chr1', 0, 20, '+'
    )
    assert_series_equal(test_series, expect_series)