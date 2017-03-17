#!/usr/env python

import os
import pytest
import pybedtools
import pandas as pd
from density import intervals

### Fixtures ###

@pytest.fixture()
def pos_chr1_0_10():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1','0','10','current','0','+']
    )

@pytest.fixture()
def pos_chr1_15_20():
    """ interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1','15','20','current','0','+']
    )

@pytest.fixture()
def pos_chr1_25_30():
    """ downstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1','25','30','current','0','+']
    )

@pytest.fixture()
def neg_chr1_25_30():
    """ upstream interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1','25','30','current','0','-']
    )

@pytest.fixture()
def neg_chr1_15_20():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '15', '20', 'current', '0', '-']
    )

@pytest.fixture()
def neg_chr1_0_10():
    """ downstream interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1','0','10','current','0','-']
    )

@pytest.fixture()
def n():
    return 0

@pytest.fixture()
def some_small_series_not_divisible():
    return pd.Series(1)

@pytest.fixture()
def some_med_series_not_divisible():
    return pd.Series([range(0,10)])

@pytest.fixture()
def some_large_series_not_divisible():
    return pd.Series([range(0,1023)])

@pytest.fixture()
def some_large_series_divisible():
    return pd.Series([range(0,1000)])

@pytest.mark.parametrize("upstream_interval", "interval", "expected", [
    (
        ["upstream_interval"],[pos_chr1_15_20()],
        ["current_interval"],[pos_chr1_25_30()]
    )
])
### Tests ###


def test_too_far_1():
    direction = 1
    anchor = 10
    offset = 5
    boundary = 16
    ref = 0
    test = intervals._too_far(anchor, offset, boundary, direction)
    assert test == ref


def test_too_far_2():
    direction = 1
    anchor = 10
    offset = 5
    boundary = 15
    ref = 0
    test = intervals._too_far(anchor, offset, boundary, direction)
    assert test == ref


def test_too_far_3():
    direction = 1
    anchor = 10
    offset = 5
    boundary = 14
    ref = 1
    test = intervals._too_far(anchor, offset, boundary, direction)
    assert test == ref


def test_too_far_4():
    direction = 1
    anchor = 10
    offset = 5
    boundary = 10
    ref = 5
    test = intervals._too_far(anchor, offset, boundary, direction)
    assert test == ref


def test_too_far_5():
    direction = -1
    anchor = 10
    offset = 5
    boundary = 9
    ref = 4
    test = intervals._too_far(anchor, offset, boundary, direction)
    assert test == ref


def test_too_far_6():
    direction = -1
    anchor = 10
    offset = 5
    boundary = 11
    ref = 6
    test = intervals._too_far(anchor, offset, boundary, direction)
    assert test == ref


def test_get_boundaries_1(pos_chr1_0_10, pos_chr1_15_20):
    left_offset = 0
    right_offset = 0
    exon_junction_site = '3p'
    stop_at_mid = False

    test_anchor, test_upper_boundary, test_upper_offset, \
    test_lower_boundary, test_lower_offset = intervals._get_boundaries(
        pos_chr1_15_20,
        pos_chr1_0_10,
        left_offset,
        right_offset,
        exon_junction_site,
        stop_at_mid
    )

    ref_anchor = 10
    ref_upper_boundary = 15
    ref_upper_offset = 0
    ref_lower_boundary = 0
    ref_lower_offset = 0

    assert test_anchor == ref_anchor
    assert test_upper_offset == ref_upper_offset
    assert test_lower_offset == ref_lower_offset
    assert test_upper_boundary == ref_upper_boundary
    assert test_lower_boundary == ref_lower_boundary

def test_get_boundaries_2(upstream_interval, interval, expected):
    print(upstream_interval, interval)
    assert intervals._get_boundaries(input) == expected


### multiply by 100 tests ###

def test_multiply_by_100(n):
    assert len(intervals.multiply_by_100(n)) == 100

### get scale tests ###
def test_get_scale_1(some_small_series_not_divisible):
    assert len(intervals.get_scale(some_small_series_not_divisible)) % 100 == 0

def test_get_scale_2(some_med_series_not_divisible):
    assert len(intervals.get_scale(some_med_series_not_divisible)) % 100 == 0

def test_get_scale_3(some_large_series_not_divisible):
    assert len(intervals.get_scale(some_large_series_not_divisible)) % 100 == 0

def test_get_scale_4(some_large_series_divisible):
    assert len(intervals.get_scale(some_large_series_divisible)) % 100 == 0


