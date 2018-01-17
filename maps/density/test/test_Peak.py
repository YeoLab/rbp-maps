#!/usr/env python

import os
import pytest
import pybedtools
import pandas as pd
from density import intervals
from density import ReadDensity
from pandas.testing import assert_series_equal

### Fixtures ###


### These are for testing the interval/boundary regions.


@pytest.fixture()
def pos_chr1_0_10():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '0', '10', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_0_10():
    """ downstream interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '0', '10', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_3_5():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '3', '5', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_3_5():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '3', '5', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_5_10():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '5', '10', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_5_10():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '5', '10', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_5_15():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '5', '15', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_5_15():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '5', '15', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_0_5():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '0', '5', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_0_5():
    """ upstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '0', '5', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_15_20():
    """ interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '15', '20', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_15_20():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '15', '20', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_25_30():
    """ downstream interval (positive) """
    return pybedtools.create_interval_from_list(
        ['chr1', '25', '30', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_25_30():
    """ upstream interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '25', '30', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_90_95():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '90', '95', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_90_95():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '90', '95', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_95_105():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '95', '105', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_95_105():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '95', '105', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_105_115():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '105', '115', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_105_115():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '105', '115', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_95_100():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '95', '100', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_95_100():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '95', '100', 'current', '0', '-']
    )

@pytest.fixture()
def pos_chr1_105_115():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '105', '115', 'current', '0', '+']
    )
@pytest.fixture()
def neg_chr1_105_115():
    """ interval (negative) """
    return pybedtools.create_interval_from_list(
        ['chr1', '105', '115', 'current', '0', '-']
    )

### These are for metagenes, the 'name' field is important!

@pytest.fixture()
def non_overlapping_gene_1():
    return 'test_intervals/gene1_nonoverlapping.bed'

@pytest.fixture()
def overlapping_gene_1():
    return 'test_intervals/gene1_overlapping.bed'

@pytest.fixture()
def overlapping_gene_1_gene_2():
    return 'test_intervals/gene1_gene2_overlapping.bed'

### These are for scale tests

@pytest.fixture()
def n():
    return 0


@pytest.fixture()
def some_small_series_not_divisible():
    return pd.Series(1)


@pytest.fixture()
def some_med_series_not_divisible():
    return pd.Series([range(0, 10)])


@pytest.fixture()
def some_large_series_not_divisible():
    return pd.Series([range(0, 1023)])


@pytest.fixture()
def some_large_series_divisible():
    return pd.Series([range(0, 1000)])

@pytest.fixture()
def get_test_rbp():
    curdir = os.path.dirname(__file__)
    pos = os.path.join(curdir,'test_intervals/test_2000bp.pos.bw')
    neg = os.path.join(curdir,'test_intervals/test_2000bp.neg.bw')
    test_rbp = ReadDensity.ReadDensity(pos=pos, neg=neg, name=None, bam=None)
    return test_rbp

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

### test get upper boundaries ###

@pytest.mark.parametrize(
    "current_interval, next_interval, strand_or_5p, stop_at_midpoint, \
    expect_boundary",
    [
        (pos_chr1_0_10(), pos_chr1_15_20(), '+', False, 15),
        (pos_chr1_15_20(), pos_chr1_0_10(), '-', False, 20),
    ]
)
def test_get_upper_boundary(
        current_interval, next_interval, strand_or_5p, stop_at_midpoint,
        expect_boundary
):
    test_boundary = intervals._get_upper_boundary(
        current_interval, next_interval, strand_or_5p, stop_at_midpoint
    )
    assert expect_boundary == test_boundary

### test get lower boundaries ###

@pytest.mark.parametrize(
    "current_interval, next_interval, strand_or_5p, stop_at_midpoint, \
    expect_boundary",
    [
        (pos_chr1_0_10(), pos_chr1_15_20(), '+', False, 0),
        (pos_chr1_15_20(), pos_chr1_0_10(), '-', False, 10),
    ]
)
def test_get_lower_boundary(
        current_interval, next_interval, strand_or_5p, stop_at_midpoint,
        expect_boundary
):
    test_boundary = intervals._get_lower_boundary(
        current_interval, next_interval, strand_or_5p, stop_at_midpoint
    )
    assert expect_boundary == test_boundary

### test get boundaries ###

@pytest.mark.parametrize(
    "downstream_interval, interval, expected_anchor, \
    upstream_offset, downstream_offset, \
    expected_upper_boundary, expected_lower_boundary, \
    expected_upper_offset, expected_lower_offset",
    [
        # downstream + upstream + anchor
        (pos_chr1_15_20(), pos_chr1_0_10(), 10,
         0, 0,  # upstream (exon) offset, downstream (intron) offset
         15, 0,  # upper genomic boundary, lower genomic boundary
         0, 0),  # upper genomic offset, lower genomic offset
        (neg_chr1_0_10(), neg_chr1_15_20(), 15,
         0, 0,
         20, 10,
         0, 0),
        (pos_chr1_15_20(), pos_chr1_0_10(), 10,
         3, 4,
         15, 0,
         4, 3),
        (neg_chr1_0_10(), neg_chr1_15_20(), 15,
         3, 4,
         20, 10,
         3, 4),
        (pos_chr1_15_20(), pos_chr1_0_10(), 10,
         11, 0,
         15, 0,
         0, 11),
        (neg_chr1_0_10(), neg_chr1_15_20(), 15,
         11, 0,
         20, 10,
         11, 0),
    ]
)
def test_get_boundaries_3p(
        downstream_interval, interval, expected_anchor,
        upstream_offset, downstream_offset,
        expected_upper_boundary, expected_upper_offset,
        expected_lower_boundary, expected_lower_offset,
        region='3p',
):
    exon_junction_site = region
    stop_at_mid = False

    test_anchor, test_upper_boundary, test_upper_offset, \
    test_lower_boundary, test_lower_offset = intervals._get_boundaries(
        downstream_interval,
        interval,
        upstream_offset,
        downstream_offset,
        exon_junction_site,
        stop_at_mid
    )
    assert test_anchor == expected_anchor
    assert test_upper_boundary == expected_upper_boundary
    assert test_lower_boundary == expected_lower_boundary
    assert test_upper_offset == expected_upper_offset
    assert test_lower_offset == expected_lower_offset


@pytest.mark.parametrize(
    "upstream_interval, interval, expected_anchor, \
    exon_offset, intron_offset, \
    expected_upper_boundary, expected_lower_boundary, \
    expected_upper_offset, expected_lower_offset",
    [
        # upstream + interval + anchor
        (pos_chr1_0_10(), pos_chr1_15_20(), 15,
         0, 0,  # exon offset, intron offset
         20, 10,  # upper genomic boundary, lower genomic boundary
         0, 0),  # upper genomic offset, lower genomic offset
        (neg_chr1_15_20(), neg_chr1_0_10(), 10,
         0, 0,
         15, 0,
         0, 0),
        (pos_chr1_0_10(), pos_chr1_15_20(), 15,
         3, 4,
         20, 10,
         3, 4),
        (neg_chr1_15_20(), neg_chr1_0_10(), 10,
         3, 4,
         15, 0,
         4, 3),
        (pos_chr1_0_10(), pos_chr1_15_20(), 15,
         6, 0,
         20, 10,
         6, 0),
        (neg_chr1_15_20(), neg_chr1_0_10(), 10,
         11, 0,
         15, 0,
         0, 11),
    ]
)
def test_get_boundaries_5p(
        upstream_interval, interval, expected_anchor,
        exon_offset, intron_offset,
        expected_upper_boundary, expected_upper_offset,
        expected_lower_boundary, expected_lower_offset,
        region='5p'
):
    exon_junction_site = region
    stop_at_mid = False

    test_anchor, test_upper_boundary, test_upper_offset, \
    test_lower_boundary, test_lower_offset = intervals._get_boundaries(
        upstream_interval,
        interval,
        exon_offset,
        intron_offset,
        exon_junction_site,
        stop_at_mid
    )
    assert test_anchor == expected_anchor
    assert test_upper_boundary == expected_upper_boundary
    assert test_lower_boundary == expected_lower_boundary
    assert test_upper_offset == expected_upper_offset
    assert test_lower_offset == expected_lower_offset


### get absolute coords and pad ###

@pytest.mark.parametrize(
    "anchor, upper_boundary, upper_offset, lower_boundary, lower_offset, \
    expect_left_pad, expect_abs_start, expect_abs_end, expect_right_pad",
    [
        # zero padding, just expects to return anchor position
        (10, 15, 0, 0, 0,
         0, 10, 10, 0),
        # 1 upper+lower padding, return +/-1 positions from anchor
        (10, 15, 1, 0, 1,
         0, 9, 11, 0),
        # 6 upper padding + 0 lower, returns 1 upper/right pad
        (10, 15, 6, 0, 0,
         0, 10, 15, 1),
        # 6 lower padding + 0 higher, returns 1 lower/left pad
        (10, 15, 0, 0, 11,
         1, 0, 10, 0),
        # 6 upper padding + 11 lower, returns 1 left+right pad
        (10, 15, 7, 0, 11,
         1, 0, 15, 2),
    ]
)

def test_get_absolute_coords_and_pad_1(
        anchor, upper_boundary, upper_offset, lower_boundary, lower_offset,
        expect_left_pad, expect_abs_start, expect_abs_end, expect_right_pad
):
    test_left_pad, test_abs_start, \
    test_abs_end, test_right_pad = intervals._get_absolute_coords_and_pad(
        anchor, upper_boundary, upper_offset, lower_boundary, lower_offset
    )
    assert test_left_pad == expect_left_pad
    assert test_right_pad == expect_right_pad
    assert test_abs_start == expect_abs_start
    assert test_abs_end == expect_abs_end


### test flip_strand ###
def test_flip_strand_pos():
    assert intervals.flip_strand('+') == '-'
    assert intervals.flip_strand('-') == '+'

### junction_site ###

def test_junction_site_1(pos_chr1_0_10, pos_chr1_15_20, get_test_rbp):
    left, wiggle, right = intervals._junction_site(
        # gets the wiggle range from 10 to 20
        get_test_rbp, pos_chr1_0_10, pos_chr1_15_20, 5, 5, '5p', False
    )
    assert len(wiggle) == 10
    assert left == 0
    assert right == 0
    assert wiggle == [3,3,3,3,3,3,3,3,3,3]


def test_junction_site_2(neg_chr1_0_10, neg_chr1_15_20, get_test_rbp):
    left, wiggle, right = intervals._junction_site(
        # gets the wiggle range from 10 to 20
        get_test_rbp, neg_chr1_0_10, neg_chr1_15_20, 5, 5, '3p', False
    )
    assert len(wiggle) == 10
    assert left == 0
    assert right == 0
    assert wiggle == [-3, -3, -3, -3, -3, -3, -3, -3, -3, -3]


def test_junction_site_3(neg_chr1_95_100, neg_chr1_105_115, get_test_rbp):
    left, wiggle, right = intervals._junction_site(
        # gets the wiggle range from 10 to 20
        get_test_rbp, neg_chr1_95_100, neg_chr1_105_115, 2, 8, '3p', False
    )
    assert len(wiggle) == 7
    assert left == 0
    assert right == 3
    assert wiggle == [-4, -4, -4, -4, -4, -4, -4]




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


### test merging (for metagene) function

def test_merge_1():
    print("Tests merge using a single interval. "
          "Should return the same interval as there "
          "are no others to merge")
    df = intervals.merge(non_overlapping_gene_1())
    print(df)
    assert df.shape[0] == 3

def test_merge_2():
    print("Tests three intervals, two of which overlap. "
          "Should return two intervals, with one being "
          "merged result from the two overlapping ones.")
    df = intervals.merge(overlapping_gene_1())
    print(df)
    assert df.shape[0] == 2

def test_merge_3():
    print("Tests three intervals, two of which overlap "
          "but belong to different genes. Should return "
          "three intervals.")
    df = intervals.merge2(overlapping_gene_1_gene_2())
    print(df)
    assert df.shape[0] == 3

### test peak overlap intervals:


def test_get_overlap_p1():
    print("tests condition when (+) newPeak partially overlaps region"
          " but shares end coords")
    peak = pos_chr1_5_10()
    region = pos_chr1_0_10()
    expect_overlap = pd.Series(
        data=[0,0,0,0,0,1,1,1,1,1],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_n1():
    print("tests condition when (-) newPeak partially overlaps region"
          " but shares end coords")
    peak = neg_chr1_5_10()
    region = neg_chr1_0_10()
    expect_overlap = pd.Series(
        data=[1,1,1,1,1,0,0,0,0,0],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_p2():
    print("tests condition when (+) newPeak partially overlaps region")
    peak = pos_chr1_5_15()
    region = pos_chr1_0_10()
    expect_overlap = pd.Series(
        data=[0,0,0,0,0,1,1,1,1,1],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_n2():
    print("tests condition when (+) newPeak partially overlaps region")
    peak = neg_chr1_5_15()
    region = neg_chr1_0_10()
    expect_overlap = pd.Series(
        data=[1,1,1,1,1,0,0,0,0,0],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_p3():
    print("tests condition when (+) newPeak completely overlaps region")
    peak = pos_chr1_5_15()
    region = pos_chr1_5_15()
    expect_overlap = pd.Series(
        data=[1,1,1,1,1,1,1,1,1,1],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_n3():
    print("tests condition when (-) newPeak completely overlaps region")
    peak = neg_chr1_5_15()
    region = neg_chr1_5_15()
    expect_overlap = pd.Series(
        data=[1,1,1,1,1,1,1,1,1,1],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_p4():
    print("tests condition when (+) newPeak doesn't overlap region")
    peak = pos_chr1_5_15()
    region = pos_chr1_0_5()
    expect_overlap = pd.Series(
        data=[0,0,0,0,0],
        index=range(5)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_n4():
    print("tests condition when (+) newPeak doesn't overlap region")
    peak = neg_chr1_5_15()
    region = neg_chr1_0_5()
    expect_overlap = pd.Series(
        data=[0,0,0,0,0],
        index=range(5)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_p5():
    print("tests condition when (+) newPeak is completely "
          "contained within the region")
    peak = pos_chr1_3_5()
    region = pos_chr1_0_10()
    expect_overlap = pd.Series(
        data=[0,0,0,1,1,0,0,0,0,0],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_n5():
    print("tests condition when (-) newPeak is completely "
          "contained within the region")
    peak = neg_chr1_3_5()
    region = neg_chr1_0_10()
    expect_overlap = pd.Series(
        data=[0,0,0,0,0,1,1,0,0,0],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region)
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)
    assert 1==
def test_get_overlap_p6():
    print("tests condition when (+) newPeak partially overlaps region"
          " but shares end coords using the region frac scoring method.")
    peak = pos_chr1_5_10()
    region = pos_chr1_0_10()
    expect_overlap = pd.Series(
        data=[0,0,0,0,0,.1,.1,.1,.1,.1],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region, 'fraction_region')
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_n6():
    print("tests condition when (-) newPeak partially overlaps region"
          " but shares end coords using the region frac scoring method.")
    peak = neg_chr1_5_10()
    region = neg_chr1_0_10()
    expect_overlap = pd.Series(
        data=[.1,.1,.1,.1,.1,0,0,0,0,0],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region, 'fraction_region')
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_p7():
    print("tests condition when (+) newPeak completely overlaps region")
    peak = pos_chr1_5_15()
    region = pos_chr1_5_15()
    expect_overlap = pd.Series(
        data=[.1,.1,.1,.1,.1,.1,.1,.1,.1,.1],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region, 'fraction_region')
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)

def test_get_overlap_n7():
    print("tests condition when (-) newPeak completely overlaps region")
    peak = neg_chr1_5_15()
    region = neg_chr1_5_15()
    expect_overlap = pd.Series(
        data=[.1,.1,.1,.1,.1,.1,.1,.1,.1,.1],
        index=range(10)
    )
    test_overlap = intervals.get_overlap(peak, region, 'fraction_region')
    print('test', test_overlap)
    print('expect', expect_overlap)
    assert_series_equal(test_overlap, expect_overlap)
