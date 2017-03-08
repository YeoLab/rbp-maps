#!/usr/env python

import os
import pytest
import pybedtools
from density import intervals

### Fixtures ###

@pytest.fixture()
def current_pos_interval():
    return pybedtools.create_interval_from_list(
        ['chr1','0','10','current','0','+']
    )

@pytest.fixture()
def downstream_pos_interval():
    return pybedtools.create_interval_from_list(
        ['chr1','15','20','current','0','+']
    )

@pytest.mark.parametrize(("input","expected"), [
    ([current_pos_interval, downstream_pos_interval, 0, 0, '3p', False],[10, 15, 0, 0, 0]),
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


def test_get_boundaries_1(downstream_pos_interval, current_pos_interval):
    left_offset = 0
    right_offset = 0
    exon_junction_site = '3p'
    stop_at_mid = False

    test_anchor, test_upper_boundary, test_upper_offset, \
    test_lower_boundary, test_lower_offset = intervals._get_boundaries(
        downstream_pos_interval,
        current_pos_interval,
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

def test_get_boundaries_2(input, expected):
    print(input)
    # assert intervals._get_boundaries(input) == expected