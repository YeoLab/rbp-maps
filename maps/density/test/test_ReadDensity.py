import os

import numpy as np
import pytest

from density import ReadDensity

curdir = os.path.dirname(__file__)


@pytest.fixture()
def rbp():
    """test_2000bp.pos.bg: chr1 0-100 = 3, 100-200 = 4, 300-400 = 5."""
    return ReadDensity.ReadDensity(
        pos=os.path.join(curdir, 'test_intervals/test_2000bp.pos.bw'),
        neg=os.path.join(curdir, 'test_intervals/test_2000bp.neg.bw'),
    )


def test_values_positive_strand(rbp):
    assert rbp.values('chr1', 98, 102, '+') == [3, 3, 4, 4]


def test_values_negative_strand_are_reversed(make_bigwig):
    rbp = ReadDensity.ReadDensity(
        pos=make_bigwig('pos.bw', [('chr1', 0, 10, 1)]),
        neg=make_bigwig('neg.bw', [('chr1', 0, 5, -1), ('chr1', 5, 10, -2)]),
    )
    assert rbp.values('chr1', 3, 7, '-') == [-2, -2, -1, -1]


def test_values_without_coverage_are_nan(rbp):
    assert np.isnan(rbp.values('chr1', 200, 203, '+')).all()


def test_values_missing_chromosome_returns_nan_per_position(rbp):
    values = rbp.values('chrMissing', 10, 15, '+')
    assert len(values) == 5
    assert np.isnan(values).all()


def test_values_invalid_strand_returns_1(rbp):
    assert rbp.values('chr1', 0, 5, '.') == 1


def test_default_name_masks_strand(rbp):
    assert rbp.name.endswith('test_2000bp.*.bw')


def test_explicit_name(make_bigwig):
    bw = make_bigwig('a.bw', [('chr1', 0, 10, 1)])
    assert ReadDensity.ReadDensity(pos=bw, neg=bw, name='AGO2').name == 'AGO2'


def test_pseudocount_and_total_mapped_come_from_bam(make_bigwig, make_bam):
    bw = make_bigwig('a.bw', [('chr1', 0, 10, 1)])
    bam = make_bam('a.bam', [('chr1', 10 * i, 20, False) for i in range(4)])
    rbp = ReadDensity.ReadDensity(pos=bw, neg=bw, bam=bam)
    assert rbp.total_mapped() == 4
    assert rbp.pseudocount() == 250000.0


def test_unopenable_bigwig_is_reported_not_raised(tmp_path, capsys):
    ReadDensity.ReadDensity(pos=str(tmp_path / 'no.bw'), neg=str(tmp_path / 'no.bw'))
    assert "couldn't open the bigwig files!" in capsys.readouterr().out
