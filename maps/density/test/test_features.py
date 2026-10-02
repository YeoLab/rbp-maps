#!/usr/env python

import os
import pandas as pd
import pybedtools
import pytest
from density import Feature

### Fixtures ###

curdir = os.path.dirname(__file__)

se_file_rmats = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-SE.MATS.JunctionCountOnly.txt'
)
se_file_miso = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-SE.MISO.JunctionCountOnly.txt'
)

ri_file_rmats = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-RI.MATS.JunctionCountOnly.txt'
)
ri_file_miso = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-RI.MISO.JunctionCountOnly.txt'
)

a3ss_file_rmats = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-A3SS.MATS.JunctionCountOnly.txt'
)
a3ss_file_miso = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-A3SS.MISO.JunctionCountOnly.txt'
)

a5ss_file_rmats = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-A5SS.MATS.JunctionCountOnly.txt'
)
a5ss_file_miso = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-A5SS.MISO.JunctionCountOnly.txt'
)

mxe_file_rmats = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-MXE.MATS.JunctionCountOnly.txt'
)
mxe_file_miso = os.path.join(
    curdir,
    'test_features/RBFOX2-BGHLV26-HepG2-MXE.MISO.JunctionCountOnly.txt'
)

### Tests ###



def shares_a5ss_boundary(longer, shorter):
    """
    Returns True if the longer/shorter a5ss exons share the correct
    5' boundary

    Parameters
    ----------
    longer
    shorter

    Returns
    -------

    """
    if longer.strand == '+':
        return True if longer.start == shorter.start else False
    elif longer.strand == '-':
        return True if longer.end == shorter.end else False
    else:
        return False


def shares_a3ss_boundary(longer, shorter):
    """
    Returns True if the longer/shorter a3ss exons share the correct
    3' boundary

    Parameters
    ----------
    longer
    shorter

    Returns
    -------

    """
    if longer.strand == '+':
        return True if longer.end == shorter.end else False
    elif longer.strand == '-':
        return True if longer.start == shorter.start else False
    return False


def shares_strand(longer, shorter):
    """
    Returns True if the strands are shared between bedtools

    Parameters
    ----------
    longer
    shorter

    Returns
    -------

    """
    return True if longer.strand == shorter.strand else False


def is_longer(longer, shorter):
    """
    Returns true if the longer interval is indeed longer.

    Parameters
    ----------
    longer
    shorter

    Returns
    -------

    """
    return True if len(longer) > len(shorter) else False


def test_a5ss_rmats_feature():
    """
    Tests the A5SS Feature.

    Returns
    -------

    """
    annotation_type = 'rmats'
    with open(a5ss_file_rmats, 'r') as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID'):
                event = line.rstrip()
                alt1, alt2, downstream = Feature.Alt_5p_splice_site(
                    event, annotation_type
                ).get_bedtools()
                assert is_longer(alt1, alt2)
                assert shares_strand(alt1, alt2)
                assert shares_a5ss_boundary(alt1, alt2)


def test_a3ss_rmats_feature():
    """
        Tests the A5SS Feature.

        Returns
        -------

        """
    annotation_type = 'rmats'
    with open(a3ss_file_rmats, 'r') as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID'):
                event = line.rstrip()
                upstream, alt1, alt2 = Feature.Alt_3p_splice_site(
                    event, annotation_type
                ).get_bedtools()
                assert is_longer(alt1, alt2)
                assert shares_strand(alt1, alt2)
                assert shares_a3ss_boundary(alt1, alt2)


### Tests on single synthetic events ###

def coords(interval):
    return interval.chrom, interval.start, interval.end, interval.strand


def rmats(chrom, strand, positions):
    fields = ['0', 'ENSG', 'GENE', chrom, strand] + [str(p) for p in positions]
    fields += ['0', '10,10', '5,5', '20,20', '5,5', '100', '50',
               '0.01', '0.05', '0.5,0.5', '0.8,0.8', '-0.3']
    return '\t'.join(fields)


def test_bed_feature_ignores_columns_past_the_sixth():
    interval = Feature.Feature('chr1\t10\t20\tname\t5\t-\textra\n', 'bed').get_bedtool()
    assert coords(interval) == ('chr1', 10, 20, '-')
    assert interval.name == 'name'


@pytest.mark.parametrize("strand,expected_up,expected_down", [
    ('+', ('chr1', 100, 200, '+'), ('chr1', 700, 800, '+')),
    ('-', ('chr1', 700, 800, '-'), ('chr1', 100, 200, '-')),
])
def test_se_rmats_upstream_follows_strand(strand, expected_up, expected_down):
    """rMATS 'upstream' is always the genomically lower exon."""
    up, se, down = Feature.Skipped_exon(
        rmats('chr1', strand, [400, 500, 100, 200, 700, 800]), 'rmats'
    ).get_bedtools()
    assert coords(up) == expected_up
    assert coords(se) == ('chr1', 400, 500, strand)
    assert coords(down) == expected_down


def test_se_rmats_file_is_ordered_5p_to_3p():
    with open(se_file_rmats) as f:
        for line in f:
            if line.startswith('ID'):
                continue
            up, se, down = Feature.Skipped_exon(line, 'rmats').get_bedtools()
            if se.strand == '+':
                assert up.end <= se.start and se.end <= down.start
            else:
                assert up.start >= se.end and se.start >= down.end


def test_se_miso_is_converted_to_zero_based():
    up, se, down = Feature.Skipped_exon(
        'chr1:1230097:1230196:-@chr1:1229781:1230008:-@chr1:1229469:1229579:-\tENSG',
        'miso'
    ).get_bedtools()
    assert coords(up) == ('chr1', 1230096, 1230196, '-')
    assert coords(se) == ('chr1', 1229780, 1230008, '-')
    assert coords(down) == ('chr1', 1229468, 1229579, '-')


@pytest.mark.parametrize("strand,expected_up,expected_down", [
    ('+', (100, 200), (700, 800)),
    ('-', (700, 800), (100, 200)),
])
def test_se_tab_swaps_flanks_on_negative_strand(strand, expected_up, expected_down):
    line = 'chr1|{}|200-400|500-700|400-500\t100-200\t400-500\t700-800\t10,10\t5,5'.format(strand)
    up, se, down = Feature.Skipped_exon(line, 'tab').get_bedtools()
    assert (up.start, up.end) == expected_up
    assert (se.start, se.end) == (400, 500)
    assert (down.start, down.end) == expected_down


@pytest.mark.parametrize("strand,expected_up,expected_down", [
    ('+', (100, 200), (700, 800)),
    ('-', (700, 800), (100, 200)),
])
def test_se_bed12(strand, expected_up, expected_down):
    line = 'chr1\t100\t800\tname\t0\t{}\t100\t800\t0\t3\t100,100,100\t0,300,600'.format(strand)
    up, se, down = Feature.Skipped_exon(line, 'bed12').get_bedtools()
    assert (up.start, up.end) == expected_up
    assert (se.start, se.end) == (400, 500)
    assert (down.start, down.end) == expected_down


def test_se_unknown_format_returns_no_intervals():
    assert Feature.Skipped_exon('anything', 'unknown').get_bedtools() == (None, None, None)


@pytest.mark.parametrize("strand,expected_up,expected_down", [
    ('+', (100, 200), (700, 800)),
    ('-', (700, 800), (100, 200)),
])
def test_ri_rmats(strand, expected_up, expected_down):
    up, down = Feature.Retained_intron(
        rmats('chr1', strand, [100, 800, 100, 200, 700, 800]), 'rmats'
    ).get_bedtools()
    assert (up.start, up.end) == expected_up
    assert (down.start, down.end) == expected_down


@pytest.mark.parametrize("strand,expected_up,expected_down", [
    ('+', (100, 200), (700, 800)),
    ('-', (700, 800), (100, 200)),
])
def test_ri_twobed(strand, expected_up, expected_down):
    line = 'chr1\t100\t200\tlow\t0\t{0}\tchr1\t700\t800\thigh\t0\t{0}'.format(strand)
    up, down = Feature.Retained_intron(line, 'twobed').get_bedtools()
    assert (up.start, up.end) == expected_up
    assert (down.start, down.end) == expected_down


def test_ri_twobed_mixed_strands_is_rejected():
    line = 'chr1\t100\t200\tlow\t0\t+\tchr1\t700\t800\thigh\t0\t-'
    assert Feature.Retained_intron(line, 'twobed').get_bedtools() == -1


@pytest.mark.parametrize("line,expected_up,expected_down", [
    ('EXOSC8_ENSG00000120699.8;RI:chr13:37577071:37577144-37578614:37578698:+',
     (37577071, 37577144), (37578614, 37578698)),
    ('CCT8_ENSG00000156261.8;RI:chr21:30434649:30434736-30434811:30434896:-',
     (30434811, 30434896), (30434649, 30434736)),
])
def test_ri_xintao(line, expected_up, expected_down):
    up, down = Feature.Retained_intron(line, 'xintao').get_bedtools()
    assert (up.start, up.end) == expected_up
    assert (down.start, down.end) == expected_down


@pytest.mark.parametrize("strand,expected_up,expected_down", [
    ('+', (100, 200), (700, 800)),
    ('-', (700, 800), (100, 200)),
])
def test_ri_tab(strand, expected_up, expected_down):
    line = 'chr1|{}|x:y:100-200|700-800:z'.format(strand)
    up, down = Feature.Retained_intron(line, 'tab').get_bedtools()
    assert (up.start, up.end) == expected_up
    assert (down.start, down.end) == expected_down


def test_mxe_rmats_file_is_ordered_5p_to_3p():
    with open(mxe_file_rmats) as f:
        for line in f:
            if line.startswith('ID'):
                continue
            up, up_mxe, down_mxe, down = Feature.Mutually_exclusive_exon(
                line, 'rmats'
            ).get_bedtools()
            if up.strand == '+':
                assert up.end <= up_mxe.start < down_mxe.start < down.start
            else:
                assert up.start >= up_mxe.end > down_mxe.end > down.end


@pytest.mark.parametrize("strand,long_exon,short_exon", [
    ('+', (100, 250), (100, 200)),
    ('-', (650, 800), (700, 800)),
])
def test_a5ss_tab(strand, long_exon, short_exon):
    flank = '700-800' if strand == '+' else '100-200'
    line = 'chr1|{}|a|b|c\t{}-{}\t{}-{}\t{}\t10,10\t5,5'.format(
        strand, short_exon[0], short_exon[1], long_exon[0], long_exon[1], flank
    )
    alt1, alt2, downstream = Feature.Alt_5p_splice_site(line, 'tab').get_bedtools()
    assert (alt1.start, alt1.end) == long_exon
    assert (alt2.start, alt2.end) == short_exon
    assert is_longer(alt1, alt2) and shares_a5ss_boundary(alt1, alt2)


def test_a5ss_miso_positive_strand():
    alt1, alt2, downstream = Feature.Alt_5p_splice_site(
        'chr1:101:250|200:+@chr1:701:800:+\tENSG', 'miso'
    ).get_bedtools()
    assert coords(alt1) == ('chr1', 100, 250, '+')
    assert coords(alt2) == ('chr1', 100, 200, '+')
    assert coords(downstream) == ('chr1', 700, 800, '+')


def test_a3ss_tab():
    line = 'chr1|+|a|b|c\t100-200\t650-800\t700-800\t10,10\t5,5'
    upstream, alt1, alt2 = Feature.Alt_3p_splice_site(line, 'tab').get_bedtools()
    assert (upstream.start, upstream.end) == (100, 200)
    assert is_longer(alt1, alt2) and shares_a3ss_boundary(alt1, alt2)


def test_a3ss_miso_positive_strand():
    upstream, alt1, alt2 = Feature.Alt_3p_splice_site(
        'chr1:101:200:+@chr1:651|701:800:+\tENSG', 'miso'
    ).get_bedtools()
    assert coords(upstream) == ('chr1', 100, 200, '+')
    assert coords(alt1) == ('chr1', 650, 800, '+')
    assert coords(alt2) == ('chr1', 700, 800, '+')


def test_metafeature_returns_sorted_intervals():
    lines = ['chr1\t300\t400\ttx\t0\t+', 'chr1\t100\t200\ttx\t0\t+']
    feature = Feature.MetaFeature(lines, 'bed').get_bedtools()
    assert [(i.start, i.end) for i in feature] == [(100, 200), (300, 400)]


def test_get_random_sample_draws_n_rows_with_replacement():
    df = pd.DataFrame({'a': range(5)}, index=list('vwxyz'))
    sample = Feature.get_random_sample(df, 20)
    assert sample.shape == (20, 1)
    assert list(sample.index) == list(range(20))
    assert set(sample['a']) <= set(range(5))
