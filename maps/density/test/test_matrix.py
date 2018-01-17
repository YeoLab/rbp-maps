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
