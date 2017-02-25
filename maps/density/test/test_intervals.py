#!/usr/env python

import os
import pytest
import pybedtools as bt

import ReadDensity
import intervals
import Feature

# Fixtures


wd = '/home/bay001/projects/codebase/rbfox2/'

regions = os.path.join(
    wd,
    'RBFOX2-BGHLV26-HepG2-SE.MATS.JunctionCountOnly.txt'
)

ip_bam = os.path.join(
    wd,
    '204_01_RBFOX2.merged.r2.bam'
)
input_bam = os.path.join(
    wd,
    'RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.bam'
)

ip_pos_bw = os.path.join(
    wd,
    '204_01_RBFOX2.merged.r2.norm.neg.bw'
) # flipped?
ip_neg_bw = os.path.join(
    wd,
    '204_01_RBFOX2.merged.r2.norm.pos.bw'
) # flipped?
input_pos_bw = os.path.join(
    wd,
    'RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.neg.bw'
) # flipped?
input_neg_bw = os.path.join(
    wd,
    'RBFOX2-204-INPUT_S2_R1.unassigned.adapterTrim.round2.rmRep.rmDup.sorted.r2.norm.pos.bw'
) # flipped?

ip_read_density = ReadDensity.ReadDensity(
    pos = ip_pos_bw,
    neg = ip_neg_bw,
    bam = ip_bam
)

input_read_density = ReadDensity.ReadDensity(
    pos = input_pos_bw,
    neg = input_neg_bw,
    bam = input_bam
)

@pytest.fixture()
def get_full_upstream_region(regions=regions):

    with open(regions,'r') as f:
        next(f)
        for line in f: # skip first line due to comment/header
            line = line.strip()
            up, skip, down = Feature.SkippedExonFeature(
                line, 'rmats'
            ).get_bedtools()
            if up.end - up.start > 100:
                print up
                return up

### Tests ###

def test__five_prime_site_full_length():
    """
    five_prime_site(
        rbp,
        upstream_interval,
        interval,
        exon_offset,
        intron_offset,
        trunc=True,
        middle_stop=False)

    Returns
    -------

    """
    expect_fivep_pad = 0
    with open('test_intervals/chr11-70266505-70266616-0-0-+_chr11-70266505-70266616-0-0-+.densities.txt','r') as f:
        for line in f:
            expect_wiggle_track.append(float(line.replace('nan',np.NaN)))
    expect_threep_pad = 0
    neg_strand_list = [u'chr5', u'140998364', u'140998566', u'0', u'0', u'-']
    negative = bt.create_interval_from_list(neg_strand_list)

    pos_strand_list = [u'chr11', u'70266505', u'70266616', u'0', u'0', u'+']
    pos_strand_upstream_list = [u'chr11', u'70266105', u'70266205', u'0', u'0', u'+']
    positive = bt.create_interval_from_list(pos_strand_list)
    positive_upstream = bt.create_interval_from_list(pos_strand_upstream_list)

    obs_fivep_pad, obs_wiggle_track, obs_threep_pad = intervals.five_prime_site(
        rbp=ip_read_density,
        upstream_interval=positive_upstream,
        interval=positive,
        exon_offset=50,
        intron_offset=300,
        trunc=True,
        middle_stop=False
    )

def test__five_prime_site_short_intron():
    pass

def test__five_prime_site_short_exon():
    pass


five_prime_site(
        rbp,
        upstream_interval,
        interval,
        exon_offset,
        intron_offset,
        trunc=True,
        middle_stop=False