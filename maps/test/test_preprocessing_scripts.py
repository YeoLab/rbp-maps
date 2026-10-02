import shutil
import sys

import pandas as pd
import pyBigWig
import pytest

from maps.conftest import rmats_line
from preprocessing_scripts import bed2bigbed
from preprocessing_scripts import subset_rmats_junctioncountonly as subset

SE_HEADER = '\t'.join([
    'ID', 'GeneID', 'geneSymbol', 'chr', 'strand', 'exonStart_0base', 'exonEnd',
    'upstreamES', 'upstreamEE', 'downstreamES', 'downstreamEE', 'ID.1',
    'IJC_SAMPLE_1', 'SJC_SAMPLE_1', 'IJC_SAMPLE_2', 'SJC_SAMPLE_2', 'IncFormLen',
    'SkipFormLen', 'PValue', 'FDR', 'IncLevel1', 'IncLevel2', 'IncLevelDifference',
])


def se_event(event_id, coords, inclusion_counts):
    """An SE row whose IJC_SAMPLE_1 and IJC_SAMPLE_2 are both `inclusion_counts`."""
    fields = rmats_line('chr1', '+', coords, event_id).split('\t')
    fields[12] = fields[14] = inclusion_counts
    return '\t'.join(fields)


### subset_jxc ###

def test_avg_inclusion_count_averages_all_replicates_of_both_samples():
    row = pd.Series({'IJC_SAMPLE_1': '10,20', 'IJC_SAMPLE_2': '30,40'})
    assert subset.get_avg_inclusion_count(row) == 25


def test_avg_inclusion_count_of_background_format():
    assert subset.get_avg_inclusion_count(pd.Series({'incl': '10,30'})) == 20


UPSTREAM_FLANK = {'flankingES': 100, 'flankingEE': 200, 'shortES': 650, 'shortEE': 700}
DOWNSTREAM_FLANK = {'flankingES': 700, 'flankingEE': 800, 'shortES': 100, 'shortEE': 150}


@pytest.mark.parametrize('event,strand,exons,expected', [
    ('se', '+', {}, (200, 700)),
    ('a3ss', '+', UPSTREAM_FLANK, (200, 650)),
    ('a3ss', '-', DOWNSTREAM_FLANK, (150, 700)),
    ('a5ss', '+', DOWNSTREAM_FLANK, (150, 700)),
    ('a5ss', '-', UPSTREAM_FLANK, (200, 650)),
])
def test_jx_region_spans_the_intron_next_to_the_flanking_exon(event, strand, exons, expected):
    """Coordinates are genomic: the flanking exon is at lower coordinates for a3ss (+) and a5ss (-)."""
    row = pd.Series(dict({'chr': 'chr1', 'strand': strand, 'upstreamEE': 200, 'downstreamES': 700}, **exons))
    interval = subset.get_jx_region_as_interval(row, 'name', event)
    assert (interval.start, interval.end) == expected


def test_subset_keeps_the_overlapping_event_with_the_most_inclusion_reads(write_lines, tmp_path):
    """Events 0 and 1 share the intron span 200-700; event 2 is alone at 1200-1700."""
    annotation = write_lines('se.txt', [
        SE_HEADER,
        se_event(0, [400, 500, 100, 200, 700, 800], '1,1'),
        se_event(1, [420, 480, 100, 200, 700, 800], '50,50'),
        se_event(2, [1400, 1500, 1100, 1200, 1700, 1800], '2,2'),
    ])
    out = str(tmp_path / 'se.nr.txt')
    merged = subset.run_subset_rmats_junctioncountonly(annotation, out, 'se')

    kept = pd.read_csv(out, sep='\t')
    assert list(kept['ID']) == [1, 2]
    assert list(kept.columns) == SE_HEADER.split('\t')
    assert list(merged['name'].astype(str)) == ['0,1', '2']


def test_subset_main_writes_the_output_file(write_lines, tmp_path, monkeypatch):
    annotation = write_lines('se.txt', [SE_HEADER, se_event(0, [400, 500, 100, 200, 700, 800], '1,1')])
    out = str(tmp_path / 'se.nr.txt')
    monkeypatch.setattr(sys, 'argv', ['subset_jxc', '-i', annotation, '-o', out, '-e', 'se'])
    subset.main()
    assert pd.read_csv(out, sep='\t').shape == (1, 23)


### bed2bigbed-eclip ###

def test_stringify_wraps_numeric_names():
    assert bed2bigbed.stringify('3.5') == '_3.5_'
    assert bed2bigbed.stringify('peak') == 'peak'


def test_filter_bed_keeps_peaks_passing_both_thresholds(write_lines, tmp_path):
    bed = write_lines('peaks.bed', [
        'chr1\t10\t20\t5.0\t4.0\t+',
        'chr1\t30\t40\t5.0\t1.0\t+',
        'chr1\t50\t60\t1.0\t4.0\t+',
    ])
    out = str(tmp_path / 'filtered.bed')
    bed2bigbed.filter_bed(bed, 3, 3, out)
    assert open(out).read().split('\t')[:3] == ['chr1', '10', '20']
    assert sum(1 for _ in open(out)) == 1


@pytest.mark.integration
@pytest.mark.parametrize('bed_type,expected', [
    ('bed6inputnorm', [(100, 200)]),
    ('bed6', [(10, 20), (100, 200)]),
])
def test_main_writes_a_sorted_bigbed(write_lines, tmp_path, monkeypatch, bed_type, expected):
    """Input-normalized peaks are filtered at -log10(p) >= 3 and log2 fold change >= 3."""
    if shutil.which('bedToBigBed') is None:
        pytest.skip('Missing required external tool: bedToBigBed')
    bed = write_lines('peaks.bed', ['chr1\t100\t200\t5.0\t4.0\t+', 'chr1\t10\t20\t1.0\t1.0\t+'])
    genome = write_lines('chr1.sizes', ['chr1\t2000'])
    out = str(tmp_path / 'peaks.bb')
    monkeypatch.setattr(sys, 'argv', [
        'bed2bigbed-eclip', '--beds', bed, '--genome', genome, '--outbbs', out, '--bedtype', bed_type,
    ])
    bed2bigbed.main()
    entries = pyBigWig.open(out).entries('chr1', 0, 2000)
    assert [(start, end) for start, end, _ in entries] == expected
