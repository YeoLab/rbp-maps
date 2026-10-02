"""Shared fixtures: small bigWig/BAM/annotation builders on a 2 kb 'chr1'."""

import pyBigWig
import pysam
import pytest

from density import ReadDensity

CHROMS = (('chr1', 2000),)


class FakeDensity(object):
    """Density whose value at every position is the position itself."""

    def values(self, chrom, start, end, strand):
        values = [float(v) for v in range(start, end)]
        return values if strand == '+' else values[::-1]

    def pseudocount(self):
        return 1.0


@pytest.fixture()
def fake_density():
    return FakeDensity()


@pytest.fixture()
def make_bigwig(tmp_path):
    """Returns a function writing (chrom, start, end, value) entries to a bigWig."""

    def _make(name, entries, chroms=CHROMS):
        path = str(tmp_path / name)
        bw = pyBigWig.open(path, 'w')
        bw.addHeader(list(chroms))
        bw.addEntries(
            [e[0] for e in entries],
            [e[1] for e in entries],
            ends=[e[2] for e in entries],
            values=[float(e[3]) for e in entries],
        )
        bw.close()
        return path

    return _make


@pytest.fixture()
def make_bam(tmp_path):
    """Returns a function writing (chrom, start, length, is_reverse) reads to an indexed BAM."""

    def _make(name, reads, chroms=CHROMS):
        path = str(tmp_path / name)
        names = [c for c, _ in chroms]
        header = {
            'HD': {'VN': '1.0', 'SO': 'coordinate'},
            'SQ': [{'SN': c, 'LN': n} for c, n in chroms],
        }
        with pysam.AlignmentFile(path, 'wb', header=header) as out:
            for i, (chrom, start, length, is_reverse) in enumerate(sorted(reads)):
                read = pysam.AlignedSegment()
                read.query_name = 'read{}'.format(i)
                read.query_sequence = 'A' * length
                read.flag = 16 if is_reverse else 0
                read.reference_id = names.index(chrom)
                read.reference_start = start
                read.mapping_quality = 255
                read.cigar = ((0, length),)
                read.query_qualities = pysam.qualitystring_to_array('I' * length)
                out.write(read)
        pysam.index(path)
        return path

    return _make


@pytest.fixture()
def density_pair(make_bigwig, make_bam):
    """(ip, input) ReadDensity objects: IP is 2 RPM and input 1 RPM across chr1."""
    reads = [('chr1', 10 * i, 20, False) for i in range(4)]
    ip = ReadDensity.ReadDensity(
        pos=make_bigwig('ip.pos.bw', [('chr1', 0, 2000, 2)]),
        neg=make_bigwig('ip.neg.bw', [('chr1', 0, 2000, -2)]),
        bam=make_bam('ip.bam', reads),
    )
    inp = ReadDensity.ReadDensity(
        pos=make_bigwig('input.pos.bw', [('chr1', 0, 2000, 1)]),
        neg=make_bigwig('input.neg.bw', [('chr1', 0, 2000, -1)]),
        bam=make_bam('input.bam', reads),
    )
    return ip, inp


def rmats_line(chrom, strand, coords, event_id=0):
    """An rMATS JunctionCountOnly row; `coords` are the event's start/end columns."""
    fields = [str(event_id), 'ENSG{}'.format(event_id), 'GENE{}'.format(event_id), chrom, strand]
    fields += [str(c) for c in coords]
    fields += [str(event_id), '10,10', '5,5', '20,20', '5,5', '100', '50',
               '0.01', '0.05', '0.5,0.5', '0.8,0.8', '-0.3']
    return '\t'.join(fields)


@pytest.fixture()
def write_lines(tmp_path):
    """Returns a function writing lines to a file in tmp_path and returning its path."""

    def _write(name, lines):
        path = str(tmp_path / name)
        with open(path, 'w') as handle:
            handle.write('\n'.join(lines) + '\n')
        return path

    return _write


@pytest.fixture()
def bed_annotation(write_lines):
    """Returns a function writing n single-base BED6 sites, 10 nt apart, on one strand."""

    def _make(name, n=20, first=100, strand='+'):
        return write_lines(name, [
            'chr1\t{}\t{}\tsite{}\t0\t{}'.format(first + 10 * i, first + 10 * i + 1, i, strand)
            for i in range(n)
        ])

    return _make
