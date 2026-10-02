"""
Stranded maps around transcription start sites (TSS).

Windows are built from TSS calls (CAGE/RAMPAGE peaks, or any BED6 file), read
density is split into sense and antisense relative to the TSS strand, and IP is
normalized over input with one denominator per window, so that sense and
antisense stay on one scale and add up to the unstranded map.

Main Functions
--------------
read_tss_calls : one single-nucleotide TSS per call, from one or more files
build_windows : fixed-size, same-strand non-overlapping windows around each TSS
shift_windows : the same windows moved downstream (a null without a TSS)
stranded_matrices : sense and antisense density per window, 5' to 3'
normalize : input-subtracted, per-window normalized sense/antisense/total
profiles : mean of each normalized matrix per position
summarize : sense and antisense sums in windows centered on the TSS
"""

import os

import numpy as np
import pandas as pd

ORIENTATIONS = ('sense', 'antisense', 'total')
OPPOSITE = {'+': '-', '-': '+'}


def read_tss_calls(filenames):
    """
    Reads TSS calls and reduces each to a single nucleotide.

    Files with 11 or more columns are read as ENCODE CAGE/RAMPAGE peaks, whose
    11th column holds the comma-separated signal at each base of the peak: the
    TSS is the leftmost base with the highest signal. Any other file is read
    as BED6: the TSS is the 5' end of the interval and the score is its height.

    Parameters
    ----------
    filenames : list
        BED files (optionally gzipped)

    Returns
    -------
    calls : pandas.DataFrame
        chrom, x (0-based TSS position), strand, h (signal at the TSS),
        c (mean signal of the call), file
    """
    rows = []
    for filename in filenames:
        df = pd.read_csv(filename, sep='\t', header=None, comment='#')
        name = os.path.basename(filename)
        for r in df.itertuples(index=False):
            chrom, start, end, strand = r[0], int(r[1]), int(r[2]), r[5]
            if len(r) >= 11:
                signal = np.array([float(v) for v in str(r[10]).split(',') if v != ''])
                if signal.size == 0:
                    continue
                rows.append((chrom, start + int(np.argmax(signal)), strand,
                             float(signal.max()), float(signal.mean()), name))
            else:
                x = start if strand == '+' else end - 1
                rows.append((chrom, x, strand, float(r[4]), float(r[4]), name))
    return pd.DataFrame(rows, columns=['chrom', 'x', 'strand', 'h', 'c', 'file'])


def build_windows(calls, slop, chrom_sizes):
    """
    Builds one window of 2 * slop + 1 nt around each TSS.

    Windows on the same strand that overlap or touch are clustered (as
    bedtools merge does) and only the window of the best call is kept:
    highest h, ties broken by highest c. Windows on opposite strands are not
    compared; overlaps_opposite marks the windows that overlap one.

    Parameters
    ----------
    calls : pandas.DataFrame
        see: read_tss_calls()
    slop : int
        bases on each side of the TSS
    chrom_sizes : dict
        {chrom: size}; calls on other chromosomes, and windows that run off
        the end of a chromosome, are dropped.

    Returns
    -------
    windows : pandas.DataFrame
        the kept calls with start, end, n_files (number of files with a call
        in the cluster), overlaps_opposite and region_id
    """
    width = 2 * slop + 1
    calls = calls[calls.chrom.isin(chrom_sizes)].reset_index(drop=True)
    calls['start'] = (calls.x - slop).clip(lower=0)
    calls['end'] = np.minimum(calls.x + slop + 1, calls.chrom.map(chrom_sizes))
    calls = calls[calls.end - calls.start == width].reset_index(drop=True)

    keep, n_files = [], []
    for strand in ['+', '-']:
        sub = calls[calls.strand == strand].sort_values(['chrom', 'start'], kind='mergesort')
        new_cluster = (sub.chrom != sub.chrom.shift()) | \
                      (sub.start > sub.groupby('chrom').end.cummax().shift())
        for _, cluster in sub.groupby(new_cluster.cumsum(), sort=False):
            keep.append(max(cluster.index, key=lambda j: (calls.at[j, 'h'], calls.at[j, 'c'])))
            n_files.append(cluster.file.nunique())
    windows = calls.loc[keep].reset_index(drop=True)
    windows['n_files'] = n_files

    overlaps = np.zeros(len(windows), dtype=bool)
    for _, chrom in windows.groupby('chrom'):
        for strand in ['+', '-']:
            own = chrom[chrom.strand == strand]
            other = np.sort(chrom[chrom.strand == OPPOSITE[strand]].x.values)
            # two windows overlap when their TSS are at most 2 * slop apart
            overlaps[own.index] = np.searchsorted(other, own.x.values + 2 * slop, 'right') > \
                                  np.searchsorted(other, own.x.values - 2 * slop, 'left')
    windows['overlaps_opposite'] = overlaps
    windows['region_id'] = windows.chrom + ':' + windows.start.astype(str) + '-' + \
                           windows.end.astype(str) + ':' + windows.strand
    return windows


def shift_windows(windows, shift, chrom_sizes):
    """
    Moves every window downstream (in the direction of transcription).
    Windows that leave the chromosome are dropped.

    Parameters
    ----------
    windows : pandas.DataFrame
        see: build_windows()
    shift : int
    chrom_sizes : dict

    Returns
    -------
    shifted : pandas.DataFrame
    """
    shifted = windows.copy()
    offset = np.where(shifted.strand == '+', shift, -shift)
    for col in ['x', 'start', 'end']:
        shifted[col] = shifted[col] + offset
    in_bounds = (shifted.start >= 0) & (shifted.end <= shifted.chrom.map(chrom_sizes))
    shifted = shifted[in_bounds].reset_index(drop=True)
    shifted['region_id'] = shifted.region_id + '_shift{}'.format(shift)
    return shifted


def stranded_matrices(density, windows):
    """
    Reads the density of every window on its own strand (sense) and on the
    opposite strand (antisense). Both are oriented 5' to 3' relative to the
    TSS, positive, and zero where the track has no value.

    Parameters
    ----------
    density : density.ReadDensity
    windows : pandas.DataFrame
        see: build_windows()

    Returns
    -------
    sense : numpy.ndarray
        windows x positions
    antisense : numpy.ndarray
    """
    def values(chrom, start, end, strand):
        return np.abs(np.nan_to_num(np.array(density.values(chrom, start, end, strand), dtype=float)))

    width = int(windows.end.iloc[0] - windows.start.iloc[0]) if len(windows) else 0
    sense = np.zeros((len(windows), width))
    antisense = np.zeros((len(windows), width))
    for i, r in enumerate(windows.itertuples(index=False)):
        sense[i] = values(r.chrom, r.start, r.end, r.strand)
        # values() returns the other strand in its own 5' to 3' direction
        antisense[i] = values(r.chrom, r.start, r.end, OPPOSITE[r.strand])[::-1]
    return sense, antisense


def normalize(ip, inp, pseudocount):
    """
    Subtracts input from IP for each orientation, then divides both
    orientations of a window by the same number: the sum of the absolute
    unstranded (sense + antisense) difference plus one pseudocount per
    position. This is per_region_subtract_and_normalize() of the unstranded
    signal, so sense + antisense equals the unstranded map.

    Parameters
    ----------
    ip : tuple
        (sense, antisense) matrices of the IP, see: stranded_matrices()
    inp : tuple
        (sense, antisense) matrices of the input
    pseudocount : float
        density of a single IP read

    Returns
    -------
    normalized : dict
        {'sense', 'antisense', 'total'} matrices
    """
    sense = ip[0] - inp[0]
    antisense = ip[1] - inp[1]
    total = sense + antisense
    denominator = np.abs(total).sum(axis=1) + pseudocount * total.shape[1]
    denominator = denominator[:, None]
    return {'sense': sense / denominator, 'antisense': antisense / denominator,
            'total': total / denominator}


def profiles(normalized, keep=None):
    """
    Mean normalized signal per position, indexed by position relative to the TSS.

    Parameters
    ----------
    normalized : dict
        see: normalize()
    keep : numpy.ndarray
        boolean mask of the windows to average (default: all)

    Returns
    -------
    profiles : pandas.DataFrame
        one column per orientation
    """
    slop = normalized['total'].shape[1] // 2
    index = pd.Index(np.arange(-slop, slop + 1), name='position')
    return pd.DataFrame({
        o: (normalized[o] if keep is None else normalized[o][keep]).mean(axis=0)
        for o in ORIENTATIONS
    }, index=index)


def summarize(profile, halves=(50, 100, 250)):
    """
    Sums the sense and antisense profiles over windows centered on the TSS.

    Parameters
    ----------
    profile : pandas.DataFrame
        see: profiles()
    halves : tuple
        half-widths of the windows to report (the full window is always added)

    Returns
    -------
    summary : pandas.DataFrame
        window, sense, antisense and sense_frac = sense / (sense + antisense),
        which is given only when both sums are positive.
    """
    slop = int(profile.index.max())
    rows = []
    for half in sorted(set(h for h in halves if h < slop) | {slop}):
        sense = profile.sense.loc[-half:half].sum()
        antisense = profile.antisense.loc[-half:half].sum()
        both_positive = sense > 0 and antisense > 0
        rows.append({
            'window': '+/-{}'.format(half), 'sense': sense, 'antisense': antisense,
            'sense_frac': sense / (sense + antisense) if both_positive else np.nan,
        })
    return pd.DataFrame(rows)
