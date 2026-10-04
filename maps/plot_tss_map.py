#!/usr/bin/env python
# encoding: utf-8
"""CLI entrypoint for stranded (sense/antisense) RBP maps around TSS calls."""

import argparse
import os

import pandas as pd

import density.ReadDensity
from density import tss
from maps import plot_map
from plotter import tss as tss_plotter


def run_make_tss_map(
        outfile, ip, inp, tss_files, slop=1000, min_files=1, min_height=0,
        exclude_opposite_overlaps=False, shift=10000
):
    """Build, normalize and plot a stranded TSS map.

    Args:
        outfile (str): Output image filename. Tables are written next to it.
        ip (density.ReadDensity.ReadDensity): IP read density.
        inp (density.ReadDensity.ReadDensity): Input read density.
        tss_files (list[str]): TSS call files (see ``density.tss.read_tss_calls``).
        slop (int): Bases plotted on each side of the TSS.
        min_files (int): Keep windows with calls from at least this many files.
        min_height (float): Keep windows whose TSS call is at least this high.
        exclude_opposite_overlaps (bool): Drop windows that overlap a window
            on the opposite strand.
        shift (int): Distance downstream of the shifted control (0: no control).

    Returns:
        pandas.DataFrame: sense/antisense sums around the TSS (and the control).
    """
    sizes = dict(zip(ip.bam.references, ip.bam.lengths))
    windows = tss.build_windows(tss.read_tss_calls(tss_files), slop, sizes)
    windows['kept'] = (windows.n_files >= min_files) & (windows.h >= min_height)
    if exclude_opposite_overlaps:
        windows['kept'] &= ~windows.overlaps_opposite
    kept = windows[windows.kept].reset_index(drop=True)
    if len(kept) == 0:
        raise ValueError("No TSS windows are left after filtering.")

    def profile(regions):
        normalized = tss.normalize(
            tss.stranded_matrices(ip, regions), tss.stranded_matrices(inp, regions),
            ip.pseudocount()
        )
        return tss.profiles(normalized)

    tss_profile = profile(kept)
    table = tss_profile.copy()
    summary = tss.summarize(tss_profile).assign(regions='TSS')
    control, control_label = None, None
    if shift:
        control_label = 'shifted {} nt'.format(shift)
        control = profile(tss.shift_windows(kept, shift, sizes))
        for orientation in tss.ORIENTATIONS:
            table['shifted_' + orientation] = control[orientation]
        summary = pd.concat(
            [summary, tss.summarize(control).assign(regions=control_label)], ignore_index=True
        )

    base = os.path.splitext(outfile)[0]
    windows.to_csv(base + '.windows.tsv', sep='\t', index=False)
    table.to_csv(base + '.profiles.csv')
    summary.to_csv(base + '.summary.csv', index=False)
    tss_plotter.plot_stranded_metagene(
        tss_profile, outfile,
        'Sense vs antisense around TSS (+/-{} nt, n={})'.format(slop, len(kept)),
        control, control_label
    )
    return summary


def main():
    """Parse CLI arguments and make a stranded TSS map."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ip", "--ipbam", dest="ipbam", required=True, help="IP BAM file")
    parser.add_argument("--input", "--inputbam", dest="inputbam", required=True, help="INPUT BAM file")
    parser.add_argument("--ip_pos_bw", default=None, help="IP positive bigwig file")
    parser.add_argument("--ip_neg_bw", default=None, help="IP negative bigwig file")
    parser.add_argument("--input_pos_bw", default=None, help="INPUT positive bigwig file")
    parser.add_argument("--input_neg_bw", default=None, help="INPUT negative bigwig file")
    parser.add_argument("--output", required=True, help="output figure (png, svg or pdf)")
    parser.add_argument(
        "--tss", nargs='+', required=True,
        help="TSS calls: ENCODE CAGE/RAMPAGE peak files (the 11th column holds "
             "the per-base signal) or BED6 files (the TSS is the 5' end)."
    )
    parser.add_argument("--slop", type=int, default=1000,
                        help="bases plotted on each side of the TSS (default: 1000)")
    parser.add_argument("--min_files", type=int, default=1,
                        help="keep TSS with calls from at least this many --tss files (default: 1)")
    parser.add_argument("--min_height", type=float, default=0,
                        help="keep TSS whose call is at least this high (default: 0)")
    parser.add_argument(
        "--exclude_opposite_overlaps", default=False, action='store_true',
        help="drop windows that overlap a window on the opposite strand "
             "(for example bidirectional promoters)"
    )
    parser.add_argument(
        "--shift", type=int, default=10000,
        help="also plot the same windows shifted this far downstream, as a "
             "control without a TSS (default: 10000; 0 to skip)"
    )
    parser.add_argument(
        "--genome", default=None,
        help="Tab-separated chrom.sizes file, needed when bigWigs have to be generated from BAM."
    )
    parser.add_argument(
        "--makebigwigfiles_direction", "--make_bigwig_files_direction", default=None,
        help="Read direction for BAM->signal conversion ('r' for reverse-stranded, 'f' for forward-stranded)."
    )
    parser.add_argument(
        "--generated_signal_dir", default=None,
        help="Directory for auto-generated .norm.pos.bw/.norm.neg.bw files (default: next to each BAM)."
    )
    parser.add_argument(
        "--makebigwigfiles_workdir", "--make_bigwig_files_workdir", default=None,
        help="Directory for the bedGraph files written while generating bigWigs."
    )
    args = parser.parse_args()

    ip_bam, input_bam, ip_pos_bw, ip_neg_bw, input_pos_bw, input_neg_bw = \
        plot_map.resolve_density_input_paths(args)
    bigwigs = {ip_bam: (ip_pos_bw, ip_neg_bw), input_bam: (input_pos_bw, input_neg_bw)}
    missing = [bw for pair in bigwigs.values() for bw in pair if not os.path.isfile(bw)]
    if missing and args.genome is None:
        parser.error(
            "Missing bigWigs detected ({}) and --genome was not provided.".format(', '.join(missing))
        )
    densities = []
    for bam, (pos_bw, neg_bw) in bigwigs.items():
        plot_map.check_for_index(bam)
        plot_map.ensure_density_bigwigs(
            bam=bam, pos_bw=pos_bw, neg_bw=neg_bw, genome_file=args.genome,
            direction=args.makebigwigfiles_direction, makebigwigfiles_cmd=None,
            makebigwigfiles_extra_args='', makebigwigfiles_workdir=args.makebigwigfiles_workdir
        )
        densities.append(density.ReadDensity.ReadDensity(pos=pos_bw, neg=neg_bw, bam=bam))

    summary = run_make_tss_map(
        args.output, densities[0], densities[1], args.tss, slop=args.slop,
        min_files=args.min_files, min_height=args.min_height,
        exclude_opposite_overlaps=args.exclude_opposite_overlaps, shift=args.shift
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
