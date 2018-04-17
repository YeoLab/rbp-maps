#!/usr/bin/env python
# encoding: utf-8
'''

Converts a list of BED files into bigbed files

@author:     brian

@copyright:  2016 organization_name. All rights reserved.

@license:    license

@contact:    bay001@ucsd.edu
@deffield    updated: 2018-01-17
'''

import argparse
import pybedtools
import subprocess
import os
import pandas as pd

def stringify(name):
    try:
        _ = float(name)
        return '_' + name + '_'
    except ValueError:
        return name

def filter_bed(in_bed, log10p, log2fc, filtered_bed):
    df = pd.read_table(in_bed, names=[
        'chrom','start','end','log10p','log2fc','strand'
    ])
    df = df[(df['log10p']>=log10p) & (df['log2fc']>=log2fc)]
    df.to_csv(filtered_bed, sep='\t', index=False, header=False)

def convert_to_bigbed(
        in_bed, genome, bed_type, out_bb, log10p, log2fc,
        blockSize=256, itemsPerSlot=512
):
    sorted_bed = os.path.join(
        os.path.dirname(out_bb),
        os.path.basename(in_bed) + '.sorted.bed'
    )

    # filter the bed file if necessary
    if bed_type == 'bed6inputnorm' and log10p is not None and log2fc is not None:
        filtered_bed = sorted_bed + '.pv{}fc{}.bed'.format(
            log10p, log2fc
        )
        filter_bed(sorted_bed, log10p, log2fc, filtered_bed)
        bed_type = 'bed6'
    else:
        filtered_bed = sorted_bed

    # must convert names (could be floats) to strings, otherwise bigbed will complain
    bedtool = pybedtools.BedTool(filtered_bed)
    out_bedtool = []
    for interval in bedtool:  # need to make sure that the name is a string
        interval.name = stringify(interval.name)
        interval.score = int(float(interval.score))
        out_bedtool.append(interval)
    out_bedtool = pybedtools.BedTool(out_bedtool)
    # sort the bed file otherwise bigbed will complain
    out_bedtool = out_bedtool.sort()
    out_bedtool.saveas(filtered_bed)


    cmds = [
        'bedToBigBed',
        filtered_bed,
        genome,
        out_bb,
        '-blockSize={}'.format(blockSize),
        '-itemsPerSlot={}'.format(itemsPerSlot),
        '-type={}'.format(bed_type)
    ]
    p = subprocess.Popen(cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    if p.returncode:
        raise ValueError("cmds: %s\nstderr:%s\nstdout:%s"
                         % (" ".join(cmds), stderr, stdout))

    return 0

def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--beds",
        required=True,
        nargs='+',
        help='bed or list of bed files to convert'
    )
    parser.add_argument(
        "--genome",
        required=True,
        help='genome sizes (chrom.sizes) file'
    )
    parser.add_argument(
        "--outbbs",
        help="output bigbed files (default: ${beds}.bb)",
        default=None,
        nargs='+'
    )
    parser.add_argument(
        "--bedtype",
        default='bed6inputnorm',
        help="BEDtype (ie. bed6). Default: bed6inputnorm (will only perform "
             "log10p/log2fc filtering if this is an input normed bed6file)",
    )
    parser.add_argument(
        "--log10p",
        help="-log10p values to filter (column 4 in inputnorm bed). Default: 3",
        default=3,
        type=float,
    )
    parser.add_argument(
        "--log2fc",
        help="log2 fold change values to filter (column 5 in inputnorm bed). "
             "Default: 3",
        default=3,
        type=float,
    )

    # Process arguments
    args = parser.parse_args()
    beds = args.beds
    genome = args.genome
    outbbs = args.outbbs if args.outbbs is not None else [bed + '.bb' for bed in beds]
    bed_type = args.bedtype

    log10p = args.log10p
    log2fc = args.log2fc
    for i in range(len(beds)):
        convert_to_bigbed(beds[i], genome, bed_type, outbbs[i], log10p, log2fc)
if __name__ == "__main__":
    main()
