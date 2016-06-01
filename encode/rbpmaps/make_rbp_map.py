'''
Created on May 3, 2016

@author: brianyee
'''

import ReadDensity
import plot
import pybedtools as bt
import seaborn as sns
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import sys
import os
# from IPython.core.display import HTML

__all__ = []
__version__ = 0.1
__date__ = '2016-5-5'
__updated__ = '2016-5-5'


def main(argv=None): # IGNORE:C0111
    '''Command line options.'''

    if argv is None:
        argv = sys.argv
    else:
        sys.argv.extend(argv)

    program_name = os.path.basename(sys.argv[0])
    program_version = "v%s" % __version__
    program_build_date = str(__updated__)
    program_version_message = '%%(prog)s %s (%s)' % (program_version, program_build_date)
    program_shortdesc = __import__('__main__').__doc__.split("\n")[1]
    program_license = '''%s

  Created by user_name on %s.
  Copyright 2016 organization_name. All rights reserved.

  Licensed under the Apache License 2.0
  http://www.apache.org/licenses/LICENSE-2.0

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE
''' % (program_shortdesc, str(__date__))

    parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        
    parser.add_argument("-o", "--output", dest="output", help="output file", required = True )
    parser.add_argument("-p", "--positive", dest="positive", help="positive *.bw file", required = True )
    parser.add_argument("-n", "--negative", dest="negative", help="negative *.bw file", required = True )
    parser.add_argument("-a", "--annotation", dest="annotation", help="bedfile or miso file containing region of interest", required = True )
    parser.add_argument("-l", "--left", dest="left", help="left margins. For a given single region, how many nt upstream of the pointsource/region boundary should we extend?", required = False, default = 300, type = int)
    parser.add_argument("-r", "--right", dest="right", help="right margins. For a given single region, how many nt downstream of the pointsource/region boundary should we extend?", required = False, default = 300, type = int)
    parser.add_argument("-e", "--exon_offset", dest="exonoffset", help="exon offset margins. For a given region, how much should we extend into exon feature?", required = False, default = 50, type = int)
    parser.add_argument("-i", "--intron_offset", dest="intronoffset", help="intron offset margins. For a given region, how much should we extend outside the exon feature?", required = False, default = 300, type = int)
    parser.add_argument("-c", "--color", dest="color", help="line color (integer from 0 to 7)", required = False, default = 4, type = int)
    parser.add_argument("-lbl", "--label", dest="label", help="label or feature", required = False, default = "feature")
    parser.add_argument("-d", "--dist", dest="dist", help="specifiy this flag if trying to plot regions of varying length", action='store_true')
    parser.add_argument("-nu", "--nucl", dest="dist", help="if regions are of same length, we can plot nucleotide resolution (conflicts with -d)", action='store_false')
    parser.add_argument("-f", "--flipped", dest="flipped", help="if positive is negative (pos.bw really means neg.bw)", action='store_true')
    parser.add_argument("-m", "--min", dest="minthreshold", help="minimum density read threshold", default=0, type = float)
    parser.add_argument("-t", "--title", dest="title", help="plot title", default=None)
    parser.add_argument("-ty", "--type", dest="type", help="--type [se, txstarts, txends, cdsstarts, cdsstops] or NONE if plotting something generic.")
    parser.add_argument("-csv", "--csv", dest="csv", help="If true, also output intermediate csv files", action='store_true', default = False)
    
    args = parser.parse_args()
 
    
    if args.flipped == True:
        print("flipped pos={}, neg={}.".format(args.negative,args.positive))
        rbp = ReadDensity.ReadDensity(pos=args.negative,
                                      neg=args.positive)
    else:
        rbp = ReadDensity.ReadDensity(pos=args.positive,
                                      neg=args.negative)
    
    if(args.type == 'se'):
        plot.plot_se(rbp, 
                     args.annotation, 
                     args.output, 
                     exon_offset = args.exonoffset, 
                     intron_offset = args.intronoffset, 
                     title = args.title, 
                     color = sns.color_palette("hls", 8)[args.color],
                     min_density_threshold = args.minthreshold,
                     csv = args.csv)
    elif(args.type == 'txstarts'):
        annotations = bt.BedTool(args.annotation)
        plot.plot_txstarts(rbp, 'txstarts', 
                           args.output, 
                           sns.color_palette("hls", 8)[args.color], 
                           args.label,
                           args.left, args.right,
                           csv = args.csv)
    elif(args.type == 'txends'):
        annotations = bt.BedTool(args.annotation)
        plot.plot_txends(rbp, 'txends', 
                           args.output, 
                           sns.color_palette("hls", 8)[args.color], 
                           args.label,
                           args.left, args.right,
                           csv = args.csv)
    elif(args.type == 'cdsstarts'):
        annotations = bt.BedTool(args.annotation)
        plot.plot_cdsstarts(rbp, 'cdsstarts', 
                           args.output, 
                           sns.color_palette("hls", 8)[args.color], 
                           args.label,
                           args.left, args.right,
                           csv = args.csv)
    elif(args.type == 'cdsends'):
        annotations = bt.BedTool(args.annotation)
        plot.plot_cdsends(rbp, 'cdsends', 
                           args.output, 
                           sns.color_palette("hls", 8)[args.color],
                           args.label,
                           args.left, args.right,
                           csv = args.csv)
    else:
        annotations = bt.BedTool(args.annotation)
        plot.plot_single_frame(rbp,
                               annotations,
                               args.output,
                               color = sns.color_palette("hls", 8)[args.color],
                               label = args.label,
                               left = args.left,
                               right = args.right,
                               distribution = args.dist,
                               title = args.title,
                               min_read_density_sum = args.minthreshold,
                               csv = args.csv)
if __name__ == "__main__":
    main()