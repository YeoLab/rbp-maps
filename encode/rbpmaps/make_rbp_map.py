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
    parser.add_argument("-m", "--min", dest="minthreshold", help="minimum density read threshold", default=0, type = int)
    parser.add_argument("-t", "--title", dest="title", help="plot title", default=None)
    parser.add_argument("-ty", "--type", dest="type", help="--type [se, txstarts, txends, cdsstarts, cdsstops] or NONE if plotting something generic.")
    
    args = parser.parse_args()
    outfile = args.output
    positive_bw = args.positive
    negative_bw = args.negative
    annotation = args.annotation
    min_read_threshold = args.minthreshold
    
    col = sns.color_palette("hls", 8)[args.color]
    lab = args.label
    mytitle = args.title
    
    if args.flipped == True:
        print("flipped pos={}, neg={}.".format(negative_bw,positive_bw))
        rbp = ReadDensity.ReadDensity(pos=negative_bw,
                                      neg=positive_bw)
    else:
        rbp = ReadDensity.ReadDensity(pos=positive_bw,
                                      neg=negative_bw)
    
    if(args.type == 'se'):
        plot.plot_se(rbp, 
                     annotation, 
                     outfile, 
                     exon_offset = args.exonoffset, 
                     intron_offset = args.intronoffset, 
                     title = mytitle, 
                     color = col)
    elif(args.type == 'txstarts'):
        annotations = bt.BedTool(annotation)
        plot.plot_txstarts(rbp, 'txstarts', 
                           outfile, col, lab,
                           args.left, args.right)
    elif(args.type == 'txends'):
        annotations = bt.BedTool(annotation)
        plot.plot_txends(rbp, 'txends', 
                           outfile, col, lab,
                           args.left, args.right)
    elif(args.type == 'cdsstarts'):
        annotations = bt.BedTool(annotation)
        plot.plot_cdsstarts(rbp, 'cdsstarts', 
                           outfile, col, lab,
                           args.left, args.right)
    elif(args.type == 'cdsends'):
        annotations = bt.BedTool(annotation)
        plot.plot_cdsends(rbp, 'cdsends', 
                           outfile, col, lab,
                           args.left, args.right)
    else:
        annotations = bt.BedTool(annotation)
        plot.plot_single_frame(rbp,
                               annotations,
                               outfile,
                               color = col,
                               label = lab,
                               left = args.left,
                               right = args.right,
                               distribution = args.dist,
                               title = mytitle,
                               min_read_density_sum = min_read_threshold)
if __name__ == "__main__":
    main()