#!/usr/local/bin/python2.7
# encoding: utf-8
'''
     up_ex       ex_up     ex_dn       dn_ex
====[=----]-----[----=]===[=----]-----[----=]====

@author:     user_name

@copyright:  2015 organization_name. All rights reserved.

@license:    license

@contact:    user_email
@deffield    updated: Updated
'''

import sys
import os
# from cluster import map_peaks
# from cluster import map_density
# from cluster import annotations
import map_peaks
import annotations
import clust

from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
from __builtin__ import True

__all__ = []
__version__ = 0.1
__date__ = '2015-12-19'
__updated__ = '2015-12-19'

DEBUG = 0
TESTRUN = 1
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg


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
    
    try:
        # Setup argument parser
        parser = ArgumentParser(description=program_license, formatter_class=RawDescriptionHelpFormatter)
        
        parser.add_argument("-o", "--output", dest="output", help="output folder", required = True )
        parser.add_argument("-i", "--input", dest="input", help="input manifest file (list of BED/BEDgraph files)", required = True )
        parser.add_argument("-i2", "--input2", dest="input2", help="input manifest file (list of BED/BEDgraph files)", required = False )
        parser.add_argument("-m", "--miso", dest="miso", help="miso annotation file", required = True )
        parser.add_argument("-c", "--celltype", dest="celltype", help="[HepG2] or [K562] or [all (Default)]", default="ALL", required = False )
        parser.add_argument('-V', "--version", action='version', version=program_version_message)
        parser.add_argument('-f', "--foldchange", dest="fc", help="log2 fold change cutoff (default = 0)", type=float, default=0, required = False)
        parser.add_argument('-p', "--pvalue", dest="pv", help= "-log10(p-value) cutoff (default = 0.05)", type=float, default=0.05, required = False)
        parser.add_argument('-s', "--hashval", dest="hash", help = "hash value (default = 100000)", type=int, default=100000, required = False)
        parser.add_argument('-eo', "--exonoverhang", dest="exonoverhang", help = "exon offset overhang (default = 50)", type=int, default=50, required = False)
        parser.add_argument('-io', "--intronoverhang", dest="intronoverhang", help="intron offset overhange (default = 500)", type=int, default=500, required = False)
        parser.add_argument('-t', "--filetype", dest="filetype", help="file type is either (b)ed or bed(g)raph (default 'b')", default="BED", required=False)
        parser.add_argument('-map', "--profile", dest="maptype", help="type of rbp map. Either: se, a3ss, a5ss (default se)", default="se", required=False)
        args = parser.parse_args()
        
        infiles = args.input
        outfolder = args.output
        
        names = [] # pretty names for clustering purposes
        outfiles = [] # array of distributions 
        
        miso = args.miso
        l2fc_cutoff = args.fc
        l10p_cutoff = args.pv
        exon_overhang = args.exonoverhang
        intron_overhang = args.intronoverhang
        hashing_val = args.hash
        celltype = args.celltype
        filetype = args.filetype
        maptype = args.maptype.lower()
        
        print("Begin building the dictionary.")
        if maptype == 'se' or maptype == 'a3ss' or maptype == 'a5ss':
            all_exons = annotations.read_four_region_miso(miso, 
                                                 hashing_val,
                                                 maptype, 
                                                 exon_overhang, 
                                                 intron_overhang)
            
        print("End dictionary build.")
        fin = open(infiles,'r')
        lines = fin.readlines()
        # generate individual maps for each line in manifest
        for line in lines:
            line = line.strip().split('\t')
            
            """
                    
            name = line[0]
            fi = line[1]
            print("output: {}".format(os.path.splitext(fi)[0]))
            outfile = outfolder+"/{}.pv_{}.fc_{}.txt".format(os.path.splitext(fi)[0],
                                                       l10p_cutoff,
                                                       l2fc_cutoff)
            names.append(name)
            outfiles.append(outfile)
            """
            code = line[0]
            name = line[1]
            cell = line[2]
            rep1bam = line[3]
            rep2bam = line[4]
            inputbam = line[5]
            rep1bed = line[6]
            rep2bed = line[7].strip()
            
            files = [rep1bed, rep2bed]
            if ((celltype.upper() == "ALL") or (celltype.upper() == cell.upper())):
                for i in range(0,len(files)):
                    print("Processing: {0}_0{1}-{2}".format(name,i+1,celltype))
                    outfile = outfolder+"/{}.pv_{}.fc_{}.txt".format(os.path.splitext(os.path.basename(files[i]))[0],
                                                           l10p_cutoff,
                                                           l2fc_cutoff)
                
                    if maptype == 'se':
                        map_peaks.make_hist_se(files[i],outfile,hashing_val,l10p_cutoff,l2fc_cutoff,all_exons,exon_overhang,intron_overhang)
                        map_peaks.rbp_plot_se(outfile,"{}.se.jpg".format(os.path.splitext(outfile)[0]),erl = exon_overhang, irl = intron_overhang)
                    elif maptype == 'a3ss':
                        map_peaks.make_hist_a3ss(files[i],outfile,hashing_val,l10p_cutoff,l2fc_cutoff,all_exons,exon_overhang,intron_overhang)
                        map_peaks.rbp_plot_a3ss(outfile,"{}.a3ss.jpg".format(os.path.splitext(outfile)[0]),erl = exon_overhang, irl = intron_overhang)
                    elif maptype == 'a5ss':
                        map_peaks.make_hist_a5ss(files[i],outfile,hashing_val,l10p_cutoff,l2fc_cutoff,all_exons,exon_overhang,intron_overhang)
                        map_peaks.rbp_plot_a5ss(outfile,"{}.a5ss.jpg".format(os.path.splitext(outfile)[0]),erl = exon_overhang, irl = intron_overhang)
                    print("makeing names")
                    names.append("{0}_0{1}-{2}-{3}".format(name,i+1,celltype,maptype))
                    outfiles.append(outfile)
            
        # cluster stuff
        """
        # generate individual maps for each line in manifest
        for line in lines:
            
            line = line.strip().split('\t')
            
            
            code = line[0]
            name = line[1]
            cell = line[2]
            rep1bam = line[3]
            rep2bam = line[4]
            inputbam = line[5]
            rep1bed = line[6]
            rep2bed = line[7].strip()
            
            files = [rep1bed, rep2bed]
            if ((celltype.upper() == "ALL") or (celltype.upper() == cell.upper())):
                for i in range(0,len(files)):
                    print("Processing: {0}_0{1}-{2}".format(name,i+1,celltype))
                    if(filetype == 'b'):
                        outfile = outfolder+"/{}.pv_{}.fc_{}.txt".format(os.path.splitext(os.path.basename(files[i]))[0],
                                                           l10p_cutoff,
                                                           l2fc_cutoff)
                        map_peaks.make_hist(files[i], 
                                  outfile, 
                                  hashing_val, 
                                  l10p_cutoff, 
                                  l2fc_cutoff, 
                                  all_exons, 
                                  exon_overhang, 
                                  intron_overhang)
                    
                        map_peaks.rbp_plot(outfile, 
                                "{}.jpg".format(os.path.splitext(outfile)[0]), 
                                erl = exon_overhang, 
                                irl = intron_overhang)
                    
                    else:
                        outfile = outfolder+"/{}.txt".format(os.path.splitext(os.path.basename(files[i]))[0])
                        outfile2 = outfolder+"/{}.2.txt".format(os.path.splitext(os.path.basename(files[i]))[0])
                        map_density.make_density(infile2, 
                                     infile, 
                                     outfile, 
                                     outfile2, 
                                     hashing_val, 
                                     all_exons, 
                                     exon_overhang, 
                                     intron_overhang)
        
        
                        map_density.concat_and_save(outfile,outfile2,"U2AF2", "output_concatenated.txt")
                        rbp_plot(outfile,
                                 outfile2, 
                                 "{}.jpg".format(os.path.splitext(outfile)[0]), 
                                 erl = exon_overhang, 
                                 irl = intron_overhang)
                    
                    names.append("{0}_0{1}-{2}".format(name,i+1,celltype))
                    outfiles.append(outfile)
         """   
        # cluster stuff
        import pandas as pd
        na = pd.Series(names)
        ot = pd.Series(outfiles)
        na.to_csv(outfolder+"/names.txt")
        ot.to_csv(outfolder+"/outfiles.txt")
        clust.heatmap(names, outfiles, outfolder+"/{}_heatmap.png".format(maptype))
        clust.get_peak_counts_table(names, outfiles, outfolder+"/counts.txt")
        return 0
    except KeyboardInterrupt:
        ### handle keyboard interrupt ###
        return 0
    except Exception, e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2

if __name__ == "__main__":
    if DEBUG:
        sys.argv.append("-h")
        sys.argv.append("-v")
        sys.argv.append("-r")
    if TESTRUN:
        import doctest
        doctest.testmod()
    if PROFILE:
        import cProfile
        import pstats
        profile_filename = 'cluster.overlap_peak_with_annot_profile.txt'
        cProfile.run('main()', profile_filename)
        statsfile = open("profile_stats.txt", "wb")
        p = pstats.Stats(profile_filename, stream=statsfile)
        stats = p.strip_dirs().sort_stats('cumulative')
        stats.print_stats()
        statsfile.close()
        sys.exit(0)
    sys.exit(main())