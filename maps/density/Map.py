#!/usr/bin/env python

"""
Created on Jun 27, 2016

Modules that help containerize the rbp-map information, including:

ip : CLIP IP density
input : CLIP Input density
output filename : name of the output file and base
norm function : normalization to use (see normalization_functions.py)
annotation : annotation file(s) to use and their type. This is a dictionary of 
  {filename:filetype(rmats, miso, etc.)} where the filename is used to 
  associated which density values belong to which condition.
offsets : either describing upstream/downstream or exon/intron bases to plot
is_scaled : whether or not we want the map to 'scale' or keep nucleotide dims
conf : percentage of scores to keep (outlier removal).

@author: brianyee
"""
import matplotlib
import matplotlib.patches as patches
matplotlib.use('Agg')
from matplotlib import rc
import matplotlib.gridspec as gridspec
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})
import seaborn as sns
import os
import pandas as pd
import matplotlib.pyplot as plt
import gzip
import colors
from collections import defaultdict, OrderedDict
import RDPlotter
import matrix as mtx
import misc
import normalization_functions as norm
import LineObject

SEP = '.'  # file delimiter

MIN_EVENT_THRESHOLD=100 # number of events required to not grey the line out

COLOR_PALETTE = sns.color_palette("hls", 8)

BG1_COLOR = 'black' # COLOR_PALETTE['black']
BG2_COLOR = COLOR_PALETTE[6]
BG3_COLOR = COLOR_PALETTE[7]
BG4_COLOR = COLOR_PALETTE[4]
POS_COLOR = COLOR_PALETTE[0]
NEG_COLOR = COLOR_PALETTE[5]

COLORS = [POS_COLOR, NEG_COLOR, BG1_COLOR, BG2_COLOR, BG3_COLOR, BG4_COLOR]

RED = COLOR_PALETTE[0]
ORANGE = COLOR_PALETTE[1]
BLUE = COLOR_PALETTE[5]
GREEN = COLOR_PALETTE[3]

class Map:
    def __init__(self, ip, output_filename, norm_function,
                 annotation, upstream_offset=0, downstream_offset=0,
                 min_density_threshold=0, is_scaled=False, conf=0.95):
        """

        Parameters
        ----------
        ip : density.ReadDensity
            ReadDensity file containing RBP density signals
        output_filename : str
            Output filename for map (ie. plot.svg).
        norm_function : function
            Function
            See: normalization_functions
        annotation : collections.OrderedDict
            {file:type} representing an annotation file and format
            See: Feature
        upstream_offset : int
            number of bases upstream of described feature to plot
        downstream_offset : int
            number of bases downstream of described feature to plot
        min_density_threshold : int
            Minimum density sum across an event to use
        is_scaled : Boolean
            True if we will re-scale features of varied lengths.
        conf : float
            Representing percentage of events to keep (default .95)
        """
        self.ip = ip
        self.output_filename = output_filename
        self.output_base = os.path.splitext(output_filename)[0]
        self.norm_function = norm_function
        self.annotation = annotation
        self.min_density_threshold = min_density_threshold
        self.upstream_offset = upstream_offset
        self.downstream_offset = downstream_offset
        self.is_scaled = is_scaled
        self.conf = conf

        self.raw_matrices = defaultdict(dict)
        self.norm_matrices = defaultdict(dict)
        self.means = defaultdict() # TODO delete when safe
        self.sems = defaultdict() # TODO delete when safe
        self.lines = []

    def create_matrix(self):
        """
        Creates a stacked density matrix for each event in each annotation file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """

        matrices = defaultdict(dict)
        for filename, filetype in self.annotation.iteritems():
            if self.is_scaled:
                matrices['ip'][filename] = mtx.scaled_region(
                    filename, self.ip, filetype,
                    self.upstream_offset, self.downstream_offset
                )
            else:
                matrices['ip'][filename] = mtx.unscaled_region(
                    filename, self.ip, filetype,
                    self.upstream_offset, self.downstream_offset
                )
        self.raw_matrices = matrices

    def create_lines(self):
        """
        For each position, remove the outlying 5% (or as specified in
        self.conf) and compute the average and mean

        Returns
        -------

        """

        c = 0  # color for iterating over default COLOR scheme
        for filename, filetype in self.annotation.iteritems():
            self.lines.append(
                LineObject.LineObject(
                    self.raw_matrices[filename],
                    filename,
                    self.conf,
                    COLORS[c],
                    MIN_EVENT_THRESHOLD
                )
            )
            c += 1

    def write_intermediate_raw_matrices_to_csv(self):
        """
        Writes to output_base + ip/input + annotation + raw_density.csv
        the raw density values for each annotation/rbp matrix created.
        """
        for filename, filetype in self.annotation.iteritems():
            output_file_ip = self.output_base + SEP + \
                             (os.path.basename(filename)) + '.ip.raw_density.txt'
            output_file_input = self.output_base + SEP + \
                                (os.path.basename(filename)) + '.input.raw_density.txt'
            self.raw_matrices['ip'][filename].to_csv(
                output_file_ip
            )
            self.raw_matrices['input'][filename].to_csv(
                output_file_input
            )

    def write_intermediate_means_to_csv(self):
        """
        Writes to output_base + annotation + norm_func + outlier + mean.txt
        the mean normalized density values over input for each
        annotation/rbp matrix created
        """
        for line in self.lines:
            output_file = self.output_base + SEP + \
                          line.file_label + '.means.txt'

            with open(output_file, 'w') as o:
                for mean in line.means:
                    o.write("{}\n".format(mean))

    def write_intermediate_sems_to_csv(self):
        """
        Writes to output_base + annotation_base + norm_func + sem.txt
        the standard error values for each position in the map.
        """
        for line in self.lines:
            output_file = self.output_base + SEP + \
                          line.file_label + '.sems.txt'

            with open(output_file, 'w') as o:
                for error in line.sems:
                    o.write("{}\n".format(error))

    def write_intermediates_to_csv(self):
        """
        Writes all intermediate files.
        """
        self.write_intermediate_raw_matrices_to_csv()
        self.write_intermediate_means_to_csv()
        self.write_intermediate_sems_to_csv()


class PhastconMap(Map):
    def __init__(self, phastcon, output_filename,
                 annotation=None, upstream_offset=0, downstream_offset=0,
                 min_density_threshold=0):
        Map.__init__(self, ip=phastcon, output_filename=output_filename,
                     norm_function=norm.get_density, annotation=annotation,
                     upstream_offset=upstream_offset,
                     downstream_offset=downstream_offset,
                     min_density_threshold=min_density_threshold,
                     is_scaled=False, conf=1)

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats
        """
        matrices = defaultdict(dict)
        for filename, filetype in self.annotation.iteritems():
            matrices['phastcon'][filename] = mtx.unscaled_cds(
                annotation=filename, density=self.ip,
                upstream_offset=self.upstream_offset,
                downstream_offset=self.downstream_offset,
                annotation_type=filetype
            )
        self.raw_matrices = matrices


class WithInput(Map):
    def __init__(self, ip, inp, output_filename, norm_function,
                 annotation=None, upstream_offset=0, downstream_offset=0,
                 min_density_threshold=0, is_scaled=False, conf=0.95):
        Map.__init__(self, ip=ip, output_filename=output_filename,
                     norm_function=norm_function, annotation=annotation,
                     upstream_offset=upstream_offset,
                     downstream_offset=downstream_offset,
                     min_density_threshold=min_density_threshold,
                     is_scaled=is_scaled, conf=conf)

        self.inp = inp
        self.lines = []

    def set_background_and_calculate_ks(self, cond_file_name, bg_file_name):
        """
        AFTER creation of all LineObjects, we can specify a condition
        matrix to run a ks-test against its distribution of values
        for each position compared to a specified background control.
        
        Parameters
        ----------
        cond_file_name : basestring
        bg_file_name : basestring
        """

        print("background: {}".format(bg_file_name))
        for line in self.lines:
            if line.annotation == cond_file_name:
                line.calculate_and_set_significance(self.norm_matrices[bg_file_name])


    def set_background_and_calculate_zscore(self, cond_file_name, bg_file_name):
        """
        AFTER creation of all LineObjects, we can specify a condition
        (bg_file_name) to run z-score comparison for each position
        of our test (cond_file_name).

        Parameters
        ----------
        cond_file_name
        bg_file_name
        """
        print("background: {}".format(bg_file_name))
        bg_line = None
        for line in self.lines:
            if line.annotation == bg_file_name:
                bg_line = line

        if bg_line is not None:
            for line in self.lines:
                if line.annotation == cond_file_name:
                    print("now testing: {} against {}".format(cond_file_name,
                                                              bg_file_name))
                    line.calculate_zscore(bg_line)

    def write_intermediate_norm_matrices_to_csv(self):
        """
        Writes intermediate normalized matrices to csv.
        """
        for line in self.lines:
            output_file = self.output_base + SEP + \
                          line.file_label + '.normed_matrix.txt'

            line.event_matrix.to_csv(output_file)

    def write_intermediate_pvalues_to_csv(self):
        """
        Writes zscores to csv.
        """
        for line in self.lines:
            output_file = self.output_base + SEP + \
                          line.file_label + '.zscores.txt'

            o = open(output_file, 'w')
            for z in line.z_scores:
                o.write("{}\n".format(z))
            o.close()
            
    def write_intermediates_to_csv(self):
        """
        Writes all intermediate matrices and values to csv.
        """
        self.write_intermediate_raw_matrices_to_csv()
        self.write_intermediate_norm_matrices_to_csv()
        self.write_intermediate_means_to_csv()
        self.write_intermediate_sems_to_csv()
        self.write_intermediate_pvalues_to_csv()
        # self.export_as_deeptool_matrix()

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats
        """
        matrices = defaultdict(dict)

        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.region(
                filename, self.ip, filetype, self.is_scaled,
                self.upstream_offset, self.downstream_offset
            )
            matrices['input'][filename] = mtx.region(
                filename, self.inp, filetype, self.is_scaled,
                self.upstream_offset, self.downstream_offset
            )
        self.raw_matrices = matrices

    def normalize_matrix(self):
        """
        For each annotation we have to parse, normalize the ip
        density values above its matched input density values.
        """
        matrices = defaultdict()
        for filename, filetype in self.annotation.iteritems():
            matrices[filename] = self.norm_function(
                self.raw_matrices['ip'][filename],
                self.raw_matrices['input'][filename],
                self.ip.pseudocount(), self.inp.pseudocount(),
                self.min_density_threshold
            )
        self.norm_matrices = matrices

    def create_lines(self):
        """
        Creates LineObject for each annotation we have to 
        parse. LineObjects contain density matrices, means, 
        errorbar positions, basically everything we need
        to construct and plot a 'line'. 
        
        Requires
        --------
        self.annotation : OrderedDict
        self.norm_matrices : dict
        self.conf : float
        """
        c = 0 # color for iterating over default COLOR scheme
        for filename, filetype in self.annotation.iteritems():
            self.lines.append(
                LineObject.LineObject(
                    self.norm_matrices[filename],
                    filename,
                    self.conf,
                    COLORS[c],
                    MIN_EVENT_THRESHOLD
                )
            )
            c+=1

    def plot(self):
        """
        Determines whether or not to plot a 1 or 2-window map.
        """
        self.plot_as_bed()
        # TODO check and fix
        """
        if self.is_scaled:
            self.plot_as_bed()
        else:
            self.plot_as_exon()
        """
    def plot_as_bed(self):
        """
        Plots the entire region in one window. This function is mostly for
        single-base regions such as CDS start sites that are slopped (flanked)
        by a fixed length on either side.
        """
        f, ax = plt.subplots(figsize=(10, 5))

        RDPlotter.plot_bed(self.lines, [ax])
        f.suptitle(misc.sane(self.output_filename))
        f.savefig(self.output_filename)

    def plot_as_exon(self):
        """
        Plots the region as two windows. This function is mostly for
        regions of unequal length (such as exons) that are not scaled
        and are plotting 3' and 5' ends as two distinct regions.
        """
        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10, 5))
        axs = [ax1, ax2]
        RDPlotter.plot_exon(self.lines, axs)
        exon1 = patches.Rectangle(
            (self.downstream_offset,-1),
            self.upstream_offset,
            2,
            alpha = 0.2
        )
        exon2 = patches.Rectangle(
            (0, -1),
            self.upstream_offset,
            2,
            alpha=0.2
        )
        ax1.add_patch(exon1)
        ax2.add_patch(exon2)
        f.suptitle(misc.sane(self.output_filename))
        f.savefig(self.output_filename)

    def export_as_deeptool_matrix(self):
        """
        Just exports as a zipped file and can be used to import into deeptools
        if you wanted to visualize using their tools.
        """
        for filename, matrix in self.norm_matrices.iteritems():
            df = misc.deeptoolify(norm.mask(matrix), self.annotation[filename])
            header = misc.create_deeptool_header(
                sample_labels=[self.output_base],
                downstream=matrix.shape[0]/2,
                upstream=matrix.shape[0]/2,
                group_boundaries=[0, matrix.shape[0]],
                sample_boundaries=[0, matrix.shape[1]],
                ref_point="BED",
                group_labels=[misc.sane(filename)],
                min_threshold=self.min_density_threshold,
            )
            o = gzip.open(self.output_base + SEP +
                os.path.splitext(os.path.basename(filename))[0] + SEP +
                '{}.deeptools.tsv.gz'.format(self.norm_function.__name__), 'wb')
            o.write(header+'\n')
            df.to_csv(o, sep='\t', index=None, header=None, compression='gzip')
            o.close()


class SkippedExon(WithInput):
    def __init__(self, ip, inp, output_filename,
                 norm_function, annotation=None, exon_offset=50,
                 intron_offset=300, min_density_threshold=0, conf=0.95):
        """

        Parameters
        ----------
        ip : density.ReadDensity
        inp : density.ReadDensity
        output_filename :
        norm_function
        annotation
        exon_offset
        intron_offset
        min_density_threshold
        conf
        """
        WithInput.__init__(self, ip=ip, inp=inp,
                           output_filename=output_filename,
                           norm_function=norm_function, annotation=annotation,
                           upstream_offset=0, downstream_offset=0,
                           min_density_threshold=min_density_threshold,
                           is_scaled=False, conf=conf)

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats
        """
        matrices = defaultdict(dict)
        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.skipped_exon(
                annotation=filename, density=self.ip,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
            matrices['input'][filename] = mtx.skipped_exon(
                annotation=filename, density=self.inp,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
        self.raw_matrices = matrices

    def plot(self):
        """
        Plots the rbp map and heatmap of z-value scores for the first 
        two items in the self.lines list.
        """

        ### Set up axes and format of plot/subplots
        map_gridspec = gridspec.GridSpec(
            ncols=4, nrows=3, width_ratios=[1, 1, 1, 1], height_ratios=[9, 1, 1]
        )
        map_gridspec.update(hspace=0.6)
        gs = gridspec.GridSpec(
            ncols=4, nrows=3, width_ratios=[1, 1, 1, 1], height_ratios=[12, 1, 1]
        )
        gs.update(hspace=0)

        f = plt.figure(figsize=(20, 10))
        plot_axs = []
        heatmap_axs = []

        plot_axs.append(f.add_subplot(map_gridspec[0]))
        for i in range(1, 4):
            plot_axs.append(f.add_subplot(map_gridspec[i], sharey=plot_axs[0]))
        for j in range(4, 8):
            heatmap_axs.append(f.add_subplot(gs[j]))
        for j in range(8, 12):
            heatmap_axs.append(f.add_subplot(gs[j]))

        ### Plot the plot stuff
        RDPlotter.plot_se(self.lines, plot_axs)

        ### Plot the heatmap stuff

        cmap_1 = colors.diverge_map(
            high=RED,  # red
            low=ORANGE  # orange/yellow
        )
        cmap_2 = colors.diverge_map(
            high=BLUE,
            low=GREEN
        )
        RDPlotter.plot_heatmap(
            self.lines[0:1], heatmap_axs[:4], cmap_1, ylabel='left',
            vmax=2, vmin=-2
        )
        RDPlotter.plot_heatmap(
            self.lines[1:2], heatmap_axs[4:], cmap_2, ylabel='right',
            vmax=2, vmin=-2
        )

        plt.tight_layout(pad=1.5 * 5.5, w_pad=0.8)
        f.savefig(self.output_filename)


class MutuallyExclusiveExon(WithInput):
    def __init__(
            self, ip, inp, output_filename,
            norm_function, annotation=None, exon_offset=50,
            intron_offset=300, min_density_threshold=0, conf=0.95
    ):
        """

        Parameters
        ----------
        ip : density.ReadDensity
        inp : density.ReadDensity
        output_filename : basestring
        norm_function : normalization_functions norm function
        annotation : dict
        exon_offset : int
        intron_offset : int
        min_density_threshold : float
        conf : float
        """
        WithInput.__init__(
            self, ip=ip, inp=inp,
            output_filename=output_filename,
            norm_function=norm_function, annotation=annotation,
            upstream_offset=0, downstream_offset=0,
            min_density_threshold=min_density_threshold,
            is_scaled=False, conf=conf
        )

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.mutually_exc_exon(
                annotation=filename, density=self.ip,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
            matrices['input'][filename] = mtx.mutually_exc_exon(
                annotation=filename, density=self.inp,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
        self.raw_matrices = matrices

    def plot(self):
        """
        Plots all lines
        Returns
        -------

        """
        # TODO fix if necessary.
        map_gridspec = gridspec.GridSpec(
            ncols=6, nrows=3, width_ratios=[1, 1, 1, 1, 1, 1], height_ratios=[9, 1, 1]
        )
        map_gridspec.update(hspace=0.6)
        gs = gridspec.GridSpec(
            ncols=6, nrows=3, width_ratios=[1, 1, 1, 1, 1, 1], height_ratios=[12, 1, 1]
        )
        gs.update(hspace=0)

        f = plt.figure(figsize=(20, 10))
        plot_axs = []
        heatmap_axs = []

        plot_axs.append(f.add_subplot(map_gridspec[0]))
        for i in range(1, 6):
            plot_axs.append(f.add_subplot(map_gridspec[i], sharey=plot_axs[0]))
        for j in range(6, 12):
            heatmap_axs.append(f.add_subplot(gs[j])) # 'left' heatmap axes
        for j in range(12, 18):
            heatmap_axs.append(f.add_subplot(gs[j])) # 'right' heatmap axes

        RDPlotter.plot_se(self.lines, plot_axs)

        cmap_1 = colors.diverge_map(
            high=RED,  # red
            low=ORANGE  # orange/yellow
        )
        cmap_2 = colors.diverge_map(
            high=BLUE,
            low=GREEN
        )
        RDPlotter.plot_heatmap(self.lines[0:1], heatmap_axs[:6], cmap_1, ylabel='left')
        RDPlotter.plot_heatmap(self.lines[1:2], heatmap_axs[6:], cmap_2, ylabel='right')

        plt.tight_layout(pad=1.5 * 5.5, w_pad=0.8)
        f.savefig(self.output_filename)

class Alt3PSpliceSite(WithInput):
    def __init__(self, ip, inp, output_filename,
                 norm_function, annotation=None, exon_offset=50,
                 intron_offset=300, min_density_threshold=0, conf=0.95):
        """

        Parameters
        ----------
        ip : density.ReadDensity
        inp : density.ReadDensity
        output_filename :
        norm_function
        annotation
        exon_offset
        intron_offset
        min_density_threshold
        conf
        """
        WithInput.__init__(self, ip=ip, inp=inp,
                           output_filename=output_filename,
                           norm_function=norm_function, annotation=annotation,
                           upstream_offset=0, downstream_offset=0,
                           min_density_threshold=min_density_threshold,
                           is_scaled=False, conf=conf)

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset
        self.lines = []

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.alt_3p_splice_site(
                annotation=filename, density=self.ip,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
            matrices['input'][filename] = mtx.alt_3p_splice_site(
                annotation=filename, density=self.inp,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
        self.raw_matrices = matrices

    def plot(self):
        map_gridspec = gridspec.GridSpec(
            ncols=3, nrows=3, width_ratios=[1, 1, 1], height_ratios=[9, 1, 1]
        )
        map_gridspec.update(hspace=0.6)
        gs = gridspec.GridSpec(
            ncols=3, nrows=3, width_ratios=[1, 1, 1], height_ratios=[12, 1, 1]
        )
        gs.update(hspace=0)
        f = plt.figure(figsize=(20, 10))
        plot_axs = []
        heatmap_axs = []

        plot_axs.append(f.add_subplot(map_gridspec[0]))
        for i in range(1, 3):
            plot_axs.append(f.add_subplot(map_gridspec[i], sharey=plot_axs[0]))
        for j in range(3, 6):
            heatmap_axs.append(f.add_subplot(gs[j]))
        for j in range(6, 9):
            heatmap_axs.append(f.add_subplot(gs[j]))

        RDPlotter.plot_a3ss(self.lines, plot_axs)

        cmap_1 = colors.diverge_map(
            high=RED,  # red
            low=ORANGE  # orange/yellow
        )
        cmap_2 = colors.diverge_map(
            high=BLUE,
            low=GREEN
        )

        RDPlotter.plot_heatmap(self.lines[0:1], heatmap_axs[:3], cmap_1, ylabel='left')
        RDPlotter.plot_heatmap(self.lines[1:2], heatmap_axs[3:], cmap_2, ylabel='right')

        plt.tight_layout(pad=1.5 * 5.5, w_pad=1)
        f.savefig(self.output_filename)


class Alt5PSpliceSite(WithInput):
    def __init__(self, ip, inp, output_filename,
                 norm_function, annotation=None, exon_offset=50,
                 intron_offset=300, min_density_threshold=0, conf=0.95):
        """

        Parameters
        ----------
        ip : density.ReadDensity
        inp : density.ReadDensity
        output_filename :
        norm_function
        annotation
        exon_offset
        intron_offset
        min_density_threshold
        conf
        """
        WithInput.__init__(self, ip=ip, inp=inp,
                           output_filename=output_filename,
                           norm_function=norm_function, annotation=annotation,
                           upstream_offset=0, downstream_offset=0,
                           min_density_threshold=min_density_threshold,
                           is_scaled=False, conf=conf)

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.alt_5p_splice_site(
                annotation=filename, density=self.ip,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
            matrices['input'][filename] = mtx.alt_5p_splice_site(
                annotation=filename, density=self.inp,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
        self.raw_matrices = matrices

    def plot(self):
        map_gridspec = gridspec.GridSpec(
            ncols=3, nrows=3, width_ratios=[1, 1, 1], height_ratios=[9, 1, 1]
        )
        map_gridspec.update(hspace=0.6)
        gs = gridspec.GridSpec(
            ncols=3, nrows=3, width_ratios=[1, 1, 1], height_ratios=[12, 1, 1]
        )
        gs.update(hspace=0)

        f = plt.figure(figsize=(20, 10))
        plot_axs = []
        heatmap_axs = []

        plot_axs.append(f.add_subplot(map_gridspec[0]))
        for i in range(1, 3):
            plot_axs.append(f.add_subplot(map_gridspec[i], sharey=plot_axs[0]))
        for j in range(3, 6):
            heatmap_axs.append(f.add_subplot(gs[j]))
        for j in range(6, 9):
            heatmap_axs.append(f.add_subplot(gs[j]))

        RDPlotter.plot_a5ss(self.lines, plot_axs)

        cmap_1 = colors.diverge_map(
            high=RED,  # red
            low=ORANGE  # orange/yellow
        )
        cmap_2 = colors.diverge_map(
            high=BLUE,
            low=GREEN
        )

        RDPlotter.plot_heatmap(self.lines[0:1], heatmap_axs[:3], cmap_1,
                               ylabel='left')
        RDPlotter.plot_heatmap(self.lines[1:2], heatmap_axs[3:], cmap_2,
                               ylabel='right')

        plt.tight_layout(pad=1.5 * 5.5, w_pad=1)
        f.savefig(self.output_filename)


class RetainedIntron(WithInput):
    def __init__(self, ip, inp, output_filename,
                 norm_function, annotation=None, exon_offset=50,
                 intron_offset=300, min_density_threshold=0, conf=0.95):
        """

        Parameters
        ----------
        ip : density.ReadDensity
        inp : density.ReadDensity
        output_filename :
        norm_function
        annotation
        exon_offset
        intron_offset
        min_density_threshold
        conf
        """
        WithInput.__init__(self, ip=ip, inp=inp,
                           output_filename=output_filename,
                           norm_function=norm_function, annotation=annotation,
                           upstream_offset=0, downstream_offset=0,
                           min_density_threshold=min_density_threshold,
                           is_scaled=False, conf=conf)

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.retained_intron(
                annotation=filename, density=self.ip,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
            matrices['input'][filename] = mtx.retained_intron(
                annotation=filename, density=self.inp,
                exon_offset=self.exon_offset, intron_offset=self.intron_offset,
                annotation_type=filetype
            )
        self.raw_matrices = matrices

    def plot(self):
        """
        f, (ax1, ax2, ) = plt.subplots(
            1, 2, sharey=True, figsize=(16, 8)
        )
        axs = [ax1, ax2]
        RDPlotter.plot_ri(self.means, self.sems, axs)
        plt.tight_layout()
        f.suptitle(misc.sane(self.output_filename))
        f.savefig(self.output_filename)
        """
        map_gridspec = gridspec.GridSpec(
            ncols=2, nrows=3, width_ratios=[1, 1], height_ratios=[9, 1, 1]
        )
        map_gridspec.update(hspace=0.6)
        gs = gridspec.GridSpec(
            ncols=2, nrows=3, width_ratios=[1, 1], height_ratios=[12, 1, 1]
        )
        gs.update(hspace=0)

        f = plt.figure(figsize=(20, 10))
        plot_axs = []
        heatmap_axs = []

        plot_axs.append(f.add_subplot(map_gridspec[0]))
        for i in range(1, 2):
            plot_axs.append(f.add_subplot(map_gridspec[i], sharey=plot_axs[0]))
        for j in range(2, 4):
            heatmap_axs.append(f.add_subplot(gs[j]))
        for j in range(4, 6):
            heatmap_axs.append(f.add_subplot(gs[j]))

        RDPlotter.plot_ri(self.lines, plot_axs)

        cmap_1 = colors.diverge_map(
            high=RED,  # red
            low=ORANGE  # orange/yellow
        )
        cmap_2 = colors.diverge_map(
            high=BLUE,
            low=GREEN
        )

        RDPlotter.plot_heatmap(self.lines[0:1], heatmap_axs[:2], cmap_1,
                               ylabel='left')
        RDPlotter.plot_heatmap(self.lines[1:2], heatmap_axs[2:], cmap_2,
                               ylabel='right')

        plt.tight_layout(pad=1.5 * 5.5, w_pad=0.5)
        f.savefig(self.output_filename)


class ATACIntron(RetainedIntron):
    def __init__(self, ip, inp, output_filename,
                 norm_function, annotation=None, exon_offset=50,
                 intron_offset=300, min_density_threshold=0, conf=0.95):
        RetainedIntron.__init__(
            self, ip, inp, output_filename,
            norm_function, annotation, exon_offset,
            intron_offset, min_density_threshold, conf
        )

class UnscaledCDS(WithInput):
    def __init__(self, ip, inp, output_filename,
                 norm_function, annotation=None, upstream_offset=50,
                 downstream_offset=50, min_density_threshold=0, conf=0.95):
        """

        Parameters
        ----------
        ip : density.ReadDensity
        inp : density.ReadDensity
        output_filename :
        norm_function
        annotation
        exon_offset
        intron_offset
        min_density_threshold
        conf
        """
        WithInput.__init__(
            self, ip=ip, inp=inp, output_filename=output_filename,
            norm_function=norm_function, annotation=annotation,
            upstream_offset=upstream_offset,
            downstream_offset=downstream_offset,
            min_density_threshold=min_density_threshold,
            is_scaled=False, conf=conf
        )

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.unscaled_cds(
                annotation=filename, density=self.ip,
                upstream_offset=self.upstream_offset,
                downstream_offset=self.downstream_offset,
                annotation_type=filetype
            )
            matrices['input'][filename] = mtx.unscaled_cds(
                annotation=filename, density=self.inp,
                upstream_offset=self.upstream_offset,
                downstream_offset=self.downstream_offset,
                annotation_type=filetype
            )
        self.raw_matrices = matrices

    def plot(self):
        """
        f, (ax1, ax2,) = plt.subplots(
            1, 2, sharey=True, figsize=(16, 8)
        )
        axs = [ax1, ax2]
        RDPlotter.plot_unscaled_cds(
            self.means, self.sems, axs,
            self.upstream_offset, self.downstream_offset
        )
        plt.tight_layout(pad=8, w_pad=3, h_pad=5)
        f.subplots_adjust(wspace=0)
        f.suptitle(misc.sane(self.output_filename))
        f.savefig(self.output_filename)
        """
        map_gridspec = gridspec.GridSpec(
            ncols=2, nrows=3, width_ratios=[1, 1], height_ratios=[9, 1, 1]
        )
        map_gridspec.update(hspace=0.6)
        gs = gridspec.GridSpec(
            ncols=2, nrows=3, width_ratios=[1, 1], height_ratios=[12, 1, 1]
        )
        gs.update(hspace=0)

        f = plt.figure(figsize=(20, 10))
        plot_axs = []
        heatmap_axs = []

        plot_axs.append(f.add_subplot(map_gridspec[0]))
        for i in range(1, 2):
            plot_axs.append(f.add_subplot(map_gridspec[i], sharey=plot_axs[0]))
        for j in range(2, 4):
            heatmap_axs.append(f.add_subplot(gs[j]))
        for j in range(4, 6):
            heatmap_axs.append(f.add_subplot(gs[j]))

        RDPlotter.plot_unscaled_cds(self.lines, plot_axs)

        cmap_1 = colors.diverge_map(
            high=RED,  # red
            low=ORANGE  # orange/yellow
        )
        cmap_2 = colors.diverge_map(
            high=BLUE,
            low=GREEN
        )

        RDPlotter.plot_heatmap(self.lines[0:1], heatmap_axs[:2], cmap_1,
                               ylabel='left')
        RDPlotter.plot_heatmap(self.lines[1:2], heatmap_axs[2:], cmap_2,
                               ylabel='right')

        plt.tight_layout(pad=1.5 * 5.5, w_pad=0.5)
        f.savefig(self.output_filename)
