#!/usr/bin/env python

"""
Created on Jun 27, 2016

Modules that help containerize the rbp-map information, including:

ip : CLIP IP density
input : CLIP Input density
output filename : name of the output file and base
norm function : normalization to use (see normalization_functions.py)
annotation_src_file : annotation_src_file file(s) to use and their type. This is a dictionary of
  {filename:filetype(rmats, miso, etc.)} where the filename is used to 
  associated which density values belong to which condition.
offsets : either describing upstream/downstream or exon/intron bases to plotter
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
import matplotlib.pyplot as plt
import gzip
from collections import defaultdict, OrderedDict
from plotter import Plotter
import matrix as mtx
import misc
import normalization_functions as norm
import LineObject
import Peak
import ReadDensity
import pandas as pd
import Feature
from tqdm import trange
import numpy as np

SEP = '.'  # file delimiter

MIN_EVENT_THRESHOLD=100 # number of events required to not grey the line out

COLOR_PALETTE = sns.color_palette("hls", 8)

BG1_COLOR = 'black' # COLOR_PALETTE['black']
BG2_COLOR = COLOR_PALETTE[6]
BG3_COLOR = COLOR_PALETTE[7]
BG4_COLOR = COLOR_PALETTE[4]
BG5_COLOR = COLOR_PALETTE[2]
BG6_COLOR = COLOR_PALETTE[1]
POS_COLOR = COLOR_PALETTE[0]
NEG_COLOR = COLOR_PALETTE[5]

COLORS = [POS_COLOR, NEG_COLOR, BG1_COLOR, BG2_COLOR, BG3_COLOR, BG4_COLOR, BG5_COLOR, BG6_COLOR]

RED = COLOR_PALETTE[0]
ORANGE = COLOR_PALETTE[1]
BLUE = COLOR_PALETTE[5]
GREEN = COLOR_PALETTE[3]

MAX_VAL = 100000000
MIN_VAL = -100000000

class Map:
    def __init__(self, ip, output_filename, norm_function,
                 annotation, upstream_offset=0, downstream_offset=0,
                 min_density_threshold=0, conf=0.95, scale=False):
        """

        Parameters
        ----------
        ip : density.ReadDensity
            ReadDensity file containing RBP density signals
        output_filename : str
            Output filename for map (ie. plotter.svg).
        norm_function : function
            Function
            See: normalization_functions
        annotation_src_file : collections.OrderedDict
            {file:type} representing an annotation_src_file file and format
            See: Feature
        upstream_offset : int
            number of bases upstream of described feature to plotter
        downstream_offset : int
            number of bases downstream of described feature to plotter
        min_density_threshold : int
            Minimum density sum across an event to use
        conf : float
            Representing percentage of events to keep (default .95)
        raw_matrices : collections.defaultdict(dict)
            raw_matrices are structured as: [clip][filename] whose keys refer to:
                clip : 'ip' or 'input' depending on the experimental variable
                filename : the filename of the annotation which was overlapped
                    by the clip.
            And whose values refer to:
                The raw density values overlapped at each position stored
                in a pandas.DataFrame.

        norm_matrices : collections.defaultdict(dict)
            norm_matrices are structured as: [filename] whose keys refer to:
                filename : the filename of the annotation which was overlapped
                    by the clip.

            And whose values refer to:
                The normalized densities resulting from normalization of
                ip and input raw matrices for each key in both dictionaries.

        num_events : collections.defaultdict(dict)
            This is a LIST containing the number of events at each position.
            I made this a list due to some positions on a map being associated
            with different numbers of events. For example, the meta maps will
            calculate the first 10 positions as having being from the 5'UTR,
            of which 100 genes are considered. The next 50 positions refer to
            CDS regions, of which 1000 genes may be considered. Finally, the
            last 40 positions are dedicated to 3'UTR densities, of which 800
            genes are considered, therefore the num_events will be:
            [100, 100, ... 100(10), 1000, 1000, ..., 1000(50), 800, ..., 800(40)]
            for a total length of 100.

            Num events are structured as: [clip][filename] whose keys refer to:
                clip : 'ip' or 'input' depending on the experimental variable
                filename : the filename of the annotation which was overlapped
                    by the clip.
            And whose values refer to:
                A list containing the number of events that are considered
                when performing certain normalizations (see above list example).

        """
        self.ip = ip
        self.map_type = self.get_map_type()
        self.output_filename = output_filename
        self.output_base = os.path.splitext(output_filename)[0]
        self.norm_function = norm_function
        self.annotation = annotation
        self.min_density_threshold = min_density_threshold
        self.upstream_offset = upstream_offset
        self.downstream_offset = downstream_offset
        self.conf = conf

        self.raw_matrices = defaultdict(dict)
        self.norm_matrices = defaultdict(dict)
        self.num_events = defaultdict(dict)

        self.lines = []

        self.scale = scale

    def get_map_type(self):
        """
        Returns the string-ified version of the type of map this is.

        Returns
        -------
        map_type: basestring
        """
        if(isinstance(self.ip, Peak.Peak)):
            return 'peak'
        elif(isinstance(self.ip, ReadDensity.ReadDensity)):
            return 'density'
        elif(isinstance(self.ip, ReadDensity.Phastcon)):
            return 'phastcon'

    def create_matrix(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """

        matrices = defaultdict(dict)
        num_events = defaultdict(dict)

        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.same_length_region(
                filename, self.ip, filetype,
                self.upstream_offset, self.downstream_offset,
                self.scale
            )
            num_events['ip'][filename] = [matrices['ip'][filename].shape[0]] * matrices['ip'][filename].shape[1]

        self.raw_matrices = matrices

        self.num_events = num_events

    def create_lines(self):
        """
        For each position, remove the outlying 5% (or as specified in
        self.conf) and compute the average and mean

        Returns
        -------

        """
        c = 0  # plotter for iterating over default COLOR scheme
        for filename, _ in self.annotation.iteritems():
            self.lines.append(
                LineObject.create_line(
                    event_matrix=self.raw_matrices['ip'][filename],
                    annotation_src_file=filename,
                    conf=self.conf,
                    color=COLORS[c],
                    min_event_threshold=MIN_EVENT_THRESHOLD,
                    map_type=self.map_type,
                    num_events=self.num_events['ip'][filename]
                )
            )
            c += 1

    def write_intermediate_raw_matrices_to_csv(self):
        """
        Writes to output_base + ip/input + annotation_src_file + raw_density.csv
        the raw density values for each annotation_src_file/rbp matrix created.
        """
        for filename, filetype in self.annotation.iteritems():
            for key in self.raw_matrices.keys():
                output_file_ip = self.output_base + SEP + \
                                 (os.path.basename(filename)) + '.{}.raw_density.txt'.format(key)
                # output_file_input = self.output_base + SEP + \
                #                     (os.path.basename(filename)) + '.input.raw_density.txt'
                self.raw_matrices[key][filename].to_csv(
                    output_file_ip
                )

                # self.raw_matrices['input'][filename].to_csv(
                #     output_file_input
                # )

    def write_intermediate_means_to_csv(self):
        """
        Writes to output_base + annotation_src_file + norm_func + outlier + mean.txt
        the mean normalized density values over input for each
        annotation_src_file/rbp matrix created
        """
        for line in self.lines:
            if line.has_mean():
                # print("label: {}".format(line.file_label))
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
        # self.write_intermediate_sems_to_csv()

    def plot(self, condition_list):
        Plotter.plot_bed(self.lines, self.output_filename, self.map_type, condition_list)


class WithInput(Map):
    def __init__(self, ip, inp, output_filename, norm_function,
                 annotation=None, upstream_offset=0, downstream_offset=0,
                 min_density_threshold=0, conf=0.95, scale=False):
        Map.__init__(self, ip=ip, output_filename=output_filename,
                     norm_function=norm_function, annotation=annotation,
                     upstream_offset=upstream_offset,
                     downstream_offset=downstream_offset,
                     min_density_threshold=min_density_threshold,
                     conf=conf, scale=scale)

        self.inp = inp
        self.lines = []

    def set_background_and_calculate_significance(
            self, cond_file_names, bg_file_name, test='mannwhitneyu', num_permutations=1000, boundary_percent=0.5
    ):
        """
        AFTER creation of all LineObjects, we can specify a condition
        matrix to run a ks-test against its distribution of values
        for each position compared to a specified background control.

        NOTE: these calculations are done BEFORE outlier removal masking.
        Well all except for zscores, which require the mean of the outlier
        removed matrices.

        We can change this later, but I don't want to deal with performing
        statistics on masked NaNs and things. To change, we need to
        call: norm.get_means_and_sems() on this background matrix, and
        also change the LineObject tests to use self.outlier_removed_matrix
        instead of self.event_matrix as test case.

        Parameters
        ----------
        cond_file_names : list
        bg_file_name : basestring
        num_permutations : int
            number of iterations we choose to randomly sample bg
        boundary_percent : float
            the extreme % value from which to (out of 1000 values, take the top and bottom 0.5%, or 5)
        """

        if test == 'permutation':
            # choose starting values for top and bottom median values
            top_values = [MIN_VAL]*len(self.num_events['ip'][bg_file_name])
            bottom_values = [MAX_VAL]*len(self.num_events['ip'][bg_file_name])

            for condition in cond_file_names:
                # select output filename TODO: move out
                tsv = os.path.join(
                    os.path.dirname(self.output_filename),
                    '{}.{}.{}.randsample.tsv'.format(
                        os.path.basename(self.output_filename),
                        os.path.basename(bg_file_name),
                        os.path.basename(condition)
                    )
                )

                subset_iterations = []
                iterations = num_permutations  # how many random samplings to take.
                percentile = boundary_percent  # the extreme % value from which to (out of 1000 values, take the top and bottom 0.5%, or 5)
                progress = trange(iterations)
                for i in range(0, iterations):
                    # since event_num is reported as a list across all positions, simply get the average for now.
                    mean_event_num = int(
                        sum(self.num_events['ip'][condition]) / float(len(self.num_events['ip'][condition]))
                    )
                    # get n random events where n is the number of events in incl/excl
                    rand_subset = Feature.get_random_sample(self.norm_matrices[bg_file_name], mean_event_num)
                    ## print("random subset: {}".format(rand_subset[:5]))
                    # remove outliers
                    means, _, _, _ = norm.get_means_and_sems(rand_subset, conf=self.conf)
                    # print("means after outlier removal: {}".format(means[:5]))
                    subset_iterations.append(pd.Series(means))

                    progress.update(1)
                # concatenate all "lines" (means of outlier-removed normalized data)
                df = pd.concat(subset_iterations, axis=1).T
                df.to_csv(tsv, sep='\t')
                bottom_values_condition, top_values_condition = norm.bottom_top_values_from_dataframe(df, percentile, percentile)

                # get the min "bottom values" and max "top_values" among each condition
                for position in range(0, len(bottom_values_condition)):
                    if bottom_values_condition[position] < bottom_values[position]:
                        bottom_values[position] = bottom_values_condition[position]

                for position in range(0, len(top_values_condition)):
                    if top_values_condition[position] > top_values[position]:
                        top_values[position] = top_values_condition[position]
                # print("top values: {}".format(top_values))

                # replace max/min defaults if nonetypes were found in any position
                bottom_values = [np.nan if (x == MAX_VAL or x == MIN_VAL) else x for x in bottom_values]
                top_values = [np.nan if (x == MAX_VAL or x == MIN_VAL) else x for x in top_values]

            # change the std error boundaries in background, and remove error boundaries in conditions.
            for line in self.lines:
                if line.annotation_src_file == bg_file_name:
                    line._set_std_error_boundaries(bottom_values, top_values)
                elif line.annotation_src_file != bg_file_name:
                    line._set_std_error_boundaries(line.means, line.means)
        else:
            for line in self.lines:
                for condition in cond_file_names:
                    # print('condition: [{}], bg: [{}]'.format(condition, bg_file_name))
                    if line.annotation_src_file == condition:
                        line.calculate_and_set_significance(
                            self.norm_matrices[bg_file_name], test
                        )

        # mark the condition which we are using for test
        for line in self.lines:
            if line.annotation_src_file == bg_file_name:
                if not line.label.endswith('*'):  # we probably just want one asterisk
                    line.label = line.label + '*'

    def write_intermediate_norm_matrices_to_csv(self):
        """
        Writes intermediate normalized matrices to csv.
        """
        for line in self.lines:
            output_file = self.output_base + SEP + \
                          line.file_label + '.normed_matrix.txt'

            line.event_matrix.to_csv(output_file)

    def write_intermediate_hist_to_csv(self):
        """
        Writes intermediate normalized matrices to csv.
        """
        for line in self.lines:
            if line.has_hist():
                output_file = self.output_base + SEP + \
                              line.file_label + '.hist.txt'
                o = open(output_file, 'w')
                for l in line.values:
                    o.write('{}\n'.format(l))
                o.close()

    def write_intermediate_pvalues_to_csv(self):
        """
        Writes zscores to csv.
        """
        for line in self.lines:
            if line.has_pvalues():
                output_file = self.output_base + SEP + \
                              line.file_label + '.pvalues.txt'

                o = open(output_file, 'w')
                pos = 0
                for p in line.p_values:
                    o.write("{}\t{}\n".format(pos,p))
                    pos+=1
                o.close()
            
    def write_intermediates_to_csv(self):
        """
        Writes all intermediate matrices and values to csv.
        Uncomment some of these if we want all of the intermediates.
        Unless something is really weird though, we probably don't need to
        be creating these every time.

        """
        self.write_intermediate_raw_matrices_to_csv()
        self.write_intermediate_norm_matrices_to_csv()
        self.write_intermediate_means_to_csv()
        # self.write_intermediate_sems_to_csv()
        self.write_intermediate_pvalues_to_csv()
        self.write_intermediate_hist_to_csv()
        # self.export_as_deeptool_matrix()
        self.write_sum_coverage_to_csv()

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats
        """
        matrices = defaultdict(dict)
        num_events = defaultdict(dict)

        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.same_length_region(
                filename, self.ip, filetype,
                self.upstream_offset, self.downstream_offset, self.scale,
            )
            matrices['input'][filename] = mtx.same_length_region(
                filename, self.inp, filetype,
                self.upstream_offset, self.downstream_offset, self.scale,
            )
            # TODO: maybe do this by position? Currently we just multiple the positions by the matrix shape, but it would be 'cleaner' if we calculated number of events at every position instead.
            num_events['ip'][filename] = [matrices['ip'][filename].shape[0]] * matrices['ip'][filename].shape[1]

        self.raw_matrices = matrices
        self.num_events = num_events

    def normalize_matrix(self):
        """
        For each annotation_src_file we have to parse, scale the ip
        density values above its matched input density values.
        """
        matrices = defaultdict()
        for filename, _ in self.annotation.iteritems():
            matrices[filename] = self.norm_function(
                self.raw_matrices['ip'][filename],
                self.raw_matrices['input'][filename],
                self.ip.pseudocount(), self.inp.pseudocount(),
                self.min_density_threshold
            )
        self.norm_matrices = matrices

    def create_lines(self):
        """
        Creates LineObject for each annotation_src_file we have to
        parse. LineObjects contain density matrices, means, 
        errorbar positions, basically everything we need
        to construct and plotter a 'line'. 
        
        Requires
        --------
        self.annotation_src_file : OrderedDict
        self.norm_matrices : dict
        self.conf : float
        """
        c = 0 # plotter for iterating over default COLOR scheme


        for filename, filetype in self.annotation.iteritems():
            self.lines.append(
                LineObject.create_line(
                    event_matrix=self.norm_matrices[filename],
                    annotation_src_file=filename,
                    conf=self.conf,
                    color=COLORS[c],
                    min_event_threshold=MIN_EVENT_THRESHOLD,
                    map_type=self.map_type,
                    num_events=self.num_events['ip'][filename]
                )
            )
            c+=1



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


    def write_sum_coverage_to_csv(self):
        """
        Writes sum coverages across each event to csv.
        """
        for filename, filetype in self.annotation.iteritems():
            for key in self.raw_matrices.keys():
                output_file_cov = self.output_base + SEP + \
                                 (os.path.basename(filename)) + '.{}.sum_coverage.txt'.format(key)
                # output_file_input = self.output_base + SEP + \
                #                     (os.path.basename(filename)) + '.input.raw_density.txt'
                raw_matrix = self.raw_matrices[key][filename]
                raw_matrix = raw_matrix.replace(-1, 0)
                sum_cov = pd.Series(raw_matrix.sum(axis=1))
                sum_cov.to_csv(
                    output_file_cov
                )


class Bed(WithInput):
    def __init__(
        self, ip, inp, output_filename, norm_function,
        annotation=None, upstream_offset=0, downstream_offset=0,
        min_density_threshold=0, conf=0.95, scale=False
    ):
        WithInput.__init__(
            self, ip, inp, output_filename, norm_function,
            annotation=annotation, upstream_offset=upstream_offset,
            downstream_offset=downstream_offset,
            min_density_threshold=min_density_threshold, conf=conf,
            scale=scale
        )

class MultiLengthBed(Bed):
    def __init__(
        self, ip, inp, output_filename, norm_function,
        annotation=None, upstream_offset=50, downstream_offset=50,
        min_density_threshold=0, conf=0.95
    ):
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
        Bed.__init__(
            self, ip, inp, output_filename, norm_function,
            annotation=annotation, upstream_offset=upstream_offset,
            downstream_offset=downstream_offset,
            min_density_threshold=min_density_threshold, conf=conf,
            scale=False
        )

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        num_events = defaultdict(dict)

        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.multi_length_regions(
                annotation=filename, density=self.ip,
                upstream_offset=self.upstream_offset,
                downstream_offset=self.downstream_offset,
                annotation_type=filetype
            )
            matrices['input'][filename] = mtx.multi_length_regions(
                annotation=filename, density=self.inp,
                upstream_offset=self.upstream_offset,
                downstream_offset=self.downstream_offset,
                annotation_type=filetype
            )
            num_events['ip'][filename] = [matrices['ip'][filename].shape[0]] * \
                                         matrices['ip'][filename].shape[1]

        self.raw_matrices = matrices
        self.num_events = num_events

    def plot(self, condition_list):
        Plotter.plot_multi_length_bed(self.lines, self.output_filename, self.map_type, condition_list)

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
                           conf=conf)

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats
        """
        matrices = defaultdict(dict)
        num_events = defaultdict(dict)

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
            num_events['ip'][filename] = [matrices['ip'][filename].shape[0]] * \
                                         matrices['ip'][filename].shape[1]

        self.raw_matrices = matrices
        self.num_events = num_events

    def plot(self, condition_list):
        """
        Plots the rbp map and heatmap of z-value scores for the first 
        two items in the self.lines list.
        """

        ### Set up axes and format of plotter/subplots

        Plotter.plot_se(self.lines, self.output_filename, self.map_type, condition_list)

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
        min_density_threshold : int
        conf : float
        """
        WithInput.__init__(
            self, ip=ip, inp=inp,
            output_filename=output_filename,
            norm_function=norm_function, annotation=annotation,
            upstream_offset=0, downstream_offset=0,
            min_density_threshold=min_density_threshold,
            conf=conf
        )

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        num_events = defaultdict(dict)
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
            num_events['ip'][filename] = [matrices['ip'][filename].shape[0]] * \
                                         matrices['ip'][filename].shape[1]

        self.raw_matrices = matrices
        self.num_events = num_events

    def plot(self, condition_list):
        """
        Plots all lines
        Returns
        -------

        """
        Plotter.plot_mxe(self.lines, self.output_filename, self.map_type, condition_list)

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
                           conf=conf)

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset
        self.lines = []

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        num_events = defaultdict(dict)

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
            num_events['ip'][filename] = [matrices['ip'][filename].shape[0]] * \
                                         matrices['ip'][filename].shape[1]

        self.raw_matrices = matrices
        self.num_events = num_events

    def plot(self, condition_list):
        Plotter.plot_a3ss(self.lines, self.output_filename, self.map_type, condition_list)


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
                           conf=conf)

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        num_events = defaultdict(dict)

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
            num_events['ip'][filename] = [matrices['ip'][filename].shape[0]] * \
                                         matrices['ip'][filename].shape[1]

        self.raw_matrices = matrices
        self.num_events = num_events

    def plot(self, condition_list):
        Plotter.plot_a5ss(self.lines, self.output_filename, self.map_type, condition_list)

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
                           conf=conf)

        self.exon_offset = exon_offset
        self.intron_offset = intron_offset

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats

        Returns
        -------

        """
        matrices = defaultdict(dict)
        num_events = defaultdict(dict)
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
            num_events['ip'][filename] = [matrices['ip'][filename].shape[0]] * \
                                         matrices['ip'][filename].shape[1]

        self.raw_matrices = matrices
        self.num_events = num_events

    def plot(self, condition_list):
        Plotter.plot_ri(self.lines, self.output_filename, self.map_type, condition_list)

### The maps below are either deprecated or experimental

class Metagene(WithInput):
    def __init__(self, ip, inp, output_filename,
                 norm_function, annotation=None, upstream_offset=0,
                 downstream_offset=0, min_density_threshold=0, conf=0.95):
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
                           upstream_offset=upstream_offset,
                           downstream_offset=downstream_offset,
                           min_density_threshold=min_density_threshold,
                           conf=conf)


    def create_matrices(self):
        three_utr_ratio = 62 # 44  # TODO: remove hardcoded percentage files
        five_utr_ratio = 16 # 17
        cds_ratio = 100

        matrices = defaultdict()
        num_events = defaultdict()

        for filename, filetype in self.annotation.iteritems():
            if filetype == '3utr' or filetype == 'utr3':
                matrices["three_prime_utr_ip"] = mtx.meta(
                    annotation=filename, density=self.ip,
                    upstream_offset=self.upstream_offset,
                    downstream_offset=self.downstream_offset,
                    annotation_type='bed', scale_to=three_utr_ratio
                )

                matrices["three_prime_utr_input"] = mtx.meta(
                    annotation=filename, density=self.inp,
                    upstream_offset=self.upstream_offset,
                    downstream_offset=self.downstream_offset,
                    annotation_type='bed', scale_to=three_utr_ratio
                )

                num_events["three_prime_utr_ip"] = [
                    matrices["three_prime_utr_ip"].shape[0]
                ] * matrices["three_prime_utr_ip"].shape[1]

            elif filetype == '5utr' or filetype == 'utr5':
                matrices["five_prime_utr_ip"] = mtx.meta(
                    annotation=filename, density=self.ip,
                    upstream_offset=self.upstream_offset,
                    downstream_offset=self.downstream_offset,
                    annotation_type='bed', scale_to=five_utr_ratio
                )

                matrices["five_prime_utr_input"] = mtx.meta(
                    annotation=filename, density=self.inp,
                    upstream_offset=self.upstream_offset,
                    downstream_offset=self.downstream_offset,
                    annotation_type='bed', scale_to=five_utr_ratio
                )

                num_events["five_prime_utr_ip"] = [
                   matrices["five_prime_utr_ip"].shape[0]
                ] * matrices["five_prime_utr_ip"].shape[1]

            elif filetype == 'cds':
                matrices["cds_ip"] = mtx.meta(
                    annotation=filename, density=self.ip,
                    upstream_offset=self.upstream_offset,
                    downstream_offset=self.downstream_offset,
                    annotation_type='bed', scale_to=cds_ratio
                )

                matrices["cds_input"] = mtx.meta(
                    annotation=filename, density=self.inp,
                    upstream_offset=self.upstream_offset,
                    downstream_offset=self.downstream_offset,
                    annotation_type='bed', scale_to=cds_ratio
                )

                num_events["cds_ip"] = [
                    matrices["cds_ip"].shape[0]
                ] * matrices["cds_ip"].shape[1]

            else:
                print('unknown filetype! for metagene!')

        # combine/merge all regions
        self.raw_matrices['ip']['meta'] = pd.merge(
            matrices['five_prime_utr_ip'],
            matrices['cds_ip'],
            how='outer', left_index=True, right_index=True
        ).merge(
            matrices['three_prime_utr_ip'],
            how='outer', left_index=True, right_index=True
        )

        self.raw_matrices['input']['meta'] = pd.merge(
            matrices['five_prime_utr_input'],
            matrices['cds_input'],
            how='outer', left_index=True, right_index=True
        ).merge(
            matrices['three_prime_utr_input'],
            how='outer', left_index=True, right_index=True
        )

        self.raw_matrices['ip']['meta'].columns = range(
            five_utr_ratio + cds_ratio + three_utr_ratio
        )
        self.raw_matrices['input']['meta'].columns = range(
            five_utr_ratio + cds_ratio + three_utr_ratio
        )
        self.annotation = {'meta':'metagene'}
        self.num_events['ip']['meta'] = num_events['five_prime_utr_ip'] + \
                                        num_events['cds_ip'] + \
                                        num_events['three_prime_utr_ip']

        # self.num_events['input']['meta'] = num_events['five_prime_utr_input'] + \
        #                                 num_events['cds_input'] + \
        #                                 num_events['three_prime_utr_input']

    def plot(self, condition_list):
        """
        Plots the rbp map and heatmap of z-value scores for the first
        two items in the self.lines list.
        """

        ### Set up axes and format of plotter/subplots
        Plotter.plot_meta(self.lines, self.output_filename, self.map_type, condition_list)

class CDS(WithInput):
    def __init__(self, ip, inp, output_filename,
                 norm_function, annotation=None, upstream_offset=0,
                 downstream_offset=0, min_density_threshold=0, conf=0.95):
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
                           upstream_offset=upstream_offset,
                           downstream_offset=downstream_offset,
                           min_density_threshold=min_density_threshold,
                           conf=conf)


    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats
        """
        matrices = defaultdict(dict)
        for filename, filetype in self.annotation.iteritems():
            matrices['ip'][filename] = mtx.meta(
                annotation=filename, density=self.ip,
                upstream_offset=self.upstream_offset,
                downstream_offset=self.downstream_offset,
                annotation_type=filetype
            )
            matrices['input'][filename] = mtx.meta(
                annotation=filename, density=self.inp,
                upstream_offset=self.upstream_offset,
                downstream_offset=self.downstream_offset,
                annotation_type=filetype
            )
        self.raw_matrices = matrices

    def plot(self, condition_list):
        """
        Plots the rbp map and heatmap of z-value scores for the first 
        two items in the self.lines list.
        """

        ### Set up axes and format of plotter/subplots
        Plotter.plot_bed(self.lines, self.output_filename, self.map_type, condition_list)

class ATACIntron(RetainedIntron):
    def __init__(self, ip, inp, output_filename,
                 norm_function, annotation=None, exon_offset=50,
                 intron_offset=300, min_density_threshold=0, conf=0.95):
        RetainedIntron.__init__(
            self, ip, inp, output_filename,
            norm_function, annotation, exon_offset,
            intron_offset, min_density_threshold, conf
        )


class PhastconMap(Map):
    def __init__(self, phastcon, peak, output_filename,
                 annotation, upstream_offset, downstream_offset,
                 min_density_threshold,
                 masked_file):
        Map.__init__(self, phastcon, output_filename, norm_function = norm.get_density,
                 annotation=annotation, upstream_offset=upstream_offset, downstream_offset=downstream_offset,
                 min_density_threshold=min_density_threshold, conf=1, scale=False)
        self.peak = peak # for masking (see: create_matrices)
        self.masked_file = masked_file # for masking (see: create_matrices)
        
    def create_matrices_meta_DEPRECATED(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats
        """
        matrices = defaultdict(dict)
        num_events = defaultdict(dict)

        for filename, filetype in self.annotation.iteritems():
            matrices['phastcon'][filename] = mtx.meta(
                annotation=filename, density=self.ip,
                upstream_offset=self.upstream_offset,
                downstream_offset=self.downstream_offset,
                annotation_type=filetype
            )
            num_events['phastcon'][filename] = [
                matrices['phastcon'][filename].shape[0]
            ] * matrices['phastcon'][filename].shape[1]
        self.raw_matrices = matrices
        self.num_events = num_events

    def create_matrices(self):
        """
        Creates a stacked density matrix for each event in each annotation_src_file file
        and sets self.raw_matrices variable.

        For Example:

        raw_matrices['ip']['condition1.rmats'] = matrix of density values for
        this RBP intersected with the events described by condition1.rmats
        """
        matrices = defaultdict(dict)
        num_events = defaultdict(dict)
        ct = 0

        for filename, filetype in self.annotation.iteritems():
            if filename == self.masked_file:
                matrices['phastcon'][filename] = mtx.phastcon_region(
                    annotation=filename, density=self.ip,
                    exon_offset=self.upstream_offset,
                    intron_offset=self.downstream_offset,
                    annotation_type=filetype,
                    peak=self.peak,
                    mask_df=True
                )
            else:
                # Don't mask the background file, leave for all

                matrices['phastcon'][filename] = mtx.phastcon_region(
                    annotation=filename, density=self.ip,
                    exon_offset=self.upstream_offset,
                    intron_offset=self.downstream_offset,
                    annotation_type=filetype,
                    peak=self.peak,
                    mask_df=False
                )
            num_events['phastcon'][filename] = [
                matrices['phastcon'][filename].shape[0]
            ] * matrices['phastcon'][filename].shape[1]
            ct += 1

        self.raw_matrices = matrices
        self.num_events = num_events

    def create_lines(self):
        """
        For each position, remove the outlying 5% (or as specified in
        self.conf) and compute the average and mean

        Returns
        -------

        """

        c = 0  # plotter for iterating over default COLOR scheme
        for filename, _ in self.annotation.iteritems():
            self.lines.append(
                LineObject.create_line(
                    event_matrix=self.raw_matrices['phastcon'][filename],
                    annotation_src_file=filename,
                    conf=self.conf,
                    color=COLORS[c],
                    min_event_threshold=MIN_EVENT_THRESHOLD,
                    map_type=self.map_type,
                    num_events=self.num_events['phastcon'][filename]
                )
            )
            c += 1
        for line in self.lines:
            if line.annotation_src_file == self.masked_file:
                if not line.label.endswith('*'):  # we probably just want one asterisk
                    line.label = line.label + '*'

    def plot(self, condition_list):
        Plotter.plot_phastcon(self.lines, self.output_filename, self.map_type, condition_list)

