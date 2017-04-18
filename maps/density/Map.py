#!/usr/bin/env python

"""
Created on Jun 27, 2016

Modules that help containerize the rbp-map information, including:

ip : CLIP IP density
input : CLIP Input density
output filename : name of the output file and base
norm function : normalization to use (see normalization_functions.py)
annotation : annotation file(s) to use
offsets : either describing upstream/downstream or exon/intron bases to plot
is_scaled : whether or not we want the map to 'scale' or keep nucleotide dims
conf : percentage of scores to keep (outlier removal).

@author: brianyee
"""
import matplotlib
import matplotlib.patches as patches
matplotlib.use('Agg')
from matplotlib import rc

rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
rc('font', **{'family': 'sans-serif', 'sans-serif': ['Helvetica']})

import os
import logging
import pandas as pd
import matplotlib.pyplot as plt
import gzip

from collections import defaultdict, OrderedDict
import RDPlotter
import matrix as mtx
import misc
import normalization_functions as norm

SEP = '.'  # file delimiter


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
        annotation : dict
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
        self.means = defaultdict()
        self.sems = defaultdict()
        self.num_events = defaultdict()

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

    def set_means_and_sems(self):
        """
        For each position, remove the outlying 5% (or as specified in
        self.conf) and compute the average and mean

        Returns
        -------

        """
        means = defaultdict()
        sems = defaultdict()

        for filename, filetype in self.annotation.iteritems():
            means[filename], sems[filename] = norm.get_means_and_sems(
                self.raw_matrices[filename], self.conf
            )
        self.means = means
        self.sems = sems

    def write_intermediate_raw_matrices_to_csv(self):
        """
        Writes to output_base + ip/input + annotation + raw_density.csv
        the raw density values for each annotation/rbp matrix created.

        Returns
        -------

        """
        for condition, annotations in self.raw_matrices.iteritems():
            for filename, matrix in annotations.iteritems():
                matrix.to_csv(
                    self.output_base + SEP +
                    os.path.splitext(os.path.basename(filename))[0] + SEP +
                    condition + SEP +
                    'raw_density.csv'
                )

    def write_intermediate_means_to_csv(self):
        """
        Writes to output_base + annotation + norm_func + outlier + mean.csv
        the mean normalized density values over input for each
        annotation/rbp matrix created

        Returns
        -------

        """
        for filename, mean in self.means.iteritems():
            pd.Series(mean['means']).to_csv(
                self.output_base + SEP +
                os.path.splitext(os.path.basename(filename))[0] + SEP +
                '{}'.format(self.norm_function.__name__) + SEP +
                '{}.mean.csv'.format(self.conf)
            )

    def write_intermediate_sems_to_csv(self):
        """
        Writes to output_base + annotation_base + norm_func + sem.csv
        the standard error values for each position in the map.

        Returns
        -------

        """
        for filename, sem in self.sems.iteritems():
            pd.Series(sem).to_csv(
                self.output_base + SEP +
                os.path.splitext(os.path.basename(filename))[0] + SEP +
                '{}'.format(self.norm_function.__name__) + SEP +
                '{}.sem.csv'.format(self.conf)
            )

    def write_intermediates_to_csv(self):
        """
        Writes all intermediate files.

        Returns
        -------

        """
        self.write_intermediate_raw_matrices_to_csv()
        self.write_intermediate_means_to_csv()
        self.write_intermediate_sems_to_csv()


    def compute_num_events(self):
        """
        Prints the number of events for each annotation used.

        Returns
        -------

        """
        num_events = dict()
        for filename, matrix in self.raw_matrices['ip'].iteritems():
            num_events[misc.sane(filename)] = matrix.shape[0]
        print(num_events)


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

    def write_intermediate_norm_matrices_to_csv(self):
        for filename, matrix in self.norm_matrices.iteritems():
            matrix.to_csv(
                self.output_base + SEP +
                os.path.splitext(os.path.basename(filename))[0] + SEP +
                '{}.csv'.format(self.norm_function.__name__)
            )

    def write_intermediates_to_csv(self):
        self.write_intermediate_raw_matrices_to_csv()
        self.write_intermediate_norm_matrices_to_csv()
        self.write_intermediate_means_to_csv()
        self.write_intermediate_sems_to_csv()
        self.export_as_deeptool_matrix()

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
        matrices = defaultdict()
        for filename, filetype in self.annotation.iteritems():
            matrices[filename] = self.norm_function(
                self.raw_matrices['ip'][filename],
                self.raw_matrices['input'][filename],
                self.ip.pseudocount(), self.inp.pseudocount(),
                self.min_density_threshold
            )
        self.norm_matrices = matrices

    def set_means_and_sems(self):
        means = OrderedDict()
        sems = OrderedDict()
        i = 0

        for filename, filetype in self.annotation.iteritems():
            i+=1
            name = "{}".format(misc.sane(filename))
            means[name] = defaultdict()

            means[name]['means'], sems[name] = norm.get_means_and_sems(
                self.norm_matrices[filename],
                self.conf
            )
            means[name]['nums'] = self.norm_matrices[filename].shape[0]
        self.means = means
        self.sems = sems

    def plot(self):
        """
        Determines whether or not to plot a 1 or 2-window map.

        Returns
        -------

        """
        if self.is_scaled:
            self.plot_as_bed()
        else:
            self.plot_as_exon()

    def plot_as_bed(self):
        """
        Plots the entire region in one window. This function is mostly for
        single-base regions such as CDS start sites that are slopped (flanked)
        by a fixed length on either side.

        :return:
        """
        f, ax = plt.subplots(figsize=(10, 5))
        RDPlotter.plot_bed(self.means, self.sems, ax)
        f.suptitle(misc.sane(self.output_filename))
        f.savefig(self.output_filename)

    def plot_as_exon(self):
        """
        Plots the region as two windows. This function is mostly for
        regions of unequal length (such as exons) that are not scaled
        and are plotting 3' and 5' ends as two distinct regions.

        :return:
        """
        f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10, 5))
        axs = [ax1, ax2]
        RDPlotter.plot_exon(self.means, self.sems, axs)
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

        Returns
        -------

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
        f, (ax1, ax2, ax3, ax4) = plt.subplots(
            1, 4, sharey=True, figsize=(16, 8)
        )
        axs = [ax1, ax2, ax3, ax4]
        RDPlotter.plot_se(self.means, self.sems, axs)

        plt.tight_layout(pad=1.5 * len(self.annotation.keys()), w_pad=1)
        # plt.tight_layout(pad=5.5, w_pad=3, h_pad=6)
        f.suptitle(
            misc.sane(self.output_filename).replace(
                '.merged.r2',''
            ).split('_')[-1]
        )

        f.savefig(self.output_filename)


class MutuallyExclusiveExon(WithInput):
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
        f, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(
            1, 6, sharey=True, figsize=(16, 8)
        )
        axs = [ax1, ax2, ax3, ax4, ax5, ax6]
        RDPlotter.plot_splice(self.means, self.sems, axs)
        plt.tight_layout(pad=1.5 * len(self.annotation.keys()), w_pad=1)
        f.suptitle(misc.sane(self.output_filename))
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
        f, (ax1, ax2, ax3) = plt.subplots(
            1, 3, sharey=True, figsize=(16, 8)
        )
        axs = [ax1, ax2, ax3]
        RDPlotter.plot_a3ss(self.means, self.sems, axs)
        plt.tight_layout(pad=1.5 * len(self.annotation.keys()), w_pad=1)
        f.suptitle(misc.sane(self.output_filename))
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
        f, (ax1, ax2, ax3) = plt.subplots(
            1, 3, sharey=True, figsize=(16, 8)
        )
        axs = [ax1, ax2, ax3]
        RDPlotter.plot_a5ss(self.means, self.sems, axs)
        plt.tight_layout(pad=1.5 * len(self.annotation.keys()), w_pad=1)
        f.suptitle(misc.sane(self.output_filename))
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
        f, (ax1, ax2, ) = plt.subplots(
            1, 2, sharey=True, figsize=(16, 8)
        )
        axs = [ax1, ax2]
        RDPlotter.plot_ri(self.means, self.sems, axs)
        plt.tight_layout(pad=1.5 * len(self.annotation.keys()), w_pad=1)
        f.suptitle(misc.sane(self.output_filename))
        f.savefig(self.output_filename)


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
