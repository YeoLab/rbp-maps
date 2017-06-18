import normalization_functions as norm
import intervals
import matrix as mtx
from scipy import stats
import os
import numpy as np


class LineObject():
    def __init__(
            self, infile, out_hist, annotation,
            l10p_cutoff, l2fc_cutoff, hashing_val,
            event_type, exon_overhang, intron_overhang,
            color, min_event_threshold
    ):
        self.event_dict = intervals.read_four_region_miso(
            annotation, hashing_val, event_type,
            exon_overhang, intron_overhang
        )
        self.peak_hist = mtx.make_hist_se(
            infile, out_hist, hashing_val,
            l10p_cutoff, l2fc_cutoff, self.event_dict,
            exon_overhang, intron_overhang
        )
        self.num_events = sum(1 for line in open(annotation))

        self.label = self._parse_filename(annotation)
        self.means = norm.norm(self.peak_hist, self.num_events)
        self.error_pos, self.error_neg = self.set_std_error_boundaries()
        self.dim = False if self.num_events > min_event_threshold else True
        self.fisher_pvalues = []
        self.color = color

    def _parse_filename(self, annotation):
        """
        Removes any extra stuff that we might not want in the name.
        
        Parameters
        ----------
        annotation : basestring

        Returns
        -------
        nice_label : basestring
        """
        return os.path.splitext(
            os.path.basename(annotation)
        )[0].replace(
            '-', ' ').replace(
            '_', ' ').replace(
            'HepG2','').replace(
            'K562','') + " ({} events)".format(
            self.num_events
        )

    def set_fisher(self, bg_peak):
        """
        Calculates the fisher exact test for each position compared
        to a background

        Parameters
        ----------
        bg_peak: Peak.Peak()

        Returns
        -------
        p-values for each position.
        """
        p_values = []
        for i in range(0, len(self.peak_hist)): # for each i, return number of peaks at position i
            total_peak = self.peak_hist[i] # total peaks at position i
            total_without_peak = self.num_events - total_peak
            bg_total_peak = bg_peak.peak_hist[i]
            bg_total_without_peak = bg_peak.num_events - bg_total_peak
            contingency_table = [
                [total_peak, total_without_peak],
                [bg_total_peak, bg_total_without_peak]
            ]

            odds, p = stats.fisher_exact(contingency_table)
            p_values.append(p)
            # if 300 < i and 307 > i:
            #     print(contingency_table, self.label, i, p)
        self.fisher_pvalues = -1*np.log10(p_values)

    def set_std_error_boundaries(self):
        """
        Sets the std error upper/lower boundaries given a mean and standard error 
        
        Parameters
        ----------
        hist
        n

        Returns
        -------

        """
        plus = [x + y for x, y in zip(
            norm.norm(self.peak_hist, self.num_events), norm.std_error(
                self.peak_hist, self.num_events)
        )]
        minus = [x - y for x, y in zip(
            norm.norm(self.peak_hist, self.num_events), norm.std_error(
                self.peak_hist, self.num_events)
        )]
        return plus, minus