import normalization_functions as norm
import intervals
import matrix as mtx
from scipy import stats

import os
import numpy as np


class LineObject():
    def __init__(
            self, event_matrix, annotation, conf,
            color, min_event_threshold
    ):
        self.event_matrix = event_matrix
        self.conf = conf
        self.annotation = annotation  # annotation file
        self.num_events = self.event_matrix.shape[0]
        self.label = self._parse_filename_for_plot()
        self.file_label = self._parse_filename()
        self.means, self.sems, self.std = self._get_means_and_sems()  # 2 lists
        self.error_pos, self.error_neg = self._get_std_error_boundaries()
        self.dim = False if self.num_events > min_event_threshold else True
        self.ks_pvalues = []
        self.z_scores = []
        self.color = color

        print(self.event_matrix.shape, self.label)


    def _parse_filename_for_plot(self):
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
            os.path.basename(self.annotation)
        )[0].replace(
            '-', ' ').replace(
            '_', ' ').replace(
            'HepG2', '').replace(
            'K562', '') + " ({} events)".format(
            self.num_events
        )

    def _parse_filename(self):
        return os.path.splitext(
            os.path.basename(self.annotation)
        )[0]

    def _get_means_and_sems(self):
        """
        Sets the means and standard error values after outlier
        removal. Replaces remove_outliers.

        Parameters
        ----------
        df : pandas.DataFrame
            table of densities or values
        conf : float
            keep {conf}% of densities present at every given position

        Returns
        -------

        means : list
            mean value for each position in the dataframe df
        sems : list
            standard error of the mean
        """

        means = list()
        sems = list()
        std_deviation = list()

        for key, value in self.event_matrix.iteritems():
            single_col = self.event_matrix[key].dropna()
            single_col = single_col.sort_values()
            nums = len(single_col)
            droppercent = (1 - self.conf) / 2.0
            dropnum = int(nums * droppercent)
            if dropnum > 0:
                single_col = single_col[dropnum:-dropnum]

            means.append(single_col.mean())
            sems.append(single_col.sem())
            std_deviation.append(single_col.std())
        return means, sems, std_deviation

    def _get_std_error_boundaries(self):
        """
        Returns the +/- error boundaries given a list of values (means)

        Parameters
        ----------
        means
        error

        Returns
        -------

        """
        pos = [x + y for x, y in zip(self.means, self.sems)]
        neg = [x - y for x, y in zip(self.means, self.sems)]
        return pos, neg

    def calculate_and_set_significance(self, bg_matrix):
        """
        Given a background event matrix, compute distribution
        and calculate 2-sample KS test
        
        Parameters
        ----------
        bg_matrix : pandas.DataFrame()
            a position matrix (event = row, positon = col)

        Returns
        -------
        list of -log10 p-values for each position
        
        """
        for position in self.event_matrix.columns:
            _, p = stats.ks_2samp(
                self.event_matrix[position], bg_matrix[position]
            )
            print('p value at position {}: {}'.format(position, p))

            self.ks_pvalues.append(-1*np.log10(p))

    def calculate_zscore(self, bg_line):
        """
        This is really messy but we are calculating the zscore
        of the outlier-removed values for both test and bg

        Parameters
        ----------
        bg_matrix

        Returns
        -------

        """

        # make sure the bg matrix is computing z scores
        # for as many positions as we have in our Line
        for i in range(0, len(self.means)):
            z_score = (self.means[i] - bg_line.means[i])/bg_line.std[i]
            # print('zscore: {}'.format(z_score))
            self.z_scores.append(z_score)

