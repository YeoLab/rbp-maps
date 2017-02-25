class DensityPlotter():
    def __init__(self, clip_with_input_map):
        """

        Parameters
        ----------
        map : density.ClipWithInput
        """
        self.map = clip_with_input_map
        self.line = self.get_mean()
        self.plus_error_bar = self.line + self.get_sems()
        self.minus_error_bar = self.line - self.get_sems()

    def get_mean(self):
        return self.map.get_means()

    def get_sems(self):
        return self.map.get_sems()

    def _matplotlib(self, ax=None):
        if ax is None:
            ax = plt.gca()
        # plot the line
        ax.plot(self.line, color = 'red')
        # upper error bar
        ax.plot(self.plus_error_bar, color = 'red')
        # lower error bar
        ax.plot(self.minus_error_bar, color = 'red')

class SEPlotter():
    def __init__(self, clip_with_input_map):
        """

        Parameters
        ----------
        map : density.ClipWithInput
        """
        self.map = clip_with_input_map
        self.full_line = self.get_mean()
        self.plus_error_bar = self.line + self.get_sems()
        self.minus_error_bar = self.line - self.get_sems()
        self._split(4)

        def _split(num):
            """
            Splits a pandas.Series into num equal parts.

            Parameters
            ----------
            num : pandas.Series

            Returns
            -------

            """

            parts = defaultdict()
            for i in range(0, num):
                parts[i]
def plot_se():
    pass