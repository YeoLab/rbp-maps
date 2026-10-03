import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

STYLES = {
    'total': dict(color='k', lw=1.2),
    'sense': dict(color='tab:blue', lw=1),
    'antisense': dict(color='tab:red', lw=1),
}


def plot_stranded_metagene(profile, output_filename, title, control=None, control_label='control'):
    """
    Plots the mean sense, antisense and total signal around the TSS.

    Total and sense share the left axis. Antisense is drawn on its own
    right-hand axis, whose zero is aligned with the zero of the left axis.

    Parameters
    ----------
    profile : pandas.DataFrame
        columns sense, antisense and total, indexed by position relative to
        the TSS (see: density.tss.profiles)
    output_filename : str
    title : str
    control : pandas.DataFrame
        profile of a control window set, drawn dashed (optional)
    control_label : str
    """
    with plt.style.context('default'):  # not the large fonts that plotter.Plotter sets globally
        _plot(profile, output_filename, title, control, control_label)


def _plot(profile, output_filename, title, control, control_label):
    fig, ax = plt.subplots(figsize=(8, 4.5))
    right = ax.twinx()
    right.tick_params(axis='y', colors='tab:red')
    for orientation in ['total', 'sense', 'antisense']:
        target = right if orientation == 'antisense' else ax
        target.plot(profile.index, profile[orientation], label=orientation, **STYLES[orientation])
        if control is not None:
            target.plot(
                control.index, control[orientation], ls='--', alpha=0.7,
                label='{} {}'.format(orientation, control_label), **STYLES[orientation]
            )
    ax.axhline(0, color='grey', lw=0.6)
    ax.axvline(0, ls=':', color='k')

    # align 0 of the antisense axis with 0 of the left axis
    low, high = ax.get_ylim()
    below = -low / (high - low)
    anti_low, anti_high = right.get_ylim()
    if 0 < below < 1:
        span = max(-anti_low / below, anti_high / (1 - below))
        right.set_ylim(-below * span, (1 - below) * span)

    handles, labels = ax.get_legend_handles_labels()
    right_handles, right_labels = right.get_legend_handles_labels()
    ax.legend(handles + right_handles, labels + right_labels, fontsize=8)
    ax.set_xlabel('Position relative to TSS (nt)')
    ax.set_ylabel('Mean normalized IP/input signal (total, sense)')
    right.set_ylabel('Mean normalized IP/input signal (antisense)', color='tab:red')
    ax.set_title(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(output_filename, dpi=150)
    plt.close(fig)
