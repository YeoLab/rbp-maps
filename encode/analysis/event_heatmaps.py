'''
Created on Jan 9, 2017

@author: brianyee
'''

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def heatmap_all_events(infile, outfile, xlabel=False, ylabel=False):
    """
    Plots a heatmap of all events of the dataframe
    """
    df = pd.read_table(infile,index_col=0)
    ax = sns.heatmap(df, xticklabels=xlabel, yticklabels=ylabel)
    plt.title("Heatmap All Events: {}")
    plt.savefig(outfile)

def main():
    pass
if __name__ == '__main__':
    main()