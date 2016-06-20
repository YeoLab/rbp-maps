'''
Created on Jun 20, 2016

@author: brianyee
'''
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from gscripts.general import dataviz
import os

class Mplot(object):
    '''
    classdocs
    '''

    def __init__(self, Map, output_file, line_color):
        '''
        Constructor
        '''
        self.map = Map
        self.output_file = output_file
        self.color = line_color
        
    def single_frame(self):
        
        ax = plt.gca()
        raw, normed, means = self.map.get_matrices()
        normed.to_csv('{}.normed_density_matrix.csv'.format(os.path.splitext(self.output_file)[0]))
        means.to_csv('{}.allmeans.txt'.format(os.path.splitext(self.output_file)[0]))
        raw.to_csv('{}.raw_density_matrix.csv'.format(os.path.splitext(self.output_file)[0]))
        
        ax.plot(means, color = self.color)
        plt.savefig(self.output_file)
    
    def four_frame(self):
        num_rows = 1
        num_cols = 4
        color = 'blue'
        matrices = self.map.get_matrices()
        
        region1_raw, region2_raw, region3_raw, region4_raw = matrices['raw']
        region1_normed, region2_normed, region3_normed, region4_normed = matrices['normed']
        region1, region2, region3, region4 = matrices['means']
        
        all_regions = pd.concat([region1_normed,region2_normed,region3_normed,region4_normed])
        
        region1_raw.to_csv("{}_region1_raw_density_matrix.csv".format(os.path.splitext(self.output_file)[0]))
        region2_raw.to_csv("{}_region2_raw_density_matrix.csv".format(os.path.splitext(self.output_file)[0]))
        region3_raw.to_csv("{}_region3_raw_density_matrix.csv".format(os.path.splitext(self.output_file)[0]))
        region4_raw.to_csv("{}_region4_raw_density_matrix.csv".format(os.path.splitext(self.output_file)[0]))
            
        region1_normed.to_csv("{}_region1_normed_pdf.csv".format(os.path.splitext(self.output_file)[0]))
        region2_normed.to_csv("{}_region2_normed_pdf.csv".format(os.path.splitext(self.output_file)[0]))
        region3_normed.to_csv("{}_region3_normed_pdf.csv".format(os.path.splitext(self.output_file)[0]))
        region4_normed.to_csv("{}_region4_normed_pdf.csv".format(os.path.splitext(self.output_file)[0]))
        
        all_regions.to_csv(os.path.splitext(self.output_file)[0]+'.allmeans.txt')
        
        with dataviz.Figure(self.output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
            
            min_height = min(min(region1),min(region2),min(region3),min(region4))
            max_height = max(max(region1),max(region2),max(region3),max(region4))
            
            linewidth = 2.5
            ax = fig.add_subplot(1,4,1)
            ax.plot(region1, linewidth=linewidth, alpha=.7, color = color)
            # ax.plot(three_upstream_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
            sns.despine(ax=ax)
            ax.set_ylim(min_height, max_height)
            # ax.set_xticklabels(np.arange(-exon_offset, intron_offset+1, 50))
            ax.set_ylabel("Mean Read Density")
            
            ax = fig.add_subplot(1,4,2)
            ax.plot(region2, linewidth=linewidth, alpha=.7, color = color)
            # ax.plot(five_skipped_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
            
            sns.despine(ax=ax, left=True)
            ax.set_ylim(min_height, max_height)
            # ax.set_xticklabels(np.arange(-intron_offset, exon_offset+1, 50))
            ax.set_yticklabels([])
            
            ax = fig.add_subplot(1,4,3)
            ax.plot(region3, linewidth=linewidth, alpha=.7, color = color)
            # ax.plot(three_skipped_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
            
            sns.despine(ax=ax, left=True)
            ax.set_ylim(min_height, max_height)
            # ax.set_xticklabels(np.arange(-exon_offset, intron_offset+1, 50))
            ax.set_yticklabels([])
            
            ax = fig.add_subplot(1,4,4)
            ax.plot(region4, linewidth=linewidth, alpha=.7, color = color)
            # ax.plot(five_downstream_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
            
            sns.despine(ax=ax, left=True)
            ax.set_ylim(min_height, max_height)
            # ax.set_xticklabels(np.arange(-intron_offset, exon_offset+1, 50))
            ax.set_yticklabels([])
            plt.suptitle(self.map.get_name(),y=1.03)