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

class Plot(object):
    '''
    classdocs
    '''

    def __init__(self, RBP, output_file, line_color, map_type):
        '''
        Constructor
        '''
        self.rbp = RBP
        self.output_file = output_file
        self.color = line_color
        self.map_type = map_type
        
        if map_type == 'se':
            self.rbp.create_se_matrices()
        else:
            self.rbp.create_matrices()

    def single_frame(self):
        
        ax = plt.gca()
        means = self.rbp.normalize()
        
        ax.plot(means, color = self.color, label = 'Mean Read Density')
        ax.legend()

        ax.set_ylabel('Read Density')
        ax.set_title(self.map.get_name(),y=1.03)
        
        ymax = max(means) * 1.1
        ymin = min(means) * 0.9
        
        ax.set_ylim([ymin,ymax])
        plt.savefig(self.output_file)
        ax.clear()
    
    def single_frame_with_error(self):
        
        ax = plt.gca()
        pdf = self.rbp.normalize()
        
        means = pdf.mean()
        error = pdf.sem()
        
        ax.plot((means+error), color = sns.color_palette("hls", 8)[0], alpha = 0.3, label = 'Standard Error')
        ax.plot(means, color = self.color, label = 'Mean Read Density')
        ax.plot((means-error), color = sns.color_palette("hls", 8)[0], alpha = 0.3)
        ax.legend()

        ax.set_ylabel('Read Density')
        ax.set_title(self.map.get_name(),y=1.03)
        
        ymax = max(means) * 1.1
        ymin = min(means) * 0.9
        
        ax.set_ylim([ymin,ymax])
        plt.savefig(self.output_file)
        ax.clear()
           
    def four_frame(self):
        num_rows = 1
        num_cols = 4
        color = 'blue'
        
        region1, region2, region3, region4 = 0
        
        
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