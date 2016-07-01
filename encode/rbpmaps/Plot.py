'''
Created on Jun 20, 2016

@author: brianyee
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from gscripts.general import dataviz


def single_frame(means, title, output_file, color='red'):
        
    ax = plt.gca()
        
    ax.plot(means, color = color, label = 'Mean Read Density')
    ax.legend()

    ax.set_ylabel('Read Density')
    ax.set_title(title,y=1.03)
        
    ymax = max(means) * 1.1
    ymin = min(means) * 0.9
        
    ax.set_ylim([ymin,ymax])
    plt.savefig(output_file)
    ax.clear()
    
def single_frame_with_error(means, error, title, output_file, color='red'):
        
    ax = plt.gca()
        
    ax.plot((means+error), color = sns.color_palette("hls", 8)[0], alpha = 0.3, label = 'Standard Error')
    ax.plot(means, color = color, label = 'Mean Read Density')
    ax.plot((means-error), color = sns.color_palette("hls", 8)[0], alpha = 0.3)
    ax.legend()

    ax.set_ylabel('Read Density')
    ax.set_title(title,y=1.03)
        
    ymax = max(means) * 1.1
    ymin = min(means) * 0.9 if min(means) > 0 else min(means)*1.1 # in case of negatives for subtraction
        
    ax.set_ylim([ymin,ymax])
    plt.savefig(output_file)
    ax.clear()

def three_frame(region1, region2, region3, 
               title, output_file, color='red'):
    num_rows = 1
    num_cols = 4
    color = 'blue'
    
        
    with dataviz.Figure(output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
            
        min_height = min(min(region1),min(region2),min(region3))
        max_height = max(max(region1),max(region2),max(region3))
        
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

    ax.clear()
    
def four_frame(region1, region2, region3, region4, 
               title, output_file, color='red'):
    num_rows = 1
    num_cols = 4
    color = 'blue'
    
        
    with dataviz.Figure(output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
            
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
        plt.suptitle(title,y=1.03)
    ax.clear()
    
def five_frame(region1, region2, region3, region4, region5,
               title, output_file, color='red'):
    num_rows = 1
    num_cols = 5
    color = 'blue'
    
        
    with dataviz.Figure(output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
            
        min_height = min(min(region1),min(region2),min(region3),min(region4))
        max_height = max(max(region1),max(region2),max(region3),max(region4))
            
        linewidth = 2.5
        ax = fig.add_subplot(1,5,1)
        ax.plot(region1, linewidth=linewidth, alpha=.7, color = color)
        # ax.plot(three_upstream_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        # ax.set_xticklabels(np.arange(-exon_offset, intron_offset+1, 50))
        ax.set_ylabel("Mean Read Density")
            
        ax = fig.add_subplot(1,5,2)
        ax.plot(region2, linewidth=linewidth, alpha=.7, color = color)
        # ax.plot(five_skipped_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
            
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        # ax.set_xticklabels(np.arange(-intron_offset, exon_offset+1, 50))
        ax.set_yticklabels([])
            
        ax = fig.add_subplot(1,5,3)
        ax.plot(region3, linewidth=linewidth, alpha=.7, color = color)
        # ax.plot(three_skipped_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
            
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        # ax.set_xticklabels(np.arange(-exon_offset, intron_offset+1, 50))
        ax.set_yticklabels([])
            
        ax = fig.add_subplot(1,5,4)
        ax.plot(region4, linewidth=linewidth, alpha=.7, color = color)
        # ax.plot(five_downstream_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
            
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        # ax.set_xticklabels(np.arange(-intron_offset, exon_offset+1, 50))
        ax.set_yticklabels([])
        
        ax = fig.add_subplot(1,5,5)
        ax.plot(region5, linewidth=linewidth, alpha=.7, color = color)
        # ax.plot(five_downstream_normed_nt, linewidth=linewidth, alpha=.7, color = 'blue')
            
        sns.despine(ax=ax, left=True)
        # ax.set_ylim(min_height, max_height)
        # ax.set_xticklabels(np.arange(-intron_offset, exon_offset+1, 50))
        ax.set_yticklabels([])
        
        plt.suptitle(title,y=1.03)
    ax.clear()
    