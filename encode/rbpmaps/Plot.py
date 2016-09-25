'''
Created on Jun 20, 2016

@author: brianyee
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from gscripts.general import dataviz
from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
import numpy as np

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

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
    plt.clf()
    plt.cla()
    plt.close()

def single_frame_with_inclusion_exclusion_events(inclusion, exclusion, both, 
                                                 title, output_file):
    ax = plt.gca()
    ax.plot(inclusion['region1'], color = sns.color_palette("hls", 8)[0], label = 'Inclusion')
    ax.plot(exclusion['region1'], color = sns.color_palette("hls", 8)[1], label = 'Exclusion')
    ax.plot(both['region1'], color = sns.color_palette("hls", 8)[2], label = 'All events')
    ax.legend()

    ax.set_ylabel('Read Density')
    ax.set_title(title,y=1.03)
    # plt.xticks([0,300,400,699],['upstream (300bp)','feature (0%)','feature (100%)','downstream (300bp)'])
    # ymax = max(means) * 1.1
    # ymin = min(means) * 0.9 if min(means) > 0 else min(means)*1.1 # in case of negatives for subtraction
        
    # ax.set_ylim([ymin,ymax])
    plt.savefig(output_file)
    plt.clf()
    plt.cla()
    plt.close()
    
def single_frame_with_up_down_events_error(up, down, both,
                                           uperr, downerr, 
                                           title, output_file,
                                           upcolor = sns.color_palette("hls", 8)[0],
                                           downcolor = sns.color_palette("hls", 8)[1],
                                           bgcolor = 'black',
                                           uplab = 'Upregulated',
                                           downlab = 'Downregulated',
                                           bglab = 'All significant events'):
    ax = plt.gca()
    
    errorbar_linewidth = 0.7
    
    ax.plot(up['region1'], color = upcolor, label = uplab)
    ax.plot((up['region1']+uperr['region1']), linewidth=errorbar_linewidth, alpha=.5, color = upcolor, linestyle = ':')
    ax.plot((up['region1']-uperr['region1']), linewidth=errorbar_linewidth, alpha=.5, color = upcolor, linestyle = ':')
    
    ax.plot(down['region1'], color = downcolor, label = downlab)
    ax.plot((down['region1']+downerr['region1']), linewidth=errorbar_linewidth, alpha=.5, color = downcolor, linestyle = ':')
    ax.plot((down['region1']-downerr['region1']), linewidth=errorbar_linewidth, alpha=.5, color = downcolor, linestyle = ':')
    
    ax.plot(both['region1'], color = bgcolor, label = bglab)
    
    ax.legend()

    ax.set_ylabel('Read Density')
    ax.set_title(title,y=1.03)
    # plt.xticks([0,300,400,699],['upstream (300bp)','feature (0%)','feature (100%)','downstream (300bp)'])
    # ymax = max(means) * 1.1
    # ymin = min(means) * 0.9 if min(means) > 0 else min(means)*1.1 # in case of negatives for subtraction
        
    # ax.set_ylim([ymin,ymax])
    plt.savefig(output_file)
    plt.clf()
    plt.cla()
    plt.close()

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
    plt.clf()
    plt.cla()
    plt.close()

def two_frame_with_inclusion_exclusion_events_with_error(inclusion, exclusion, both,
                                                          inclusion_err, exclusion_err,
                                                          title, output_file, 
                                                          color1=sns.color_palette("hls", 8)[0], 
                                                          color2=sns.color_palette("hls", 8)[5],
                                                          color3='black'):
    
    """
    Special plot:
    plots a 4-region map that contains three separate plots for inclusion, 
    exclusion, and all spliced events. 
    
    Args:
        inclusion: {region1, region2, region3, region4}
        exclusion: {region1, region2, region3, region4}
        both: {region1, region2, region3, region4}
        
    """

    num_rows = 1
    num_cols = 2
    ax = plt.gca()
    min_height = min(min(inclusion['region1']),min(exclusion['region1']),min(both['region1']),
                     min(inclusion['region2']),min(exclusion['region2']),min(both['region2']))
    max_height = max(max(inclusion['region1']),max(exclusion['region1']),max(both['region1']),
                     max(inclusion['region2']),max(exclusion['region2']),max(both['region2']))
        
    with dataviz.Figure(output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
        
        sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})
        
        linewidth = 2
        errorbar_linewidth = 0.7
        
        ax = fig.add_subplot(1,2,1)
        ax.plot(inclusion['region1'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(exclusion['region1'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(both['region1'], linewidth=linewidth, alpha=.3, color = color3)
        
        ax.plot((inclusion['region1']+inclusion_err['region1']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((inclusion['region1']-inclusion_err['region1']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((exclusion['region1']+exclusion_err['region1']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((exclusion['region1']-exclusion_err['region1']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_ylabel("Mean Read Density")
        ax.set_xticklabels(range(-50,351,50),rotation=90)
        ax.axvline(x=50,linestyle=':',alpha=0.5)
        
        sns.set_style({'ytick.major.size':0})
        
        ax = fig.add_subplot(1,2,2)
        ax.plot(inclusion['region2'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(exclusion['region2'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(both['region2'], linewidth=linewidth, alpha=.3, color = color3)
        
        ax.plot((inclusion['region2']+inclusion_err['region2']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((inclusion['region2']-inclusion_err['region2']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((exclusion['region2']+exclusion_err['region2']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((exclusion['region2']-exclusion_err['region2']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.set_xticklabels(range(-300,51,50),rotation=90)
        ax.axvline(x=300,linestyle=':',alpha=0.5)

        ax.legend()
        plt.suptitle(title,y=1.03)
    plt.clf()
    plt.cla()
    plt.close()
def three_frame(region1, region2, region3, 
               title, output_file, color='red'):
    num_rows = 1
    num_cols = 4
    
    ax = plt.gca()
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

    plt.clf()
    plt.cla()
    plt.close()
    
def four_frame(region1, region2, region3, region4, 
               title, output_file, color='red'):
    num_rows = 1
    num_cols = 4
    
    ax = plt.gca()
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
    plt.clf()
    plt.cla()
    plt.close()


def four_frame_with_inclusion_exclusion_events_with_error(inclusion, exclusion, both,
                                                          inclusion_err, exclusion_err,
                                                          title, output_file, 
                                                          color1=sns.color_palette("hls", 8)[0], 
                                                          color2=sns.color_palette("hls", 8)[5],
                                                          color3='black'):
    
    """
    Special plot:
    plots a 4-region map that contains three separate plots for inclusion, 
    exclusion, and all spliced events. 
    
    Args:
        inclusion: {region1, region2, region3, region4}
        exclusion: {region1, region2, region3, region4}
        both: {region1, region2, region3, region4}
        
    """
    num_rows = 1
    num_cols = 4
    ax = plt.gca()
    min_height = min(min(inclusion['region1']),min(exclusion['region1']),min(both['region1']),
                     min(inclusion['region2']),min(exclusion['region2']),min(both['region2']),
                     min(inclusion['region3']),min(exclusion['region3']),min(both['region3']),
                     min(inclusion['region4']),min(exclusion['region4']),min(both['region4']))
    max_height = max(max(inclusion['region1']),max(exclusion['region1']),max(both['region1']),
                     max(inclusion['region2']),max(exclusion['region2']),max(both['region2']),
                     max(inclusion['region3']),max(exclusion['region3']),max(both['region3']),
                     max(inclusion['region4']),max(exclusion['region4']),max(both['region4']))
        
    with dataviz.Figure(output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
        
        sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})
        
        linewidth = 2
        errorbar_linewidth = 0.7
        
        ax = fig.add_subplot(1,4,1)
        ax.plot(inclusion['region1'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(exclusion['region1'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(both['region1'], linewidth=linewidth, alpha=.3, color = color3)
        
        ax.plot((inclusion['region1']+inclusion_err['region1']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((inclusion['region1']-inclusion_err['region1']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((exclusion['region1']+exclusion_err['region1']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((exclusion['region1']-exclusion_err['region1']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_ylabel("Mean Read Density")
        ax.set_xticklabels(range(-50,351,50),rotation=90)
        ax.axvline(x=50,linestyle=':',alpha=0.5)
        
        sns.set_style({'ytick.major.size':0})
        
        ax = fig.add_subplot(1,4,2)
        ax.plot(inclusion['region2'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(exclusion['region2'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(both['region2'], linewidth=linewidth, alpha=.3, color = color3)
        
        ax.plot((inclusion['region2']+inclusion_err['region2']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((inclusion['region2']-inclusion_err['region2']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((exclusion['region2']+exclusion_err['region2']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((exclusion['region2']-exclusion_err['region2']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.set_xticklabels(range(-300,51,50),rotation=90)
        ax.axvline(x=300,linestyle=':',alpha=0.5)
            
        ax = fig.add_subplot(1,4,3)
        ax.plot(inclusion['region3'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(exclusion['region3'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(both['region3'], linewidth=linewidth, alpha=.3, color = color3)
        ax.plot((inclusion['region3']+inclusion_err['region3']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((inclusion['region3']-inclusion_err['region3']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((exclusion['region3']+exclusion_err['region3']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((exclusion['region3']-exclusion_err['region3']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.axvline(x=50,linestyle=':',alpha=0.5)
        ax.set_xticklabels(range(-50,351,50),rotation=90)
        
        ax = fig.add_subplot(1,4,4)
        ax.plot(inclusion['region4'], linewidth=linewidth, alpha=.8, color = color1, label="incl in kd")
        ax.plot(exclusion['region4'], linewidth=linewidth, alpha=.8, color = color2, label="excl in kd")
        ax.plot(both['region4'], linewidth=linewidth, alpha=.3, color = color3, label="background")
        ax.plot((inclusion['region4']+inclusion_err['region4']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((inclusion['region4']-inclusion_err['region4']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((exclusion['region4']+exclusion_err['region4']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((exclusion['region4']-exclusion_err['region4']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        ax.axvline(x=300,linestyle=':',alpha=0.5)
        ax.set_xticklabels(range(-300,51,50),rotation=90)
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.legend()
        plt.suptitle(title,y=1.03)
    plt.clf()
    plt.cla()
    plt.close()

def four_frame_with_inclusion_exclusion_events(inclusion, exclusion, both,
                                               title, output_file, 
                                               color1=sns.color_palette("hls", 8)[0], 
                                               color2=sns.color_palette("hls", 8)[5],
                                               color3='black'):
    """
    Special plot:
    plots a 4-region map that contains three separate plots for inclusion, 
    exclusion, and all spliced events. 
    
    Args:
        inclusion: {region1, region2, region3, region4}
        exclusion: {region1, region2, region3, region4}
        both: {region1, region2, region3, region4}
        
    """
    num_rows = 1
    num_cols = 4
    ax = plt.gca()
    min_height = min(min(inclusion['region1']),min(exclusion['region1']),min(both['region1']),
                     min(inclusion['region2']),min(exclusion['region2']),min(both['region2']),
                     min(inclusion['region3']),min(exclusion['region3']),min(both['region3']),
                     min(inclusion['region4']),min(exclusion['region4']),min(both['region4']))
    max_height = max(max(inclusion['region1']),max(exclusion['region1']),max(both['region1']),
                     max(inclusion['region2']),max(exclusion['region2']),max(both['region2']),
                     max(inclusion['region3']),max(exclusion['region3']),max(both['region3']),
                     max(inclusion['region4']),max(exclusion['region4']),max(both['region4']))
        
    with dataviz.Figure(output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
        
        sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})
        
        linewidth = 2.5
        ax = fig.add_subplot(1,4,1)
        ax.plot(inclusion['region1'], linewidth=linewidth, alpha=.7, color = color1)
        ax.plot(exclusion['region1'], linewidth=linewidth, alpha=.7, color = color2)
        ax.plot(both['region1'], linewidth=linewidth, alpha=.3, color = color3)
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_ylabel("Mean Read Density")
        ax.set_xticklabels(range(-50,351,50),rotation=90)
        ax.axvline(x=50,linestyle=':',alpha=0.5)
        
        sns.set_style({'ytick.major.size':0})
        
        ax = fig.add_subplot(1,4,2)
        ax.plot(inclusion['region2'], linewidth=linewidth, alpha=.7, color = color1)
        ax.plot(exclusion['region2'], linewidth=linewidth, alpha=.7, color = color2)
        ax.plot(both['region2'], linewidth=linewidth, alpha=.3, color = color3)
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.set_xticklabels(range(-300,51,50),rotation=90)
        ax.axvline(x=300,linestyle=':',alpha=0.5)
            
        ax = fig.add_subplot(1,4,3)
        ax.plot(inclusion['region3'], linewidth=linewidth, alpha=.7, color = color1)
        ax.plot(exclusion['region3'], linewidth=linewidth, alpha=.7, color = color2)
        ax.plot(both['region3'], linewidth=linewidth, alpha=.3, color = color3)
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.axvline(x=50,linestyle=':',alpha=0.5)
        ax.set_xticklabels(range(-50,351,50),rotation=90)
        
        ax = fig.add_subplot(1,4,4)
        ax.plot(inclusion['region4'], linewidth=linewidth, alpha=.7, color = color1, label="incl in kd")
        ax.plot(exclusion['region4'], linewidth=linewidth, alpha=.7, color = color2, label="excl in kd")
        ax.plot(both['region4'], linewidth=linewidth, alpha=.3, color = color3, label="background")
        ax.axvline(x=300,linestyle=':',alpha=0.5)
        ax.set_xticklabels(range(-300,51,50),rotation=90)
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.legend()
        plt.suptitle(title,y=1.03)
    plt.clf()
    plt.cla()
    plt.close()

def four_frame_with_inclusion_exclusion_events_from_one_region(inclusion, exclusion, both,
                                               title, output_file, 
                                               color1=sns.color_palette("hls", 8)[0], 
                                               color2=sns.color_palette("hls", 8)[5],
                                               color3='black'):
    """
    Special plot:
    plots a 4-region map that contains three separate plots for inclusion, 
    exclusion, and all spliced events. 
    
    Args:
        inclusion: {region1, region2, region3, region4}
        exclusion: {region1, region2, region3, region4}
        both: {region1, region2, region3, region4}
        
    """
    
    i = {}
    e = {}
    b = {}
    
    i['region1'], i['region2'], i['region3'], i['region4'] = np.array_split(inclusion['region1'],4)
    e['region1'], e['region2'], e['region3'], e['region4'] = np.array_split(exclusion['region1'],4)
    b['region1'], b['region2'], b['region3'], b['region4'] = np.array_split(both['region1'],4)
    
    four_frame_with_inclusion_exclusion_events(i,e,b,title,output_file,color1,color2,color3)

def plot_ri(inclusion, exclusion, both, inclusion_err, exclusion_err, title, output_file, 
            color1=sns.color_palette("hls", 8)[0], color2=sns.color_palette("hls", 8)[5], color3='black'):
    """
    Formerly: two_frame_with_inclusion_exclusion_events_from_one_region_with_error
    Special plot:
    
    plots a 2-region map that contains three separate plots for inclusion, 
    exclusion, and all spliced events. 
    
    Args:
        inclusion: {region1, region2, region3, region4}
        exclusion: {region1, region2, region3, region4}
        both: {region1, region2, region3, region4}
        
    """
    
    i = {}
    e = {}
    b = {}
    ie = {}
    ee = {}
    
    """
    Split the single region into four. This will break if the regions are of different dimensions!!!
    """
    i['region1'], i['region2'] = np.array_split(inclusion['region1'],2)
    e['region1'], e['region2'] = np.array_split(exclusion['region1'],2)
    b['region1'], b['region2'] = np.array_split(both['region1'],2)
    ie['region1'], ie['region2'] = np.array_split(inclusion_err['region1'],2)
    ee['region1'], ee['region2'] = np.array_split(exclusion_err['region1'],2)
    
    two_frame_with_inclusion_exclusion_events_with_error(i,e,b,ie,ee,title,output_file,color1,color2,color3)
    
def plot_se(inclusion, exclusion, both, inclusion_err, exclusion_err, title, output_file, 
            color1=sns.color_palette("hls", 8)[0], color2=sns.color_palette("hls", 8)[5], color3='black'):
    """
    Formerly: four_frame_with_inclusion_exclusion_events_from_one_region_with_error
    Special plot:
    
    plots a 4-region map that contains three separate plots for inclusion, 
    exclusion, and all spliced events. 
    
    Args:
        inclusion: {region1, region2, region3, region4}
        exclusion: {region1, region2, region3, region4}
        both: {region1, region2, region3, region4}
        
    """
    
    i = {}
    e = {}
    b = {}
    ie = {}
    ee = {}
    
    """
    Split the single region into four. This will break if the regions are of different dimensions!!!
    """
    i['region1'], i['region2'], i['region3'], i['region4'] = np.array_split(inclusion['region1'],4)
    e['region1'], e['region2'], e['region3'], e['region4'] = np.array_split(exclusion['region1'],4)
    b['region1'], b['region2'], b['region3'], b['region4'] = np.array_split(both['region1'],4)
    ie['region1'], ie['region2'], ie['region3'], ie['region4'] = np.array_split(inclusion_err['region1'],4)
    ee['region1'], ee['region2'], ee['region3'], ee['region4'] = np.array_split(exclusion_err['region1'],4)
    
    four_frame_with_inclusion_exclusion_events_with_error(i,e,b,ie,ee,title,output_file,color1,color2,color3)

    
def plot_a3ss(inclusion, exclusion, both, inclusion_err, exclusion_err, title, output_file, 
              color1=sns.color_palette("hls", 8)[0], color2=sns.color_palette("hls", 8)[5], color3='black'):
    
    i = {}
    e = {}
    b = {}
    ie = {}
    ee = {}
    i['three_upstream'] = np.array(inclusion['region1'][:350])
    i['five_alt1'] = np.array(inclusion['region1'][350:700])
    i['five_alt2'] = np.array(inclusion['region1'][700:])
    
    e['three_upstream'] = np.array(exclusion['region1'][:350])
    e['five_alt1'] = np.array(exclusion['region1'][350:700])
    e['five_alt2'] = np.array(exclusion['region1'][700:])
    
    b['three_upstream'] = np.array(both['region1'][:350])
    b['five_alt1'] = np.array(both['region1'][350:700])
    b['five_alt2'] = np.array(both['region1'][700:])
    
    ie['three_upstream'] = np.array(inclusion_err['region1'][:350])
    ie['five_alt1'] = np.array(inclusion_err['region1'][350:700])
    ie['five_alt2'] = np.array(inclusion_err['region1'][700:])
    
    ee['three_upstream'] = np.array(exclusion_err['region1'][:350])
    ee['five_alt1'] = np.array(exclusion_err['region1'][350:700])
    ee['five_alt2'] = np.array(exclusion_err['region1'][700:])
    
    """
    Special plot:
    plots a 4-region map that contains three separate plots for inclusion, 
    exclusion, and all spliced events. 
    
    Args:
        inclusion: {region1, region2, region3, region4}
        exclusion: {region1, region2, region3, region4}
        both: {region1, region2, region3, region4}
        
    """
    num_rows = 1
    num_cols = 3
    ax = plt.gca()
    min_height = min(min(i['three_upstream']),min(e['three_upstream']),min(b['three_upstream']),
                     min(i['five_alt1']),min(e['five_alt1']),min(b['five_alt1']),
                     min(i['five_alt2']),min(e['five_alt2']),min(b['five_alt2']))
    max_height = max(max(i['three_upstream']),max(e['three_upstream']),max(b['three_upstream']),
                     max(i['five_alt1']),max(e['five_alt1']),max(b['five_alt1']),
                     max(i['five_alt2']),max(e['five_alt2']),max(b['five_alt2']))
    
    with dataviz.Figure(output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
        
        sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})
        
        linewidth = 2
        errorbar_linewidth = 0.7
        
        ax = fig.add_subplot(1,4,1)
        ax.plot(i['three_upstream'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(e['three_upstream'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(b['three_upstream'], linewidth=linewidth, alpha=.3, color = color3)
        
        ax.plot((i['three_upstream']+ie['three_upstream']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((i['three_upstream']-ie['three_upstream']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((e['three_upstream']+ee['three_upstream']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((e['three_upstream']-ee['three_upstream']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_ylabel("Mean Read Density")
        ax.set_xticklabels(range(-50,351,50),rotation=90)
        ax.axvline(x=50,linestyle=':',alpha=0.5)
        
        sns.set_style({'ytick.major.size':0})
        
        ax = fig.add_subplot(1,4,2)
        ax.plot(i['five_alt1'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(e['five_alt1'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(b['five_alt1'], linewidth=linewidth, alpha=.3, color = color3)
        
        ax.plot((i['five_alt1']+ie['five_alt1']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((i['five_alt1']-ie['five_alt1']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((e['five_alt1']+ee['five_alt1']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((e['five_alt1']-ee['five_alt1']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.set_xticklabels(range(-300,51,50),rotation=90)
        ax.axvline(x=300,linestyle=':',alpha=0.5)
            
        ax = fig.add_subplot(1,4,3)
        ax.plot(i['five_alt2'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(e['five_alt2'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(b['five_alt2'], linewidth=linewidth, alpha=.3, color = color3)
        ax.plot((i['five_alt2']+ie['five_alt2']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((i['five_alt2']-ie['five_alt2']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((e['five_alt2']+ee['five_alt2']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((e['five_alt2']-ee['five_alt2']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.axvline(x=300,linestyle=':',alpha=0.5)
        ax.set_xticklabels(range(-300,51,50),rotation=90)
        
        ax.legend()
        plt.suptitle(title,y=1.03)
    plt.clf()
    plt.cla()
    plt.close()

def plot_a5ss(inclusion, exclusion, both, inclusion_err, exclusion_err, title, output_file, 
              color1=sns.color_palette("hls", 8)[0], color2=sns.color_palette("hls", 8)[5], color3='black'):
    
    i = {}
    e = {}
    b = {}
    ie = {}
    ee = {}
    i['three_alt1'] = np.array(inclusion['region1'][:350])
    i['three_alt2'] = np.array(inclusion['region1'][350:700])
    i['five_downstream'] = np.array(inclusion['region1'][700:])
    
    e['three_alt1'] = np.array(exclusion['region1'][:350])
    e['three_alt2'] = np.array(exclusion['region1'][350:700])
    e['five_downstream'] = np.array(exclusion['region1'][700:])
    
    b['three_alt1'] = np.array(both['region1'][:350])
    b['three_alt2'] = np.array(both['region1'][350:700])
    b['five_downstream'] = np.array(both['region1'][700:])
    
    ie['three_alt1'] = np.array(inclusion_err['region1'][:350])
    ie['three_alt2'] = np.array(inclusion_err['region1'][350:700])
    ie['five_downstream'] = np.array(inclusion_err['region1'][700:])
    
    ee['three_alt1'] = np.array(exclusion_err['region1'][:350])
    ee['three_alt2'] = np.array(exclusion_err['region1'][350:700])
    ee['five_downstream'] = np.array(exclusion_err['region1'][700:])
    """
    Special plot:
    plots a 4-region map that contains three separate plots for inclusion, 
    exclusion, and all spliced events. 
    
    Args:
        inclusion: {region1, region2, region3, region4}
        exclusion: {region1, region2, region3, region4}
        both: {region1, region2, region3, region4}
        
    """
    num_rows = 1
    num_cols = 3
    ax = plt.gca()
    min_height = min(min(i['three_alt1']),min(e['three_alt1']),min(b['three_alt1']),
                     min(i['three_alt2']),min(e['three_alt2']),min(b['three_alt2']),
                     min(i['five_downstream']),min(e['five_downstream']),min(b['five_downstream']))
    max_height = max(max(i['three_alt1']),max(e['three_alt1']),max(b['three_alt1']),
                     max(i['three_alt2']),max(e['three_alt2']),max(b['three_alt2']),
                     max(i['five_downstream']),max(e['five_downstream']),max(b['five_downstream']))
        
    with dataviz.Figure(output_file, figsize=(num_cols * 2.5,num_rows * 2.5)) as fig:
        
        sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})
        
        linewidth = 2
        errorbar_linewidth = 0.7
        
        ax = fig.add_subplot(1,4,1)
        ax.plot(i['three_alt1'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(e['three_alt1'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(b['three_alt1'], linewidth=linewidth, alpha=.3, color = color3)
        
        ax.plot((i['three_alt1']+ie['three_alt1']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((i['three_alt1']-ie['three_alt1']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((e['three_alt1']+ee['three_alt1']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((e['three_alt1']-ee['three_alt1']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_ylabel("Mean Read Density")
        ax.set_xticklabels(range(-50,351,50),rotation=90)
        ax.axvline(x=50,linestyle=':',alpha=0.5)
        
        sns.set_style({'ytick.major.size':0})
        
        ax = fig.add_subplot(1,4,2)
        ax.plot(i['three_alt2'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(e['three_alt2'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(b['three_alt2'], linewidth=linewidth, alpha=.3, color = color3)
        
        ax.plot((i['three_alt2']+ie['three_alt2']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((i['three_alt2']-ie['three_alt2']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((e['three_alt2']+ee['three_alt2']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((e['three_alt2']-ee['three_alt2']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.set_xticklabels(range(-50,351,50),rotation=90)
        ax.axvline(x=50,linestyle=':',alpha=0.5)
            
        ax = fig.add_subplot(1,4,3)
        ax.plot(i['five_downstream'], linewidth=linewidth, alpha=.8, color = color1)
        ax.plot(e['five_downstream'], linewidth=linewidth, alpha=.8, color = color2)
        ax.plot(b['five_downstream'], linewidth=linewidth, alpha=.3, color = color3)
        ax.plot((i['five_downstream']+ie['five_downstream']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((i['five_downstream']-ie['five_downstream']), linewidth=errorbar_linewidth, alpha=.5, color = color1, linestyle = ':')
        ax.plot((e['five_downstream']+ee['five_downstream']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        ax.plot((e['five_downstream']-ee['five_downstream']), linewidth=errorbar_linewidth, alpha=.5, color = color2, linestyle = ':')
        
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        ax.axvline(x=300,linestyle=':',alpha=0.5)
        ax.set_xticklabels(range(-300,51,50),rotation=90)
        
        ax.legend()
        plt.suptitle(title,y=1.03)
    plt.clf()
    plt.cla()
    plt.close()
    
def five_frame(region1, region2, region3, region4, region5,
               title, output_file, color='red'):
    num_rows = 1
    num_cols = 5
    
    ax = plt.gca()
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
    plt.clf()
    plt.cla()
    plt.close()
    