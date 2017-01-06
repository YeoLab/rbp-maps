'''
Created on Jun 20, 2016

@author: brianyee
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from gscripts.general import dataviz
import seaborn as sns

from matplotlib import rc
rc('text', usetex=False)
matplotlib.rcParams['svg.fonttype'] = 'none'
import numpy as np
import misc

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

def plot_err(ax, region, inclusion, exclusion, both, 
             inclusion_err, exclusion_err,
             inclusion_color, exclusion_color, both_color,
             inclusion_label = 'included in KD', 
             exclusion_label = 'excluded in KD', 
             both_label = 'all events', 
             linewidth=2, 
             errorbar_linewidth=0.7):
    """Helps plot a region for three events, usually inclusion, 
    exclusion, and both (background) densities, plus errorlines.
    
    Parameters
    ----------
    ax : matplotlib.axes._subplots.AxesSubplot
    region : string
    inclusion : dictionary
    exclusion : dictionary
    both : dictionary
    inclusion_err : dictionary
    exclusion_err : dictionary
    inclusion_color : string
    exclusion_color : string
    both_color : string
    inclusion_label : string
    exclusion_label : string
    both_label : string
    linewidth : int
    errorbar_linewidth : int
    """
    ax.plot(both[region], linewidth=linewidth, 
            alpha=.5, color = both_color, label = both_label)
    ax.plot(inclusion[region], linewidth=linewidth, 
            alpha=.8, color = inclusion_color, label = inclusion_label)
    ax.plot(exclusion[region], linewidth=linewidth, 
            alpha=.8, color = exclusion_color, label = exclusion_label)
    ax.plot((inclusion[region]+inclusion_err[region]), 
            linewidth=errorbar_linewidth, alpha=.5, color = inclusion_color, linestyle = ':')
    ax.plot((inclusion[region]-inclusion_err[region]), 
            linewidth=errorbar_linewidth, alpha=.5, color = inclusion_color, linestyle = ':')
    ax.plot((exclusion[region]+exclusion_err[region]), 
            linewidth=errorbar_linewidth, alpha=.5, color = exclusion_color, linestyle = ':')
    ax.plot((exclusion[region]-exclusion_err[region]), 
            linewidth=errorbar_linewidth, alpha=.5, color = exclusion_color, linestyle = ':')
    return min(min(inclusion[region]),min(exclusion[region])), max(max(inclusion[region]),max(exclusion[region]))
def single_frame(dic, title, output_file):
    """Plots a single frame given a set of points, 
    in this case the means of densities across a predescribed
    event.
    
    Parameters
    ----------
    dict (dictionary of lists) : list of mean density values to plot
    title (string) : plot title
    output_file (string) : output file including extension (.svg)
    color (string) : color of the means line (eg. 'red')

    """
    colors = sns.color_palette("hls", len(dic))
    ax = plt.gca()
    i = 0
    for key, value in dic.iteritems():
        ax.plot(value, color = colors[i], label = key)
        i = i + 1
    ax.legend()

    ax.set_ylabel('Read Density')
    ax.set_title(title,y=1.10)
    
    plt.savefig(output_file)
    plt.clf()
    plt.cla()
    plt.close()

def plot_bed(up, down, both, uperr, downerr, title, output_file, 
             cond1lab = 'cond1', cond2lab = 'cond2', bglab = 'background'):
    """Plots a single frame feature of events. 
    Just array-ifies the dictionary of lists and calls
    single_frame_with_up_down_events_error
    
    Parameters
    ----------
    up : dictionary of lists
    down : dictionary of lists
    both : dictionary of lists
    uperr : dictionary of lists
    downerr : dictionary of lists
    title : string
    output_file : string
    """
    single_frame_with_up_down_events_error(misc.toarray(up),
                                           misc.toarray(down),
                                           misc.toarray(both),
                                           misc.toarray(uperr),
                                           misc.toarray(downerr), 
                                           title, 
                                           output_file,
                                           uplab = cond1lab,
                                           downlab = cond2lab,
                                           bglab = bglab)
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
    
    min1, max1 = plot_err(ax, 'region1', up, down, both, uperr, downerr, upcolor, downcolor, bgcolor, uplab, downlab, bglab)
    
    ax.legend()

    ax.set_ylabel('Read Density')
    ax.set_title(title,y=1.10)
    # plt.xticks([0,300,599],['upstream (300bp)','feature','downstream (300bp)'])
    plt.ylim(min1+min1*0.1, max1+max1*0.1)
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
   
    with dataviz.Figure(output_file, figsize=(num_cols * 4,num_rows * 4)) as fig:
        
        sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})

        ax1 = fig.add_subplot(1,2,1)
        min1, max1 = plot_err(ax1, 'region1', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
        ax1.set_ylabel("Normalized Signal")
        ax1.set_xticklabels(range(-50,351,50),rotation=90)
        ax1.axvline(x=50,linestyle=':',alpha=0.5)
        sns.despine(ax=ax1)
        sns.set_style({'ytick.major.size':0})
        
        ax2 = fig.add_subplot(1,2,2, sharey=ax1)
        min2, max2 = plot_err(ax2, 'region2', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
        sns.despine(ax=ax2, left=True)
        ax2.set_yticklabels([])
        ax2.set_xticklabels(range(-300,51,50),rotation=90)
        ax2.axvline(x=300,linestyle=':',alpha=0.5)

        ax2.legend()
        
    
        print(max1, max2)
        print(min1, min2)
        mx = max(max1,max2)
        mi = min(min1,min2)
        # print(mi,mx)
        plt.ylim(mi-mi*0.1,mx+mx*0.1)
        increment = (mx)/6
        try:
            ax1.set_yticks(range(mi,mx,increment))
            ax1.set_yticklabels(range(mi,mx,increment))
        except Exception as e:
            title = title + "\ny-lim: {} - {}".format(mi,mx)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.suptitle(title,y=1.10)
        plt.tight_layout()
    plt.clf()
    plt.cla()
    plt.close()
def three_frame(region1, region2, region3, 
               title, output_file, color='red'):
    num_rows = 1
    num_cols = 4
    
    with dataviz.Figure(output_file, figsize=(num_cols * 4,num_rows * 4)) as fig:
            
        min_height = min(min(region1),min(region2),min(region3))
        max_height = max(max(region1),max(region2),max(region3))
        
        linewidth = 2.5
        ax = fig.add_subplot(1,4,1)
        ax.plot(region1, linewidth=linewidth, alpha=.7, color = color)
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_ylabel("Normalized Signal")
            
        ax = fig.add_subplot(1,4,2)
        ax.plot(region2, linewidth=linewidth, alpha=.7, color = color)
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
            
        ax = fig.add_subplot(1,4,3)
        ax.plot(region3, linewidth=linewidth, alpha=.7, color = color)
            
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])

    plt.clf()
    plt.cla()
    plt.close()
    
def four_frame(region1, region2, region3, region4, 
               title, output_file, color='red'):
    num_rows = 1
    num_cols = 4
    
    with dataviz.Figure(output_file, figsize=(num_cols * 4,num_rows * 4)) as fig:
            
        min_height = min(min(region1),min(region2),min(region3),min(region4))
        max_height = max(max(region1),max(region2),max(region3),max(region4))
            
        linewidth = 2.5
        ax = fig.add_subplot(1,4,1)
        ax.plot(region1, linewidth=linewidth, alpha=.7, color = color)
        sns.despine(ax=ax)
        ax.set_ylim(min_height, max_height)
        ax.set_ylabel("Normalized Signal")
            
        ax = fig.add_subplot(1,4,2)
        ax.plot(region2, linewidth=linewidth, alpha=.7, color = color)
            
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
            
        ax = fig.add_subplot(1,4,3)
        ax.plot(region3, linewidth=linewidth, alpha=.7, color = color)
            
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
            
        ax = fig.add_subplot(1,4,4)
        ax.plot(region4, linewidth=linewidth, alpha=.7, color = color)
            
        sns.despine(ax=ax, left=True)
        ax.set_ylim(min_height, max_height)
        ax.set_yticklabels([])
        plt.suptitle(title,y=1.10)
        plt.tight_layout()
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
    # ax = plt.gca()
    # with dataviz.Figure(output_file, figsize=(num_cols * 4,num_rows * 4)) as fig:
    fig = plt.figure(figsize=(2.5*num_cols, 2.5*num_rows)) 
    
    sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})
      
    # sns.set(font_scale=1)
    ax1 = fig.add_subplot(1,4,1)
    min1, max1 = plot_err(ax1, 'region1', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
        
    ax1.set_ylabel("Normalized signal")
    ax1.set_xticklabels(labels=range(-50,351,50),rotation=90)
    ax1.axvline(x=50,linestyle=':',alpha=0.5)
    sns.despine(ax=ax1)
    sns.set_style({'ytick.major.size':0})

    # sns.set(font_scale=0)

    ax2 = fig.add_subplot(1,4,2,sharey=ax1)
    min2, max2 = plot_err(ax2, 'region2', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
    ax2.set_yticklabels([])
    ax2.set_xticklabels(range(-300,51,50),rotation=90, visible=True)
    # ax2.set_xticklabels(range(-50,351,50),visible=True,color="black",rotation=90)
    ax2.axvline(x=300,linestyle=':',alpha=0.5)
    sns.despine(ax=ax2, left=True)
    ax3 = fig.add_subplot(1,4,3,sharey=ax1)
    min3, max3 = plot_err(ax3, 'region3', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
    ax3.set_yticklabels([])
    ax3.set_xticklabels(labels=range(-50,351,50),rotation=90)
    ax3.axvline(x=50,linestyle=':',alpha=0.5)
    sns.despine(ax=ax3, left=True)
        
    ax4 = fig.add_subplot(1,4,4,sharey=ax1)
    min4, max4 = plot_err(ax4, 'region4', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
    ax4.axvline(x=300,linestyle=':',alpha=0.5)
    ax4.set_xticklabels(labels=range(-300,51,50),rotation=90)
    
    sns.despine(ax=ax4, left=True)

    ax4.legend()
    
    """
    Set min/max for y-axis, hide y-axis labels for all but the first
    """
    mx = max(max1,max2,max3,max4)
    mi = min(min1,min2,min3,min4)
    plt.ylim(mi-mi*0.1,mx+mx*0.1)
    increment = (mx)/6
    try:
        ax1.set_yticks(range(mi,mx,increment))
        ax1.set_yticklabels(range(mi,mx,increment))
    except Exception as e:
        title = title + "\ny-lim: {} - {}".format(mi,mx)
    plt.setp(ax2.get_yticklabels(), visible=False)
    plt.setp(ax3.get_yticklabels(), visible=False)
    plt.setp(ax4.get_yticklabels(), visible=False)
    
    plt.suptitle(title,y=1.10)
    plt.tight_layout()
    plt.savefig(output_file)
    plt.clf()
    plt.cla()
    plt.close()

def six_frame_with_inclusion_exclusion_events_with_error(inclusion, exclusion, both,
                                                         inclusion_err, exclusion_err,
                                                         title, output_file, 
                                                         color1=sns.color_palette("hls", 8)[0], 
                                                         color2=sns.color_palette("hls", 8)[5],
                                                         color3='black'):
    
    """
    Special plot:
    plots a 6-region map that contains three separate plots for inclusion, 
    exclusion, and all spliced events. 
    
    Args:
        inclusion: {region1, region2, region3, region4}
        exclusion: {region1, region2, region3, region4}
        both: {region1, region2, region3, region4}
        
    """
    num_rows = 1
    num_cols = 6
    
    with dataviz.Figure(output_file, figsize=(num_cols * 4,num_rows * 4)) as fig:
        
        sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})
        
        ax1 = fig.add_subplot(1,6,1)
        min1, max1 = plot_err(ax1, 'region1', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
        ax1.set_ylabel("Normalized Signal")
        ax1.set_xticklabels(range(-50,351,50),rotation=90)
        ax1.axvline(x=50,linestyle=':',alpha=0.5)
        sns.despine(ax=ax1)
        sns.set_style({'ytick.major.size':0})
        
        ax2 = fig.add_subplot(1,6,2,sharey=ax1)
        min2, max2 = plot_err(ax2, 'region2', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
        ax2.set_yticklabels([])
        ax2.set_xticklabels(range(-300,51,50),rotation=90)
        ax2.axvline(x=300,linestyle=':',alpha=0.5)
        sns.despine(ax=ax2, left=True)
        
        ax3 = fig.add_subplot(1,6,3,sharey=ax1)
        min3, max3 = plot_err(ax3, 'region3', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
        ax3.set_yticklabels([])
        ax3.set_xticklabels(range(-50,351,50),rotation=90)
        ax3.axvline(x=50,linestyle=':',alpha=0.5)
        sns.despine(ax=ax3, left=True)
        
        ax4 = fig.add_subplot(1,6,4,sharey=ax1)
        min4, max4 = plot_err(ax4, 'region4', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
        ax4.set_yticklabels([])
        ax4.axvline(x=300,linestyle=':',alpha=0.5)
        ax4.set_xticklabels(range(-300,51,50),rotation=90)
        sns.despine(ax=ax4, left=True)
        
        ax5 = fig.add_subplot(1,6,5,sharey=ax1)
        min5, max5 = plot_err(ax5, 'region5', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
        ax5.set_yticklabels([])
        ax5.set_xticklabels(range(-50,351,50),rotation=90)
        ax5.axvline(x=50,linestyle=':',alpha=0.5)
        sns.despine(ax=ax5, left=True)
        
        ax6 = fig.add_subplot(1,6,6,sharey=ax1)
        min6, max6 = plot_err(ax6, 'region6', inclusion, exclusion, both, inclusion_err, exclusion_err, color1, color2, color3)
        ax6.set_yticklabels([])
        ax6.axvline(x=300,linestyle=':',alpha=0.5)
        ax6.set_xticklabels(range(-300,51,50),rotation=90)
        sns.despine(ax=ax6, left=True)
        ax6.legend()
        
        
        mx = max(max1,max2,max3,max4,max5,max6)
        mi = min(min1,min2,min3,min4,min5,min6)
        plt.ylim(mi-mi*0.1,mx+mx*0.1)
        increment = (mx)/6
        try:
            ax1.set_yticks(range(mi,mx,increment))
            ax1.set_yticklabels(range(mi,mx,increment))
        except Exception as e:
            title = title + "\ny-lim: {} - {}".format(mi,mx)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_yticklabels(), visible=False)
        plt.setp(ax4.get_yticklabels(), visible=False)
        plt.setp(ax5.get_yticklabels(), visible=False)
        plt.setp(ax6.get_yticklabels(), visible=False)
        plt.tight_layout()
        plt.suptitle(title,y=1.10)
    plt.clf()
    plt.cla()
    plt.close()
    
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
        inclusion (dictionary) : dictionary of means
        exclusion (dictionary) : dictionary of means
        both: dictionary of means
        
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

def plot_mxe(inclusion, exclusion, both, inclusion_err, exclusion_err, title, output_file,
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
    i['region1'], i['region2'], i['region3'], i['region4'], i['region5'], i['region6'] = np.array_split(inclusion['region1'],6)
    e['region1'], e['region2'], e['region3'], e['region4'], e['region5'], e['region6'] = np.array_split(exclusion['region1'],6)
    b['region1'], b['region2'], b['region3'], b['region4'], b['region5'], b['region6'] = np.array_split(both['region1'],6)
    ie['region1'], ie['region2'], ie['region3'], ie['region4'], ie['region5'], ie['region6'] = np.array_split(inclusion_err['region1'],6)
    ee['region1'], ee['region2'], ee['region3'], ee['region4'], ee['region5'], ee['region6'] = np.array_split(exclusion_err['region1'],6)
    
    six_frame_with_inclusion_exclusion_events_with_error(i,e,b,ie,ee,title,output_file,color1,color2,color3)
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

    with dataviz.Figure(output_file, figsize=(num_cols * 4,num_rows * 4)) as fig:
        
        sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})
        
        ax1 = fig.add_subplot(1,4,1)
        min1, max1 = plot_err(ax1, 'three_upstream', i, e, b, ie, ee, color1, color2, color3)
        sns.despine(ax=ax1)
        ax1.set_ylabel("Normalized Signal")
        ax1.set_xticklabels(range(-50,351,50),rotation=90)
        ax1.axvline(x=50,linestyle=':',alpha=0.5)
        
        sns.set_style({'ytick.major.size':0})
        
        ax2 = fig.add_subplot(1,4,2, sharey=ax1)
        min2, max2 = plot_err(ax2, 'five_alt1', i, e, b, ie, ee, color1, color2, color3)
        
        sns.despine(ax=ax2, left=True)
        ax2.set_yticklabels([])
        ax2.set_xticklabels(range(-300,51,50),rotation=90)
        ax2.axvline(x=300,linestyle=':',alpha=0.5)
            
        ax3 = fig.add_subplot(1,4,3, sharey=ax1)
        min3, max3 = plot_err(ax3, 'five_alt2', i, e, b, ie, ee, color1, color2, color3)
        
        sns.despine(ax=ax3, left=True)
        ax3.set_yticklabels([])
        ax3.axvline(x=300,linestyle=':',alpha=0.5)
        ax3.set_xticklabels(range(-300,51,50),rotation=90)
        
        ax3.legend()
        
    
        mx = max(max1,max2,max3)
        mi = min(min1,min2,min3)
        plt.ylim(mi-mi*0.1,mx+mx*0.1)
        increment = (mx)/6
        try:
            ax1.set_yticks(range(mi,mx,increment))
            ax1.set_yticklabels(range(mi,mx,increment))
        except Exception as e:
            title = title + "\ny-lim: {} - {}".format(mi,mx)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_yticklabels(), visible=False)
        plt.suptitle(title,y=1.10)
        plt.tight_layout()
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

    with dataviz.Figure(output_file, figsize=(num_cols * 4,num_rows * 4)) as fig:
        
        sns.set_style({'xtick.major.size':5,
                   'ytick.major.size':5,
                   'xtick.color':'.15'})
        
        ax1 = fig.add_subplot(1,4,1)
        min1, max1 = plot_err(ax1, 'three_alt1', i, e, b, ie, ee, color1, color2, color3)
        sns.despine(ax=ax1)
        ax1.set_ylabel("Normalized Signal")
        ax1.set_xticklabels(range(-50,351,50),rotation=90)
        ax1.axvline(x=50,linestyle=':',alpha=0.5)
        sns.set_style({'ytick.major.size':0})
        
        ax2 = fig.add_subplot(1,4,2)
        min2, max2 = plot_err(ax2, 'three_alt2', i, e, b, ie, ee, color1, color2, color3)
        sns.despine(ax=ax2, left=True)
        ax2.set_yticklabels([])
        ax2.set_xticklabels(range(-50,351,50),rotation=90)
        ax2.axvline(x=50,linestyle=':',alpha=0.5)
            
        ax3 = fig.add_subplot(1,4,3)
        min3, max3 = plot_err(ax3, 'five_downstream', i, e, b, ie, ee, color1, color2, color3)
        sns.despine(ax=ax3, left=True)
        ax3.set_yticklabels([])
        ax3.axvline(x=300,linestyle=':',alpha=0.5)
        ax3.set_xticklabels(range(-300,51,50),rotation=90)
        
        ax3.legend()
    
        mx = max(max1,max2,max3)
        mi = min(min1,min2,min3)
        plt.ylim(mi-mi*0.1,mx+mx*0.1)
        increment = (mx)/6
        try: # for very small ranges, increments are too tiny?
            ax1.set_yticks(range(mi,mx,increment))
            ax1.set_yticklabels(range(mi,mx,increment))
        except Exception as e:
            title = title + "\ny-lim: {} - {}".format(mi,mx)
        plt.setp(ax2.get_yticklabels(), visible=False)
        plt.setp(ax3.get_yticklabels(), visible=False)
        plt.suptitle(title,y=1.10)
        plt.tight_layout()
    plt.clf()
    plt.cla()
    plt.close()