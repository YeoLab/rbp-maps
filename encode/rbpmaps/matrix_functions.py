'''
Created on Jun 18, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import intervals
import misc
import Feature
import logging
logger = logging.getLogger('plot_features')

def create_matrix(annotation, density, 
                  upstream_offset, downstream_offset, 
                  is_scaled = True, annotation_type = 'bed'):
    """
    Creates an r x c pandas dataframe of r events for a feature of length c.
    
    Args:
        annotation (string) : path of file containing the annotation
        density (ReadDensity) : object containing positive and negative BigWig files
        upstream_offset (integer) : number of bases upstream of feature to plot
        downstream_offset (integer) : number of bases downstream of feature to plot
        is_scaled (boolean) : if all features are of different length, this must be true
            to resize all features to fit on a 0-100% scale.
        annotation_type (string) : may be bed format or any additional defined format in Feature
    
    Returns:
        pandas.DataFrame : a dataframe of r events for a feature of length c.
    """
    count = 0
    densities = {}
    logger.info("Start matrix creation [ANNOTATION:{},DENSITY:{},UP:{},DOWN:{},SCALED:{},TYPE:{}".format(
                                                                                                            annotation,
                                                                                                            density.name,
                                                                                                            upstream_offset,
                                                                                                            downstream_offset,
                                                                                                            is_scaled,
                                                                                                            annotation_type))
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID'):
                count = count + 1
                if count % 50000 == 0:
                    logger.info('Processed {} features'.format(count))
                event = line.rstrip() # .split('\t')[0]
                interval = Feature.Feature(event, annotation_type).get_bedtool()
                wiggle = pd.Series(intervals.some_range(density, interval, 0, 0))
                wiggle = wiggle.fillna(0) # convert all nans to 0
                wiggle = abs(wiggle) # convert all values to positive
                if(is_scaled == True):
                    wiggle = intervals.get_scale(wiggle)
                densities[intervals.rename_index(interval)] = wiggle
                """ Leave out for later discussion, but even if there is nothing there, 
                    we want to keep the region. When we add pseudocounts, it will
                    contribute to the map as a 'flat' region.
                if not all(np.isnan(wiggle)):
                    wiggle = wiggle.fillna(0) # convert all nans to 0
                    wiggle = abs(wiggle) # convert all values to positive
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    densities[intervals.rename_index(interval)] = wiggle
                """
    logger.info("Finished matrix creation.")
    return pd.DataFrame(densities).T

def create_mxe_matrix(annotation, density, 
                      exon_offset, intron_offset, 
                      is_scaled = False, combine_regions = True, 
                      annotation_type="rmats"):
    logger.info("Starting MXE matrix creation [ANNOTATION:{},DENSITY:{},UP:{},DOWN:{},SCALED:{},TYPE:{}".format(
                                                                                                            annotation,
                                                                                                            density.name,
                                                                                                            exon_offset,
                                                                                                            intron_offset,
                                                                                                            is_scaled,
                                                                                                            annotation_type))
    """
    Creates an r x c pandas dataframe of r events for a mutually exclusive
    exon feature. An MXE matrix will contain six distinct regions: 
    
    |_]----||----[__||__]----||----[__||__]----||----[_|
    
    - the [..exon_offset]--intron_offset--... 3' site of an upstream exon
    - the ...--intron_offset--[exon_offset..] 5' site of the upstream skipped exon
    - the [..exon_offset]--intron_offset--... 3' site of the upstream skipped exon
    - the ...--intron_offset--[exon_offset..] 5' site of the downstream skipped exon
    - the [..exon_offset]--intron_offset--... 3' site of the downstream skipped exon
    - the ..--intron_offset--[exon_offset..] 5' site of the downstream exon
    Args:
        annotation (string) : path of file containing the annotation
        density (ReadDensity) : object containing positive and negative BigWig files
        exon_offset (integer) : how far into the exon boundary to plot
        intron_offset (integer) : how far after the exon boundary to plot
        is_scaled (boolean) : if all features are of different length, this must be true
            to resize all features to fit on a 0-100% scale.
        combine_regions (boolean) : if False, return six DataFrames instead of one.
        annotation_type (string) : may be rmats format or any additional defined format in Feature
    
    Returns:
        pandas.DataFrame : a dataframe of r events for an MXE feature.
    """
    three_upstream = {}
    three_up_mxe = {}
    five_up_mxe = {}
    three_down_mxe = {}
    five_down_mxe = {}
    five_downstream = {}
    
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID'):
                event = line.rstrip() # .split('\t')[0]
                upstream_interval, upstream_mxe_interval, \
                downstream_mxe_interval, downstream_interval = Feature.MXEFeature(event, annotation_type).get_bedtools() 
                
                """three prime upstream region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        upstream_mxe_interval,
                                                                        upstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                    
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
        
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) 
                
                three_upstream[event] = wiggle
                
                """five prime site of mxe1 (upstream mxe) region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream_interval,
                                                                        upstream_mxe_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle)
                five_up_mxe[event] = wiggle
                
                """three prime site of mxe1 (upstream mxe) region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        downstream_mxe_interval,
                                                                        upstream_mxe_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                    
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
        
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) 
                
                three_up_mxe[event] = wiggle
                
                """five prime site of mxe2 (downstream mxe) region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream_mxe_interval,
                                                                        downstream_mxe_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle)
                five_down_mxe[event] = wiggle
                
                """three prime site of mxe2 (downstream mxe) region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        downstream_interval,
                                                                        downstream_mxe_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                    
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
        
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) 
                
                three_down_mxe[event] = wiggle
                
                
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        downstream_mxe_interval,
                                                                        downstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) # convert all nans to 0

                five_downstream[event] = wiggle
    
        three_upstream = pd.DataFrame(three_upstream).T
        five_up_mxe = pd.DataFrame(five_up_mxe).T
        three_up_mxe = pd.DataFrame(three_up_mxe).T
        five_down_mxe = pd.DataFrame(five_down_mxe).T
        three_down_mxe = pd.DataFrame(three_down_mxe).T
        five_downstream = pd.DataFrame(five_downstream).T
    logger.info("Finished matrix creation.")
    if combine_regions == False:
        return three_upstream, five_up_mxe, three_up_mxe, five_down_mxe, three_down_mxe, five_downstream
    else:
        ra = pd.concat([three_upstream, five_up_mxe, three_up_mxe, five_down_mxe, three_down_mxe, five_downstream],axis=1)
        ra.columns = range(0,ra.shape[1])
        return ra
def create_ri_matrix(annotation, density, 
                     exon_offset, intron_offset, 
                     is_scaled = False, combine_regions = True, 
                     annotation_type="rmats"):
    """
    Creates an r x c pandas dataframe of r events for a Retained Intron feature.
    A RI matrix will contain two distinct regions: 
    
    |_]----||----[_|
    
    - the [..exon_offset]--intron_offset--... 3' site of an upstream exon
    - the ..--intron_offset--[exon_offset..] 5' site of the downstream exon
    Args:
        annotation (string) : path of file containing the annotation
        density (ReadDensity) : object containing positive and negative BigWig files
        exon_offset (integer) : how far into the exon boundary to plot
        intron_offset (integer) : how far after the exon boundary to plot
        is_scaled (boolean) : if all features are of different length, this must be true
            to resize all features to fit on a 0-100% scale.
        combine_regions (boolean) : if False, return two DataFrames instead of one
        annotation_type (string) : may be rmats format or any additional defined format in Feature
    
    Returns:
        pandas.DataFrame : a dataframe of r events for an MXE feature.
    """
    logger.info("Starting RI matrix creation [ANNOTATION:{},DENSITY:{},UP:{},DOWN:{},SCALED:{},TYPE:{}".format(
                                                                                                            annotation,
                                                                                                            density.name,
                                                                                                            exon_offset,
                                                                                                            intron_offset,
                                                                                                            is_scaled,
                                                                                                            annotation_type))
    
    three_upstream = {}
    five_downstream = {}
            
    with open(annotation) as f:
        # f.next() # for title
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID'):
                event = line.rstrip() # .split('\t')[0]
                upstream_interval, downstream_interval = Feature.RIFeature(event,annotation_type).get_bedtools()
                
                """three prime upstream region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        downstream_interval,
                                                                        upstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                    
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) 
                three_upstream[event] = wiggle
                
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream_interval,
                                                                        downstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                five_downstream[event] = wiggle
        
        three_upstream = pd.DataFrame(three_upstream).T
        five_downstream = pd.DataFrame(five_downstream).T
    logger.info("Finished matrix creation.")
    if combine_regions == False:
        return three_upstream, five_downstream
    else:
        ra = pd.concat([three_upstream,five_downstream],axis=1)
        ra.columns = range(0,ra.shape[1])
        # print("TYPE OF MATRIX: {}".format(type(ra)))
        return ra
    
def create_a5ss_matrix(annotation, density, 
                       exon_offset, intron_offset, 
                       is_scaled = True, combine_regions = True, 
                       annotation_type = "rmats"):
    """
    Creates an r x c pandas dataframe of r events for an 
    alternative 5' splice site feature. An A5ss matrix will 
    contain three distinct regions: 
    
    ______|__]------||------[__|
    __|__]------|    |------[__|
    
    - the [..exon_offset]--intron_offset--... 3' site of an alt1 spliced exon
    - the [..exon_offset]--intron_offset--... 3' site of an alt2 spliced exon
    - the ..--intron_offset--[exon_offset..] 5' site of the downstream exon
    Args:
        annotation (string) : path of file containing the annotation
        density (ReadDensity) : object containing positive and negative BigWig files
        exon_offset (integer) : how far into the exon boundary to plot
        intron_offset (integer) : how far after the exon boundary to plot
        is_scaled (boolean) : if all features are of different length, this must be true
            to resize all features to fit on a 0-100% scale.
        combine_regions (boolean) : if False, return three DataFrames instead of one.
        annotation_type (string) : may be rmats format or any additional defined format in Feature
    
    Returns:
        pandas.DataFrame : a dataframe of r events for an A5SS feature.
    """
    logger.info("Starting A5SS matrix creation [ANNOTATION:{},DENSITY:{},UP:{},DOWN:{},SCALED:{},TYPE:{}".format(
                                                                                                            annotation,
                                                                                                            density.name,
                                                                                                            exon_offset,
                                                                                                            intron_offset,
                                                                                                            is_scaled,
                                                                                                            annotation_type))
    three_alt1 = {}
    three_alt2 = {}
    five_downstream = {}

    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID'):
                event = line.rstrip() # .split('\t')[0]
                alt1, alt2, downstream = Feature.A5ssFeature(event,annotation_type).get_bedtools()
                """three prime alt1  (shorter) region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        downstream,
                                                                        alt1,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
        
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) 
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_alt1[event] = wiggle
                    
                """three prime site of alt2 (longer) region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                         downstream,
                                                                         alt2,
                                                                         exon_offset,
                                                                         intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) #
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_alt2[event] = wiggle
                    # print("length of 3p skipped: {}".format(len(wiggle)))
        
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        alt2,
                                                                        downstream,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_downstream[event] = wiggle
                    # print("length of downstream: {}".format(len(wiggle)))
        three_alt1 = pd.DataFrame(three_alt1).T
        three_alt2 = pd.DataFrame(three_alt2).T
        five_downstream = pd.DataFrame(five_downstream).T
    logger.info("Finished matrix creation.")
    if combine_regions == False:
        return three_alt1, three_alt2, five_downstream
    else:
        ra = pd.concat([three_alt1,three_alt2,five_downstream],axis=1)
        ra.columns = range(0,ra.shape[1])
        return ra      
def create_a3ss_matrix(annotation, density, exon_offset, intron_offset, is_scaled, combine_regions=True, annotation_type="rmats"):
    """
    Creates an r x c pandas dataframe of r events for an 
    alternative 3' splice site feature. An A3SS matrix will 
    contain three distinct regions: 

    __|__]------||-----[__|____
    __|__]------|    |------[__|
    
    - the [..exon_offset]--intron_offset--... 3' site of an upstream exon
    - the ..--intron_offset--[exon_offset..] 5' site of the alt1 spliced exon
    - the ..--intron_offset--[exon_offset..] 5' site of the alt2 spliced exon
    
    Args:
        annotation (string) : path of file containing the annotation
        density (ReadDensity) : object containing positive and negative BigWig files
        exon_offset (integer) : how far into the exon boundary to plot
        intron_offset (integer) : how far after the exon boundary to plot
        is_scaled (boolean) : if all features are of different length, this must be true
            to resize all features to fit on a 0-100% scale.
        combine_regions (boolean) : if False, return three DataFrames instead of one.
        annotation_type (string) : may be rmats format or any additional defined format in Feature
    
    Returns:
        pandas.DataFrame : a dataframe of r events for an A3SS feature.
    """
    logger.info("Starting a3ss matrix creation [ANNOTATION:{},DENSITY:{},UP:{},DOWN:{},SCALED:{},TYPE:{}".format(
                                                                                                            annotation,
                                                                                                            density.name,
                                                                                                            exon_offset,
                                                                                                            intron_offset,
                                                                                                            is_scaled,
                                                                                                            annotation_type))
    three_upstream = {}
    five_alt1 = {}
    five_alt2 = {}
    
    
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID'):
                event = line.rstrip() # .split('\t')[0]
                
                upstream, alt1, alt2 = Feature.A3ssFeature(event,annotation_type).get_bedtools()
                # print('three prime site upstream')
                
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        alt1,
                                                                        upstream,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
        
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) 
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    three_upstream[event] = wiggle
        
                """five prime site of alt1 (longer exon)"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream,
                                                                        alt1,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle)
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_alt1[event] = wiggle
                
        
                """five prime site of alt2 (shorter exon) """
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream,
                                                                        alt2,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                if not all(np.isnan(wiggle)):
                    wiggle = abs(wiggle) # convert all values to positive
                    wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                    wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                    if(is_scaled == True):
                        wiggle = intervals.get_scale(wiggle)
                    five_alt2[event] = wiggle
                  
        three_upstream = pd.DataFrame(three_upstream).T
        five_alt1 = pd.DataFrame(five_alt1).T
        five_alt2 = pd.DataFrame(five_alt2).T
    logger.info("Finished matrix creation.")
    if combine_regions == False:
        return three_upstream, five_alt1, five_alt2
    else:
        ra = pd.concat([three_upstream,five_alt1,five_alt2],axis=1)
        ra.columns = range(0,ra.shape[1])
        return ra                           
                        
def create_se_matrix(annotation, density, exon_offset, intron_offset, is_scaled, combine_regions=True, annotation_type="rmats"):
    """
    Creates an r x c pandas dataframe of r events for a skipped
    exon feature. An SE matrix will contain four distinct regions: 
    
    |_]----||----[__||__]----||----[_|
    
    - the [..exon_offset]--intron_offset--... 3' site of an upstream exon
    - the ...--intron_offset--[exon_offset..] 5' site of the upstream skipped exon
    - the [..exon_offset]--intron_offset--... 3' site of the downstream skipped exon
    - the ..--intron_offset--[exon_offset..] 5' site of the downstream exon
    Args:
        annotation (string) : path of file containing the annotation
        density (ReadDensity) : object containing positive and negative BigWig files
        exon_offset (integer) : how far into the exon boundary to plot
        intron_offset (integer) : how far after the exon boundary to plot
        is_scaled (boolean) : if all features are of different length, this must be true
            to resize all features to fit on a 0-100% scale.
        combine_regions (boolean) : if False, return four DataFrames instead of one.
        annotation_type (string) : may be rmats format or any additional defined format in Feature
    
    Returns:
        pandas.DataFrame : a dataframe of r events for an SE feature.
    """
    logger.info("Starting SE matrix creation [ANNOTATION:{},DENSITY:{},UP:{},DOWN:{},SCALED:{},TYPE:{}".format(
                                                                                                            annotation,
                                                                                                            density.name,
                                                                                                            exon_offset,
                                                                                                            intron_offset,
                                                                                                            is_scaled,
                                                                                                            annotation_type))
    three_upstream = {}
    five_skipped = {}
    three_skipped = {}
    five_downstream = {}
    
    with open(annotation) as f:
        for line in f:
            if not line.startswith('event_name') and not line.startswith('ID'):
                event = line.rstrip()
                upstream_interval, interval, downstream_interval = Feature.SkippedExonFeature(event,annotation_type).get_bedtools()
                
                """three prime upstream region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                        interval,
                                                                        upstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
        
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) 

                three_upstream[event] = wiggle
                """five prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        upstream_interval,
                                                                        interval,
                                                                        exon_offset,
                                                                        intron_offset)
                
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle)
                five_skipped[event] = wiggle
                """three prime site of skipped region"""
                left_pad, wiggle, right_pad = intervals.three_prime_site(density, 
                                                                         downstream_interval,
                                                                         interval,
                                                                         exon_offset,
                                                                         intron_offset)
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) #
                three_skipped[event] = wiggle
                """five prime site of downstream region"""
                left_pad, wiggle, right_pad = intervals.five_prime_site(density, 
                                                                        interval,
                                                                        downstream_interval,
                                                                        exon_offset,
                                                                        intron_offset)
                wiggle = pd.Series(wiggle)
                wiggle = abs(wiggle) # convert all values to positive
                wiggle = np.pad(wiggle,(left_pad,right_pad),'constant',constant_values=(-1))
                wiggle = np.nan_to_num(wiggle) # convert all nans to 0
                five_downstream[event] = wiggle

        three_upstream = pd.DataFrame(three_upstream).T
        five_skipped = pd.DataFrame(five_skipped).T
        three_skipped = pd.DataFrame(three_skipped).T
        five_downstream = pd.DataFrame(five_downstream).T
    logger.info("Finished matrix creation.")
    if combine_regions == False:
        return three_upstream, five_skipped, three_skipped, five_downstream
    else:
        ra = pd.concat([three_upstream,five_skipped,three_skipped,five_downstream],axis=1)
        ra.columns = range(0,ra.shape[1])
        return ra

def make_hist_se(infile, outfile, 
                 l10p_cutoff, l2fc_cutoff, 
                 all_exons, exon_offset, 
                 intron_offset):
    """
    Creates an r x c pandas dataframe of r events for a skipped
    exon feature. An SE matrix will contain four distinct regions: 
    
    |_]----||----[__||__]----||----[_|
    
    - the [..exon_offset]--intron_offset--... 3' site of an upstream exon
    - the ...--intron_offset--[exon_offset..] 5' site of the upstream skipped exon
    - the [..exon_offset]--intron_offset--... 3' site of the downstream skipped exon
    - the ..--intron_offset--[exon_offset..]  5' site of the downstream exon
    Args:
        infile (string) : input file (input normed bedfile)
        outfile (string) : output file containing peaks per position
        l10p_cutoff (float) : l10 pvalue cutoff
        l2fc_cutoff (float) : l2 fold change cutoff
        all_exons (string) : MISO-style annotation
        exon_offset (integer) : how far into the exon boundary to plot
        intron_offset (integer) : how far after the exon boundary to plot
    Writes:
        file containing the number of peaks at each position. Line-delimited.
    """
    
    try:
        region_types = ["upstream_region_skipped_exon",
                        "upstream_region_downstream_exon",
                        "downstream_region_skipped_exon",
                        "downstream_region_upstream_exon"]
        position_sum = {}
        count = 0
        hashing_val = 10000
        with open(infile,'r') as f:
            for line in f:
                line = line.split('\t')
                chrom = line[0]
                pstart = int(line[1])
                pstop = int(line[2])
                l10p = float(line[3])
                l2fc = float(line[4])
                stra = line[5].strip()
                
                # correct bed files being 0-based, open ended
                pstart = pstart + 1
                
                if l10p < l10p_cutoff:
                    continue
                if l2fc < l2fc_cutoff:
                    continue
                
                x = int(pstart / hashing_val)
                y = int(pstop / hashing_val)

                # for each peak, find ALL regions that intersect it
                for region_type in region_types: # within a region
                    tmphash = {}
                    for i in range(x,y+1): # within a bin
                        for event in all_exons[chrom,stra,i,region_type]:
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            if pstop < int(exregstart): # pass if peak stop occurs before exon region start
                                continue
                            if pstart > int(exregstop): # pass if peak start occurs after exon region end
                                continue
                            tmphash[event] = 1 # otherwise peak falls within event region
                    for event in tmphash:
                        if stra == "+":
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart)) # peak start OR region start
                            end_val = min(int(pstop), int(exregstop)) # peak stop OR region stop
                            for j in range(start_val, end_val+1): # count intersecting positions between peak and region
                                relative_pos = j - int(exstart) # calculate relative position
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos) # count + 1 for the region
                        elif stra == '-':
                            exchr, exregstart, exstart, exregstop, exstr = event.split(':')
                            start_val = max(int(pstart), int(exregstart))
                            end_val = min(int(pstop), int(exregstop))
                            for j in range(start_val, end_val+1):
                                relative_pos = -1 * (j - int(exstart))
                                position_sum[region_type, relative_pos] = misc.ini(position_sum,
                                                                              region_type, 
                                                                              relative_pos)
                        else:
                            print("strand error\n")
                    # we have a peak that maps to every region
                    
        # count from 0 to max
        current_pos = 0
        o = open(outfile,'w')
        for j in range(-exon_offset, intron_offset+1):
            if misc.exists(position_sum,"downstream_region_upstream_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_upstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_offset, exon_offset+1):
            if misc.exists(position_sum,"upstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-exon_offset, intron_offset+1):
            if misc.exists(position_sum,"downstream_region_skipped_exon",j):
                o.write("{}\n".format(position_sum["downstream_region_skipped_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        for j in range(-intron_offset, exon_offset+1):
            if misc.exists(position_sum,"upstream_region_downstream_exon",j):
                o.write("{}\n".format(position_sum["upstream_region_downstream_exon",j]))
            else:
                o.write("{}\n".format(0))
            current_pos = current_pos + 1
        o.close()
    except Exception as e:
        print(e)