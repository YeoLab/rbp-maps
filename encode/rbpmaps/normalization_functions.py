'''
Created on Jun 27, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import logging
logger = logging.getLogger('plot_features')


def clean(density):
    """
    These functions expect a dataframe with density values (columns)
    across a number of regions (rows). These dataframes may also contain
    information regarding premature boundaries for each region (marked as -1)
    and no-density regions (marked by nan). This cleans the dataframe.
    """
    density = density.fillna(0) # NaNs are regions which contain zero density
    return density.replace(-1, np.nan) # -1 are regions which should not be counted at all

def remove_outliers(rbpdataframe, conf = 0.95):
    logger.info("Removing outliers (keep {})".format(conf))
    means = list()
    sems = list()
    for key, value in rbpdataframe.iteritems():
        df = rbpdataframe[key].dropna()
        
        nums = len(df)
        droppercent = (1-conf)/2.0
        dropnum = int(nums*(droppercent))
        df = df.sort_values()
        if(dropnum>0):
            df = df[dropnum:-dropnum]
        
        means.append(df.mean())
        sems.append(df.sem())
    logger.info("Finished removing outliers (keep {})".format(conf))
    return means, sems

def pdf_entropy(density, input_density, 
                pseudocount, ipseudocount,
                min_density_threshold = 0):
    """
    Return the entropy of pdf of each position
    Logic: 
        (in calc pdf):Fill NaNs with zero - we want to count all regions and add pseudocount
        (in calc pdf):Fill -1 with NaNs - we want to negate any -1, which signifies a premature exon boundary
        (in calc pdf):Add minimum pseudocount
        (in calc pdf):Calculate PDF
        Calculate entropy
    """
    logger.info("Starting normalization (pdf_entropy)")
    df_indices = density.index
    dfi_indices = input_density.index
    missing = set(df_indices) - set(dfi_indices)
    
    input_density = input_density.append(input_density.ix[missing])
    
    pdf = calculate_pdf(density, pseudocount, min_density_threshold)
    input_pdf = calculate_pdf(input_density, ipseudocount, min_density_threshold)
    
    en = pdf.multiply(np.log2(pdf.div(input_pdf)))
    logger.info("Finished normalization (pdf_entropy)")
    return en
def read_entropy(density, input_density, 
                 pseudocount, ipseudocount,
                 min_density_threshold = 0):
    """
    Return the entropy of each position
    Logic: 
        Fill NaNs with zero - we want to count all regions and add pseudocount
        Fill -1 with NaNs - we want to negate any -1, which signifies a premature exon boundary
        Add minimum pseudocount
        Calculate entropy
    """
    logger.info("Starting normalization (read_entropy)")
    df_indices = density.index
    dfi_indices = input_density.index
    missing = set(df_indices) - set(dfi_indices)
    
    input_density = input_density.append(input_density.ix[missing])
    
    rpm = clean(density)
    rpmi = clean(input_density)
    
    r = (rpm)/(pseudocount)
    ri = (rpmi)/(ipseudocount)
    
    r = r + 1
    ri = ri + 1
    
    en = r.multiply(np.log2(r.div(ri)))
    logger.info("Finished normalization (read_entropy)")
    return en
def pdf_read_entropy(density, input_density, 
                     pseudocount, ipseudocount,
                     min_density_threshold = 0):
    """
    Logic: 
        Calculate the read_entropy
        Take PDF
    """
    logger.info("Starting normalization (pdf_read_entropy)")
    en = read_entropy(density, input_density, 
                      pseudocount, ipseudocount,
                      min_density_threshold)
    pdf = en.div(en.sum(axis=1), axis=0)
    logger.info("Starting normalization (pdf_read_entropy)")
    return pdf

def get_density(density, input_density, 
                pseudocount, ipseudocount,
                min_density_threshold = 0):
    logger.info("Returning density")
    return density

def get_input(density, input_density, 
              pseudocount, ipseudocount,
              min_density_threshold = 0):
    logger.info("Returning input density")
    return input_density

def normalize_and_subtract(density, input_density, 
                           pseudocount, ipseudocount,
                           min_density_threshold = 0):
    logger.info("Starting normalization (pdf and subtraction)")
    pdf = calculate_pdf(density, pseudocount, min_density_threshold)
    input_pdf = calculate_pdf(input_density, ipseudocount, min_density_threshold)
        
    subtracted = pd.DataFrame(pdf.mean() - input_pdf.mean()).T
    logger.info("Starting normalization (pdf and subtraction)")
    return subtracted

def normalize_and_per_region_subtract(density, input_density, 
                                      pseudocount, ipseudocount, 
                                      min_density_threshold = 0):
    """
    Normalizes ip matrix of m x n (where m is the row of each event in a feature,
    and n is the column relating to nucleotide position). 
    """
    logger.info("Starting normalization (per region subtraction)")
    df_indices = density.index
    dfi_indices = input_density.index
    missing = set(df_indices) - set(dfi_indices)
    
    input_density = input_density.append(input_density.ix[missing])
    
    pdf = calculate_pdf(density, pseudocount, min_density_threshold)
    # pdf.to_csv('/Users/brianyee/git/encode/encode/rbpmaps/testfiles/rbfox2/outputs/ip_pdf.csv')
    pdfi = calculate_pdf(input_density, ipseudocount, min_density_threshold)
    # pdfi.to_csv('/Users/brianyee/git/encode/encode/rbpmaps/testfiles/rbfox2/outputs/input_pdf.csv')
    subtracted = pdf.sub(pdfi)
    logger.info("Starting normalization (per region subtraction)")
    return subtracted

def calculate_pdf(density, pseudocount = None, min_density_threshold = 0):
    """
    Calculates the PDF of a density matrix.
    Logic:
    
    Args: 
        density (pandas.DataFrame) : r x c matrix of densities. May contain
            NaN corresponding to values in which no density was returned.
            These values should be counted.
            May contain -1 corresponding to values in which a particular
            region is shorter than the full DataFrame length. These 
            values should not be counted.
        min_density_threshold (integer) : minimum total density across
            a row. (Deprecated - may be removed in the future)
    
    Returns:
        pdf (pandas.DataFrame) : r x c matrix of densities normalized
            across each respective (r)ow as a probability density func.
    """
    df = clean(density)
    min_read = pseudocount if pseudocount else min([item for item in df.unstack().values if item > 0])

    df = df + min_read
    pdf = df.div(df.sum(axis=1), axis=0)
    return pdf # , mean, sem