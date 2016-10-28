'''
Created on Jun 27, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import os

def remove_outliers(rbpdataframe, conf = 0.95):
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
    return means, sems

def entropy(density, input_density, min_density_threshold = 0):
    # print("TYPE OF OBJECT: {}".format(type(density)))
    # print("KLDivergence (entropy of PDF)")
    PDF_CONST = 1.0/len(density.columns)
    
    density = density.replace(-1,np.nan)
    
    pdf = calculate_pdf(density,min_density_threshold)
    input_pdf = calculate_pdf(input_density,min_density_threshold)
        
    # pdft = pd.merge(pdf,input_pdf,how='left',left_index=True,right_index=True).fillna(PDF_CONST)
        
    # pdf = pdft.filter(regex='\d+_x')
    # pdfi = pdft.filter(regex='\d+_y')
    # pdf = pdf.rename(columns=lambda x: x.replace('_x', ''))
    # pdfi = pdfi.rename(columns=lambda x: x.replace('_y', ''))

    en = pdf.multiply(np.log2(pdf.div(input_pdf)))
    # print("TYPE AFTER ENTROPY OF READS: {}".format(type(en)))
    return en
def entropy_of_reads(density, input_density, min_density_threshold = 0):
    """
    Return the entropy of each position
    Logic: 
        Fill NaNs with zero - we want to count all regions and add pseudocount
        Fill -1 with NaNs - we want to negate any -1, which signifies a premature exon boundary
        Add minimum pseudocount
        Calculate entropy
    """
    df_indices = density.index
    dfi_indices = input_density.index
    missing = pd.index(set(df_indices) - set(dfi_indices))
    
    input_density = input_density.append(input_density.ix[missing])
    
    density = density.fillna(0) # NaNs are regions which contain zero density
    df = density.replace(-1, np.nan) # -1 are regions which should not be counted at all
    
    input_density = input_density.fillna(0) # NaNs are regions which contain zero density
    dfi = input_density.replace(-1, np.nan) # -1 are regions which should not be counted at all
    
    min_ip_read_number = min([item for item in df.unstack().values if item > 0])
    min_in_read_number = min([item for item in dfi.unstack().values if item > 0])
    min_read_number = min(min_ip_read_number,min_in_read_number)
    
    df = df + min_read_number
    dfi = dfi + min_read_number
    
    en = df.multiply(np.log2(df.div(dfi)))
    # print("TYPE AFTER ENTROPY OF READS: {}".format(type(en)))
    return en
def pdf_of_entropy_of_reads(density, input_density, min_density_threshold = 0):
    """
    globally for input, add pseudocount of 1 read
    divide each position by 1,000,000
    do en
    do pdf
    
    """

    en = entropy_of_reads(density, input_density, min_density_threshold)
    # min_normalized_read_number = abs(min([item for item in en.unstack().values if abs(item) > 0]))
    # en = en + min_normalized_read_number
    pdf = en.div(en.sum(axis=1), axis=0)
    # print("TYPE AFTER ENTROPY OF READS: {}".format(type(pdf)))
    return pdf

def get_density(density, input_density, min_density_threshold = 0):
    # df = density[density.sum(axis=1) > min_density_threshold]
    # print("TYPE AFTER DENSITY: {}".format(type(density)))
    return density

def get_input(density, input_density, min_density_threshold = 0):
    # df = input_density[input_density.sum(axis=1) > min_density_threshold]
    # print("TYPE AFTER INPUT: {}".format(type(input_density)))
    return input_density
    
def normalize(density, input_density, min_density_threshold = 0):
    """
    This is identical to calculate_pdf.
    """
    # print("Normalize (calculate PDF)")
    pdf = calculate_pdf(density, min_density_threshold)

    return pdf 

def baseline_rpm_mean_subtraction(density, input_density, min_density_threshold = 0):
    pass

def normalize_and_subtract(density, input_density, min_density_threshold = 0):
    # print("normalization and subtraction")

    pdf = calculate_pdf(density,min_density_threshold)
    input_pdf = calculate_pdf(input_density,min_density_threshold)
        
    subtracted = pd.DataFrame(pdf.mean() - input_pdf.mean()).T
    # print("TYPE AFTER normalize_and_subtract: {}".format(type(input_density)))
    return subtracted

def normalize_and_per_region_subtract(density, input_density, min_density_threshold = 0):
    """
    Normalizes ip matrix of m x n (where m is the row of each event in a feature,
    and n is the column relating to nucleotide position). 
    """
    # PDF_CONST = 1.0/len(density.columns)
    
    dft = pd.merge(density,input_density, how='outer',left_index=True,right_index=True)
    
    dfx = dft.filter(regex='\d+_x')
    dfy = dft.filter(regex='\d+_y')
    
    pdf = calculate_pdf(dfx, min_density_threshold)
    pdfi = calculate_pdf(dfy, min_density_threshold)
    
    pdf = pdf.rename(columns=lambda x: x.replace('_x', ''))
    pdfi = pdfi.rename(columns=lambda x: x.replace('_y', ''))
    
    # pdfi = fill_all_nans_with_minpdf(pdfi, PDF_CONST)
    # pdf = fill_all_nans_with_minpdf(pdf, PDF_CONST)
    
    subtracted = pdf.sub(pdfi)
    # print("TYPE AFTER PER REGION SUBTRACT: {}".format(type(subtracted)))
    return subtracted

def calculate_pdf(density, min_density_threshold = 0):
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
    density = density.fillna(0) # NaNs are regions which contain zero density
    df = density.replace(-1, np.nan) # -1 are regions which should not be counted at all
    
    # df = df[df.sum(axis=1) > min_density_threshold]
    min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    df = df + min_normalized_read_number
    pdf = df.div(df.sum(axis=1), axis=0)
    return pdf # , mean, sem