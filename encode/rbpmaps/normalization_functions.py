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

def KLDivergence(density, input_density, min_density_threshold = 0):
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
    globally for input, add pseudocount of 1 read
    divide each position by 1,000,000
    do en
    plot mean
    """
    # print('Entropy of reads')
    ipdf = density.replace(-1, np.nan)
    indf = input_density.replace(-1, np.nan)
    
    # ipdf = density[density.sum(axis=1) > min_density_threshold]
    # indf = input_density[input_density.sum(axis=1) > min_density_threshold]
    
    min_ip_read_number = min([item for item in ipdf.unstack().values if item > 0])
    min_in_read_number = min([item for item in indf.unstack().values if item > 0])
    min_read_number = min(min_ip_read_number,min_in_read_number)
    
    ipdf = ipdf + min_read_number
    indf = indf + min_read_number
    
    ipdfdiv = ipdf # .div(1000000)
    indfdiv = indf # .div(1000000)
    
    dft = pd.merge(ipdfdiv,indfdiv,how='left',left_index=True,right_index=True)
    
    ipdfdiv = dft.filter(regex='\d+_x')
    indfdiv = dft.filter(regex='\d+_y')
    
    ipdfdiv = ipdfdiv.rename(columns=lambda x: x.replace('_x', ''))
    indfdiv = indfdiv.rename(columns=lambda x: x.replace('_y', '')).fillna(min_read_number)
    
    en = ipdfdiv.multiply(np.log2(ipdfdiv.div(indfdiv)))
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

def calculate_pdf(density, min_density_threshold = 0):
    # print("calculate PDF")
    df = density.replace(-1, np.nan)
    
    # df = densities[densities.sum(axis=1) > min_density_threshold]
    min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    df = df + min_normalized_read_number
    pdf = df.div(df.sum(axis=1), axis=0)
    return pdf # , mean, sem
    
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
    # print("normalization and per region subtraction")
    PDF_CONST = 1.0/len(density.columns)
        
    pdf = calculate_pdf(density, min_density_threshold)
    input_pdf = calculate_pdf(input_density, min_density_threshold)
        
    pdft = pd.merge(pdf,input_pdf, how='left',left_index=True,right_index=True)
        
    pdf = pdft.filter(regex='\d+_x')
    pdfi = pdft.filter(regex='\d+_y')
    
    
    pdf = pdf.rename(columns=lambda x: x.replace('_x', ''))
    pdfi = pdfi.rename(columns=lambda x: x.replace('_y', '')).fillna(PDF_CONST)
    
    subtracted = pdf.sub(pdfi)
    # print("TYPE AFTER PER REGION SUBTRACT: {}".format(type(subtracted)))
    return subtracted