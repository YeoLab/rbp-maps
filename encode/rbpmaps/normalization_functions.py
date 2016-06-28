'''
Created on Jun 27, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import os

def KLDivergence(density, input_density):
        PDF_CONST = 1.0/len(density.columns)
        
        pdf = calculate_pdf(density,'ip')
        input_pdf = calculate_pdf(input_density,'input')
        
        pdft = pd.merge(pdf,input_pdf,how='left',left_index=True,right_index=True).fillna(PDF_CONST)
        
        pdf = pdft.filter(regex='\d+_x')
        pdfi = pdft.filter(regex='\d+_y')
        pdf = pdf.rename(columns=lambda x: x.replace('_x', ''))
        pdfi = pdfi.rename(columns=lambda x: x.replace('_y', ''))

        en = pdf.multiply(np.log2(pdf.div(pdfi)))
        
        return pdf, pdfi, en, en.mean()
    
def calculate_pdf(density, min_density_threshold, prefix, base = None):
    densities = density.replace(-1, np.nan)   
    df = densities[densities.sum(axis=1) > min_density_threshold]
    min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    df = df + min_normalized_read_number
    pdf = df.div(df.sum(axis=1), axis=0)
        
    if(base):
        pdf.to_csv("{}.{}.pdf.csv".format(base,prefix))
            
    return pdf # , mean, sem
    
def normalize(density, min_density_threshold, prefix, base = None):
    """
    This is identical to calculate_pdf.
    """
    pdf = calculate_pdf(density, min_density_threshold)
        
    if(base):
        pdf.to_csv("{}.{}.pdf.csv".format(base,prefix))
    return pdf 
    
def normalize_and_subtract(density, input_density, min_density_threshold, base = None):
        
    pdf = calculate_pdf(density,'ip')
    input_pdf = calculate_pdf(input_density,'input')
        
    subtracted = pdf.mean() - input_pdf.mean()
        
    if(base):
            
        subtracted.to_csv("{}.means_subtracted.csv".format(base))
            
    return subtracted
    
def normalize_and_per_region_subtract(density, input_density, min_density_threshold, base = None):
        
    PDF_CONST = 1.0/len(density.columns)
        
    pdf = calculate_pdf(density, 'ip')
    input_pdf = calculate_pdf(input_density, 'input')
        
    pdft = pd.merge(pdf,input_pdf, how='left',left_index=True,right_index=True).fillna(PDF_CONST)
        
    pdf = pdft.filter(regex='\d+_x')
    pdfi = pdft.filter(regex='\d+_y')
    
    
    pdf = pdf.rename(columns=lambda x: x.replace('_x', ''))
    pdfi = pdfi.rename(columns=lambda x: x.replace('_y', ''))
        
    subtracted = pdf.sub(pdfi)
        
    if(base):
        base = os.path.join(base)
            
        subtracted.to_csv("{}.means_per_region_subtracted.csv".format(base))
            
    return subtracted