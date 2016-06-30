'''
Created on Jun 27, 2016

@author: brianyee
'''
import pandas as pd
import numpy as np
import os

def KLDivergence(density, input_density, min_density_threshold = 0):
        PDF_CONST = 1.0/len(density.columns)
        
        pdf = calculate_pdf(density,min_density_threshold)
        input_pdf = calculate_pdf(input_density,min_density_threshold)
        
        pdft = pd.merge(pdf,input_pdf,how='left',left_index=True,right_index=True).fillna(PDF_CONST)
        
        pdf = pdft.filter(regex='\d+_x')
        pdfi = pdft.filter(regex='\d+_y')
        pdf = pdf.rename(columns=lambda x: x.replace('_x', ''))
        pdfi = pdfi.rename(columns=lambda x: x.replace('_y', ''))

        en = pdf.multiply(np.log2(pdf.div(pdfi)))
        
        return en
    
def calculate_pdf(density, min_density_threshold = 0):
    densities = density.replace(-1, np.nan)   
    df = densities[densities.sum(axis=1) > min_density_threshold]
    min_normalized_read_number = min([item for item in df.unstack().values if item > 0])
    df = df + min_normalized_read_number
    pdf = df.div(df.sum(axis=1), axis=0)
  
    return pdf # , mean, sem
    
def normalize(density, min_density_threshold = 0):
    """
    This is identical to calculate_pdf.
    """
    pdf = calculate_pdf(density, min_density_threshold)

    return pdf 
    
def normalize_and_subtract(density, input_density, min_density_threshold = 0):
        
    pdf = calculate_pdf(density,min_density_threshold)
    input_pdf = calculate_pdf(input_density,min_density_threshold)
        
    subtracted = pd.DataFrame(pdf.mean() - input_pdf.mean()).T
      
    return subtracted
    
def normalize_and_per_region_subtract(density, input_density, min_density_threshold = 0):
        
    PDF_CONST = 1.0/len(density.columns)
        
    pdf = calculate_pdf(density, min_density_threshold)
    input_pdf = calculate_pdf(input_density, min_density_threshold)
        
    pdft = pd.merge(pdf,input_pdf, how='left',left_index=True,right_index=True).fillna(PDF_CONST)
        
    pdf = pdft.filter(regex='\d+_x')
    pdfi = pdft.filter(regex='\d+_y')
    
    
    pdf = pdf.rename(columns=lambda x: x.replace('_x', ''))
    pdfi = pdfi.rename(columns=lambda x: x.replace('_y', ''))
        
    subtracted = pdf.sub(pdfi)
  
    return subtracted