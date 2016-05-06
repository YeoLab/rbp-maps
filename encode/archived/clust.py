'''
Created on Jan 22, 2016

@author: brianyee
'''
import matplotlib
matplotlib.use('Agg')
import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import networkx as nx
from matplotlib import pyplot as plt

try:
    from networkx import graphviz_layout
except ImportError:
    raise ImportError("This example needs Graphviz and either PyGraphviz or Pydot")

# concatenates each RBP histogram into a clusterable file
def stack(names, dist_list, binned=False, percent=False):
    fs = pd.DataFrame()
    x = []
    for i in range(0, len(dist_list)):
        f = pd.read_table(dist_list[i],names=[names[i]])
        x.append(names[i])
        if(binned==True):
            f = pd.Series(f[names[i]])
            b = binit(f)
            if(percent==True):
                fs = pd.concat([fs,pd.DataFrame(b)/sum(b)],axis=1)
            else:
                fs = pd.concat([fs,pd.DataFrame(b)],axis=1)
            
        else:
            fs = pd.concat([fs,f],axis=1)
    fs.columns = x
    return fs

# rough bin func
def binit(f, n=100):
    a = 0
    b = []
    every = len(f)/n
    for i in range(1,len(f)+1):
        a = a + f[i-1]
        if i % every == 0 and i != 0:
            b.append(a)
            a = 0

    return b

# generates a heatmap of positions (rows = rbp, cols = position)
def heatmap(names,dist_list,output_file):
    cmap = sns.cubehelix_palette(as_cmap=True, rot=-.3, light=1)
    percentages = stack(names,dist_list,True,True).T
    sns.clustermap(percentages,col_cluster = False, cmap=cmap, linewidths=0.5, xticklabels=False)
    plt.savefig(output_file)

# saves file for R import
def get_peak_counts_table(names,dist_list,output_file):
    dat = stack(names, dist_list, False, False)
    dat.to_csv(output_file,sep="\t",index=False)
    print("finished!")
    
def main():
    o = []
    names = open('names.txt','r')
    outfiles = open('outfiles.txt','r')
    lines = outfiles.readlines()
    for i in range(0,len(lines)):
        o.append(lines[i])
    print o
        
    """names = ['U2AF1', 'FMR1','3','4','5','6','7','8','9','10']
    dist_list = ['/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/242_01.pv_0.05.fc_0.txt', '/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/FMR1_01.pv_0.05.fc_0.txt','/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/242_01.pv_0.05.fc_0.txt', '/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/FMR1_01.pv_0.05.fc_0.txt','/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/242_01.pv_0.05.fc_0.txt', '/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/FMR1_01.pv_0.05.fc_0.txt','/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/242_01.pv_0.05.fc_0.txt', '/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/FMR1_01.pv_0.05.fc_0.txt','/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/242_01.pv_0.05.fc_0.txt', '/Users/brianyee/Documents/workspace/encode_clip/rbpmaps/FMR1_01.pv_0.05.fc_0.txt']
    
    heatmap(names,dist_list)
    # stack(names,dist_list,True,True)"""
if __name__ == '__main__':
    sys.exit(main())