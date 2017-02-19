'''
Created on Jun 20, 2016

@author: brianyee
'''

import numpy as np


def toarray(dic):
    tmp = {}
    for key in dic.keys():
        tmp[key] = np.asarray(dic[key])
    return tmp
