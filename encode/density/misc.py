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

def isin(row,lst):
    for g in row['gene']:
        if g in lst:
            return True
    return False

def ensembl_from_gencode(gencode_id):
    return gencode_id.split('.')[0]

# returns True if key combinations exist in a dictionary, False otherwise
def exists(dictionary, *args):
    if args in dictionary:
        return True
    else:
        return False
# auto initializes a dictionary with key to 0 value otherwise increments
def ini(dictionary, *args):
    if args in dictionary:
        # if 499 in args and 'upex' in args:
        #    print("incrementing position by 1")
        return dictionary[args]+1
    else:
        # if 499 in args and 'upex' in args:
        #    print("initializing position")
        return 1