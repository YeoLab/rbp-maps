#!/bin/env python
# encoding: utf-8
"""

This module contains methods for selecting regions given a single feature.

Main Functions
--------------

- five_prime_site : given a ReadDensity object and an interval (optionally
    its neighboring upstream interval), returns the read density across
    the 5' (exon) splice site
- three_prime_site : given a ReadDensity object and an interval (optionally
    its neighboring downstream interval), returns the read density across
    the 3' (exon) splice site
- generic_site : given a ReadDensity object and an interval, returns the
    density across a given interval as described by a single BedTool.

Created on May 3, 2016

@author: brianyee
"""

import itertools
import numpy as np
import pandas as pd
import sys


def split(lst, n):
    """
    Splits list (lst) into n equal parts.

    Parameters
    ----------
    lst : list
    n : int

    Returns
    -------
    newlist : list
        a list of equally portioned n sublists
    """
    newlist = []
    division = len(lst) / float(n)
    for i in xrange(n):
        newlist.append(
            lst[int(round(division * i)):int(round(division * (i + 1)))])
    return newlist
