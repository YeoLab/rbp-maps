#!/usr/env python

import os
import pytest
import pybedtools
import pandas as pd
from density import Feature

### Fixtures ###

se_file = 'test_features/RBFOX2-BGHLV26-HepG2-SE.MATS.JunctionCountOnly.txt'
ri_file = 'test_features/RBFOX2-BGHLV26-HepG2-RI.MATS.JunctionCountOnly.txt'
a3ss_file = 'test_features/RBFOX2-BGHLV26-HepG2-A3SS.MATS.JunctionCountOnly.txt'
a5ss_file = 'test_features/RBFOX2-BGHLV26-HepG2-A5SS.MATS.JunctionCountOnly.txt'
mxe_file = 'test_features/RBFOX2-BGHLV26-HepG2-MXE.MATS.JunctionCountOnly.txt'


