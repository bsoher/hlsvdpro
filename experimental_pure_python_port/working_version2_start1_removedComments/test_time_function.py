#!/usr/bin/env python

# Python modules
from __future__ import division
import cProfile
import os
import math

# Third party modules
import numpy as np



def dlapy1(a, b):
    return math.sqrt((a ** 2) + (b ** 2))


def dlapy2(a, b):
    return np.hypot(a,b)


def test(reps):
    
    a = 12.5
    b = 31.8
    
    for i in range(reps):
#         r = math.sqrt((a ** 2) + (b ** 2))
#         r = np.hypot(a,b)
#         r = dlapy1(a,b)
        r = dlapy2(a,b)
        
        
        



if __name__ == "__main__":

    if os.path.exists("profile2.data"):
        os.remove("profile2.data")

    #test(filename)
     
    reps = 100000
     
    cProfile.run('test(reps)', 'profile2.data')
    import pstats as ps
    p = ps.Stats('profile2.data')
    p.strip_dirs().sort_stats('cumulative').print_stats()