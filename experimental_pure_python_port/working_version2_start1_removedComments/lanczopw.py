# Python modules
from __future__ import division
import ctypes
import math
import pdb

# 3rd party modules
import numpy as np

# Our modules
import zlansvdw

# Variable name translations 
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# kuser      nsv_sought         # of singular values requested by caller
# kfit       nsv_found          # of singular values found
# lsinval    singular_values    Array containing the singular values (floats)
# ndp        n_data_points      
# kmax       MAX_SINGULAR_VALUES  Hardcoded to 50.
# lrow       n_rows             # of rows in the very important uuu matrix.
#                               Once inside zlansvdw(), this is called 'm'.
# mcol       n_columns          # of rows? in the very important vvv matrix.
#                               Once inside zlansvdw(), this is called 'n'.
#                               Since it counts rows, I'm not sure why it was
#                               called 'mcol'.


def lanczopw(signals, n_rows, n_columns, nsv_sought):

    n_data_points = len(signals)

    # computation of the first column of the circulant matrix
    
    # PS - There's a couple of odd loops here that set up fvect. I think the 'f'
    # might stand for 'fased' (phased)? These loops basically shove all the
    # points to the right by (ndp/2) + 1.

    fvect   = np.roll(signals.copy(), (n_data_points // 2) + 1)
    lambda_ = np.fft.fft(fvect)
    fvect   = np.conjugate(fvect)

    # PS - Reverse elements of fvect, except for first element which isn't changed.

    fvect = np.concatenate( (np.array([fvect[0]]), fvect[1:].copy()[::-1]))

    trlambda = np.fft.fft(fvect)

    # PS - lambda and trlambda never change from this point on.
    lambda_  = lambda_  / n_data_points
    trlambda = trlambda / n_data_points

    uuu  = None
    done = False

    while not done:
        
        res = zlansvdw.zlansvdw(n_data_points, n_rows, n_columns, nsv_sought, lambda_, trlambda, uuu)
        uuu, singular_values, info, nsv_found = res
            
        if info == -1:
            if nsv_found > 0:
                # Need to call zlansvdw() again.
                nsv_sought = nsv_found
                # print "Python: calling zlansvdw again, nsv_sought = {}".format(nsv_sought)
                pass
            else:
                # FIXME - raise an error; PROPACK didn't converge.
                done = True
        else:
            done = True

    return uuu, singular_values, nsv_found

