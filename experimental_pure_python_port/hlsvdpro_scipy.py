# Python modules
from __future__ import division
import math
import pprint
pp=pprint.pprint

# 3rd party modules
import numpy as np
import numpy.linalg
import scipy.linalg
import scipy.sparse.linalg
import scipy.linalg.lapack as lapack

# Our modules


MAX_SINGULAR_VALUES = 50

# Complex LAPACK functions
zgelss, = scipy.linalg.lapack.get_lapack_funcs( ['gelss'], 
                                                np.array([1j])
                                              )

# Variable name translations 
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# kuser      nsv_sought         # of singular values requested by caller
# kfit       nsv_found          # of singular values found
# lsinval    singular_values    Array containing the singular values (floats)
# ndp        n_data_points      
# kmax       MAX_SINGULAR_VALUES  Hardcoded to 50.
# lrow       n_rows             # of rows in the critical work matrix 
# mcol       n_columns          # of columns in the critical work matrix 


def hlsvdpro_scipy(data, nsv_sought, m):

    xx   = data
    k    = nsv_sought
    mode = 'f'          # or 'b'

    n = len(xx)
    L = n - m - 1

    if mode == "f":
        x = scipy.linalg.hankel(xx[:L + 1], xx[L:])
    else:
        # for backward LP we need to make the hankel matrix:
        # x_N-1 x_N-2 ... x_N-m-1
        # x_N-2 x_N-3 ... x_N-m-2
        # ...
        # x_M   x_M-1 ... x_0
        x = scipy.linalg.hankel(xx[:m - 1:-1], xx[m::-1])    # SVD of data matrix and truncation of U to form Uk

    #----------------------------------------------------------------
    # bjs - trying a bunch of scipy options here for SVD

    u, s, vh = scipy.linalg.svd(x, full_matrices=False)  # cp5 100reps -> 7.9sec (best) vs 1.16 for hlsvdpro (fortran call)
#    u, s, vh = scipy.linalg.svd(x, full_matrices=False, lapack_driver='gesvd')   # same 10.2sec
#    u, s, vh = scipy.sparse.linalg.svds(x, k=k)         # slow x2 18.9sec
#    u, s, vh = scipy.linalg.svd(x, full_matrices=True)   # slow x4 30.8

    uk = np.mat(u[:, :k])   # trucated U matrix of rank k
    ub = uk[:-1]            # Uk with bottom row removed
    ut = uk[1:]             # Uk with top row removed

    zp, resid, rank, ss = scipy.linalg.lstsq(ub, ut)

    #----------------------------------------------------------------
    #  Diagonalization of Z' (=hx), yields 'signal'poles' aka 'roots'

    roots = scipy.linalg.eigvals(zp)

    # Eigenvalues are returned unordered. I sort them here to make it easier
    # to compare them with output from the Fortran code.

    roots = np.array(sorted(roots))

    #----------------------------------------------------------------
    #  Calculation of dampings (damp) and frequencies (freq) from roots

    dampings    = np.log(np.abs(roots))
    frequencies2 = np.arctan2(roots.imag, roots.real) / (math.pi * 2)

    #----------------------------------------------------------------
    #  Calculation of complex-valued amplitudes , using the
    #  pseudoinverse of the Lrow*kfit Vandermonde matrix zeta.

    zeta = np.vander(roots, N=len(data), increasing=True).T
    v2, x1, s22, rank, _, info = lapack.zgelss(zeta, data)

    # Discard the uneeded values of x.

    x1 = x1[:k]
    amplitudes = np.abs(x1)
    phases = np.arctan2(x1.imag, x1.real)

    return k, s[0:k], frequencies2, dampings, amplitudes, phases, u, vh







