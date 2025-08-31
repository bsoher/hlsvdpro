# Python modules
from __future__ import division
import math
import pdb
import pprint
from unittest import signals
pp=pprint.pprint

# 3rd party modules
import numpy as np
import scipy.linalg
import scipy.sparse.linalg
import scipy.linalg.lapack as lapack

# Our modules
import pypropack



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


# NOTE NOTE NOTE
# this module is basically the same as an import/call to hlsvdpro (just older/messier)

def hlsvdpro_propack(signals, nsv_sought, M=None):

    xx   = signals
    K    = nsv_sought
    mode = 'f'          # or 'b'
    N = len(xx)

    if M is None:
        M = int(N/2)

    L = N - M - 1

    # SVD of data matrix and truncation of U to form Uk
    if mode == "f":
        X = scipy.linalg.hankel(xx[:L + 1], xx[L:])
    else:
        # for backward LP we need to make the hankel matrix:
        X = scipy.linalg.hankel(xx[:M - 1:-1], xx[M::-1])

        # x_N-1 x_N-2 ... x_N-M-1
        # x_N-2 x_N-3 ... x_N-M-2
        # ...
        # x_M   x_M-1 ... x_0

    U, s, Vh = pypropack.svdp_aprod(X, K)
    Uk = np.mat(U[:, :K])   # trucated U matrix of rank K
    Ub = Uk[:-1]            # Uk with bottom row removed
    Ut = Uk[1:]             # Uk with top row removed

    Zp, resid, rank, ss = scipy.linalg.lstsq(Ub, Ut)

    #===============================================================
    # Diagonalization of Z' (=hx), yields 'signal'poles' aka 'roots'
    #   Eigenvalues are returned unordered. I sort them here to make
    #   it easier to compare them with output from the Fortran code.

    roots = scipy.linalg.eigvals(Zp)
    roots = np.array(sorted(roots))
    
    #===============================================================
    # Calculation of dampings (damp) and frequencies (freq) from roots

    dampings    = np.log(np.abs(roots))
    frequencies = np.arctan2(roots.imag, roots.real) / (math.pi * 2)

    #===============================================================
    #  Calculation of complex-valued amplitudes , using the
    #  pseudoinverse of the Lrow*kfit Vandermonde matrix zeta.

    zeta = np.vander(roots, N=len(signals), increasing=True).T
    v2, x2, s2, rank2, _, info2 = lapack.zgelss(zeta, signals)

    #===============================================================
    # Discard the uneeded values of x.

    x2 = x2[:K]
    amplitudes = np.abs(x2)
    phases = np.arctan2(x2.imag, x2.real)

    return K, s, frequencies, dampings, amplitudes, phases, U, Vh







def test():


    A = [0.54881350 + 0.31179588j, 0.71518937 + 0.69634349j, 0.60276338 + 0.37775184j, 0.54488318 + 0.17960368j, 0.42365480 + 0.02467873j, 0.64589411 + 0.06724963j, 0.43758721 + 0.67939277j, 0.89177300 + 0.45369684j, 0.96366276 + 0.53657921j, 0.38344152 + 0.89667129j]
    A = np.array(A)
    A = A.astype(np.complex128)

    sing0 = np.array([3.844931797977525, 1.0606943021737258, 0.8539664440754421, 0.6082379328434804, 0.47619445709816216], dtype=float)

    k = 5
    M = 4

    kfound, sing, freq, damp, ampl, phas, u, vt = hlsvdpro_propack(A, k, M)
    
    if kfound == k:
        print(' Results hlsvdpro_propack, max(abs(sing-sing0)) = ', max(abs(sing-sing0)))
    else:
        print(' Warning, only '+str(kfound)+' singular values found, k requested = '+str(k))


if __name__ == "__main__":
    test()