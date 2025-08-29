# Python modules
from __future__ import division
import math
import pdb
import pprint
pp=pprint.pprint

# 3rd party modules
import numpy as np
import scipy.linalg

# Our modules
import lanczopw

MAX_SINGULAR_VALUES = 50

# Max # of data points is set to 8192 in the Fortran. That's only necessary
# because of Fortran's lack (at the time) of dynamic array allocation. In 
# Python there's no such restriction.
# NDPMAX = 8192

# Complex LAPACK functions

zgelss, = scipy.linalg.lapack.get_lapack_funcs( ['gelss'], np.array([1j]) )

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


def hlsvdpro(signals, nsv_sought, step_size):

    n_data_points = len(signals)
    n_columns = n_data_points // 2
    n_rows = n_columns + 1


    # PS Making a copy of the signal array is necessary in Fortran because 
    # zgelss() would overwrite signal otherwise. In Python's version of zgelss
    # we have the option not to overwrite the input array so copying signals
    # is not necessary.

    #========================================== 
    # Lanczos SVD of hankel datamatrix.
    #==========================================

    # PS - lanczopw() is called in a loop until it returns a satisfactory 
    # answer. I moved that loop from this code into the call to lanczopw().
    
    uuu, singular_values, nsv_found = lanczopw.lanczopw(signals, n_rows, n_columns, nsv_sought)

    #========================================== 
    # Calculation of the Z' matrix, with ...
    # =========================================

    zprime = _zcalc(nsv_found, n_rows, uuu)


    #=============================================================================
    # Diagonalization of Z' (=hx), yielding the 'signal'poles', also called roots.
    #=============================================================================

    # zgeev - procedure, which calculates eigvals and eigvecs
    #         of a general matrix - only eigvals are used and computed
    #         into variable root
  
    #       call zgeev('N','N',kfit,zprime,kmax,root,U,Lrow,U,Lrow, cwork,Lcwork,rwork,info)

    # This is a LAPACK routine
    # zgeev('N', 'N', kfit, zprime, kmax, root, U, lrow, U, lrow, cwork, lcwork, rwork, info)

    roots, _ = np.linalg.eig(zprime)

    # Eigenvalues are returned unordered. I sort them here to make it easier
    # to compare them with output from the Fortran code.
    
    roots = np.array(sorted(roots))



    #=================================================================
    # Calculation of dampings (damp) and frequencies (freq) from roots
    #=================================================================
    
    dampings = np.log(np.abs(roots))

    frequencies = np.arctan2(roots.imag, roots.real) / (math.pi * 2)


    #=================================================================
    # Calculation of complex-valued amplitudes , using the 
    # pseudoinverse of the Lrow*kfit Vandermonde matrix zeta.
    # 
    # First calculation of zeta:
    #=================================================================
    
    zeta = _vanmon(len(signals), roots)

    # zgells writes solution in space of vector x, but the vector 'signal' is preserved. 
    # 5-8-99: rcond was -1.0; g77 makes rank = 0!! Not with SUN.
    #         So I set rcond = +1.0d-10!!!!
    # 
    # call zgelss(ndp,kfit,1,zeta,ndp,x,ndp, sinval,1.0d-10,rank,cwork,Lcwork,rwork,info)
    #
    # v,x,s,rank,info = zgelss(a,b,cond=-1.0,lwork=2*minmn+MAX(maxmn,nrhs),overwrite_a=0,overwrite_b=0)
    
    lwork = 2*MAX_SINGULAR_VALUES+64*(len(signals)+MAX_SINGULAR_VALUES)
    v, x, s, rank, _, info = zgelss(zeta, signals, 
                                     #cond=1e-10, 
                                     cond=-1.0,
                                     lwork=lwork,
                                     overwrite_a=False, overwrite_b=False)

    # FIXME this code, like the Fortran, ignores possible errors reported in the "info" return value.

    # Discard the uneeded values of x.
    
    x = x[:nsv_found]

    amplitudes = np.abs(x)
    phases = np.arctan2(x.imag, x.real)


    return nsv_found, singular_values, frequencies, dampings, amplitudes, phases




def _vanmon(n_data_points, roots):
    ''' vanmon(.) calculates the ndp*kfit Vandermonde matrix zeta '''

    nsv_found = len(roots)

    zeta = np.zeros( (n_data_points, nsv_found), np.complex128)

    zeta[1, :nsv_found] = (1+0j)

    for j in range(nsv_found):
        root = roots[j]
        temp = (1+0j)
        for i in range(1, n_data_points):
            temp *= root
            zeta[i, j] = temp

    return zeta


def _zcalc(nsv_found, n_rows, uuu):
    
    # PS - zprime is the output of this function
    # PS - unit is only used in this function
    # PS - us is only used in this function

    uuu_sum = np.zeros( (nsv_found, nsv_found), np.complex128)

    m = n_rows

    for i in range(nsv_found):
        for j in range(nsv_found):
    
            # PS This is the Fortran equivalent, but slow --
            # sum_ = (0+0j)
            # for k in range(m - 1):
            #     sum_ += uuu[k][i].conjugate() * uuu[k + 1][j]

            # PS - This is the fast way to do the same --
            # FIXME - is there an off by one error here? Double check sums
            # against Fortran.
            
            sum_ = (uuu[:m - 1, i].conjugate() * uuu[1:m, j]).sum()

            uuu_sum[i][j] = sum_


    # PS this is the simple Python equivalent of the fortran loop --
    # sum_ = (0+0j)
    # for i in range(kfit):
    #     sum_ += uuu[m - 1, i].conjugate() * uuu[m - 1, i]

    # Here's the fast way to do it -- 

    sum_ = (uuu[m - 1, :nsv_found].conjugate() * uuu[m - 1, :nsv_found]).sum()
    uot  = 1.0 - sum_.real

    unit = np.zeros( (nsv_found, nsv_found), np.complex128)

    for i in range(nsv_found):
        for j in range(nsv_found):
            temp = ((1+0j) if j == i else (0+0j))
            unit[i, j] = temp + (uuu[m - 1][i].conjugate() * uuu[m - 1][j] / uot)

    zprime = np.zeros( (nsv_found, nsv_found), np.complex128)

    for i in range(nsv_found):
        for j in range(nsv_found):
            zprime[i, j] = (unit[i, :nsv_found] * uuu_sum[:nsv_found, j]).sum()
            #print "zprime({}, {}) = {:.17E}".format(i + 1, j + 1, zprime[i, j])

    return zprime
