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

# Our modules


MAX_SINGULAR_VALUES = 50

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


def hlsvdpro_scipy(signals, nsv_sought, M):

    # M - Hankel Matrix size in pts

    xx   = signals
    K    = nsv_sought
    mode = 'f'          # or 'b'

    N = len(xx)
    L = N - M - 1

    if mode == "f":
        X = scipy.linalg.hankel(xx[:L + 1], xx[L:])
    else:
        # for backward LP we need to make the hankel matrix:
        # x_N-1 x_N-2 ... x_N-M-1
        # x_N-2 x_N-3 ... x_N-M-2
        # ...
        # x_M   x_M-1 ... x_0
        X = scipy.linalg.hankel(xx[:M - 1:-1], xx[M::-1])    # SVD of data matrix and truncation of U to form Uk
    
#    U2,  s2,  Vh2  = scipy.linalg.svd(X, full_matrices=False)
#    Uk1 = np.mat(U1[:, :K])   # trucated U matrix of rank K
#    Ub1 = Uk1[:-1]            # Uk with bottom row removed
#    Ut1 = Uk1[1:]             # Uk with top row removed
#    
#    Zp1, resid1, rank1, ss1 = scipy.linalg.lstsq(Ub1, Ut1)
    
    U2, s2, Vh2 = scipy.sparse.linalg.svds(X, k=K)
    Uk2 = np.mat(U2[:, :K])   # trucated U matrix of rank K
    Ub2 = Uk2[:-1]            # Uk with bottom row removed
    Ut2 = Uk2[1:]             # Uk with top row removed

    Zp2, resid2, rank2, ss2 = scipy.linalg.lstsq(Ub2, Ut2)

#  Diagonalization of Z' (=hx), yields 'signal'poles' aka 'roots'
#  ===============================================================

#    roots1 = scipy.linalg.eigvals(Zp1)
    roots2 = scipy.linalg.eigvals(Zp2)

    # Eigenvalues are returned unordered. I sort them here to make it easier
    # to compare them with output from the Fortran code.
#    roots1 = np.array(sorted(roots1))
    roots2 = np.array(sorted(roots2))
    

#  Calculation of dampings (damp) and frequencies (freq) from roots
#  ===============================================================

#    dampings1    = np.log(np.abs(roots1))
#    frequencies1 = np.arctan2(roots1.imag, roots1.real) / (math.pi * 2)
    dampings2    = np.log(np.abs(roots2))
    frequencies2 = np.arctan2(roots2.imag, roots2.real) / (math.pi * 2)

#  Calculation of complex-valued amplitudes , using the 
#  pseudoinverse of the Lrow*kfit Vandermonde matrix zeta.
#  ===============================================================

#  First calculate zeta:

#    zeta1 = _vanmon(len(signals), roots1)
    zeta2 = _vanmon(len(signals), roots2)

# c     zgells writes solution in space of vector x,
# c     but the vector 'signal' is preserved. 
# c     5-8-99: rcond was -1.0; g77 makes rank = 0!! Not with SUN.
# c             So I set rcond = +1.0d-10!!!!
# c

    lwork = 2*MAX_SINGULAR_VALUES+64*(len(signals)+MAX_SINGULAR_VALUES)

#    v11, x11, s11, rank11, _, info11 = zgelss(zeta1, signals, cond=-1.0, lwork=lwork,
#                                     overwrite_a=False, overwrite_b=False)

    v22, x22, s22, rank22, _, info22 = zgelss(zeta2, signals, cond=-1.0, lwork=lwork,
                                     overwrite_a=False, overwrite_b=False)




    # FIXME this code, like the Fortran, ignores possible errors reported 
    # in the "info" return value.


    # print "------- V  ----------"
    # print v.shape, v
    # print "------- X  ----------"
    # print x.shape, x
    # print "------- S  ----------"
    # print s
    # print "------- Rank  ----------"
    # print rank
    # print "------- info  ----------"
    # print info

    # Discard the uneeded values of x.

#    x11 = x11[:K]
#    amplitudes11 = np.abs(x11)
#    phases11 = np.arctan2(x11.imag, x11.real)

    x22 = x22[:K]
    amplitudes22 = np.abs(x22)
    phases22 = np.arctan2(x22.imag, x22.real)


    #return K, s1, frequencies1, dampings1, amplitudes11, phases11
    return K, s2, frequencies2, dampings2, amplitudes22, phases22




def _vanmon(n_data_points, roots):
    """ calculates the ndp*kfit Vandermonde matrix zeta. """

    nsv_found = len(roots)

    zeta = np.zeros( (n_data_points, nsv_found), np.complex128)

    zeta[1, :nsv_found] = (1+0j)

    for j in range(nsv_found):
        root = roots[j]
        #print "vanmon, roots[{}] = {:.17E}".format(j + 1, root)
        temp = (1+0j)
        for i in range(1, n_data_points):
            temp *= root

            zeta[i, j] = temp
            
            #print "zeta[{},{}] = {:.17E}".format(i + 1, j + 1, temp)

    return zeta


def _zcalc(nsv_found, n_rows, uuu):
    # c=======================================================================
    #       subroutine Zcalc(kfit,kmax,Lrow,U,zprime,us,unit)

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

    uot = 1.0 - sum_.real

    unit = np.zeros( (nsv_found, nsv_found), np.complex128)

    for i in range(nsv_found):
        for j in range(nsv_found):
            temp = ((1+0j) if j == i else (0+0j))
            unit[i, j] = temp + (uuu[m - 1][i].conjugate() * uuu[m - 1][j] / uot)
            #print "unit({}, {}) = {:.17E}".format(i + 1, j + 1, unit[i, j])


    zprime = np.zeros( (nsv_found, nsv_found), np.complex128)

    for i in range(nsv_found):
        for j in range(nsv_found):
            zprime[i, j] = (unit[i, :nsv_found] * uuu_sum[:nsv_found, j]).sum()

    return zprime
