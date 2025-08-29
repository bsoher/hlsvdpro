#
#  (C) Brian J Soher, 2020
#

from __future__ import division, print_function

import numpy as np

from zlarnv import zlarnv
from zreorth import zreorth
from zblasext import pdznrm2

KAPPA = 0.717

USE_FORTRAN_RANDOMS = True


# def _get_predictable_randoms(size):
#     """Given an array size, reads from the file first_random.txt and turns
#     the values therein into a complex array of length size. The array is
#     returned.
#     """
#     lines = open('first_random.txt').read().split('\n')
#     lines = ['complex' + line for line in lines]
#     lines = map(eval, lines)
#     lines = lines[:size]
#
#     return np.array(lines)



def zgetu0(transa, m, n, j, ntry, u0, U, ldu, aprod, parm, icgs):
    """
    DGETU0: Attempt to generate a pseudo-random vector in SPAN(Op(A)) 
       orthogonal to span(U(:,1:j)), where Op(A) = A if transa='n' and
       Op(A) = A^H if transa='c'.
    
    # u0norm, anormest, and ierr are returned
    # u0 (a vector) is altered

    """
    transa = transa.lower()

    if transa == 'n':
        # %-------------------------%
        # | u0 is to be an m-vector |
        # %-------------------------%
        rsize = n
        usize = m
    else:
        # %-------------------------%
        # | u0 is to be an n-vector |
        # %-------------------------%
        rsize = m
        usize = n

    iseed = np.array([1, 3, 5, 7], dtype='i', order='F')
    idist = 2
    ierr  = 0

    set_ierr = True
    zworkr = np.zeros((rsize,), np.complex128)

    for itry in range(ntry):
        
        if USE_FORTRAN_RANDOMS:
            #zworkr = _get_predictable_randoms(rsize)
            zlarnv(idist, iseed, rsize, zworkr)
        else:
            reals = np.random.uniform(-1.0, 1.0, rsize)
            imags = np.random.uniform(-1.0, 1.0, rsize)
            zworkr = np.array( [complex(r, i) for r, i in zip(reals, imags)] )

        nrm = pdznrm2(rsize, zworkr, 1)

        aprod(transa,m,n,zworkr,u0,parm)

        u0norm   = pdznrm2(usize,u0,1)
        anormest = u0norm / nrm

        if j >= 0:        # from v2.1
            index = [1, j, j + 1]
            u0norm = zreorth(usize,j,U,ldu,u0,u0norm,index,KAPPA,icgs)

        if u0norm > 0:
            set_ierr = False
            break


    if set_ierr:
        ierr = -1

    return u0norm, anormest, ierr


