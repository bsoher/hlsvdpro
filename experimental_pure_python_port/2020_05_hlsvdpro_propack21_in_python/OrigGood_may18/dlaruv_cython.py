#
#  (C) Brian J Soher, 2020
#

# Python modules
from __future__ import division
import ctypes

# 3rd party modules
import numpy as np

from scipy.linalg import cython_lapack

# NB. This version assumes that the user passes in properly sized arrays!
#     I do no array size checking here. However there is some code commentd
#     out that could be added at a later time.


def dlaruv(iseed, n, x):
    """
    SUBROUTINE DLARUV( ISEED, N, X )
    
          .. Scalar Arguments ..
          INTEGER            N
          ..
          .. Array Arguments ..
          INTEGER            ISEED( 4 )
          DOUBLE PRECISION   X( N )
          ..
    
    
    Purpose:
    =============
    
    DLARUV returns a vector of n random real numbers from a uniform (0,1)
    distribution (n <= 128).
    
    This is an auxiliary routine called by DLARNV and ZLARNV.

    Arguments:
    ==========

    ISEED
              ISEED is INTEGER array, dimension (4)
              On entry, the seed of the random number generator; the array
              elements must be between 0 and 4095, and ISEED(4) must be
              odd.
              On exit, the seed is updated.
    N
            N is INTEGER
            The number of random numbers to be generated. N <= 128.
    X
            X is DOUBLE PRECISION array, dimension (N)
            The generated random numbers.
    
    Further Details:
    =====================
    
    This routine uses a multiplicative congruential method with modulus
    2**48 and multiplier 33952834046453 (see G.S.Fishman,
    'Multiplicative congruential random number generators with modulus
    2**b: an exhaustive analysis for b = 32 and a partial analysis for
    b = 48', Math. Comp. 189, pp 331-344, 1990).
    
    48-bit integers are stored in 4 integer array elements with 12 bits
    per element. Hence the routine is portable across machines with
    integers of 32 bits or more.

    """


    # prep step to access the DLARUV function in cython_lapack
    #
    # DLARUV( ISEED, N, X )

    ctypes.pythonapi.PyCapsule_GetPointer.restype  =  ctypes.c_void_p
    ctypes.pythonapi.PyCapsule_GetPointer.argtypes = [ctypes.py_object, ctypes.c_char_p]
    ctypes.pythonapi.PyCapsule_GetName.restype     =  ctypes.c_char_p
    ctypes.pythonapi.PyCapsule_GetName.argtypes    = [ctypes.py_object]
    
    ptr_type = ctypes.CFUNCTYPE(ctypes.c_void_p, 
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_double))
    dlaruv2 = ptr_type(
                 ctypes.pythonapi.PyCapsule_GetPointer(
                     cython_lapack.__pyx_capi__['dlaruv'],
                     ctypes.pythonapi.PyCapsule_GetName(
                         cython_lapack.__pyx_capi__['dlaruv'])))

    # FIXME - assumption here is that user passes in Fortran oriented arrays
    #         that are of the proper size. AND that they are Numpy arrays, so
    #         that they can be modified in place to return output values.


#     xiseed = np.ndarray([4,], dtype='i', order='F')
#     for i, value in enumerate(iseed[:4]):
#         xiseed[i] = value
# 
#     xx = np.ndarray([n,], dtype='d', order='F')
#     for i, value in enumerate(x[:n]):
#         xx[i] = value
# 

    n = ctypes.c_int(n)

    args = (iseed.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.byref(n), 
            x.ctypes.data_as(ctypes.POINTER(ctypes.c_double)) 
           )

    dlaruv2(*args)

    n = n.value
 
    # Derefernce pointers

#     x[:] = np.array([item for item in xd], dtype='d')


    return 


#--------------------------------------------------------------------

def test():


    n     = 33      # NB. Must be <= 128!!!
    iseed = np.array([1,3,5,7], 'i', order='F')
    x     = np.ndarray([n,],      'd', order='F')
    
    print(" x before =", x)

    dlaruv(iseed, n, x)

    print("random values")
    print(x)
    print("len(x) = ", len(x))
    print(" ")
    print("Current seed values ")
    print(iseed)


if __name__ == "__main__":
    
    test()
