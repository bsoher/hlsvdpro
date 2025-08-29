# Python modules
from __future__ import division
import ctypes
import math
import pprint
import re
import pdb
pp = pprint.pprint

# 3rd party modules
import numpy as np
import scipy
import scipy.linalg.blas
import scipy.linalg.lapack

from scipy.linalg import cython_lapack



def dbdsqr(uplo, n, ncvt, nru, ncc, ddd, eee, vt, ldvt, uuu, ldu, ccc, ldc):

    # prep step to access the DBDSQR function in cython_lapack
    #
    # SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )

    try:
        n = n[0]  # convenience, since only used locally
    except:
        n = n

    ctypes.pythonapi.PyCapsule_GetPointer.restype  =  ctypes.c_void_p
    ctypes.pythonapi.PyCapsule_GetPointer.argtypes = [ctypes.py_object, ctypes.c_char_p]
    ctypes.pythonapi.PyCapsule_GetName.restype     =  ctypes.c_char_p
    ctypes.pythonapi.PyCapsule_GetName.argtypes    = [ctypes.py_object]
    
    ptr_type = ctypes.CFUNCTYPE(ctypes.c_void_p, 
                                ctypes.c_char_p, 
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_double), # * n),
                                ctypes.POINTER(ctypes.c_double), # * (n-1)),
                                ctypes.POINTER(ctypes.c_double), # * (ldvt * ncvt)),
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_double), # * (ldu * n)),
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_double), # * (ldc * ncc)),
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_double), # * (n * 4)),
                                ctypes.POINTER(ctypes.c_int))
    dbdsqr2 = ptr_type(
                 ctypes.pythonapi.PyCapsule_GetPointer(
                     cython_lapack.__pyx_capi__['dbdsqr'],
                     ctypes.pythonapi.PyCapsule_GetName(
                         cython_lapack.__pyx_capi__['dbdsqr'])))

    uplo = uplo.upper()

    xd = np.ndarray([n,], dtype='d', order='F')
    for i, value in enumerate(ddd[:n]):
        xd[i] = value

    xe = np.ndarray([n-1,], dtype='d', order='F')
    for i, value in enumerate(eee[:n-1]):
        xe[i] = value

    xvt = np.zeros([ldvt*ncvt,], dtype='d', order='F')
    if ncvt:
        for i, value in enumerate(vt):
            xvt[i] = value

    xu = np.zeros([ldu*n,], dtype='d', order='F')
    if nru:
        for i, value in enumerate(uuu[:ldu * n]):
            xu[i] = value

    xc = np.ndarray([ldc * ncc,], dtype='d', order='F')
    if ncc:
        for i in range(ldc * ncc):
            value = ccc[i]
            xc[i] = value

    xwork = np.ndarray([n*4,], dtype='d', order='F')

    info = 0

    uplo = ctypes.c_char(uplo)
    n    = ctypes.c_int(n)
    ncvt = ctypes.c_int(ncvt)
    nru  = ctypes.c_int(nru)
    ncc  = ctypes.c_int(ncc)
    ldvt = ctypes.c_int(ldvt)
    ldu  = ctypes.c_int(ldu)
    ldc  = ctypes.c_int(ldc)
    info = ctypes.c_int(info)

    args = (ctypes.byref(uplo), 
            ctypes.byref(n), 
            ctypes.byref(ncvt), 
            ctypes.byref(nru), 
            ctypes.byref(ncc),
            xd.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            xe.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            xvt.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.byref(ldvt),
            xu.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.byref(ldu),
            xc.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.byref(ldc),
            xwork.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(info))

    dbdsqr2(*args)

    info = info.value
    ncvt = ncvt.value
    nru  = nru.value
    ncc  = ncc.value
    n    = n.value
 
    # Derefernce pointers

    ddd = [item for item in xd]
    eee = [item for item in xe]
    
    vt  = None if vt  is None else [item for item in xvt]
    uuu = None if uuu is None else [item for item in xu]
    ccc = None if ccc is None else [item for item in xc]

    # FIXME - probably want to turn these into numpy arrays and reshape them.

    return info, ddd, eee, vt, uuu, ccc