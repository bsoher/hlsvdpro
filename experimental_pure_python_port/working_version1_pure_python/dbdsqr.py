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

# Our modules
import util

# Crashed on Win64 Miniconda2, so commented out
#
# # dbdsqr() is available via scipy on OS X, but not elsewhere.
# lib = ctypes.CDLL(scipy.linalg.lapack.clapack.__file__)
# if hasattr(lib, "dbdsqr"):
#     dbdsqr_function = lib.dbdsqr
# else:
#     libhlsvd = util.load_the_library()
#     dbdsqr_function = libhlsvd.dbdsqr_

libhlsvd = util.load_the_library()
dbdsqr_function = libhlsvd.dbdsqr_


def dbdsqr(uplo, n, ncvt, nru, ncc, ddd, eee, vt, ldvt, uuu, ldu, ccc, ldc):
    uplo = uplo.upper()

    # FIXME - uuu, ccc and vt can be passed in as 2D, but of course they will be
    # C-style row-major ordering and they need to be passed to Fortran as
    # column-major. At present I only have examples where the second dimension
    # of the 2-D arrays is 1, so this code will need to be fixed later when
    # I have true multi-dimension arrays.

    DddArrayType = ctypes.c_double * n
    ddd_fortran = DddArrayType()
    for i, value in enumerate(ddd[:n]):
        ddd_fortran[i] = value

    EeeArrayType = ctypes.c_double * (n - 1)
    eee_fortran = EeeArrayType()
    for i, value in enumerate(eee[:n - 1]):
        eee_fortran[i] = value


    if ncvt:
        VtArrayType = ctypes.c_double * (ldvt * ncvt)
        vt_fortran = VtArrayType()
        for i, value in enumerate(vt):
            #print "vt_fortran[{}] = {}".format(i + 1, value)
            vt_fortran[i] = value

        vt_fortran = ctypes.pointer(vt_fortran)
    else:
        vt_fortran = None

    if nru:
        UuuArrayType = ctypes.c_double * (ldu * n)
        uuu_fortran = UuuArrayType()
        for i, value in enumerate(uuu[:ldu * n]):
            uuu_fortran[i] = value

        uuu_fortran = ctypes.pointer(uuu_fortran)
    else:
        uuu_fortran = None


    if ncc:
        # print "ldc * ncc = {}, len(ccc) = {}".format(ldc * ncc, len(ccc))
        CccArrayType = ctypes.c_double * (ldc * ncc)
        ccc_fortran = CccArrayType()
        # for i, value in enumerate(ccc):
        for i in range(ldc * ncc):
            value = ccc[i]
            # print "ccc_fortran[{}] = {}".format(i + 1, value)
            ccc_fortran[i] = value

        ccc_fortran = ctypes.pointer(ccc_fortran)
    else:
        ccc_fortran = None

    WorkArrayType = ctypes.c_double * (n * 4)
    work = WorkArrayType()

    info = 0


    uplo = ctypes.c_char(uplo)
    n = ctypes.c_long(n)
    ncvt = ctypes.c_long(ncvt)
    nru = ctypes.c_long(nru)
    ncc = ctypes.c_long(ncc)
    ldvt = ctypes.c_long(ldvt)
    ldu = ctypes.c_long(ldu)
    ldc = ctypes.c_long(ldc)
    info = ctypes.c_long(info)

    dbdsqr_function(ctypes.pointer(uplo),
                     ctypes.pointer(n),
                     ctypes.pointer(ncvt),
                     ctypes.pointer(nru),
                     ctypes.pointer(ncc),
                     ctypes.pointer(ddd_fortran),
                     ctypes.pointer(eee_fortran),
                     vt_fortran,
                     ctypes.pointer(ldvt),
                     uuu_fortran,
                     ctypes.pointer(ldu),
                     ccc_fortran,
                     ctypes.pointer(ldc),
                     ctypes.pointer(work),
                     ctypes.pointer(info)
                    )

    info = info.value

    ncvt = ncvt.value
    nru = nru.value
    ncc = ncc.value
    n = n.value

    # Derefernce pointers
    if vt_fortran:
        vt_fortran = vt_fortran.contents
    if uuu_fortran:
        uuu_fortran = uuu_fortran.contents
    if ccc_fortran:
        ccc_fortran = ccc_fortran.contents

    # Convert these from ctypes arrays to ordinary Python lists.
    ddd = [value for value in ddd_fortran]

    # for i, value in enumerate(ddd):
    #     print "ddd[{}] = {}".format(i + 1, value)

    eee = [value for value in eee_fortran]
    # for i, value in enumerate(eee):
    #     print "eee[{}] = {}".format(i + 1, value)


    if vt_fortran:
        vt = [value for value in vt_fortran]
    else:
        vt = None

    if uuu_fortran:
        uuu = [value for value in uuu_fortran]

        # for i, value in enumerate(uuu):
        #     print "uuu[{}] = {}".format(i + 1, value)
    else:
        uuu = None

    if ccc_fortran:
        delta = len(ccc) - len(ccc_fortran)
        tail = ccc[-delta:].tolist()
        ccc = [value for value in ccc_fortran] + tail

        # for i, value in enumerate(ccc):
        #     print "ccc[{}] = {}".format(i + 1, value)
    else:
        ccc = None

    # FIXME - probably want to turn these into numpy arrays and reshape them.

    return info, ddd, eee, vt, uuu, ccc




def dbdsqr_cython(uplo, n, ncvt, nru, ncc, ddd, eee, vt, ldvt, uuu, ldu, ccc, ldc):

    # prep step to access the DBDSQR function in cython_lapack
    #
    # SUBROUTINE DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )

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

    # FIXME - uuu, ccc and vt can be passed in as 2D, but of course they will be
    # C-style row-major ordering and they need to be passed to Fortran as
    # column-major. At present I only have examples where the second dimension
    # of the 2-D arrays is 1, so this code will need to be fixed later when
    # I have true multi-dimension arrays.

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

    print("info = ", str(info))
    print(" ")
    print("singular values")
    print(xd)
    print(" ")
    print("Right singular vectors, by row ")
    print(xvt)
    print(" ")
    print("Left singular vectors, by column ")
    print(xu)


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
    
#     if vt_fortran:
#         vt_fortran = vt_fortran.contents
#     if uuu_fortran:
#         uuu_fortran = uuu_fortran.contents
#     if ccc_fortran:
#         ccc_fortran = ccc_fortran.contents
 
#     # Convert these from ctypes arrays to ordinary Python lists.
#     ddd = [value for value in ddd_fortran]
#  
#     # for i, value in enumerate(ddd):
#     #     print "ddd[{}] = {}".format(i + 1, value)
#  
#     eee = [value for value in eee_fortran]
#     # for i, value in enumerate(eee):
#     #     print "eee[{}] = {}".format(i + 1, value)
# 
#     if vt_fortran:
#         vt = [value for value in vt_fortran]
#     else:
#         vt = None
#  
#     if uuu_fortran:
#         uuu = [value for value in uuu_fortran]
#  
#         # for i, value in enumerate(uuu):
#         #     print "uuu[{}] = {}".format(i + 1, value)
#     else:
#         uuu = None
#  
#     if ccc_fortran:
#         delta = len(ccc) - len(ccc_fortran)
#         tail = ccc[-delta:].tolist()
#         ccc = [value for value in ccc_fortran] + tail
#  
#         # for i, value in enumerate(ccc):
#         #     print "ccc[{}] = {}".format(i + 1, value)
#     else:
#         ccc = None
 
    # FIXME - probably want to turn these into numpy arrays and reshape them.
 
    return info, ddd, eee, vt, uuu, ccc