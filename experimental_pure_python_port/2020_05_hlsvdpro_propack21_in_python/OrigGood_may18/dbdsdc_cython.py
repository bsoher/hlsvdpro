#
#  (C) Brian J Soher, 2020
#

# Python modules
from __future__ import division

import ctypes
import numpy as np

from scipy.linalg import cython_lapack

# NB. This version assumes that the user passes in properly sized arrays!
#     I do no array size checking here. However there is some code commentd
#     out that could be added at a later time.


def dbdsdc(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info):

    #  DBDSDC( UPLO,COMPQ,N,D,E,U,LDU,VT,LDVT,Q,IQ,WORK,IWORK,INFO )
    #
    # Purpose:
    # =============
    #
    # DBDSDC computes the singular value decomposition (SVD) of a real
    # N-by-N (upper or lower) bidiagonal matrix B:  B = U * S * VT,
    # using a divide and conquer method, where S is a diagonal matrix
    # with non-negative diagonal elements (the singular values of B), and
    # U and VT are orthogonal matrices of left and right singular vectors,
    # respectively. DBDSDC can be used to compute all singular values,
    # and optionally, singular vectors or singular vectors in compact form.
    #
    # This code makes very mild assumptions about floating point
    # arithmetic. It will work on machines with a guard digit in
    # add/subtract, or on those binary machines without guard digits
    # which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
    # It could conceivably fail on hexadecimal or decimal machines
    # without guard digits, but we know of none.  See DLASD3 for details.
    #
    # The code currently calls DLASDQ if singular values only are desired.
    # However, it can be slightly modified to compute singular values
    # using the divide and conquer method.
    #
    #  Arguments:
    #  ==========
    #
    # UPLO
    #          UPLO is CHARACTER*1
    #          = 'U':  B is upper bidiagonal.
    #          = 'L':  B is lower bidiagonal.
    #
    # COMPQ
    #          COMPQ is CHARACTER*1
    #          Specifies whether singular vectors are to be computed
    #          as follows:
    #          = 'N':  Compute singular values only;
    #          = 'P':  Compute singular values and compute singular
    #                  vectors in compact form;
    #          = 'I':  Compute singular values and singular vectors.
    #
    # N
    #          N is INTEGER
    #          The order of the matrix B.  N >= 0.
    #
    # D
    #          D is DOUBLE PRECISION array, dimension (N)
    #          On entry, the n diagonal elements of the bidiagonal matrix B.
    #          On exit, if INFO=0, the singular values of B.
    #
    # E
    #          E is DOUBLE PRECISION array, dimension (N-1)
    #          On entry, the elements of E contain the offdiagonal
    #          elements of the bidiagonal matrix whose SVD is desired.
    #          On exit, E has been destroyed.
    #
    # U
    #          U is DOUBLE PRECISION array, dimension (LDU,N)
    #          If  COMPQ = 'I', then:
    #             On exit, if INFO = 0, U contains the left singular vectors
    #             of the bidiagonal matrix.
    #          For other values of COMPQ, U is not referenced.
    #
    # LDU
    #          LDU is INTEGER
    #          The leading dimension of the array U.  LDU >= 1.
    #          If singular vectors are desired, then LDU >= max( 1, N ).
    #
    # VT
    #          VT is DOUBLE PRECISION array, dimension (LDVT,N)
    #          If  COMPQ = 'I', then:
    #             On exit, if INFO = 0, VT**T contains the right singular
    #             vectors of the bidiagonal matrix.
    #          For other values of COMPQ, VT is not referenced.
    #
    # LDVT
    #          LDVT is INTEGER
    #          The leading dimension of the array VT.  LDVT >= 1.
    #          If singular vectors are desired, then LDVT >= max( 1, N ).
    #
    # Q
    #          Q is DOUBLE PRECISION array, dimension (LDQ)
    #          If  COMPQ = 'P', then:
    #             On exit, if INFO = 0, Q and IQ contain the left
    #             and right singular vectors in a compact form,
    #             requiring O(N log N) space instead of 2*N**2.
    #             In particular, Q contains all the DOUBLE PRECISION data in
    #             LDQ >= N*(11 + 2*SMLSIZ + 8*INT(LOG_2(N/(SMLSIZ+1))))
    #             words of memory, where SMLSIZ is returned by ILAENV and
    #             is equal to the maximum size of the subproblems at the
    #             bottom of the computation tree (usually about 25).
    #          For other values of COMPQ, Q is not referenced.
    #
    # IQ
    #          IQ is INTEGER array, dimension (LDIQ)
    #          If  COMPQ = 'P', then:
    #             On exit, if INFO = 0, Q and IQ contain the left
    #             and right singular vectors in a compact form,
    #             requiring O(N log N) space instead of 2*N**2.
    #             In particular, IQ contains all INTEGER data in
    #             LDIQ >= N*(3 + 3*INT(LOG_2(N/(SMLSIZ+1))))
    #             words of memory, where SMLSIZ is returned by ILAENV and
    #             is equal to the maximum size of the subproblems at the
    #             bottom of the computation tree (usually about 25).
    #          For other values of COMPQ, IQ is not referenced.
    #
    # WORK
    #          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
    #          If COMPQ = 'N' then LWORK >= (4 * N).
    #          If COMPQ = 'P' then LWORK >= (6 * N).
    #          If COMPQ = 'I' then LWORK >= (3 * N**2 + 4 * N).
    #
    # IWORK
    #          IWORK is INTEGER array, dimension (8*N)
    #
    # INFO
    #          INFO is INTEGER
    #          = 0:  successful exit.
    #          < 0:  if INFO = -i, the i-th argument had an illegal value.
    #          > 0:  The algorithm failed to compute a singular value.
    #                The update process of divide and conquer failed.

    # prep step to access the DBDSQR function in cython_lapack
    #
    #  DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, WORK, IWORK, INFO )

    ctypes.pythonapi.PyCapsule_GetPointer.restype  =  ctypes.c_void_p
    ctypes.pythonapi.PyCapsule_GetPointer.argtypes = [ctypes.py_object, ctypes.c_char_p]
    ctypes.pythonapi.PyCapsule_GetName.restype     =  ctypes.c_char_p
    ctypes.pythonapi.PyCapsule_GetName.argtypes    = [ctypes.py_object]
    
    ptr_type = ctypes.CFUNCTYPE(ctypes.c_void_p, 
                                ctypes.c_char_p, 
                                ctypes.c_char_p, 
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_int), 
                                ctypes.POINTER(ctypes.c_int))
    
    dbdsdc2 = ptr_type(
                 ctypes.pythonapi.PyCapsule_GetPointer(
                     cython_lapack.__pyx_capi__['dbdsdc'],
                     ctypes.pythonapi.PyCapsule_GetName(
                         cython_lapack.__pyx_capi__['dbdsdc'])))

    uplo  = uplo.upper()
    compq = compq.upper()

#     xd = np.ndarray([n,], dtype='d', order='F')
#     for i, value in enumerate(d[:n]):
#         xd[i] = value
# 
#     xe = np.ndarray([n-1,], dtype='d', order='F')
#     for i, value in enumerate(e[:n-1]):
#         xe[i] = value
# 
#     xu = np.zeros([ldu*n,], dtype='d', order='F')
#     for i, value in enumerate(u.flatten()[:ldu*n]):
#         xu[i] = value
# 
#     xvt = np.zeros([ldvt*n,], dtype='d', order='F')
#     for i, value in enumerate(vt.flatten()[:ldvt*n]):
#         xvt[i] = value
# 
#     xq  = np.ndarray([2*n*n,],   'd', order='F')   # may not be used depend on INFO=0
#     xiq = np.ndarray([2*n*n,],   'i', order='F')   # may not be used depend on INFO=0
#     
#     xwork = np.ndarray([(3*n*n+4*n),], 'd', order='F')
#     for i, value in enumerate(work.flatten()[:(3*n*n+4*n)]):
#         xwork[i] = value
#         
#     xiwork = np.ndarray([n*8,],         'i', order='F')
#     for i, value in enumerate(iwork.flatten()[:n*8]):
#         xiwork[i] = value
# 
#     info = 0

    uplo  = ctypes.c_char(uplo)
    compq = ctypes.c_char(compq)
    n     = ctypes.c_int(n)
    ldu   = ctypes.c_int(ldu)
    ldvt  = ctypes.c_int(ldvt)
    info  = ctypes.c_int(info)

    args = (ctypes.byref(uplo), 
            ctypes.byref(compq), 
            ctypes.byref(n), 
            d.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            e.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            u.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(ldu),
            vt.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.byref(ldvt),
            q.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            iq.ctypes.data_as(ctypes.POINTER(ctypes.c_int)), 
            work.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            iwork.ctypes.data_as(ctypes.POINTER(ctypes.c_int)),
            ctypes.byref(info))
    
    dbdsdc2(*args)

    #     print("info = ", str(info))
    #     print(" ")
    #     print("singular values")
    #     print(xd)
    #     print(" ")
    #     print("Right singular vectors, by row ")
    #     print(xvt)
    #     print(" ")
    #     print("Left singular vectors, by column ")
    #     print(xu)

    info = info.value
    n    = n.value
 
#     # Dereference pointers - if we did array size checking
# 
#     d[:] = np.array([item for item in xd], dtype='d')
#     e[:] = np.array([item for item in xe], dtype='d')
#     
#     vt[:] = None if vt is None else np.array([item for item in xvt], dtype='d').reshape(ldvt,n)
#     u[:]  = None if u  is None else np.array([item for item in xu],  dtype='d').reshape(ldu,n)
#     q[:]  = None if q  is None else np.array([item for item in xq],  dtype='d')
#     iq[:] = None if iq is None else np.array([item for item in xiq], dtype='i')
 
    return 


#------------------------------------------------------------------------------

def test():
    
    uplo  = b'U'
    compq = b'I'
    n     = 4
    ldu   = n
    ldvt  = n
    info  = 0
    
    d  = np.array([3.62, -2.41, 1.92, -1.43], 'd', order='F')
    e  = np.array([1.26, -1.53, 1.19],        'd', order='F')
    u  = np.ndarray([ldu,n],                  'd', order='F')
    vt = np.ndarray([ldvt,n],                 'd', order='F')
    q  = np.ndarray([1,],                     'd', order='F')   # hard coded to 1 in example
    iq = np.ndarray([1,],                     'i', order='F')   # hard coded to 1 in example
    
    work  = np.ndarray([n*(3*n+4),], 'd', order='F')
    iwork = np.ndarray([n*(3*n+4),], 'i', order='F')    
    
    dbdsdc(uplo, compq, n, d, e, u, ldu, vt, ldvt, q, iq, work, iwork, info)

    #info, d, e, vt, u, q, iq = r

    print("DBDSDC info = ", str(info))
    print(" ")
    print("singular values")
    print(d)
    print(" ")
    print("Right singular vectors, by row ")
    print(vt)
    print(" ")
    print("Left singular vectors, by column ")
    print(u)

    # check results

    d0  = np.array([4.0001, 3.0006, 1.9960, 0.9998])
    vt0 = np.array([[ 0.8261,  0.5246,  0.2024,  0.0369],
                    [ 0.4512, -0.4056, -0.7350, -0.3030],
                    [ 0.2823, -0.5644,  0.1731,  0.7561],
                    [ 0.1852, -0.4916,  0.6236, -0.5789]])
    u0  = np.array([[ 0.9129,  0.3740,  0.1556,  0.0512],
                    [-0.3935,  0.7005,  0.5489,  0.2307],
                    [ 0.1081, -0.5904,  0.6173,  0.5086],
                    [-0.0132,  0.1444, -0.5417,  0.8280]])

    print('    SingVals all < 1e-4 - '+str(np.all( np.abs(d-d0) < 1e-4)))     # np.allclose() was giving odd answers
    print('Lft SingVals all < 1e-4 - '+str(np.all( np.abs(u-u0) < 1e-4)))
    print('Rgt SingVals all < 1e-4 - '+str(np.all( np.abs(vt-vt0) < 1e-4)))



if __name__ == "__main__":
    
    test()