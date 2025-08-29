# Python modules
from __future__ import division
import ctypes

# 3rd party modules
import numpy as np

from scipy.linalg import cython_lapack

# NB. This version assumes that the user passes in properly sized arrays!
#     I do no array size checking here. However there is some code commentd
#     out that could be added at a later time.


def dbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info):

    #  DBDSQR( UPLO,N,NCVT,NRU,NCC,D,E,VT,LDVT,U,LDU,C,LDC,WORK,INFO )
    #
    #       .. Scalar Arguments ..
    #       CHARACTER          UPLO
    #       INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
    #       ..
    #       .. Array Arguments ..
    #       DOUBLE PRECISION   C( LDC, * ), D( * ), E( * ), U( LDU, * ),
    #                          VT( LDVT, * ), WORK( * )
    #
    #  Purpose:
    #  =============
    # 
    #  DBDSQR computes the singular values and, optionally, the right and/or
    #  left singular vectors from the singular value decomposition (SVD) of
    #  a real N-by-N (upper or lower) bidiagonal matrix B using the implicit
    #  zero-shift QR algorithm.  The SVD of B has the form
    # 
    #     B = Q * S * P**T
    # 
    #  where S is the diagonal matrix of singular values, Q is an orthogonal
    #  matrix of left singular vectors, and P is an orthogonal matrix of
    #  right singular vectors.  If left singular vectors are requested, this
    #  subroutine actually returns U*Q instead of Q, and, if right singular
    #  vectors are requested, this subroutine returns P**T*VT instead of
    #  P**T, for given real input matrices U and VT.  When U and VT are the
    #  orthogonal matrices that reduce a general matrix A to bidiagonal
    #  form:  A = U*B*VT, as computed by DGEBRD, then
    # 
    #     A = (U*Q) * S * (P**T*VT)
    # 
    #  is the SVD of A.  Optionally, the subroutine may also compute Q**T*C
    #  for a given real input matrix C.
    # 
    #  See "Computing  Small Singular Values of Bidiagonal Matrices With
    #  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
    #  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
    #  no. 5, pp. 873-912, Sept 1990) and
    #  "Accurate singular values and differential qd algorithms," by
    #  B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
    #  Department, University of California at Berkeley, July 1992
    #  for a detailed description of the algorithm.
    #  \endverbatim
    #
    #  Arguments:
    #  ==========
    #
    # UPLO
    #           UPLO is CHARACTER*1
    #           = 'U':  B is upper bidiagonal;
    #           = 'L':  B is lower bidiagonal.
    # N
    #           N is INTEGER
    #           The order of the matrix B.  N >= 0.
    # NCVT
    #           NCVT is INTEGER
    #           The number of columns of the matrix VT. NCVT >= 0.
    # NRU
    #           NRU is INTEGER
    #           The number of rows of the matrix U. NRU >= 0.
    # NCC
    #           NCC is INTEGER
    #           The number of columns of the matrix C. NCC >= 0.
    # D
    #           D is DOUBLE PRECISION array, dimension (N)
    #           On entry, the n diagonal elements of the bidiagonal matrix B.
    #           On exit, if INFO=0, the singular values of B in decreasing
    #           order.
    # E
    #           E is DOUBLE PRECISION array, dimension (N-1)
    #           On entry, the N-1 offdiagonal elements of the bidiagonal
    #           matrix B.
    #           On exit, if INFO = 0, E is destroyed; if INFO > 0, D and E
    #           will contain the diagonal and superdiagonal elements of a
    #           bidiagonal matrix orthogonally equivalent to the one given
    #           as input.
    # VT
    #           VT is DOUBLE PRECISION array, dimension (LDVT, NCVT)
    #           On entry, an N-by-NCVT matrix VT.
    #           On exit, VT is overwritten by P**T * VT.
    #           Not referenced if NCVT = 0.
    # LDVT
    #           LDVT is INTEGER
    #           The leading dimension of the array VT.
    #           LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
    # U
    #           U is DOUBLE PRECISION array, dimension (LDU, N)
    #           On entry, an NRU-by-N matrix U.
    #           On exit, U is overwritten by U * Q.
    #           Not referenced if NRU = 0.
    # LDU
    #           LDU is INTEGER
    #           The leading dimension of the array U.  LDU >= max(1,NRU).
    # C
    #           C is DOUBLE PRECISION array, dimension (LDC, NCC)
    #           On entry, an N-by-NCC matrix C.
    #           On exit, C is overwritten by Q**T * C.
    #           Not referenced if NCC = 0.
    # LDC
    #           LDC is INTEGER
    #           The leading dimension of the array C.
    #           LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
    # WORK
    #           WORK is DOUBLE PRECISION array, dimension (4*(N-1))
    # INFO
    #           INFO is INTEGER
    #           = 0:  successful exit
    #           < 0:  If INFO = -i, the i-th argument had an illegal value
    #           > 0:
    #              if NCVT = NRU = NCC = 0,
    #                 = 1, a split was marked by a positive value in E
    #                 = 2, current block of Z not diagonalized after 30*N
    #                      iterations (in inner while loop)
    #                 = 3, termination criterion of outer while loop not met
    #                      (program created more than N unreduced blocks)
    #              else NCVT = NRU = NCC = 0,
    #                    the algorithm did not converge; D and E contain the
    #                    elements of a bidiagonal matrix which is orthogonally
    #                    similar to the input matrix B;  if INFO = i, i
    #                    elements of E have not converged to zero.
    #
    # Internal Parameters:
    # =========================
    # 
    #   TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))
    #           TOLMUL controls the convergence criterion of the QR loop.
    #           If it is positive, TOLMUL*EPS is the desired relative
    #              precision in the computed singular values.
    #           If it is negative, abs(TOLMUL*EPS*sigma_max) is the
    #              desired absolute accuracy in the computed singular
    #              values (corresponds to relative accuracy
    #              abs(TOLMUL*EPS) in the largest singular value.
    #           abs(TOLMUL) should be between 1 and 1/EPS, and preferably
    #              between 10 (for fast convergence) and .1/EPS
    #              (for there to be some accuracy in the results).
    #           Default is to lose at either one eighth or 2 of the
    #              available decimal digits in each computed singular value
    #              (whichever is smaller).
    # 
    #   MAXITR  INTEGER, default = 6
    #           MAXITR controls the maximum number of passes of the
    #           algorithm through its inner loop. The algorithms stops
    #           (and so fails to converge) if the number of passes
    #           through the inner loop exceeds MAXITR*N**2.


    # prep step to access the DBDSQR function in cython_lapack
    #
    # DBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U, LDU, C, LDC, WORK, INFO )

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

    # FIXME - assumption here is that user passes in Fortran oriented arrays
    #         that are of the proper size. AND that they are Numpy arrays, so
    #         that they can be modified in place to return output values.

    uplo = uplo.upper()

#     xd = np.ndarray([n,], dtype='d', order='F')
#     for i, value in enumerate(ddd[:n]):
#         xd[i] = value
# 
#     xe = np.ndarray([n-1,], dtype='d', order='F')
#     for i, value in enumerate(eee[:n-1]):
#         xe[i] = value
# 
# 
#     xvt = np.zeros([ldvt*ncvt,], dtype='d', order='F')
#     if ncvt:
#         for i, value in enumerate(vt):
#             xvt[i] = value
# 
# 
#     xu = np.zeros([ldu*n,], dtype='d', order='F')
#     if nru:
#         for i, value in enumerate(uuu[:ldu * n]):
#             xu[i] = value
# 
# 
#     xc = np.ndarray([ldc * ncc,], dtype='d', order='F')
#     if ncc:
#         for i in range(ldc * ncc):
#             value = ccc[i]
#             xc[i] = value
# 
#     xwork = np.ndarray([n*4,], dtype='d', order='F')
# 
#     info = 0

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
            d.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            e.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            vt.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.byref(ldvt),
            u.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.byref(ldu),
            c.ctypes.data_as(ctypes.POINTER(ctypes.c_double)), 
            ctypes.byref(ldc),
            work.ctypes.data_as(ctypes.POINTER(ctypes.c_double)),
            ctypes.byref(info))

    dbdsqr2(*args)

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
    ncvt = ncvt.value
    nru  = nru.value
    ncc  = ncc.value
    n    = n.value
 
    # Derefernce pointers

#     d[:] = np.array([item for item in xd], dtype='d')
#     e[:] = np.array([item for item in xe], dtype='d')
#     
#     vt[:] = None if vt is None else np.array([item for item in xvt], dtype='d').reshape(ldvt,ncvt)
#     u[:]  = None if u  is None else np.array([item for item in xu],  dtype='d').reshape(ldu,n)
#     c[:]  = None if c  is None else np.array([item for item in xc],  dtype='d').reshape(ldc,ncc)


    return 


def test():

    uplo = b'U'
    n    = 4
    ldc  = 1
    ldu  = n
    ldvt = n
    ncc  = 0
    ncvt = n
    nru  = n    # from input to dbdsqr()
    info = 0
    
    d  = np.array([3.62, -2.41, 1.92, -1.43], 'd', order='F')
    e  = np.array([1.26, -1.53, 1.19],        'd', order='F')
    vt = np.eye(  n,                    dtype='d', order='F')
    u  = np.eye(  n,                    dtype='d', order='F')
    c  = np.zeros([ldc, 1],                   'd', order='F')   # hard coded to ldc,1 in example
    
    work = np.zeros([4*n,], 'd', order='F')

    dbdsqr(uplo, n, ncvt, nru, ncc, d, e, vt, ldvt, u, ldu, c, ldc, work, info)

    print("DBDSQR info = ", str(info))
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
