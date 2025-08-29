import numpy as np
import ctypes as ct
from scipy.linalg import cython_lapack

# This is a working example of how to access DBDSDC() from the Cython LAPACK
# library, since it is not accessible from scipy.linalg.lapack module.
#
# It was built of an example from:
#
#  https://gist.github.com/insertinterestingnamehere/b50390fc720e1e554af0
#
# that accesses the DGEMM function from cython_blas module.  I've just modified
# it to access the DBDSDC function in cython_lapack.  And suprisingly it works. 
#
#
# EXAMPLE from:  https://www.nag.com/numeric/fl/nagdoc_fl26.2/html/f08/f08mdf.html
# to test if the code is working properly ...
#
#  This example computes the singular value decomposition of the upper bidiagonal matrix B, where
#
#       /                           \
#       | 3.62   1.26   0.00   0.00 |
#   B = | 0.00  -2.41  -1.53   0.00 |
#       | 0.00   0.00   1.92   1.19 |
#       | 0.00   0.00   0.00  -1.42 |
#       \                           /
#
# INPUTS from https://www.nag.com/numeric/fl/nagdoc_fl26.2/examples/data/f08mefe.d.html
#        and  https://www.nag.com/numeric/fl/nagdoc_fl26.2/examples/source/f08mefe.f90.html
#
#
# F08MEF Example Program Data
#   4                           :Value of N
#   3.62  -2.41   1.92  -1.43                         this is 'd' matrix
#   1.26  -1.53   1.19          :End of matrix B      this is 'e' matrix
#   'U'                         :Value of UPLO
#
# RESULTS from https://www.nag.com/numeric/fl/nagdoc_fl26.2/examples/baseresults/f08mefe.r.html
#
# Singular values
#     4.0001  3.0006  1.9960  0.9998
#
# Right singular vectors, by row
#          1       2       3       4
# 1   0.8261  0.5246  0.2024  0.0369
# 2   0.4512 -0.4056 -0.7350 -0.3030
# 3   0.2823 -0.5644  0.1731  0.7561
# 4   0.1852 -0.4916  0.6236 -0.5789
#
# Left singular vectors, by column
#          1       2       3       4
# 1   0.9129  0.3740  0.1556  0.0512
# 2  -0.3935  0.7005  0.5489  0.2307
# 3   0.1081 -0.5904  0.6173  0.5086
# 4  -0.0132  0.1444 -0.5417  0.8280


ct.pythonapi.PyCapsule_GetPointer.restype  =  ct.c_void_p
ct.pythonapi.PyCapsule_GetPointer.argtypes = [ct.py_object, ct.c_char_p]
ct.pythonapi.PyCapsule_GetName.restype     =  ct.c_char_p
ct.pythonapi.PyCapsule_GetName.argtypes    = [ct.py_object]

# SUBROUTINE DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, WORK, IWORK, INFO )

ptr_type = ct.CFUNCTYPE(ct.c_void_p, 
                        ct.c_char_p, 
                        ct.c_char_p, 
                        ct.POINTER(ct.c_int), 
                        ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_int), 
                        ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_int), 
                        ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_int), 
                        ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_int), 
                        ct.POINTER(ct.c_int))

dbdsdc2 = ptr_type(
             ct.pythonapi.PyCapsule_GetPointer(
                 cython_lapack.__pyx_capi__['dbdsdc'],
                 ct.pythonapi.PyCapsule_GetName(
                     cython_lapack.__pyx_capi__['dbdsdc'])))




n    = 4
ldu  = n
ldvt = n

d  = np.array([3.62, -2.41, 1.92, -1.43], 'd', order='F')
e  = np.array([1.26, -1.53, 1.19],        'd', order='F')
vt = np.ndarray([ldvt,n],                 'd', order='F')
u  = np.ndarray([ldu,n],                  'd', order='F')
q  = np.ndarray([1,],                     'd', order='F')   # hard coded to 1 in example
iq = np.ndarray([1,],                     'i', order='F')   # hard coded to 1 in example

work  = np.ndarray([n*(3*n+4),], 'd', order='F')
iwork = np.ndarray([n*(3*n+4),], 'i', order='F')

uplo  = ct.c_char(b'U')
compq = ct.c_char(b'I')
n     = ct.c_int(n)
ldu  = ct.c_int(ldu)
ldvt = ct.c_int(ldvt)
info = ct.c_int(0)

# DBDSDC( UPLO, COMPQ, N, D, E, U, LDU, VT, LDVT, Q, IQ, WORK, IWORK, INFO )


args = (ct.byref(uplo), 
        ct.byref(compq), 
        ct.byref(n), 
        d.ctypes.data_as(ct.POINTER(ct.c_double)),
        e.ctypes.data_as(ct.POINTER(ct.c_double)),
        u.ctypes.data_as(ct.POINTER(ct.c_double)),
        ct.byref(ldu),
        vt.ctypes.data_as(ct.POINTER(ct.c_double)), 
        ct.byref(ldvt),
        q.ctypes.data_as(ct.POINTER(ct.c_double)), 
        iq.ctypes.data_as(ct.POINTER(ct.c_int)), 
        work.ctypes.data_as(ct.POINTER(ct.c_double)),
        iwork.ctypes.data_as(ct.POINTER(ct.c_int)),
        ct.byref(info))

dbdsdc2(*args)


print("info = ", str(info))
print(" ")
print("singular values")
print(d)
print(" ")
print("Right singular vectors, by row ")
print(vt)
print(" ")
print("Left singular vectors, by column ")
print(u)

bob = 10
bob += 1
#print a.dot(b)