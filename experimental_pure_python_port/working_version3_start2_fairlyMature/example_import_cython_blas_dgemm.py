import numpy as np
import ctypes as ct
from scipy.linalg import cython_blas

#ct.cdll.LoadLibrary("libopenblas.dll")
#openblas = ct.CDLL("libopenblas.dll")
#dgemm = openblas.dgemm

ct.pythonapi.PyCapsule_GetPointer.restype = ct.c_void_p
ct.pythonapi.PyCapsule_GetPointer.argtypes = [ct.py_object, ct.c_char_p]
ct.pythonapi.PyCapsule_GetName.restype = ct.c_char_p
ct.pythonapi.PyCapsule_GetName.argtypes = [ct.py_object]
ptr_type = ct.CFUNCTYPE(ct.c_void_p, ct.c_char_p, ct.c_char_p,
                        ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                        ct.POINTER(ct.c_int), ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_double), ct.POINTER(ct.c_int),
                        ct.POINTER(ct.c_double), ct.POINTER(ct.c_int),
                        ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_int))
dgemm2 = ptr_type(
             ct.pythonapi.PyCapsule_GetPointer(
                 cython_blas.__pyx_capi__['dgemm'],
                 ct.pythonapi.PyCapsule_GetName(
                     cython_blas.__pyx_capi__['dgemm'])))

a = np.array([[1,2],[3,4]], 'd', order='F')
b = np.array([[5,6],[7,8]], 'd', order='F')
c = np.empty((2,2), order='F')

transa = ct.c_char(b'N')
transb = ct.c_char(b'N')
alpha = ct.c_double(1.)
beta = ct.c_double(0.)
lda = ct.c_int(2)
ldb = ct.c_int(2)
ldc = ct.c_int(2)
m = ct.c_int(2)
n = ct.c_int(2)
k = ct.c_int(2)
args = (ct.byref(transa), ct.byref(transb), ct.byref(m), ct.byref(n),
        ct.byref(k), ct.byref(alpha),
        a.ctypes.data_as(ct.POINTER(ct.c_double)), ct.byref(lda),
        b.ctypes.data_as(ct.POINTER(ct.c_double)), ct.byref(ldb),
        ct.byref(beta), c.ctypes.data_as(ct.POINTER(ct.c_double)),
        ct.byref(ldc))
#dgemm(*args)
#print( c)
c[:] = 0
dgemm2(*args)
print( c)
print( a.dot(b))