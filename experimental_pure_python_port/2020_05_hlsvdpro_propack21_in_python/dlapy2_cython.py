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


def dlapy2(x,y):
    """
    DOUBLE PRECISION FUNCTION DLAPY2( X, Y )
     
        DOUBLE PRECISION   X, Y
    
    Purpose:
    =============

        DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
        overflow.

    Arguments:
    ==========

    X
             X is DOUBLE PRECISION
    Y
             Y is DOUBLE PRECISION
             X and Y specify the values x and y.

    """

    x = np.float(x)     # a little error checking ... not doing this worked too for 3 or 3.0 coming in.
    y = np.float(y)

    # prep step to access the DLAPY2 function in cython_lapack
    #
    # DLAPY2( X, Y )

    ctypes.pythonapi.PyCapsule_GetPointer.restype = ctypes.c_void_p
    ctypes.pythonapi.PyCapsule_GetPointer.argtypes = [ctypes.py_object, ctypes.c_char_p]
    ctypes.pythonapi.PyCapsule_GetName.restype = ctypes.c_char_p
    ctypes.pythonapi.PyCapsule_GetName.argtypes = [ctypes.py_object]

    ptr_type = ctypes.CFUNCTYPE(ctypes.c_double,                    # this is val I changed to return a double from Fortran Function
                                ctypes.POINTER(ctypes.c_double),
                                ctypes.POINTER(ctypes.c_double))
    dlapy22 = ptr_type(
        ctypes.pythonapi.PyCapsule_GetPointer(
            cython_lapack.__pyx_capi__['dlapy2'],
            ctypes.pythonapi.PyCapsule_GetName(
                cython_lapack.__pyx_capi__['dlapy2'])))

    # FIXME - assumption here is that user passes in Fortran oriented arrays
    #         that are of the proper size. AND that they are Numpy arrays, so
    #         that they can be modified in place to return output values.

    x = ctypes.c_double(x)
    y = ctypes.c_double(y)

    args = (ctypes.byref(x),
            ctypes.byref(y)
           )

    r = dlapy22(*args)

    return r


#--------------------------------------------------------------------

def test():


    x = 3
    y = 4
    r0 = np.sqrt(x*x + y*y)
    
    print(" x,y before =", x, y)

    r = dlapy2(x,y)

    print("result dlayp2 = ", r)
    if r == 5:
        print("  - {0:3.1f} is the correct result".format(r0))
    else:
        print("  - Uh oh, result should be {0:3.1f} !!!!".format(r0))


if __name__ == "__main__":
    
    test()
