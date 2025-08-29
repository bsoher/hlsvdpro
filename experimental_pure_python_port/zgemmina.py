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

# Our modules
import util

def zgemmina(m, n, k, aaa, bbb, ldb):
    # subroutine zgemmina(m,n,k,A,lda,B,ldb,dwork,ldwork)
      # integer m,n,k,lda,ldb,ldwork 
      # double precision B(ldb,*)
      # double complex A(lda,*),dwork(ldwork)

    lda, oda = aaa.shape
    # pdb.set_trace()
    #lbd, odb = bbb.shape

    libhlsvd = util.load_the_library()

    aaa = np.transpose(aaa).flatten()
    AaaArrayType = ctypes.c_double * (len(aaa))
    aaa_r = AaaArrayType()
    aaa_i = AaaArrayType()
    for i, value in enumerate(aaa):
        aaa_r[i] = value.real
        aaa_i[i] = value.imag

    #bbb = np.transpose(bbb).flatten()
    bbb = bbb.flatten()
    BbbArrayType = ctypes.c_double * len(bbb)
    bbb_fortran = BbbArrayType()
    for i, value in enumerate(bbb):
        bbb_fortran[i] = value

    m = ctypes.c_long(m)
    n = ctypes.c_long(n)
    k = ctypes.c_long(k)
    lda = ctypes.c_long(lda)
    oda = ctypes.c_long(oda)
    ldb = ctypes.c_long(ldb)

    libhlsvd.zgemmina_python_(ctypes.pointer(m),
                     ctypes.pointer(n),
                     ctypes.pointer(k),
                     ctypes.pointer(aaa_r),
                     ctypes.pointer(aaa_i),
                     ctypes.pointer(lda),
                     ctypes.pointer(oda),
                     ctypes.pointer(bbb_fortran),
                     ctypes.pointer(ldb),
                    )

    lda = lda.value
    oda = oda.value

    aaa = [complex(r, i) for r, i in zip(aaa_r, aaa_i)]
    aaa = util.flat_fortran_to_2d_python(aaa, lda, oda)

    bbb = np.array([value for value in bbb])

    return aaa, bbb
