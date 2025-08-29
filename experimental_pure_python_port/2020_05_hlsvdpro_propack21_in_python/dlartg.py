# Python modules
from __future__ import division
import math


# 3rd party modules
import numpy as np
import scipy.linalg.lapack

# Our modules

# Float LAPACK functions
dlamch, = scipy.linalg.lapack.get_lapack_funcs( ['lamch'], 
                                                np.array([0.0])
                                              )
sqrt = math.sqrt


# Calculate some constants
SAFMIN = dlamch('S')
EPS = np.finfo(np.float64).eps
SAFMN2 = int(math.log(SAFMIN / EPS) / math.log(dlamch('B')) / 2.0)
SAFMN2 = dlamch('B') ** SAFMN2
SAFMX2 = 1.0 / SAFMN2



def dlartg(fff, ggg):
    """
    returns cs, sn, r

    # SUBROUTINE DLARTG( F, G, CS, SN, R )
    # *
    # DOUBLE PRECISION   CS, F, G, R, SN
    # 
    # Purpose:
    # =============
    # 
    # DLARTG generate a plane rotation so that
    # 
    #    [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
    #    [ -SN  CS  ]     [ G ]     [ 0 ]
    # 
    # This is a slower, more accurate version of the BLAS1 routine DROTG,
    # with the following other differences:
    #    F and G are unchanged on return.
    #    If G=0, then CS=1 and SN=0.
    #    If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
    #       floating point operations (saves work in DBDSQR when
    #       there are zeros on the diagonal).
    # 
    # If F exceeds G in magnitude, CS will be positive.
    # 
    # Arguments:
    # ==========
    # F
    #     F is DOUBLE PRECISION
    #     The first component of vector to be rotated.
    # G
    #     G is DOUBLE PRECISION
    #     The second component of vector to be rotated.
    # CS
    #     CS is DOUBLE PRECISION
    #     The cosine of the rotation.
    # SN
    #     SN is DOUBLE PRECISION
    #     The sine of the rotation.
    # R
    #     R is DOUBLE PRECISION
    #     The nonzero component of the rotated vector.
    # 
    #  This version has a few statements commented out for thread safety
    #  (machine parameters are computed on each entry). 10 feb 03, SJH.

    """

    # Check for pathologically trivial cases first
    if ggg == 0.0:
        cs = 1.0
        sn = 0.0
        r = fff
    elif fff == 0.0:
        cs = 0.0
        sn = 1.0
        r = ggg
    else:

        # fff and ggg are non-zero; this case requires proper math.

        f1 = fff
        g1 = ggg
        scale = max(abs(f1), abs(g1))
        if scale >= SAFMX2:
            count = 0
            while scale >= SAFMX2:
                count += 1
                f1 *= SAFMN2
                g1 *= SAFMN2
                scale = max(abs(f1), abs(g1))

            r = sqrt(f1*f1+g1+g1)
            cs = f1 / r
            sn = g1 / r

            for i in range(count):
                r *= SAFMX2

        elif scale <= SAFMN2:
            count = 0
            scale = max(abs(f1), abs(g1))
            while scale <= SAFMN2:
                count += 1
                f1 *= SAFMX2
                g1 *= SAFMX2
                scale = max(abs(f1), abs(g1))

            r = sqrt(f1*f1+g1+g1)
            cs = f1 / r
            sn = g1 / r

            for i in range(count):
                r *= SAFMN2
        else:
            r = sqrt(f1*f1 + g1*g1)
            cs = f1 / r
            sn = g1 / r

        if (abs(fff) > abs(ggg)) and (cs < 0.0):
            cs = -cs
            sn = -sn
            r = -r

    return cs, sn, r
