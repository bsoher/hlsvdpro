# Python modules
from __future__ import division
import pdb

# 3rd party modules
import numpy as np


# Our modules



def aprodw(transa, ndp, m, n, z, y, z1, y1, lambda_, trlambda):
    """
    Computes y=TH*T*z where lambda contains the eigenvalues of the circulant
    matrix used for T and trlambda contains the eigenvalues of the circulant
    matrix used for TH.

    parameter (zeroc=(0.0d0,0.0d0))
    integer         m,n,ndp
    complex*16      z(*),y(*),z1(ndp),y1(ndp),lambda(*),trlambda(*)
    character*1     transa
    
    """
    transa = transa.lower()
    if transa == 't':
        z1[:m] = z[:m]

        z1[m:ndp] *= 0j

        z1_copy, y1_copy = ztmultz(trlambda, ndp, z1, y1)

        z1[:ndp] = z1_copy[:ndp]
        y1[:ndp] = y1_copy[:ndp]
        y[:n] = y1[:n]

    else:
        z1[:n] = z[:n]

        z1[n:ndp] *= 0j

        z1_copy, y1_copy = ztmultz(lambda_, ndp, z1, y1)
        z1[:ndp] = z1_copy[:ndp]
        y1[:ndp] = y1_copy[:ndp]

        y[:m] = y1[:m]

    return z, y, z1, y1




def ztmultz(lambda_, ndp, z, y):
    """
    Multiply a Toeplitz matrix by a vector with the help of FFT

    integer         ndp
    double complex  lambda(*), z(*), y(*)

    """
    y = np.fft.fft(z)

    z = y * lambda_

    y = np.fft.ifft(z)

    # The FFTW doc says, "FFTW computes an unnormalized transform, in that 
    # there is no coefficient in front of the summation in the DFT. In other 
    # words, applying the forward and then the backward transform will multiply 
    # the input by n."
    # http://www.fftw.org/faq/section3.html#whyscaled
    # Apparently numpy's FFT is normalized, so in order to get the same 
    # results as FFTW computes, we have to multiply by n a.k.a. ndp.
    y *= ndp

    return z, y
