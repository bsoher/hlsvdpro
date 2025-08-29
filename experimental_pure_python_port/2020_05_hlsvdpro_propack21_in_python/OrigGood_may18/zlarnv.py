#
#  (C) Brian J Soher, 2020
#

import numpy as np

import dlaruv_cython as dlaruv

def zlarnv( IDIST, ISEED, N, X):
    """
     .. Scalar Arguments ..
     INTEGER            IDIST, N
     ..
     .. Array Arguments ..
     INTEGER            ISEED( 4 )
     COMPLEX*16         X( * )
     ..
    
    
    #  Purpose:
    *  =============
 
    #  ZLARNV returns a vector of n random complex numbers from a uniform or
    #  normal distribution.
    #  \endverbatim
 
    *  Arguments:
    *  ==========

    #  IDIST
    #           IDIST is INTEGER
    #           Specifies the distribution of the random numbers:
    #           = 1:  real and imaginary parts each uniform (0,1)
    #           = 2:  real and imaginary parts each uniform (-1,1)
    #           = 3:  real and imaginary parts each normal (0,1)
    #           = 4:  uniformly distributed on the disc abs(z) < 1
    #           = 5:  uniformly distributed on the circle abs(z) = 1
    #  ISEED
    #           ISEED is INTEGER array, dimension (4)
    #           On entry, the seed of the random number generator; the array
    #           elements must be between 0 and 4095, and ISEED(4) must be
    #           odd.
    #           On exit, the seed is updated.
    #  N
    #           N is INTEGER
    #           The number of random numbers to be generated.
    #  X
    #           X is COMPLEX*16 array, dimension (N)
    #           The generated random numbers.
    *
    #  Further Details:
    *  =====================
    # 
    #   This routine calls the auxiliary routine DLARUV to generate random
    #   real numbers from a uniform (0,1) distribution, in batches of up to
    #   128 using vectorisable code. The Box-Muller method is used to
    #   transform numbers from a uniform to a normal distribution.
    #  \endverbatim

    """
    TWOPI = 6.2831853071795864769252867663
    ZERO  = 0.0
    ONE   = 1.0
    TWO   = 2.0
    
    LV = 128 

    U = np.ndarray([LV,],  'd', order='F')
    

    for IV in range(0,N,int(LV/2)):
        
        IL = np.min( [int(LV/2), N-IV] )
    
        #  Call DLARUV to generate 2*IL real numbers from a uniform (0,1)
        #  distribution (2*IL <= LV)
    
        dlaruv.dlaruv(ISEED, 2*IL, U)
    
        if IDIST == 1:

            #  Copy generated numbers
    
            for I in range(IL):
                X[IV+I] = U[2*I] + (U[2*I+1] * 1j)
        
        elif IDIST == 2:
    
            # Convert generated numbers to uniform (-1,1) distribution
    
            for I in range(IL):
                X[IV+I] = TWO*U[2*I]-ONE + (( TWO*U[2*I+1]-ONE )*1j)
    
        elif IDIST == 3:
    
            #  Convert generated numbers to normal (0,1) distribution
    
            for I in range(IL):
                X[IV+I] = np.sqrt( -TWO*np.log(U[2*I]) ) * np.exp( 0.0+TWOPI*U[2*I+1]*1j )

        elif IDIST == 4:
    
            # Convert generated numbers to complex numbers uniformly
            # distributed on the unit disk
    
            for I in range(IL):
                X[IV+I] = np.sqrt( U[2*I] ) * np.exp( ZERO + TWOPI*U[2*I+1]*1j )

        elif IDIST == 5: 
    
            # Convert generated numbers to complex numbers uniformly
            # distributed on the unit circle
    
            for I in range(IL):
                X[IV+I] = np.exp( ZERO + TWOPI*U[2*I+1]*1j )
    
    return


# --------------------------------------------------------------------

def test():

    idist = 2
    n     = 250
    iseed = np.array([1, 3, 5, 7], 'i', order='F')
    x = np.ndarray([n, ], dtype = np.complex, order='F')

    print(" x before =", x)

    zlarnv(idist, iseed, n, x)

    print("random values")
    print(x)
    print(" ")
    print("Current seed values ")
    print(iseed)


if __name__ == "__main__":
    test()
