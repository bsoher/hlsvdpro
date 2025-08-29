from __future__ import division, print_function, absolute_import

import os
import numpy as np
import scipy as sp
import cProfile
from propack import svdp, svdp_aprod


def test(flavor, profiling=False):
    np.random.seed(0)
   
    m = 20      # was 50x10 orig
    n = 10
   
    # Create a random matrix 
    if flavor == 's':
        A = np.random.random((m, n)).astype(np.float32)
        flavorlabel = 'SINGLE'
    elif flavor =='d':
        A = np.random.random((m, n)).astype(np.float64)
        flavorlabel = 'DOUBLE'
    elif flavor =='c':
        A = np.random.random((m, n)) + 1j*np.random.random((m, n))
        A = A.astype(np.complex64)
        flavorlabel = 'COMPLEX8'
    elif flavor =='z':
        A = np.random.random((m, n)) + 1j*np.random.random((m, n))
        A = A.astype(np.complex128)
        flavorlabel = 'COMPLEX16'
    
    if not profiling:
        reps = 1
    else:
        reps = 1000
    
    k = 5

    for i in range(reps):
        # u, sigma, v = np.linalg.svd(A, full_matrices=False)
        u, sigma, v = sp.sparse.linalg.svds(A, k=k)
        sigma = sigma[::-1]     # sp returns singvals reverse of np or PROPACK - go figure
    
    # compute SVD via propack and lapack
    for i in range(reps):
        u1, sigma1, v1 = svdp_aprod(A, k)

    for i in range(reps):
        u2, sigma2, v2 = svdp_aprod(A, k, irl_mode=True)

    for i in range(reps):
        u3, sigma3, v3 = svdp(A, k)
 
    for i in range(reps):
        u4, sigma4, v4 = svdp(A, k, irl_mode=True)

    
    # print the results
    np.set_printoptions(suppress=True, precision=8)
    
    print( 'Flavor Label = '+flavorlabel)
    print( 'SingVals svd linalg = ', sigma[:k])
    print( 'SingVals svdp       = ', sigma1)
    print( 'SingVal Diffs (linalg vs svdp_aprod    ) = ', sigma[:k]-sigma1)
    print( 'SingVal Diffs (linalg vs svdp_irl_aprod) = ', sigma[:k]-sigma2)
    print( 'SingVal Diffs (linalg vs svdp          ) = ', sigma[:k]-sigma3)
    print( 'SingVal Diffs (linalg vs svdp_irl      ) = ', sigma[:k]-sigma4)
    print('')
    print( 'SingVal Diffs (svdp       vs svdp_aprod    )   = ', sigma3-sigma1)
    print( 'SingVal Diffs (svdp_irl   vs svdp_irl_aprod)   = ', sigma4-sigma2)
    print( 'SingVal Diffs (svdp       vs svdp_irl      )   = ', sigma3-sigma4)
    print( 'SingVal Diffs (svdp_aprod vs svdp_irl_aprod)   = ', sigma1-sigma2)
    print('')
    
    #print(np.dot(v, v1.T))


def test_profiling():
    for flavor in ['z']:  #['s','d','c','z']:
        test(flavor, profiling=True)


if __name__ == "__main__":

    for flavor in ['s','d','c','z']:
#    for flavor in ['c',]:
        print(' ')
        print('----------------------------------------------------------------')
        test(flavor, profiling=False)
        print(' ')
    
    print('----------------------------------------------------------------')
     
#    if os.path.exists("profile.data"):
#        os.remove("profile.data")
#    cProfile.run('test_profiling()', 'profile.data')
#    import pstats as ps
#    p = ps.Stats('profile.data')
#    p.strip_dirs().sort_stats('cumulative').print_stats()
