# Python modules
from __future__ import division
import ctypes
import math
import os
import pprint
import pdb
pp = pprint.pprint

# 3rd party modules
import numpy as np
import scipy
import scipy.linalg.blas

# Our modules

KMAX    = 50
NDPMAX  = 8192
GAMMA   = 0.98
NTRY    = 4
MAX_SINGULAR_VALUES = KMAX


# Complex BLAS functions
zgemv, zdotc = scipy.linalg.blas.get_blas_funcs( ['gemv', 'dotc'], np.array([0j]) )

# Float BLAS functions
dznrm2, = scipy.linalg.blas.get_blas_funcs( ['znrm2'], np.array([0.0]) )

# Variable name translations 
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# iflag      classical          > 0 or True ==> use classical Gram-Schmidt,
#                               otherwise use modified Gram-Schmidt. 



def zreorth(n, k, vvv, vnew, normvnew, index, alpha, classical):
    """

    Orthogonalize the N-vector VNEW against a subset of the columns of
    the N-by-K matrix V(1:N,1:K) using iterated classical or modified
    Gram-Schmidt. LDV is the leading dimension of the array containing
    V.
    
    Which columns to orthogonalize against is decided by the integer
    array INDEX = [s_1,e_1, s_2,e_2,..., s_k,e_l, s_{l+1}], which
    selects the columns V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]). s_{l+1}
    must be larger than k and marks the end of INDEX.

    The reorthogonalization is repeated until

      ||VNEW'|| > ALPHA * ||VNEW|| , 

    where VNEW' is the vector obtained by orthogonalizing VNEW.  If
    VNEW' fails to satisfy this after 4 tries, VNEW is deemed to lie
    numerically in the span of V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]), and
    is set to the zero vector.

    On return NORMVNEW contains ||VNEW||.

    WORK is a workspace array of length at least 

      max e_i-s_i+1, i=1...l.
    
    WORK is only used if IFLAG==1.

    If IFLAG==0 then iterated modified Gram-Schmidt is used.
    If IFLAG==1 then iterated classical Gram-Schmidt is used.


    References: 
      Aake Bjorck, "Numerical Methods for Least Squares Problems",
      SIAM, Philadelphia, 1996, pp. 68-69.
    
      J.~W. Daniel, W.~B. Gragg, L. Kaufman and G.~W. Stewart, 
      ``Reorthogonalization and Stable Algorithms Updating the
      Gram-Schmidt QR Factorization'', Math. Comp.,  30 (1976), no.
      136, pp. 772-795.

      B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
      Prentice-Hall, Englewood Cliffs, NJ, 1980. pp. 105-109

    Rasmus Munk Larsen, Stanford, 1999.
    
    """
    if (k > 0) and (n > 0):

        skip_zeroing = False

        for itry in range(NTRY):

            normvnew_0 = normvnew
            if classical:
                zcgs(n, k, vvv, ldv, vnew, index)
            else:
                zmgs2(n, k, vvv, ldv, vnew, index)

            #normvnew = dznrm2(n, vnew, 1)
            normvnew = scipy.linalg.norm(vnew)
            # BJS replace with scipy.linalg.blas.dznrm2

            if normvnew > (alpha * normvnew_0):
                skip_zeroing = True
                break

        if not skip_zeroing:
            normvnew = 0j
            vnew[:n] *= 0j

    return vnew, normvnew



#****************************************************************************
#       subroutine zCGS(n,k,V,ldv,vnew,index,work)

def zcgs(n, k, vvv, ldv, vnew, index):
    """
    Block  Gram-Schmidt orthogonalization:

    """
    i = 1
    while (index[i] < k) and (index[i] > 0):
        # Select the next block of columns from V
        p = index[i]
        q = index[i + 1]
        l = q - p + 1

        # According to the documentation for zgemv(), work must be  --
        #    ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
        #    and at least
        #    ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
        shape = ( max(m, n), 1)
        work = np.zeros(shape, np.complex128)

#          call zgemv('C',n,l,one,V(1,p),ldv,vnew(1),1,zero,work(1),1)
        zgemv('C', n, l, (1+0j), vvv[0, :p], ldv, vnew, 1, 0j, work, 1)

#          call zgemv('N',n,l,-one,V(1,p),ldv,work(1),1,one,vnew(1),1)
        zgemv('N', n, l, -(1+0j), vvv[0, :p], ldv, work, 1, (1+0j), vnew, 1)

        i += 2



#****************************************************************************
#       subroutine zMGS2(n,k,V,ldv,vnew,index)

def zmgs2(n, k, vvv, ldv, vnew, index):

    # PS - This function is never called. It has not been tested.

    # Check for quick return
    if (k >= 0) and (n >= 0):
        iblck = 1

        while (index[iblck] <= k) and (index[iblck] > 0):
            p = index[iblck]
            q = index[iblck + 1]
            # ndot = ndot + (q - p + 1)

#          do i=p,q
#             s = (0.0d0,0.0d0)
#             do j=1,n
#                s = s + conjg(V(j,i))*vnew(j)
#             enddo
# Cc             s = - zdotc(n,V(1,i),1,vnew(1),1)

#             do j=1,n
#                vnew(j) = vnew(j) - s*V(j,i)
#             enddo
# Cc             call zaxpy(n,s,V(1,i),1,vnew(1),1)     
#          enddo
            
            # FIXME check for off by one error
            for i in range(p, q):
                s = 0j
                for j in range(n):
                    s += vvv[j, i].conjugate() * vnew[j]

                for j in range(n):
                    # orig PS code changed by BJS vnew[j] = vnew[j] - (s * V[j, i])
                    vnew[j] = vnew[j] - (s * vvv[j, i])

            iblck += 2


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#       subroutine zreorth2(n,k,V,ldv,vnew,nrm,index)

def zreorth2(n, k, vvv, vnew, nrm, index):
    """
    Modified Gram-Schmidt orthogonalization:
    Orthogalizes vnew against the k vectors in V by the
    iterative process     
    
    """
     
    # PS - slight difference from the fortran: k == 0 is OK since Python uses
    # 0-based indices. However, as in the Fortran, n must be > 0 since it is
    # a counter. 
    if (k >= 0) and (n > 0):
        nrm0 = nrm * nrm
        thr = GAMMA * nrm0
        iblck = 0

#        # part of the speedup attempt below - worked OK, but did not get faster
#        vvvc = np.conjugate(vvv)

        # PS - python/fortran difference: it's OK for index[iblck] to be 0
        # since Python uses 0-based indices
        while (index[iblck] <= k) and (index[iblck] >= 0):
            p = index[iblck]
            q = index[iblck + 1]
            # ndot unused

            # FIXME PS this might should be q + 1?
            for i in range(p, q+1):
                
                # PS - here's the slow way of doing the above
                # s = 0j
                # for j in range(n):
                #     s += (vvv[j,i].conjugate() * vnew[j])
 
                # PS - Here's the fast way.
                s = (vvv[:n,i].conjugate() * vnew[:n]).sum()
                h = s

#                 # Attempted to speed up this line by moving .conjugate() 
#                 # it worked OK outside loop, but it got slower somehow!
#                 s = (vvvc[:n,i] * vnew[:n]).sum()
#                 h = s
                
                vnew[:n] -= (s * vvv[:n, i])

                if (s.conjugate() * s).real > thr:
                    # PS - With the test data I have, I wasn't able to 
                    # trigger this case so I haven't tested this code.
#                    s = (vvv[:n,i].conjugate() * vnew[:n]).sum()
                    
                    tmp = np.conjugate(vvv[:n,i])
                    s = (tmp * vnew[:n]).sum()
                    
                    h += s
                    vnew[:n] -= (s * vvv[:n, i])
 
                nrm0 -= (h.conjugate() * h).real
                thr = nrm0 * GAMMA

            iblck += 2

        nrm = dznrm2(vnew)
        #nrm  = scipy.linalg.norm(vnew)
        



    

