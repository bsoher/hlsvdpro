#
#  (C) Brian J Soher, 2020
#

import numpy as np
import scipy.linalg.blas as blas

#from zmgs import zmgs
from zblasext import pdznrm2, pzzero




GAMMA = 0.98
NTRY = 5


def zreorth(n, k, V, ldv, vnew, normvnew, index, alpha, iflag):
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
    
    """
    if (k >= 0) and (n > 0):

        skip_zeroing = False

        for itry in range(NTRY):
            normvnew_0 = normvnew

            if iflag == 1:
                zcgs(n, k, V, ldv, vnew, index)
            else:
                zmgs(n, k, V, ldv, vnew, index)

            normvnew = pdznrm2(n, vnew, 1)

            if normvnew > (alpha * normvnew_0):
                skip_zeroing = True
                break

        if not skip_zeroing:
            # vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
            normvnew = 0.0
            vnew[:n] *= 0j      # replace pzzero()
    
    return normvnew


#------------------------------------------------------------------------------

def zcgs(n, k, V, ldv, vnew, index):
    """
    Block  Gram-Schmidt orthogonalization:
    FOR i= 1:l
        vnew = vnew - V(:,[s_i:e_i])*(V(:,[s_i:e_i])'*vnew)
     
    If l=1 and s_1=1 and e_1=k then this becomes classical Gram-Schmidt.

    """
    work = np.zeros([n,], np.complex128)
    yloc = np.zeros([n,], np.complex128)

    i = 0

    # PS - python/fortran difference: it's OK for index[iblck] to be 0
    # since Python uses 0-based indices
    while index[i] <= k and index[i] >= 0:

        p = index[i]        # Select the next block of columns from V
        q = index[i+1]
        l = q-p+1

        # Classical Gram-Schmidt: vnew = vnew - V(:,p:q)*(V(:,p:q)'*vnew)
        if l >= 0:

            # call zgemv('C', n, l, 1+0j, V(1,p), ldv, vnew(1),1, 0+0j, ylocal, 1)
            yloc = blas.zgemv((1+0j), V[:n,p:p+l], vnew, (0+0j), yloc, 0, 1, 0, 1, 2, 1)

            work[:l] = yloc[:l]

            # call zgemv('N',n,l,dcmplx(-1d0,0d0),V(1,p),ldv,work,1, dcmplx(0d0,0d0),ylocal,1)
            yloc = blas.zgemv((-1-0j), V[:n,p:p+l], work, (0+0j), yloc, 0, 1, 0, 1, 0, 0)

            vnew[:] += yloc

            print( 'i,p,q,l =',i,p,q,l)

        i = i+2



#------------------------------------------------------------------------------

def zmgs(n, k, V, ldv, vnew, index):

    # NB. This code taken from v1 pure python zreorth.zreorth2() where MGS was
    #     hard-coded, thus CGS was never called.

    # PS - slight difference from the fortran: k == 0 is OK since Python uses
    # 0-based indices. However, as in the Fortran, n must be > 0 since it is
    # a counter.
    if (k >= 0) and (n > 0):

        # PS - python/fortran difference: it's OK for index[iblck] to be 0
        # since Python uses 0-based indices

        iblck = 0

        while (index[iblck] <= k) and (index[iblck] >= 0) and (index[iblck] <= index[iblck + 1]):

            p = index[iblck]        # Select the next block of columns from V
            q = index[iblck + 1]

            for i in range(p, q+1):
                s = 0+0j

                s = (V[:n, i].conjugate() * vnew[:n]).sum()
                #s = np.sum(V[:n, i].conjugate() * vnew[:n])    # slower bjs

                vnew[:n] -= (s * V[:n, i])

                print('p,q =', p, q)

            iblck += 2

    return



