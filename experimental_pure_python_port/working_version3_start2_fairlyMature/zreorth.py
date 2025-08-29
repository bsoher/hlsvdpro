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

KMAX = 50
MAX_SINGULAR_VALUES = KMAX
NDPMAX = 8192


GAMMA = 0.98
NTRY = 4


# Complex BLAS functions
zgemv, zdotc = scipy.linalg.blas.get_blas_funcs( ['gemv', 'dotc'], 
                                                 np.array([0j]) 
                                               )

# Float BLAS functions
dznrm2, = scipy.linalg.blas.get_blas_funcs( ['znrm2'], 
                                            np.array([0.0]) 
                                          )
# Variable name translations 
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# iflag      classical          > 0 or True ==> use classical Gram-Schmidt,
#                               otherwise use modified Gram-Schmidt. 



#       subroutine zreorth(n,k,V,ldv,vnew,normvnew,index,alpha,work,
#      c     iflag)
def zreorth(n, k, vvv, vnew, normvnew, index, alpha, classical):


# c     Orthogonalize the N-vector VNEW against a subset of the columns of
# c     the N-by-K matrix V(1:N,1:K) using iterated classical or modified
# c     Gram-Schmidt. LDV is the leading dimension of the array containing
# c     V.
# c     
# c     Which columns to orthogonalize against is decided by the integer
# c     array INDEX = [s_1,e_1, s_2,e_2,..., s_k,e_l, s_{l+1}], which
# c     selects the columns V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]). s_{l+1}
# c     must be larger than k and marks the end of INDEX.
# c
# c     The reorthogonalization is repeated until
# c
# c       ||VNEW'|| > ALPHA * ||VNEW|| , 
# c
# c     where VNEW' is the vector obtained by orthogonalizing VNEW.  If
# c     VNEW' fails to satisfy this after 4 tries, VNEW is deemed to lie
# c     numerically in the span of V(:,[s_1:e_1 s_2:e_2 ... s_l:e_l]), and
# c     is set to the zero vector.
# c
# c     On return NORMVNEW contains ||VNEW||.
# c
# c     WORK is a workspace array of length at least 
# c
# c       max e_i-s_i+1, i=1...l.
# c     
# c     WORK is only used if IFLAG==1.
# c
# c     If IFLAG==0 then iterated modified Gram-Schmidt is used.
# c     If IFLAG==1 then iterated classical Gram-Schmidt is used.
# c

# c     References: 
# c       Aake Bjorck, "Numerical Methods for Least Squares Problems",
# c       SIAM, Philadelphia, 1996, pp. 68-69.
# c     
# c       J.~W. Daniel, W.~B. Gragg, L. Kaufman and G.~W. Stewart, 
# c       ``Reorthogonalization and Stable Algorithms Updating the
# c       Gram-Schmidt QR Factorization'', Math. Comp.,  30 (1976), no.
# c       136, pp. 772-795.
# c
# c       B. N. Parlett, ``The Symmetric Eigenvalue Problem'', 
# c       Prentice-Hall, Englewood Cliffs, NJ, 1980. pp. 105-109

# c     Rasmus Munk Larsen, Stanford, 1999.

# c     %-----------%
# c     | Arguments |
# c     %-----------%
#       implicit none
#       include 'stat.h'
#       integer n,k,ldv,iflag,index(*)
#       double complex V(ldv,*),vnew(*),work(*)
#       double precision normvnew

# c     %------------%
# c     | Parameters |
# c     %------------%
#       integer NTRY
#       double precision one, zero
#       parameter(one = 1.0d0, zero = 0.0d0, NTRY=4)
      
# c     %-----------------%
# c     | Local variables |
# c     %-----------------%
#       integer i,itry
#       double precision alpha,normvnew_0
#       real t2,t3
      
# c     %----------------------%
# c     | External Subroutines |
# c     %----------------------%
#       external zgemv
      
# c     %--------------------%
# c     | External Functions |
# c     %--------------------%
#       double precision dznrm2
#       external dznrm2      

#       if (k.le.0 .or. n.le.0) return

    ldv, _ = vvv.shape

    if (k > 0) and (n > 0):

#       call second(t2)

        skip_zeroing = False
#       do itry=1,NTRY
        for itry in range(NTRY):

#           normvnew_0 = normvnew         
#           if (iflag.eq.1 ) then
#              call zCGS(n,k,V,ldv,vnew,index,work)
#           else  
#              call zMGS2(n,k,V,ldv,vnew,index)
#           endif

            normvnew_0 = normvnew
            if classical:
                zgcs(n, k, vvv, ldv, vnew, index)
            else:
                zmgs2(n, k, vvv, ldv, vnew, index)

# c         ndot = ndot + k
#           normvnew = dznrm2(n,vnew,1)

            #normvnew = dznrm2(n, vnew, 1)
            normvnew = scipy.linalg.norm(vnew)
            # print "normvnew = {}".format(normvnew)
#           if (normvnew.gt.alpha*normvnew_0) goto 9999
            if normvnew > (alpha * normvnew_0):
                skip_zeroing = True
                break
#       enddo
#  8888 normvnew = zero
        if not skip_zeroing:
            normvnew = 0j
    # c     vnew is numerically in span(V) => return vnew = (0,0,...,0)^T
    #       do i=1,n
    #          vnew(i) = zero
    #       enddo
            vnew[:n] *= 0j


    return vnew, normvnew

#  9999 call second(t3)
#       treorth = treorth + (t3-t2)
#       nreorth = nreorth + 1
#       return
#       end



# c
# c****************************************************************************
# c

#       subroutine zCGS(n,k,V,ldv,vnew,index,work)

def zgcs(n, k, vvv, ldv, vnew, index):
# c     Block  Gram-Schmidt orthogonalization:
# c     FOR i= 1:l
# c         vnew = vnew - V(:,[s_i:e_i])*(V(:,[s_i:e_i])'*vnew)
# c      
# c     If l=1 and s_1=1 and e_1=k then this becomes classical Gram-Schmidt.

 
# c     %-----------%
# c     | Arguments |
# c     %-----------%
#       implicit none
#       include 'stat.h'
#       integer n,k,ldv,index(*)
#       double complex V(ldv,*),vnew(*),work(*)
# c     %------------%
# c     | Parameters |
# c     %------------%
#       double complex one, zero
#       parameter(one = (1.0d0,0.0d0), zero = (0.0d0,0.0d0))
#       integer i,p,q,l

#       i=1
#       do while(index(i).le.k .and. index(i).gt.0)
    i = 1
    while (index[i] < k) and (index[i] > 0):
# c
# c     Select the next block of columns from V
# c         
#          p = index(i)
#          q = index(i+1)
#          l = q-p+1
        p = index[i]
        q = index[i + 1]
        l = q - p + 1

# c         ndot = ndot + l

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

#          i = i+2
        i += 2
#       enddo
#       end      



# c
# c****************************************************************************
# c


#       subroutine zMGS2(n,k,V,ldv,vnew,index)
def zmgs2(n, k, vvv, ldv, vnew, index):

    # PS - This function is never called. It has not been tested.


# c      subroutine MGS2(n,k,V,ldv,vnew)
#       implicit none
#       include 'stat.h'
#       integer n,k,ldv,index(*)
#       double complex V(ldv,*),vnew(*)
#       integer i,j,p,q,iblck
#       double complex s,zdotc
#       external zdotc

# c     Check for quick return
#       if ((k.le.0).or.(n.le.0)) return
    if (k >= 0) and (n >= 0):
#       iblck = 1
        iblck = 1

#       do while(index(iblck).le.k.and.index(iblck).gt.0)
        while (index[iblck] <= k) and (index[iblck] > 0):
#          p = index(iblck)
#          q = index(iblck+1)
#          ndot = ndot + (q-p+1)
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
                    vnew[j] = vnew[j] - (s * V[j, i])

#           iblck = iblck + 2

            iblck += 2

#       enddo
#       end

# c
# c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# c

#       subroutine zreorth2(n,k,V,ldv,vnew,nrm,index)




def zreorth2(n, k, vvv, vnew, nrm, index):

#       implicit none
#       integer n,k,ldv,index(*)
#       complex*16 V(ldv,*),vnew(*),h,s
#       double precision  nrm

# c     Local variables
#       include 'stat.h'
#       integer i,j,p,q,iblck
#       double precision gamma,nrm0,thr
#       parameter(gamma = 0.98d0)
#       real t2,t3

# c     External subroutines
#       double precision dznrm2
#       external dznrm2


# c     
# c     Modified Gram-Schmidt orthogonalization:
# c     Orthogalizes vnew against the k vectors in V by the
# c     iterative process     
# c     
# c     FOR i=1...k DO          
# c       vnew = vnew - DOT( V(:,i), vnew ) * V(:,i) 
# c

# c     This simple version is faster on Pentium machines.
# c     Compile with "g77 -O6 -funroll-all-loops -fomit-frame-pointer"
     

# c     Check for quick return
#       if (k.le.0 .or. n.le.0) return
#       call second(t2)

    #pdb.set_trace()

    # PS - slight difference from the fortran: k == 0 is OK since Python uses
    # 0-based indices. However, as in the Fortran, n must be > 0 since it is
    # a counter. 
    if (k >= 0) and (n > 0):
#       nrm0 = nrm*nrm   
#       thr = gamma*nrm0
#       iblck = 1
        nrm0 = nrm * nrm
        thr = GAMMA * nrm0
        iblck = 0

#       do while(index(iblck).le.k.and.index(iblck).gt.0)

        # PS - python/fortran difference: it's OK for index[iblck] to be 0
        # since Python uses 0-based indices
        while (index[iblck] <= k) and (index[iblck] >= 0):
#           p = index(iblck)
#           q = index(iblck+1)
#           ndot = ndot + (q-p+1)
            p = index[iblck]
            q = index[iblck + 1]
            # print "zreorth2: p = {}, q = {}".format(p, q)
            # ndot unused

#           do i=p,q
            # FIXME PS this might should be q + 1?
            for i in range(p, q+1):
#               s = (0.0d0,0.0d0)
#               do j=1,n
#                  s = s + conjg( V(j,i))*vnew(j)
#               enddo
                
                # PS - here's the slow way of doing the above
                # s = 0j
                # for j in range(n):
                #     s += (vvv[j,i].conjugate() * vnew[j])

                # PS - Here's the fast way.
                s = (vvv[:n,i].conjugate() * vnew[:n]).sum()

#               h = s
                h = s

                # print "zreorth2: i = {}, h = {}".format(i, h)

#               do j=1,n
#                  vnew(j) = vnew(j) - s*V(j,i)
#               enddo
                # for j in range(n):
                #     vnew[j] -= (s * vvv[j, i])
                vnew[:n] -= (s * vvv[:n, i])


   
#               if  ((dble(conjg(s)*s)).gt.thr)  then
                if (s.conjugate() * s).real > thr:
                    # PS - With the test data I have, I wasn't able to 
                    # trigger this case so I haven't tested this code.

#                   ndot = ndot+1
#                   s = (0.0d0,0.0d0)
#                   do j=1,n
#                       s = s + conjg(V(j,i))*vnew(j)
#                   enddo
                    # PS - here's the slow way of doing the above
                    # s = 0j
                    # for j in range(n):
                    #     s += vvv[j,i].conjugate() * vnew[j]

                    s = (vvv[:n,i].conjugate() * vnew[:n]).sum()

#                   h = s + h
                    h += s


#                   do j=1,n
#                      vnew(j) = vnew(j) - s*V(j,i)
#                   enddo
                    # for j in range(n):
                    #     vnew[j] -= (s * vvv[j,i])
                    vnew[:n] -= (s * vvv[:n, i])

#               endif
#               nrm0 = nrm0 - dble(conjg(h)*h)
                nrm0 -= (h.conjugate() * h).real
#               thr = nrm0*gamma
                thr = nrm0 * GAMMA
#           enddo

#           iblck = iblck + 2
            iblck += 2
#       enddo

#       nrm = dznrm2(n,vnew,1)
        nrm = scipy.linalg.norm(vnew)
#       call second(t3)
#       treorth = treorth + (t3-t2)
#       nreorth = nreorth + 1



#       end
# c
# c****************************************************************************
# c
