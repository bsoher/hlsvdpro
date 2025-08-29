#
#  (C) Brian J Soher, 2020
#

import numpy as np
import scipy.linalg.blas as blas


def dgemm_ovwr(transa,m,n,k,alpha,A,lda,beta,B,ldb, dwork,ldwork):
    """ compute B <- alpha*op(A)*B + beta*B """

    if m <= 0 or n <= 0 or k <= 0:
        print('dgemm_ovwr: returning, m<=0 or n<=0 or k<=0')
        return

    if k>m or k>n:
        print('dgemm_ovwr: returning, k>m or k>n')
        return

    # FIXME bjs, hardcoded for HLSVDPro for the moment

    boba = A.reshape(m,m)
    bobb = B.reshape(n,n)

    res = blas.dgemm(alpha, boba.T[:k, :m], bobb.T[:k, :n], beta=beta, trans_a=1, trans_b=0)
    bobb[:n,:m] = res.T


    return 


    # if transa in ['n','N']:
    #     ka = k
    #     if lda != m: raise ValueError("dgemm_ovwr: transa in ['n','N'] and lda!=m")
    # else:
    #     ka = m
    #     if lda != k: raise ValueError("dgemm_ovwr: transa not ['n','N'] and lda!=k")
    #
    # if transb in ['n','N']:
    #     kb = n
    #     if ldb != k: raise ValueError("dgemm_ovwr: transb in ['n','N'] and  ldb!=k")
    # else:
    #     kb = k
    #     if ldb != n: raise ValueError("dgemm_ovwr: transb not ['n','N'] and  ldb!=n")


#       implicit none
#       character*1 transa
#       integer m,n,k,lda,ldb,ldwork
#       double precision alpha,beta,A(lda,*),B(ldb,*),dwork(ldwork)
#       integer i,j,l,blocksize
# 
#       if((m.le.0).or.(n.le.0).or.(k.le.0)) return
#       if (ldwork.lt.m) stop 'Too little workspace in DGEMM_OVWR'
#       if (m.gt.ldb) stop 'm>ldb in DGEMM_OVWR'
#       blocksize = int(ldwork/m)
#       do i=1,n-blocksize+1,blocksize
#          call dgemm(transa,'N',m,blocksize,k,alpha,A,lda, B(1,i),ldb,0D0,dwork,m)
#          if (beta.eq.0D0) then
#             do j=0,blocksize-1
#                do l=1,m
#                   B(l,i+j)  = dwork(j*m+l)
#                enddo
#             enddo
#          else
#             do j=0,blocksize-1
#                do l=1,m
#                   B(l,i+j)  = dwork(j*m+l) + beta*B(l,i+j)
#                enddo
#             enddo
#          endif
#       enddo
#       call dgemm(transa,'N',m,n-i+1,k,alpha,A,lda,
#      c           B(1,i),ldb,0D0,dwork,m)
#       if (beta.eq.0D0) then
#          do j=0,n-i
#             do l=1,m
#                B(l,i+j)  = dwork(j*m+l)
#             enddo
#          enddo
#       else
#          do j=0,n-i
#             do l=1,m
#                B(l,i+j)  = dwork(j*m+l) + beta*B(l,i+j)
#             enddo            
#          enddo
#       endif      
#       return
#       end


def dgemm_ovwr_left(transb,m,n,k,alpha,A,lda,beta,B,ldb, dwork,ldwork):
    """  compute  A <- alpha*A*op(B) """

    # As of April 4, 2020, I do not see this being used in
    # the PROPACK v2.1 code, so I am stopping the port, other
    # than removing any broken code.

    # NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED

    if m <= 0 or n <= 0 or k <= 0:
        raise ValueError('dgemm_ovwr_left: m<=0 or n<=0 or k<=0')

    if k>m or k>n:
        raise ValueError('dgemm_ovwr_left: k>m or k>n')

    # A figure out ka and kb depending on transa
    
#    c = np.zeros([m,n]).astype(np.double)
    itransa = 0 #'N'                                # default for dgemm_ovwr_left
    itransb = 0 if transb in ['n', 'N'] else 1

    if transa in ['n','N']:
        ka = k
        if lda != m: raise ValueError("dgemm_ovwr_left: transa in ['n','N'] and lda!=m")
    else:
        ka = m
        if lda != k: raise ValueError("dgemm_ovwr_left: transa not ['n','N'] and lda!=k")
    
    if transb in ['n','N']:
        kb = n
        if ldb != k: raise ValueError("dgemm_ovwr_left: transb in ['n','N'] and  ldb!=k")
    else:
        kb = k
        if ldb != n: raise ValueError("dgemm_ovwr_left: transb not ['n','N'] and  ldb!=n")
        
    # A[:m,:n] = blas.dgemm(alpha, A[:lda,:ka], B[:ldb,:kb], beta=np.float(0.0), trans_a=itransa, trans_b=itransb)      # FIXME bjs, See dgemm_ovwr() above!!
    A[:m*n] = blas.dgemm(alpha, A[:m*n], B[:m*n], beta=np.float(0.0), trans_a=itransa, trans_b=itransb)
    
    return 

#       implicit none
#       character*1 transb
#       integer m,n,k,lda,ldb,ldwork
#       double precision alpha,beta,A(lda,*),B(ldb,*),dwork(ldwork)
#       integer i,j,l,blocksize
# 
#       if((m.le.0).or.(n.le.0).or.(k.le.0)) return
#       if (ldwork.lt.n) stop 'Too little workspace in DGEMM_OVWR_LEFT'
#       blocksize = int(ldwork/n)
#       do i=1,m-blocksize+1,blocksize
#          call dgemm('n',transb,blocksize,n,k,alpha,A(i,1),lda, B,ldb,0d0,dwork,blocksize)
#          do j=0,n-1
#             do l=0,blocksize-1
#                A(i+l,j+1) = dwork(j*blocksize+1+l)
#             enddo
#          enddo
#       enddo
#       call dgemm('n',transb,m-i+1,n,k,alpha,A(i,1),lda,
#      c           B,ldb,0d0,dwork,m-i+1)
#       do j=0,n-1
#          do l=0,m-i
#             A(i+l,j+1) = dwork(j*(m-i+1)+1+l)
#          enddo
#       enddo
#       return
#       end


