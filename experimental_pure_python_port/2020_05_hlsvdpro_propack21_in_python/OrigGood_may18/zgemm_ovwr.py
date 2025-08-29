#
#     (C) Brian J Soher, 2020
#

import numpy as np




def zdgemm_ovwr_left(transb,m,n,k,A,lda,B,ldb,zwork,lzwork):
    """compute  A <- A*op(B) """

    # NB. bjs NOT USING ANYMORE, REWROTE IN NUMPY IN ZRITZVEC.PY

    if m <= 0 or n <= 0 or k <= 0:
        return
        # raise ValueError('dgemm_ovwr: m<=0 or n<=0 or k<=0')

#     transa = 'N'                                # default for zdgemm_ovwr_left
#     
#     if transa in ['n','N']:
#         ka = k
#         if lda != m: raise ValueError("dgemm_ovwr: transa in ['n','N'] and lda!=m")
#     else:
#         ka = m
#         if lda != k: raise ValueError("dgemm_ovwr: transa not ['n','N'] and lda!=k")
#     
#     if transb in ['n','N']:
#         kb = n
#         if ldb != k: raise ValueError("dgemm_ovwr: transb in ['n','N'] and  ldb!=k")
#     else:
#         kb = k
#         if ldb != n: raise ValueError("dgemm_ovwr: transb not ['n','N'] and  ldb!=n")

    ccc = zdgemm(transb,m,n,k,A,lda,B,ldb,zwork,lzwork)

    # for j in range(n):
    #     for i in range(m):
    #         A[i,j] = zwork(j*m + i)

    A[:m,:n] = ccc

    return

#     implicit none
#     character*1 transb
#     integer m,n,k,lda,ldb,lzwork
#     complex*16 A(lda,*),zwork(lzwork)
#     double precision B(ldb,*)
#     integer i,j,l,blocksize
#     
#     if((m.le.0).or.(n.le.0).or.(k.le.0)) return
#     if (lzwork.lt.n) stop 'Too little workspace in ZDGEMM_OVWR_LEFT'
#     blocksize = int(lzwork/n)
#     i = 1
#     do i=1,m-blocksize+1,blocksize
#        call zdgemm(transb,blocksize,n,k, A(i,1),lda,
#     c              B,ldb,zwork,blocksize)
#        do j=0,n-1
#           do l=0,blocksize-1
#              A(i+l,j+1) = zwork(j*blocksize+1+l)
#           enddo
#        enddo
#     enddo
#     blocksize = m-i+1
#     call zdgemm(transb,blocksize,n,k,A(i,1),lda,
#     c           B,ldb,zwork,blocksize)
#     do j=0,n-1
#        do l=0,blocksize-1
#           A(i+l,j+1) = zwork(j*(m-i+1)+1+l)
#        enddo
#     enddo
#     return
#     end

def zdgemm(transb, m, n, k, A, lda, B, ldb, C, ldc):

    # NB. bjs  NOT USING ANYMORE, REWROTE IN NUMPY IN ZRITZVEC.PY


    # zero C, then C = A * B basically

    aaa = A.reshape(lda,len(A.flatten())/lda)
    bbb = B.reshape(ldb,len(B.flatten())/ldb)
    bbt = bbb.T
    ccc = np.zeros((m,n), dtype=np.complex128)
    for l in range(k):
        for j in range(n):
            for i in range(m):
                #C[i+j*m] += aaa[i, l] * bbt[j, l]
                ccc[i,j] += aaa[i, l] * bbt[j, l]


    # for i in range(m):
    #     for j in range(n):
    #         C[i+j*m] = 0.0+0.0j  # FIXME change to C *= 0.0+0.0j  ??
    #
    #
    # for l in range(k):
    #     for j in range(n):
    #         for i in range(m):
    #             C[i+j*m] += A[i,l] * B[j,l]

    return ccc


#     implicit none
#     character*1 transb
#     integer m,n,k,lda,ldb,ldc
#     complex*16 A(lda,*), C(ldc,*)
#     double precision B(ldb,*),btmp
#     integer i,j,l
#     
#     do i=1,m
#        do j=1,n
#           C(i,j) = dcmplx(0d0,0d0)
#        enddo
#     enddo
#     do l=1,k
#        do j=1,n
#           do i=1,m
#              C(i,j) = C(i,j) + A(i,l)*B(j,l)
#           enddo
#        enddo
#     enddo
#     end


def zdgemmblk(A,lda,B,ldb,C,ldc):

    # As of April 4, 2020, I do not see this being used in
    # the PROPACK v2.1 code, so I am stopping the port, other
    # than removing any broken code.

    # NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED

    blksz=96
    for l in range(blksz):
        for j in range(blksz):
            for i in range(blksz):
                C[i,j] = (A[i,l].real*B[j,l]+C[i,j].real) + (A[i,l].imag*B[j,l]+C[i,j].imag)*1j

    return C

#     implicit none
#     integer blksz
#     parameter (blksz=96)
#     integer lda,ldb,ldc
#     complex*16 A(lda,blksz), C(ldc,blksz)
#     double precision B(ldb,blksz)
#     integer i,j,l, i2,j2,l2
#     
#     do l=1,blksz
#        do j=1,blksz
#           do i=1,blksz
#              C(i,j) = dcmplx(dreal(A(i,l))*B(j,l)+dreal(C(i,j)),
#     c              dimag(A(i,l))*B(j,l)+dimag(C(i,j)))
#           enddo
#        enddo
#     enddo
#     end


def zdgemm1(transb,m,n,k,A,lda,B,ldb,C,ldc):
    """ compute C = A * OP(B) """

    # As of April 4, 2020, I do not see this being used in
    # the PROPACK v2.1 code, so I am stopping the port, other
    # than removing any broken code.

    # NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED NOT PORTED

    blksz=96
    
    BB = np.ndarray([blksz,blksz], dtype=np.double)
    CC = np.ndarray([blksz,blksz], dtype=np.complex)

    if transb in ['t', 'T']:
        # compute C = A*B^T
        
        for lblk in range(1,k-blksz+1,blksz):
            for jblk in range(1,n-blksz+1,blksz):
                for l in range(1,blksz):
                    for j in range(1,blksz):
                        BB[j,l] = B[jblk-1+j,lblk-1+l]

                for iblk in range(1,m-blksz+1,blksz):
                    if (lblk == 1):
                        for j in range(jblk,jblk+blksz-1):
                            for i in range(iblk,iblk+blksz-1):
                                C[i,j] = 0.0+0.0j
                    r = zdgemmblk(A[iblk,lblk],lda,BB,blksz, C[iblk,jblk],ldc)
                    C[iblk:iblk+blksz, jblk:jblk+blksz] = r                             # FIXME

                # clean up loops for i
                if (lblk == 1):
                    for j in range(jblk,jblk+blksz-1):
                        for i in range(iblk,m):
                            C[i,j] = 0.0+0.0j

                for l in range(lblk,lblk+blksz-1):
                    for j in range(jblk,jblk+blksz-1):
                        btmp = B[j,l]
                        for i in range(iblk,m):
                            C[i,j] =  A[i,l].real*btmp+C[i,j].real + A[i,l].imag*btmp+C[i,j].imag

            # clean up loops for j
            if (lblk == 1):
                for j in range(jblk,n):
                    for i in range(1,m):
                        C[i,j] = 0.0+0.0j

            for l in range(lblk,lblk+blksz-1):
                for j in range(jblk,n):
                    btmp = B[j,l]
                    for i in range(1,m):
                        C[i,j] =  A[i,l].real*btmp+C[i,j].real + A[i,l].imag*btmp+C[i,j].imag

        # clean up loop for l
        for l in range(lblk,k):
            if (l == 1):
                for j in range(1,n):
                    for i in range(1,m):
                        C[i,j] = 0.0+0.0j

            for jblk in range(1,n-blksz+1,blksz):
                for iblk in range(1,m-blksz+1,blksz):
                    for j in range(jblk,jblk+blksz-1):
                        btmp = B[j,l]
                        for i in range(iblk,iblk+blksz-1):
                            C[i,j] =  A[i,l].real*btmp+C[i,j].real + A[i,l].imag*btmp+C[i,j].imag
            
                for j in range(jblk,jblk+blksz-1):
                    btmp = B[j,l]
                    for i in range(iblk,m):
                        C[i,j] =  A[i,l].real*btmp+C[i,j].real + A[i,l].imag*btmp+C[i,j].imag
            
            for j in range(jblk,n):
                btmp = B[j,l]
                for i in range(1,m):
                    C[i,j] =  A[i,l].real*btmp+C[i,j].real + A[i,l].imag*btmp+C[i,j].imag
      
    else:
        # compute C = A*B
        
        for iblk in range(1,m-blksz+1,blksz):
            for jblk in range(1,n-blksz+1,blksz):
                for j in range(1,blksz):
                    for i in range(1,blksz):
                        CC[i,j] = 0.0+0.0j

                for j in range(1,blksz):
                    for l in range(1,k):
                        for i in range(1,blksz):
                            CC[i,j] = A[iblk-1+i,l]*B[l,jblk-1+j] + CC[i,j]

                for j in range(1,blksz):
                    for i in range(1,blksz):
                        C[iblk-1+i,jblk-1+j] = CC[i,j]

            for j in range(jblk,n):
                for i in range(iblk,iblk+blksz-1):
                    C[i,j] = 0.0+0.0j

            # clean up loop for j
            for j in range(jblk,n):
                for l in range(1,k):
                    for i in range(iblk,iblk+blksz-1):
                        C[i,j] = A[i,l]*B[l,j] + C[i,j]
        for j in range(1,n):
            for i in range(iblk,m):
                C[i,j] = 0.0+0.0j
             
        # clean up loop for i
        for j in range(1,n):
            for l in range(1,k):
                for i in range(iblk,m):
                    C[i,j] = A[i,l]*B[l,j] + C[i,j]

    return 



      

