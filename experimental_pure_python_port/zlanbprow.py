# Python modules
from __future__ import division
import ctypes
import math
import pdb

# 3rd party modules
import numpy as np
import numpy.linalg
import scipy
import scipy.linalg.blas
import scipy.linalg.lapack

# Our modules
import vespa.common.util.misc as util_misc
import zget0w
import dlartg
import dbdsqr
import util
import aprodw
import zreorth
import zsafescal


MAX_SINGULAR_VALUES = 50
FUDGE = 1.01
KAPPA = 0.717

# Machine- and compiler-dependent constants
# On my Mac, EPS in Python is 2.22044604925e-16. In Fortran
# is is half as big (1.11022302462515654E-016).
EPS = np.finfo(np.float64).eps
EPS34 = EPS ** (3.0 / 4.0)
EPS1 = 100 * EPS


# In the Fortran, the constant MGS (==> modified Gram-Schmidt) is hardcoded 
# to 1. It's not present in the original PROPACK code. Since it's always 1,
# I've opted to ignore it.

# I create a shortcut for sqrt() since I use it a lot.
sqrt = math.sqrt

# http://www.scipy.org/doc/api_docs/SciPy.lib.lapack.info.html

# Float BLAS functions
dznrm2, = scipy.linalg.blas.get_blas_funcs( ['znrm2'], 
                                            np.array([0.0]) 
                                          )

# Complex BLAS functions
zdotc, = scipy.linalg.blas.get_blas_funcs( ['dotc'], 
                                            np.array([0j]) 
                                          )


# Float LAPACK functions
dlamch, = scipy.linalg.lapack.get_lapack_funcs( ['lamch'], 
                                                np.array([np.float64(0.0)])
                                              )




def fortran_range(a, b):
    the_range = range(a, b)
    if b >= a:
        the_range.append(b)
    return the_range


# cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# c........1.........2.........3.........4.........5.........6.........7..
# c23456789012345678901234567890123456789012345678901234567890123456789012
# c
#       subroutine zlanbprow_python(ndp, m, n, k0, k, ierr,
#      c                            uuu_r, uuu_i, ldu,
#      c                            vvv_r, vvv_i, ldv,
#      c                            bbb, ldb,
#      c                            rnorm, doption, ioption, 
#      c                            lambda_r, lambda_i,
#      c                            trlambda_r, trlambda_i)


def zlanbprow(ndp, m, n, k0, k, uuu, vvv, bbb_a, bbb_b, rnorm,
              options, lambda_, trlambda):
    # FIXME complete the docstring
    """This docstring is incomplete.

    bbb_a and bbb_b are modified in-place.

    This function accepts a dict of options that exert some control over how
    it behaves. The dict keys and their meaning are as follows --
       - classical_gs: True ==> classical Gram-Schmidt reorth (see zreorth)
       - extended_reorth_iterations: # of iterations in one case of 
             extended local reorth
       - delta: used if >= 0, otherwise a reasonable default is calculated.
       - eta:   used if >= 0, otherwise a reasonable default is calculated.
       - anorm: used if >  0, otherwise a reasonable default is calculated.

    None of these options are changed by this function EXCEPT FOR "anorm"
    which is rewritten every time (presumably so it can be passed back in 
    on the next iteration).
    """

# C     params
#       integer ndp, m, n, k, k0, kmax, ierr, ldu, ldv, ldb
#       parameter (kmax=50)
#       integer ndpmax
#       parameter (ndpmax=8192)
# C       integer lwork1, lworkp, lworkq, 
#       integer lworktotal
#       integer info
#       real*8 lambda_r(ndp), lambda_i(ndp)
#       real*8 trlambda_r(ndp), trlambda_i(ndp)
#       real*8 uuu_r(ldu, kmax + 1), uuu_i(ldu, kmax + 1)
#       real*8 vvv_r(ldv, kmax + 1), vvv_i(ldv, kmax + 1)
#       real*8 rnorm
#       real*8 bbb(ldb, 2)        ! dims = (lanmax, 2)
# C       real*8 bbb1(ldb, 2)       ! dims = (lanmax, 2)
# C       real*8 work1(lwork1), workp(lworkp), workq(lworkq)


# c     Local variables.
#       integer ips, jps, kps, work_index
#       complex*16 uuu(ldu, kmax + 1)
#       complex*16 vvv(ldv, kmax)
#       complex*16 lambda(ndp)
#       complex*16 trlambda(ndp)

#       integer ioption(2)
#       real*8 doption(3)
#       real*8 a_work

#       integer, allocatable :: iwork(:)
#       real*8, allocatable :: work(:)
#       complex*16, allocatable :: zwork(:)


#       allocate(iwork(2*kmax+1))
            

# C      i = m+n+13*kmax+8*kmax**2+32*m+ndp*kmax
# C      lworktotal = lwork1 + (ldb * 2) + (ldb * 2) + lworkp + lworkq
#       lworktotal = 295115
#       allocate(work(lworktotal))


#       subroutine zlanbprow( ndp, m, n, k0, k,  U, ldu, V, ldv, B, ldb,
#      c     rnorm, doption, ioption, work, iwork,
#      c      ierr,zwork,lambda,trlambda,planF,planB)

# c     %-----------%
# c     | Arguments |
# c     %-----------%
#       implicit none
#       include 'stat.h'
#       integer m, n, k0, k, ldb, ldu, ldv, ierr
#       integer ips, jps
#       integer ioption(*), iwork(*),ndp,PLANF,PLANB
#       double precision rnorm,B(ldb,*), doption(*), work(*)
#       complex*16  U(ldu,*),V(ldv,*), zwork(*), czero,s,cone,
#      c            lambda(*),trlambda(*)
#       external aprodw
#       integer kmax
#       parameter (kmax=50)


# c     %------------%
# c     | Parameters |
# c     %------------%
#       logical DEBUG
#       integer MGS
#       double precision one, zero, FUDGE, kappa
#       parameter(one = 1.0d0, zero = 0.0d0, FUDGE = 1.01)
# c      parameter (DEBUG=.TRUE., MGS=0, kappa = 0.717)
#       parameter (DEBUG=.FALSE., MGS=1, kappa = 0.717)
#       parameter (czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0))
# c     %-----------------%
# c     | Local variables |
# c     %-----------------%
#       integer i,j,iv,iu,inu,imu,is,iidx,j0,iy1,iz1
#       double precision eps,eps34,epsn2,epsn,sfmin,delta,eta,anorm
#       double precision mumax,numax,alpha,beta,a1,b1,amax,anormest
#       logical force_reorth,full_reorth,debito
#       real t1,t2,t3

# c-------------------- Here begins executable code ---------------------
#       call second(t1)

    # print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ zlanbprow() "

# c     %---------------------------------%
# c     | Set machine dependent constants |
# c     %---------------------------------%
#       eps = dlamch('e')
#       eps34 = eps**(3d0/4d0)
#       epsn = dble(max(m,n))*eps
#       epsn2 = sqrt(dble(max(m,n)))*eps
#       sfmin = dlamch('s')
    # PS - Fortran populates variable sfmin but otherwise doesn't use it.
    
    EPSN = max(m, n) * EPS
    EPSN2 = sqrt(max(m, n)) * EPS

    ierr = 0

# c      kappa =1.0d0/dsqrt(2.0d0)
#       debito = .false.

# c     %------------------------%
# c     | Set default parameters |
# c     %------------------------%
#       if (doption(1).lt.zero) then
#          delta = sqrt(eps/k)
#       else
#          delta = doption(1)
#       endif
    delta = sqrt(EPS / k) if (options["delta"] < 0) else options["delta"]

#       if (doption(2).lt.zero) then
#          eta = eps34/sqrt(dble(k))
#       else
#          eta = doption(2)
#       endif

    eta = EPS34 / sqrt(k) if (options["eta"] < 0) else options["eta"]


#       if (delta.le.eta .or. delta.eq.zero) then
#          full_reorth = .true.
#       else
#          full_reorth = .false.
#       endif

    full_reorth = ((delta <= eta) or (delta == 0.0))

#       if (doption(3).gt.zero) then
#          anorm = doption(3)
#       else if (k0.gt.0) then
#          anorm = dlapy2(B(1,1),B(1,2))
#          if (anorm.le.zero) then
#             ierr = -1
#             goto 9999         
#          endif
#       else
#          anorm = zero
#       endif
# c               delta =3.107106895105165e-09
# c              eta = 3.792855096563922e-13

    if options["anorm"] > 0:
        anorm = options["anorm"]
    elif k0 > 0:
        anorm = util.dlapy2(bbb_a[0], bbb_b[0])
        if anorm <= 0:
            ierr = -1
            # FIXME raise an error instead?
            return
    else:
        anorm = 0.0

# c     %---------------------%
# c     | Get starting vector |
# c     %---------------------%
#       if (rnorm .eq. zero) then
# c         call zgetu0('n',m, n, k0, 3, U(1,k0+1) , rnorm, U,ldu, aprod,
# c     c         ierr, ioption(1), anormest,zwork,lambda,trlambda)
# c         anorm = max(anorm,anormest)
#          print *, "rnorm is zero; calling zgetu0w()"
#          call zgetu0w('n',ndp,m, n, k0, 3, U(1,k0+1) , rnorm, U,ldu, 
#      c         ierr, ioption(1), anormest,zwork,lambda,
#      c           trlambda,planF,planB) ! it was zwork instead of zwork(is) ! diana

# ccc         call zgetu0('n',m,n,0,1,U,rnorm,U,ldu,
# ccc     c        ierr,ioption(1),anorm,zwork(iwrk),lambda,
# ccc     c        trlambda,wsave)

# c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#          anorm = max(anorm,anormest)
#       endif

    # print "ndp, m, n = ", ndp, m, n
    # print "k0, k = ", k0, k
    # print "U(1,k0+1) = ", uuu[0, k0]
    # print "anorm = ", anorm
    # print "rnorm = ", rnorm
    # print "doption: ", doption
    # print "ioption: ", ioption

    if not rnorm:
        rnorm, anorm, ierr = zget0w.zgetu0w('n', ndp, m, n, k0, 3, uuu[:,k0:],
                                               uuu, options["classical_gs"],
                                               lambda_, trlambda)
        # print "done zgetu0w()"
        # print "rnorm, anorm = {}".format(rnorm, anorm)

# c %------------------------------%
# c | Set pointers into work array |
# c %------------------------------%
#   iu = 1
#   iv = iu+m
#   imu = iv+n
# c  iv = iu+ndp
# c  imu = iv+ndp
#   inu = imu+k+1
#   is = inu+k+1
#   iz1 = is+3*ndp
#   iy1 = iz1+ndp
#   iidx = 1
#   call dzero(m+n+max(m,n)+2*k+2,work,1)
#   call zzero(m+n+max(m,n)+2*k+2,zwork,1)
#   call izero(k,iwork,1)

    zwork  = np.zeros( (m, ), np.complex128)  # replaces zwork[iu] a.k.a. zwork
    zworkv = np.zeros( (n, ), np.complex128)  # replaces zwork[iv]
    zworks = np.zeros( (3 * ndp, ), np.complex128)  # replaces zwork[is]
    zworkz1 = np.zeros( (ndp, ), np.complex128)  # replaces zwork[iz1]
    # I'm not 100% sure about the size of zworky1 since it is the last pointer
    # in the Fortran zwork array. Since it doesn't have a next neighbor in 
    # Fortran, I can't tell where iy1 stops. Based on how it's used, I think 
    # a length of ndp is OK.
    zworky1 = np.zeros( (ndp, ), np.complex128)  # replaces zwork[iy1]

    # mus and nus track the mu and nu values. They replace the work array in
    # Fortran which appears to be much bigger than it needs to be.
    mus = np.zeros( (k + 1), np.float64)
    nus = np.zeros( (k + 1), np.float64)

    # indices replaces iwork from the Fortran code. It might not need to be this
    # big.
    indices = np.zeros( ((2 * MAX_SINGULAR_VALUES) + 1, ), np.int)
    
    indices -= 1


    # The Fortran code does something fishy in the first calls to zgetu0w()
    # below. It passes a complex
    # variable 's' in the slot that zgetu0w() calls u0norm. u0norm
    # is declared as a real in zgetu0w() despite being passed a complex
    # by zlanbprow. Whether or not this is valid Fortran I don't know.
    # In any case, the practical result is that the real part of s
    # is altered while the imaginary part remains unchanged. 
    # In Fortran, s is initially populated with stack garbage. I "declare"
    # it here to simulate Fortran's behavior.
    s = 0j

# c %---------------------------%
# c | Prepare Lanczos iteration |
# c %---------------------------%
# c  write(*,*)k0,'******************************** k0'

#   if (k0.eq.0) then
#      print *, "zlanbprow: k0 is zero"
#      amax = zero
#      alpha = zero
#      beta = rnorm
#      force_reorth = .false.
    if not k0:
        # print "zlanbprow: k0 is zero"
        amax = 0.0
        alpha = 0.0
        beta = rnorm
        force_reorth = False

# c     %-----------------------------------------------%
# c     | Compute ||A x|| / ||x|| for a random vector x |
# c     | to make it less likely that ||A|| is grossly  |
# c     | underestimated.                               |
# c     %-----------------------------------------------%
#       if (n.gt.m) then      
#           print *, "zlanbprow: n > m"
# c         call zgetu0('n',m,n,0,1,zwork(iu),s,U,ldu,aprod,
# c  c           ierr,ioption(1),anormest,zwork(is),lambda,
# c  c           trlambda,wsave)
#          call zgetu0w('n',ndp,m,n,0,1,zwork(iu),s,U,ldu,
#   c           ierr,ioption(1),anormest,zwork(is),lambda,
#   c           trlambda,planF,planB)

        if n > m:
            # print "zlanbprow: n > m"
            s_real, anormest, ierr = zget0w.zgetu0w('n', ndp, m, n, 0, 1, zwork,
                                               uuu, options["classical_gs"],
                                               lambda_, trlambda)
            s.real = s_real
#       else  
#           print *, "zlanbprow: n <= m"

# C            do i = 1, n
# C              do j = 1, kmax
# C                PRINT *, "V(", i, j, "):", V(i,j)
# C             end do
# C            end do

# c         call zgetu0('y',m,n,0,1,zwork(iv),s,V,ldv,aprod,
# c  c           ierr,ioption(1),anormest,zwork(is),lambda,
# c  c           trlambda,wsave)
#          call zgetu0w('t',ndp,m,n,0,1,zwork(iv),s,V,ldv,
#   c           ierr,ioption(1),anormest,zwork(is),lambda,
#   c           trlambda,planF,planB)

# C          do i=iv, iv + m
# C              print *, "zwork(", i, ") = ", zwork(i)
# C          end do
# C          print *, "s = ", s
#       endif

        else:
            # print "zlanbprow: n <= m"
            s_real, anormest, ierr = zget0w.zgetu0w('t', ndp, m, n, 0, 1, 
                                                zworkv,
                                                vvv, options["classical_gs"],
                                                lambda_, trlambda)
            s = complex(s_real, s.imag)
            # print "after zgetu0w 531"
            # print "s = {}, anormest = {}".format(s, anormest)

# C     print *, "ndp, m, n = ", ndp, m, n
# C     print *, "k0, k = ", k0, k
# C     print *, "U(1,k0+1) = ", U(1,k0+1)
# C     print *, "ldu, ldv, ldb = ", ldu, ldv, ldb
# C     print *, "anorm = ", anorm
# C     print *, "rnorm = ", rnorm
# C     print *, "planF, planB = ", planF, planB
# C     print *, "doption: ", doption(1), doption(2), doption(3)
# C     print *, "ioption: ", ioption(1), ioption(2)

# C     print *, "banana"

#       anorm = max(anorm,FUDGE*anormest)
        anorm = max(anorm, FUDGE * anormest)
        # print "max anorm = {}".format(anorm)
#       j0 = 1
        j0 = 0
#       if (beta.ne.zero) then
#          call donothing()
#          call zsafescal(m,beta,U(1,k0+1))
#       endif
        if beta:
            # print "beta = {}; calling zsafescal w/k0={}".format(beta, 0)
            zsafescal.zsafescal(m, beta, uuu[:, 0])


#       call zcopy(m,U(1,k0+1),1,zwork(iu),1)

        zwork = uuu[:, 0].copy()


#       work(imu) = one
#       work(inu) = one
        mus[0] = 1.0
        nus[0] = 1.0

        # I suspect that the two lines below from the Fortran code are a 
        # harmless bug. In every other case where they're used, imu and inu 
        # are indices into work, not zwork. I'm just going to ignore these
        # two lines.
#       zwork(imu) = cone
#       zwork(inu) = cone
        

#   else
    else:
        # k0 != 0

#       alpha = B(k0,1)
#       call zcopy(n,V(1,k0),1,zwork(iv),1)
#       beta = rnorm
#       force_reorth = .true.
    
        # for ips, value in enumerate(bbb_a):
        #     print "bbb_a[{}] = {}".format(ips + 1, value)
        # print "boom, k0 = {}".format(k0)

        alpha = bbb_a[k0 - 1]
        zworkv = vvv[:, k0].copy()
        # for ips, value in enumerate(zworkv):
        #     print "zworkv[{}] = {}".format(ips + 1, value)
        beta = rnorm
        force_reorth = True

#       if (k0.lt.k .and. beta*delta.lt.anorm*eps) then
#          full_reorth = .true.
#          ierr = k0
#       endif
#       iwork(iidx) = 1
#       iwork(iidx+1) = k0
#       iwork(iidx+2) = k0+1

        if (k0 < k) and (beta * delta < anorm * EPS):
            full_reorth = True
            ierr = k0

        indices[0] = 0
        indices[1] = k0 - 1 
        indices[2] = k0 

#       call second(t2)

        # PS - In the Fortran, MGS is hardcoded to 1 and the MGS == 0 case is 
        # commented out. I didn't port the dead code.
# c     if (MGS.eq.0) then
# c        call zreorth(m,k0,U,ldu,U(1,k0+1),s,iwork(iidx),kappa,
# c c           zwork(is),ioption(1))
# c     else
# c        call zreorth2(m,k0,U,ldu,U(1,k0+1),s,iwork(iidx))
#         call zreorth2(m,k0,U,ldu,U(1,k0+1),rnorm,iwork(iidx))
# c     endif
#       call second(t3)
#       treorthu = treorthu+(t3-t2)

        # print "before zreorth2 624"
        zreorth.zreorth2(m, k0, uuu, uuu[:, k0], rnorm, indices)

#       call dset_mu(k0,work(imu),iwork(iidx),epsn2)
        dset_mu(k0, mus, indices, EPSN2)


# c     %--------------------------------------%
# c     | Estimate ||B||_2^2 as ||B^T * B||_1  |   
# c     %--------------------------------------%

# c
# c     beta = s*beta
# c

# c*******************************
# c         beta = dble(s)*beta
# c*******************************
#       B(k0,2) = beta
#       amax = zero

        bbb_b[k0 - 1] = beta
        amax = 0.0

#       do j=1,k0
#          amax = max(amax,B(j,1),B(j,2))
#          if (j.eq.1) then
#             anorm = max(anorm,FUDGE*alpha)
#          else if (j.eq.2) then
#             a1 = B(1,2)/amax
#             a1 = FUDGE*amax*sqrt((B(1,1)/amax)**2 + a1**2 +
#   c              B(2,1)/amax*a1)
#             anorm = max(anorm,a1)
#          else
#             a1 = B(j-1,1)/amax
#             b1 = B(j-1,2)/amax 
#             a1 = FUDGE*amax*sqrt( a1**2 + b1**2 +
#   c              a1*B(j-2,2)/amax + B(j,1)/amax*b1)
#             anorm = max(anorm,a1)
#          endif
#          work(imu+j-1) = epsn2
#          work(inu+j-1) = epsn2
#       enddo

        for j in range(k0):
            amax = max(amax, bbb_a[j], bbb_b[j])
            if j == 0:
                anorm = max(anorm, FUDGE * alpha)
            elif j == 1:
                a1 = bbb_b[0] / amax
                a1 = FUDGE * amax * sqrt((bbb_a[0] / amax) ** 2 + \
                                         a1 ** 2 +                 \
                                         bbb_a[1] / amax * a1)
            else:
                a1 = bbb_a[j - 1] / amax
                b1 = bbb_b[j - 1] / amax
                a1 = FUDGE * amax * sqrt(a1**2 +                        \
                                         b1**2 +                        \
                                         a1 * bbb_b[j - 2] / amax +    \
                                         bbb_a[j] / amax * b1
                                        )
                anorm = max(anorm, a1)
            mus[j - 1] = EPSN2
            nus[j - 1] = EPSN2

#       j0 = k0+1
#       call zcopy(m,U(1,k0+1),1,zwork(iu),1)
        j0 = k0
        zwork = uuu[:, k0].copy()
#   endif
#   numax = zero
#   mumax = zero

    numax = 0.0
    mumax = 0.0


    
# c  %-------------------------------------------%
# c  | Start Lanczos bidiagonalization iteration |
# c  %-------------------------------------------%            

    # print "Starting Lanczos bidiagonalization iteration, j0 = {}".format(j0)
    # print "alpha = ", alpha

#   do j=j0,k

    # FIXME check for off-by-one error
    try:
        range(j0, k)
    except:
        pdb.set_trace()


    for j in range(j0, k):
        # print "\nLanczos iteration, j = " + ((str(j + 1) + ' ') * 15)
#       if (j.eq.1) then
#          call aprodw('t',ndp,m,n,zwork(iu),zwork(iv),zwork(iz1),
#   c                  zwork(iy1),lambda,trlambda,
#   c                  planF,planB )

#          alpha = dznrm2(n,zwork(iv),1)
#          anorm = max(anorm,FUDGE*alpha)
        if j == 0:
            # print "before aprodw 711"
            # for ips, value in enumerate(zwork[:10]):
            #     print "zwork[{}] = {}".format(ips, value)
            # print "-"
            # for ips, value in enumerate(zworkv[:10]):
            #     print "zworkv[{}] = {}".format(ips, value)
            # print "-"
            # for ips, value in enumerate(zworkz1[:10]):
            #     print "zworkz1[{}] = {}".format(ips, value)
            # print "-"
            # for ips, value in enumerate(zworky1[:10]):
            #     print "zworky1[{}] = {}".format(ips, value)
            # print "-"
            aprodw.aprodw('t', ndp, m, n, zwork, zworkv, 
                          zworkz1, zworky1, lambda_, trlambda)
            # f = open("1p.txt", 'w')
            # for value in zworkv:
            #     f.write("%.17f\t%.17f\n" % (value.real, value.imag))
            # f.close()
            # kjhafsdkgjadfh

            # print "after aprodw 711"
            # for ips, value in enumerate(zwork[:10]):
            #     print "zwork[{}] = {}".format(ips, value)
            # print "-"
            # for ips, value in enumerate(zworkv[:10]):
            #     print "zworkv[{}] = {}".format(ips, value)
            # print "-"
            # for ips, value in enumerate(zworkz1[:10]):
            #     print "zworkz1[{}] = {}".format(ips, value)
            # print "-"
            # for ips, value in enumerate(zworky1[:10]):
            #     print "zworky1[{}] = {}".format(ips, value)
            # print "-"

            #pdb.set_trace()
            #alpha = dznrm2(n, zworkv, 1)

            alpha = scipy.linalg.norm(zworkv)
            anorm = max(anorm, FUDGE * alpha)
            # print "alpha (e norm) = {}, anorm = {}".format(alpha, anorm)
#       else
        else:
            
# c        %---------------------------------------------%
# c        | alpha_{j} v_{j} = A'*u_{j} - beta_{j} v_{j} |
# c        %---------------------------------------------%
# c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#          call aprodw('t',ndp,m,n,zwork(iu),zwork(is),zwork(iz1),
#      c               zwork(iy1),lambda,trlambda,planF,planB )

# c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            # print "before aprodw 748"
            # for ips, value in enumerate(zwork[:10]):
            #     print "zwork[{}] = {}".format(ips, value)
            # for ips, value in enumerate(zworks[:10]):
            #     print "zworks[{}] = {}".format(ips, value)
            # for ips, value in enumerate(zworkz1[:10]):
            #     print "zworkz1[{}] = {}".format(ips, value)
            # for ips, value in enumerate(zworky1[:10]):
            #     print "zworky1[{}] = {}".format(ips, value)
            aprodw.aprodw('t', ndp, m, n, zwork, zworks, 
                          zworkz1, zworky1, lambda_, trlambda)
            # print "after aprodw 748"
            # for ips, value in enumerate(zworks[:10]):
            #     print "zworks[{}] = {}".format(ips, value)

#             do i=0,n-1
#                zwork(iv+i) = zwork(is+i) - beta*zwork(iv+i)
#             enddo

            # PS - simple (but slow) Python equivalent -- 
            # for i in range(n):
            #     zworkv[i] = zworks[i] - (beta * zworkv[i])

            # Fast Python equivalent -- 
            zworkv[:n] = zworks[:n] - (beta * zworkv[:n])            


# c         %------------------------------------%
# c         | Extended local reorthogonalization |
# c         %------------------------------------%
#           alpha = dznrm2(n,zwork(iv),1)
            alpha = scipy.linalg.norm(zworkv)


#           if (j.gt.1 .and. ioption(2).gt.0 .and.
#      c         alpha.lt.kappa*beta) then
#              do i=1,2
#                 s = zdotc(n,V(1,j-1),1,zwork(iv),1)

#                 call zaxpy(n,-s,V(1,j-1),1,zwork(iv),1)
#                 if (beta .ne. zero) then
#                    beta = beta + dble(s)
#                    B(j-1,2) = beta
#                 endif
#                 s = dznrm2(n,zwork(iv),1)
# c                if (dsqrt(dble(conjg(s)*s)) .ge. kappa*alpha) goto 10 
#                  if (dble(s) .ge. kappa*alpha) goto 10                
#              enddo
#  10          work(inu+j-2) = eps
# c             alpha =dsqrt( dble(conjg(s)*s))
#              alpha =dble(s)
#           endif

            # FIXME I think 'j > 0' is always True -- look at the surrounding
            # if statement.
            # if j == 2:
            #     pdb.set_trace()
            if (j > 0) and (options["extended_reorth_iterations"] > 0) and \
               (alpha < (KAPPA * beta)):
                for i in (0, 1):
                    # for ips, value in enumerate(vvv[:, j - 1]):
                    #     print "647 V[{}, {}] = {}".format(ips + 1, j, value)
                    # for ips, value in enumerate(zworkv):
                    #     print "647 zworkv[{}] = {}".format(ips + 1, value)

                    s = np.vdot(vvv[:, j - 1], zworkv)
                    # print "vdot 647, s = {}".format(s)

#                   call zaxpy(n,-s,V(1,j-1),1,zwork(iv),1)

                    # This is a slow inline Python equivalent of the call 
                    # to zaxpy --
                    # for k in range(n):
                    #     zworkv[k] += (-s * vvv[k, j - 1])

                    # This is a fast version --
                    zworkv[:n] += (-s * vvv[:n, j - 1])


                    if beta != 0.0:
                        beta += s
                        bbb_b[j - 1] = beta.real
                    s = scipy.linalg.norm(zworkv)

                    if s >= KAPPA * alpha:
                        # Bail out of for loop
                        break
                # print "setting nus[{}] = {}".format(j - 1, eps)
                nus[j - 1] = EPS
                alpha = s

#           B(j,1) = alpha
#           amax = max(amax,alpha)
            bbb_a[j] = alpha
            amax = max(amax, alpha)

# c     %---------------------------%
# c     | Update estimate of ||A||_2 |         
# c     %---------------------------%
#           if (j.eq.2) then
#              a1 = B(1,2)/amax
#              a1 = FUDGE*amax*sqrt((B(1,1)/amax)**2 + a1**2 +
#   c               B(2,1)/amax*a1)
#           else
#              a1 = B(j-1,1)/amax
#              b1 = B(j-1,2)/amax 
#              a1 = FUDGE*amax*sqrt( a1**2 + b1**2 +
#   c               a1*B(j-2,2)/amax + B(j,1)/amax*b1)
#           endif
            if j == 1:
                a1 = bbb_b[0] / amax
                a1 = FUDGE * amax * sqrt( (bbb_a[0] / amax) ** 2 +     \
                                          (a1 ** 2) +                   \
                                          bbb_a[1] / amax * a1
                                        )
            else:
                a1 = bbb_a[j - 1] / amax
                b1 = bbb_b[j - 1] / amax
                a1 = FUDGE * amax * sqrt(a1**2 +                        \
                                         b1**2 +                        \
                                         a1 * bbb_b[j - 2] / amax +    \
                                         bbb_a[j] / amax * b1
                                        )
#           anorm = max(anorm,a1)
            anorm = max(anorm, a1)
#       endif
        # end if j == 1/else
        # still in 'for j in range(j0, k):'

# c     %--------------------------%
# c     | Update the nu recurrence |
# c     %--------------------------%
#       if (.not.full_reorth .and. alpha.ne.zero) then
#          call dupdate_nu(numax,work(imu),work(inu),j,
#      c                   B(1,1),B(1,2),anorm)
#       endif


        if (not full_reorth) and (alpha != 0.0):
            # print "before dupdate_nu()"
            # for ips, value in enumerate(nus[:10]):
            #     print "nus[{}] = {}".format(ips, value)
            # for ips, value in enumerate(mus[:10]):
            #     print "mus[{}] = {}".format(ips, value)
            # if j == 2:
            #     pdb.set_trace()
            numax = dupdate_nu(mus, nus, j, bbb_a, bbb_b, anorm)
            # print "after dupdate_nu(), numax = {}".format(numax)
            # for ips, value in enumerate(nus[:10]):
            #     print "nus[{}] = {}".format(ips, value)

# c     %------------------------------%
# c     | Reorthogonalize if necessary |
# c     %------------------------------%
#       if ( .true. .and.
#   c        (full_reorth .or. numax.gt.delta .or. force_reorth) 
#   c        .and. alpha.ne.zero) then
        # if j == 6:
        #     pdb.set_trace()
        if (full_reorth or force_reorth or (numax > delta)) and (alpha != 0.0):
            # print "Reorthogonalizing..."
#           if (full_reorth .or. eta.eq.zero) then
#              iwork(iidx) = 1
#              iwork(iidx+1) = j-1
#              iwork(iidx+2) = j
#           else if (.not. force_reorth) then
#              call dcompute_int(work(inu),j-1,delta,eta,iwork(iidx))
#           endif
            if full_reorth or (eta == 0.0):
                indices[0] = 0
                indices[1] = j - 1
                indices[2] = j 
                indices[3] = -1
                indices[4] = -1
            elif not force_reorth:
                #pdb.set_trace()
                # if j == 1:
                #     pdb.set_trace()                
                # print "before dcompute_int 981"
                indices = dcompute_int(nus, j - 1, delta, eta)

#           if (MGS.eq.0) then
#              call zreorth(n,j-1,V,ldv,zwork(iv),alpha,iwork(iidx),
#   c               kappa,zwork(is),ioption(1))
#           else
#              call zreorth2(n,j-1,V,ldv,zwork(iv),alpha,iwork(iidx))
#           endif
            # Note: MGS is hardcoded to 1
            # if MGS == 0:
            #     zreorth.zreorth(n, j - 1, vvv, ldv, zworkv, alpha, 
            #                     indices, KAPPA, zworks, 
            #                     options["classical_gs"])
            # else:

            # print "before zreorth2 990"
            zreorth.zreorth2(n, j - 1, vvv, zworkv, alpha, indices)

#           call dset_mu(j-1,work(inu),iwork(iidx),epsn2)
            # if j == 1:
            #     pdb.set_trace()
            # print "before dset_mu 883"
            # for ips, value in enumerate(nus[:10]):
            #     print "nus[{}] = {}".format(ips, value)
            # print "indices = " + str(indices)
            dset_mu(j - 1, nus, indices, EPSN2)
            # print "after dset_mu 883"
            # for ips, value in enumerate(nus[:10]):
            #     print "nus[{}] = {}".format(ips, value)

#           numax = eta
            numax = eta

#           if (force_reorth) then
#              force_reorth = .false.
#           else
#              force_reorth = .true.
#           endif
            force_reorth = not force_reorth
#       endif
        # end if (full_reorth or force_reorth or (numax > delta)) and (alpha != 0.0)


# c     %-----------------------------------------------%
# c     | Check whether an invariant subspace was found |
# c     %-----------------------------------------------% 
        # print "Checking for invariant subspace..."
#       if (alpha .lt. anorm*epsn .and. j.lt.k) then
        if (alpha < (anorm * EPSN)) and (j < k):
            # print "alpha = {}, anorm = {}, anorm * epsn = {}".format(alpha, anorm, anorm * epsn)
#           rnorm = alpha
#           alpha = zero
            rnorm = alpha
            alpha = 0.0
# c     %------------------------------------------------%
# c     | Try to build and orthogonal subspace, starting |
# c     | with a random vector.                          |
# c     %------------------------------------------------%
#           call zgetu0w('t',ndp, m, n, j-1, 3, zwork(iv), alpha, V, 
#      c                 ldv, ierr,ioption(1),anormest,
#      c                 zwork(is),lambda,trlambda,planF,planB)
            alpha, anormest, ierr = \
                zget0w.zgetu0w('t', ndp, m, n, j - 1, 3, zworkv, vvv, 
                               options["classical_gs"], lambda_, trlambda)


#           if (alpha .eq. zero) then
#              k = j-1
#              ierr = -j
#              goto 9999
#           else
#              nrstrt = nrstrt + 1
#              call donothing()
#              call zsafescal(n,alpha,zwork(iv))
#              alpha = zero
#              force_reorth = .true.
#              if (delta.gt.zero) then
#                 full_reorth = .false.
#              endif
#           endif
            if alpha == 0.0:
                k = j - 1
                ierr = -j
                # FIXME this should merely bail out of the loop, not raise an 
                # error.
                raise ValueError
            else:
                zsafescal.zsafescal(n, alpha, zworkv)
                alpha = 0.0
                force_reorth = True
                if delta > 0.0:
                    full_reorth = False

#        else if (j.gt.1 .and. .not. full_reorth .and. j.lt.k .and.
#      c         (delta*alpha .lt. anorm*eps)) then
# c          full_reorth = .true.
# c          full_reorth = .true.
#           ierr = j
#        endif            
        elif (j > 1) and (not full_reorth) and (j < k) and \
             (delta * alpha < anorm * EPS):
             # print "setting ierr = {}".format(j)
             ierr = j



#       B(j,1) = alpha

        bbb_a[j] = alpha

#       if (alpha.ne.zero) then
#          call zsafescal(n,alpha,zwork(iv))
#       endif
        if alpha != 0.0:
            # print "alpha = {}; calling zsafescal()".format(alpha)
            zsafescal.zsafescal(n, alpha, zworkv)
            # for ips, value in enumerate(zworkv):
            #     print "zworkv[{}] = {}".format(ips + 1, value)

#       call zcopy(n,zwork(iv),1,V(1,j),1)

        # PS - This is the slow but simple Python equivalent of the fortran 
        # call to zcopy() --
        # for i, value in enumerate(zworkv):
        #     vvv[i, j] = value
            #print "V[{}, {}] = {}".format(i + 1, j, vvv[i, j - 1])

        # This is the fast version --
        vvv[:, j] = zworkv.copy()

 
# c     %------------------------------------------------%
# c     | beta_{j+1} u_{j+1} = A*v_{j} - alpha_{j} u_{j} |
# c     %------------------------------------------------%

             
#       call aprodw('n',ndp,m,n,zwork(iv),zwork(is),zwork(iz1),
#      c            zwork(iy1),lambda,trlambda, planF,planB )

        # print "before aprodw 982"
        # if j == 8:
        #     f = open("8p.txt", 'w')
        #     for value in zworkv:
        #         f.write("%.17f\t%.17f\n" % (value.real, value.imag))
        #     f.close()
        #     kjhafsdkgjadfh


        # for ips, value in enumerate(zworkv[:10]):
        #     print "zworkv[{}] = {}".format(ips, value)
        # print '-'
        # for ips, value in enumerate(zworks[:10]):
        #     print "zworks[{}] = {}".format(ips, value)
        # print '-'
        # for ips, value in enumerate(zworkz1[:10]):
        #     print "zworkz1[{}] = {}".format(ips, value)
        # print '-'
        # for ips, value in enumerate(zworky1[:10]):
        #     print "zworky1[{}] = {}".format(ips, value)
        # if j in (5, ):
        #     pdb.set_trace()
        aprodw.aprodw('n', ndp, m, n, zworkv, zworks, zworkz1,
                      zworky1, lambda_, trlambda)
        # if j == 8:
        #     f = open("8p.txt", 'w')
        #     for value in zworkv:
        #         f.write("%.17f\t%.17f\n" % (value.real, value.imag))
        #     f.close()
        #     kjhafsdkgjadfh


        # print "after aprodw 982"
        # for ips, value in enumerate(zworkv[:10]):
        #     print "zworkv[{}] = {}".format(ips, value)
        # print '-'
        # for ips, value in enumerate(zworks[:10]):
        #     print "zworks[{}] = {}".format(ips, value)
        # print '-'
        # for ips, value in enumerate(zworkz1[:10]):
        #     print "zworkz1[{}] = {}".format(ips, value)
        # print '-'
        # for ips, value in enumerate(zworky1[:10]):
        #     print "zworky1[{}] = {}".format(ips, value)
        # print '-'


#       do i=0,m-1
#          zwork(iu+i) = zwork(is+i) - alpha*zwork(iu+i)
#       enddo

        # PS - This is the naive (slow) Python port of the loop above -- 
#         for i in range(m):
# # #            print "zwork[{}] = {}".format(i, zwork[i])
# #             if i < 10:
# #                 print "zwork[{}] = {}".format(i, zwork[i])
# #                 print "zworks[{}] = {}".format(i, zworks[i])
#             zwork[i] = zworks[i] - (alpha * zwork[i])

        # PS - This is the fast version -- 
        zwork[:m] = zworks[:m] - (alpha * zwork[:m])


# c     %------------------------------------%
# c     | Extended local reorthogonalization |
# c     %------------------------------------%
#       do i=1,ioption(2)
# C        PS - call to zdotc() does not modify U
#          s = zdotc(m,U(1,j),1,zwork(iu),1)

# C        PS - I don't think this call to zaxpy() modifies U
#          call zaxpy(m,-s,U(1,j),1,zwork(iu),1)
#          if (alpha .ne. zero) then
#             alpha = alpha + dble(s)

#             B(j,1) = alpha
#          endif
#          work(imu+j-1) = eps
#       enddo

        for i in range(options["extended_reorth_iterations"]):
            #s = zdotc(m, uuu[:, j], 1, zwork, 1)
            # for ips, value in enumerate(uuu[:, j]):
            #     print "vdot input[{}] = {}".format(ips, value)
            # for ips, value in enumerate(zwork):
            #     print "vdot input[{}] = {}".format(ips, value)

            s = np.vdot(uuu[:, j], zwork)
            # print "vdot s = {}".format(s)

            # PS - This is the naive (and slow) equivalent of the fortran 
            # call to zaxpy() --
            # for ips in range(m):
            #     zwork[ips] = zwork[ips] + (-s * uuu[ips, j])

            # PS - This is the fast version --
            zwork[:m] += (-s * uuu[:m, j])

            if alpha != 0.0:
                alpha += s
                bbb_a[j] = alpha.real

            # FIXME - why is this inside a for loop? It is invariant. The
            # Fortran code does the same.
            mus[j] = EPS

#       beta= dznrm2(m,zwork(iu),1)
#       B(j,2) = beta
#       amax = max(amax,beta)

#       DZNRM2 returns the euclidean norm, hopefully this scipy call does
#       the same.
        beta = scipy.linalg.norm(zwork)
        # print "beta = {}".format(beta)
        bbb_b[j] = beta
        amax = max(amax, beta)
      
# c     %---------------------------%
# c     | Update estimate of ||A||_2 |         
# c     %---------------------------%
#       if (j.eq.1) then
#          a1 = dlapy2(B(1,1), B(1,2))
#       else
#          a1 = B(j,1)/amax
#          a1 = amax*sqrt(a1**2 + (B(j,2)/amax)**2 +
#   c           a1*B(j-1,2)/amax)
#       endif      
#       anorm=max(anorm,a1)

        if j == 0:
            #print "dlapy2 args = {}, {}".format(bbb_a[0], bbb_b[0])
            a1 = util.dlapy2(bbb_a[0], bbb_b[0])
        else:
            a1 = bbb_b[j] / amax
            a1 = amax * sqrt(a1 ** 2 + \
                             (bbb_b[j] / amax) ** 2 + 
                             a1 * bbb_b[j - 1] / amax
                            )
        # print "a1 = {}".format(a1)
        anorm = max(anorm, a1)

# c     %--------------------------%
# c     | Update the mu recurrence |
# c     %--------------------------%
#       if (.not.full_reorth .and. beta.ne.zero) then
#          call dupdate_mu(mumax,work(imu),work(inu),j,B(1,1),
#   c           B(1,2),anorm)
#       endif

        if (not full_reorth) and (beta != 0.0):
            # print "before dupdate_mu() 967, anorm = {}".format(anorm)
            # for ips, value in enumerate(nus):
            #     print "nus[{}] = {}".format(ips + 1, value)
            # for ips, value in enumerate(mus[:10]):
            #     print "mus[{}] = {}".format(ips + 1, value)
            # for ips, value in enumerate(nus[:10]):
            #     print "nus[{}] = {}".format(ips + 1, value)
            # for ips, value in enumerate(bbb_a[:10]):
            #     print "bbb_a[{}] = {}".format(ips + 1, value)
            # for ips, value in enumerate(bbb_b[:10]):
            #     print "bbb_b[{}] = {}".format(ips + 1, value)
            # if j == 1:
            #     pdb.set_trace()
            mumax = dupdate_mu(mus, nus, j, bbb_a, bbb_b, anorm)
            # print "after dupdate_mu() 967, mumax = {}".format(mumax)
            # for ips, value in enumerate(mus[:10]):
            #     print "mus[{}] = {}".format(ips + 1, value)
            # if j == 2:
            #     lkdfkgkhasdh

# c     %--------------------------------------%
# c     | Reorthogonalize u_{j+1} if necessary |
# c     %--------------------------------------%
#       if ( .true. .and.
#   c        (full_reorth .or. mumax.gt.delta .or. force_reorth) 
#   c        .and. beta.ne.zero) then
#          if (full_reorth .or. eta.eq.zero) then
#             iwork(iidx) = 1
#             iwork(iidx+1) = j
#             iwork(iidx+2) = j+1
#          else if (.not. force_reorth) then
#             call dcompute_int(work(imu),j,delta,eta,iwork(iidx))
#          else
#             do i=1,j+1
#                if (iwork(iidx+i-1).eq.j) then
#                   iwork(iidx+i-1) = j+1
#                   goto 25
#                endif
#             enddo
#          endif

        if (beta != 0.0) and (full_reorth or force_reorth or (mumax > delta)):
            # print "Reorthogonalizing u_{j+1}"
            if full_reorth or (eta == 0.0):
                # print "1109 setting indices"
                # indices[iidx]     = 1
                # indices[iidx + 1] = j
                # indices[iidx + 2] = j + 1
                indices[0] = 0
                indices[1] = j
                indices[2] = j + 1
                indices[3] = -1
                indices[4] = -1
            elif not force_reorth:
                # if j == 2:
                #     pdb.set_trace()
                # print "before dcompute_int 1278"
                indices = dcompute_int(mus, j, delta, eta)
            else:
                # ==> force_reorth = True
                # print "1118 setting indices"
                for i in range(j + 1):
                    if indices[i] == j:
                        indices[i] = j + 1
                        break

#  25         call second(t2)

# c 25         if (MGS.eq.0) then
# c 25       call zreorth(m,j,U,ldu,zwork(iu),beta,iwork(iidx),
# c     c              kappa, zwork(is),ioption(1))
# c            else
#                call zreorth2(m,j,U,ldu,zwork(iu),beta,iwork(iidx))
# c            endif

            # print "before zreorth2 1303"
            #pdb.set_trace()
            zreorth.zreorth2(m, j, uuu, zwork, beta, indices)
            
#           call dset_mu(j,work(imu),iwork(iidx),epsn2)
#           mumax = eta
#           if (force_reorth) then
#              force_reorth = .false.
#           else
#              force_reorth = .true.
#           endif
            # if j == 1:
            #     pdb.set_trace()
            # print "before dset_mu() 1146"
            # for ips, value in enumerate(mus[:10]):
            #     print "mus[{}] = {}".format(ips, value)
            # if j == 2:
            #     pdb.set_trace()
            dset_mu(j, mus, indices, EPSN2)       
            # print "after dset_mu() 1146"
            # for ips, value in enumerate(mus[:10]):
            #     print "mus[{}] = {}".format(ips, value)

            mumax = eta
            force_reorth = not force_reorth

#       endif


# c     %-----------------------------------------------%
# c     | Check whether an invariant subspace was found |
# c     %-----------------------------------------------%
#       if (beta .lt. anorm*epsn .and. j.lt.k) then
#          rnorm = beta
#          beta = zero

        if (beta < (anorm * EPSN)) and (j < k):
            rnorm = beta
            beta = 0.0
# c         %------------------------------------------------%
# c         | Try to build an orthogonal subspace, starting |
# c         | with a random vector.                          |
# c         %------------------------------------------------%
#           call zgetu0w('n',ndp, m, n, j, 3, zwork(iu), beta, U, ldu, 
#      c           ierr,ioption(1),anormest,zwork(is),
#      c           lambda,trlambda,planF,planB)

            beta, anormest, ierr = zget0w.zgetu0w('n', ndp, m, n, j, 3, zwork,
                                                  uuu, options["classical_gs"],
                                                  lambda_, trlambda)
#           if (beta .eq. zero) then
#              k = j
#              ierr = -j
#              goto 9999
#           else
#              nrstrt = nrstrt + 1
#              call donothing()
#              call zsafescal(n,beta,zwork(iu))
#              beta = zero
#              force_reorth = .true.
#              if (delta .gt. zero) then
#                 full_reorth = .false.
#              endif
#           endif
            if beta == 0.0:
                k = j
                ierr = -j
                # FIXME This should merely exit the loop, not raise an error
                raise ValueError
            else:
                zsafescal.zsafescal(n, beta, zwork)
                beta = 0.0
                force_reorth = True
                if delta > 0:
                    full_reorth = False

#        else if (.not.full_reorth .and. j.lt.k .and. 
#      c           (delta*beta .lt. anorm*eps)) then
# c         full_reorth = .true.
#           debito=.true.
#           ierr = j
#        endif            

        elif not full_reorth and (j < k) and ((delta * beta) < (anorm * EPS)):
            ierr = j 

#       B(j,2) = beta

        bbb_b[j] = beta

#       if (beta.ne.zero .and. beta.ne.one) then
#          call donothing()
#          call zsafescal(m,beta,zwork(iu))
#       endif

        if (beta != 0) and (beta != 1):
            zsafescal.zsafescal(m, beta, zwork)

#       call zcopy(m,zwork(iu),1,U(1,j+1),1)

        # PS - this is the naive equivalent of the call to zcopy() -- 
        # for ips in range(m):
        #     uuu[ips, j + 1] = zwork[ips]

        # PS - this is the fast version
        uuu[:m, j + 1] = zwork[:m].copy()

 

#       rnorm = beta
#       call second(t2)
        rnorm = beta

#      enddo
#  9999 doption(3) = anorm      
#       call second(t2)
#       tlanbpro = tlanbpro + (t2-t1)
    options["anorm"] = anorm


#  100  format(1i6,6e12.3)
#  101  format(1i6,4e12.3)

    # print "zlanbprow done"

    return ierr, rnorm

#       return
#       end

# c
# c**********************************************************************
# c

#       subroutine dset_mu(k,mu,index,val)
# c     %-----------%
# c     | Arguments |
# c     %-----------%
#       implicit none
#       integer k,index(*)
#       double precision mu(*), val


def dset_mu(k, mus, indices, val):

# c     %-----------------%
# c     | Local variables |
# c     %-----------------%
#       integer i,j,p,q

#       i=1
#       do while(index(i).le.k .and. index(i).gt.0)
#          p = index(i)
#          q = index(i+1)
#          do j=p,q
#             mu(j) = val
#          enddo

#          i = i+2
#       enddo      
#       end

    i = 0
    while (indices[i] <= k) and (indices[i] >= 0):
        p = indices[i]
        q = indices[i + 1]
        # FIXME check for off-by-one error
        for j in fortran_range(p, q):
            mus[j] = val

        i += 2



# c
# c**********************************************************************
# c
#       subroutine dcompute_int(mu,j,delta,eta,index)
#def dcompute_int(mus, j, delta, eta, indices):
def dcompute_int(mus, j, delta, eta):
    """I think "int" is short for "intervals". returns indices."""

# c     %-----------%
# c     | Arguments |
# c     %-----------%
#       implicit none
#       include 'stat.h'
#       integer j,index(*)
#       double precision mu(*)
#       double precision delta,eta

# c     %-----------------%
# c     | Local variables |
# c     %-----------------%
#       integer i,k,s,ip
#       real t1,t2

#   if (delta.lt.eta) then
#      return
#   endif

    # print "begin dcompute_int, j = {}".format(j)
    # print "delta = {}, eta = {}".format(delta, eta)

    # for ips, value in enumerate(mus):
    #     print "mus[{}] = {}".format(ips, value)
    #     if not value:
    #         break

    #indices = [0, 0, 0]
    indices = [ ]

    ip = 0

    i = -1
    # if j == 1:
    #     pdb.set_trace()
    #pdb.set_trace()
    if delta >= eta:
        while i < j: 
            # find the next mu(k), k>i where abs(mu(k)) > delta.
            # If there is no such mu(k), set k = i + 1.
            k = i + 1

            while (k < j) and (abs(mus[k]) <= delta):
                k += 1

            # FIXME the 'if' block below seemed like a good idea but
            # doesn't actually work. Que pasa???
            # if k >= j:
            #     # If the loop above completes without finding
            #     # (abs(mus[k]) <= delta), it should skip the rest. Fortran
            #     # implements this with a goto.
            #     break

            # find smallest s<k such that for all j=s,..,k, m(j) >= eta
            #the_range = [value for value in reversed(range(max(i, 0), k))]
            # the_range will be something like [3, 2, 1, 0] and .pop() will 
            # remove entries from the right-hand side.
            s = k
            done = False
            while not done:
                s -= 1
                if s < max(i, 0):
                    # Oooops, we're off the end of our search range.
                    done = True
                else:
                    if abs(mus[s]) < eta:
                        done = True
            # We actually want to record the index one to the right of 
            # the position we found.
            s += 1
            indices.append(s)

            # Now search to the right
            found = False
            for i in range(s, j):
                if abs(mus[i]) < eta:
                    found = True
                    break

            # When a loop iterator is exhausted in fortran, the iterator 
            # variable is set to the end of the range + 1. This code mimics
            # that behavior.
            if not found:
                i += 1

            indices.append(i)

        indices.append(j + 1)
        # Append end-of-list flags
        indices.append(-1)
        indices.append(-1)
    # else:
    #     delta < eta, so raise an error?

    
    # print "end dcompute_int"
    # for ips, value in enumerate(indices):
    #     print "indices({}) = {}".format(ips, value)
    

    return indices


#       ip = 0
#       index(1) = 0
#       i=0
#         ip = -1
#         indices[0] = 0
#         i = -1
# #       do while(i.lt.j)
#         done = False
#         while (i < j) and not done:
# # c     find the next mu(k), k>i where abs(mu(k)) > delta
# #           do k=i+1,j
# #              if (abs(mu(k)).gt.delta) goto 10
# #           enddo         
# #           goto 40
#             done = True
            
#             # This dance around the_range results from the fact that Fortran 
#             # handles pathological ranges differently than Python. 
#             # range(a, b) generates an empty tuple when a == b. In Fortran,
#             # do i = a,b ==> execute the loop once.
#             # if (i + 1) == j:
#             #     the_range = (j, )
#             # else:
#             #     the_range = range(i + 1, j)

#             the_range = fortran_range(i + 1, j)

#             for k in the_range:
#                 if abs(mus[k]) > delta:
#                     done = False
#                     break

# # c         find smallest i<k such that for all j=i,..,k, m(j) >= eta
# #  10       do s=k,max(i,1),-1
# #              if (abs(mu(s)).lt.eta) goto 20
# #           enddo
#             # The code below relies on a Fortran subtlety. At the end of a do
#             # loop, the loop index (in this case, 's') is one beyond the end
#             # of the loop. Since we have a negative step in this case, s will
#             # be one *less* than the end of the loop range.
#             the_range = range(k, max(i, 1), -1)
#             if not the_range:
#                 # If the_range is empty, the for loop won't execute
#                 # at all and 's' will be undefined, so I have to explicitly
#                 # give it a value here.
#                 s = max(i, 1) - 1
#             for s in the_range:
#                 if abs(mus[s]) < eta:
#                     break
#             print "dcompute_int: k = {}, the_range = {}, s = {}".format(k, the_range, s)       
# #  20       ip= ip+1
# #           index(ip) = s+1
#             ip += 1
#             indices[ip] = s
# #           do i=s+1,j            
# #              if (abs(mu(i)).lt.eta) goto 30
# #           enddo    

#             the_range = fortran_range(s, j)     

#             exited_early = False
#             for i in the_range:
#                 if abs(mus[i]) < eta:
#                     exited_early = True
#                     break

#             if not exited_early:
#                 # Account for difference between Python and Fortran loop counter
#                 # incrementation.
#                 i += 1
# #  30       ip= ip+1
# #           index(ip) = i-1
#             ip += 1
            
#             # try:
#             #     indices[ip]
#             # except:
#             #     pdb.set_trace()    
#             indices[ip] = i - 1
# #       enddo
# #  40   ip = ip+1
# #       index(ip) = j+1
#         ip += 1
#         indices[ip] = j 

#     print "Exiting dcompute_int, ip = {}".format(ip)
#     for ips, value in enumerate(indices[:10]):
#         print "indices[{}] = {}".format(ips, value)

# c      write (*,*) 'i  index(i)'
# c      do i=1,j
# c         write(*,*) i,index(i)
# c         if (index(i).gt.j) goto 50
# c      enddo
# c 50   write (*,*) 'Exit compute_int, ip=',ip
#       call second(t2)
#       tintv = tintv + (t2-t1)
#       end
# c
# c**********************************************************************
# c
#       subroutine dupdate_mu(mumax,mu,nu,j,alpha,beta,anorm)

#            mumax = dupdate_mu(mus, nus, j, bbb_a, bbb_b, anorm)

def dupdate_mu(mus, nus, j, alpha, beta, anorm):
    """Alters mus, returns mumax. I'm not sure if this is truly the max value
    in abs(mus) or just the max of the subset that this function touches."""

# c     %-----------%
# c     | Arguments |
# c     %-----------%
#       implicit none
#       include 'stat.h'
#       integer j,formula
#       parameter(formula=1)
#       double precision mumax,eps1,eps,anorm
#       double precision mu(*),nu(*),alpha(*),beta(*)

# c     %------------%
# c     | Parameters |
# c     %------------%
#       double precision one, zero, FUDGE
#       parameter(one = 1.0d0, zero = 0.0d0, FUDGE = 1.01)

# c     %-----------------%
# c     | Local variables |
# c     %-----------------%
#       double precision d
#       integer k
#       real t1,t2
      
# c     %--------------------%
# c     | External Functions |
# c     %--------------------%
#       double precision dlamch,dlapy2
#       external dlamch,dlapy2

# c      write (*,*) 'Enter update_mu'
#       call second(t1)
#       eps = dlamch('e')
#       eps1 = 100*eps

    # In the Fortran, forumla is hardcoded to 1 
#       if (formula.eq.1) then



#       if (j.eq.1) then
#          d = eps1*(dlapy2(alpha(j), beta(j)) + alpha(1)) + eps1*anorm
#          mu(1) = eps1/beta(1)
#          mumax = abs(mu(1))
    if j == 0:
        # FIXME is this a bug in the fortran? d is computed but not used.
        d = EPS1 * (util.dlapy2(alpha[j], beta[j]) + alpha[0]) + \
            EPS1 * anorm
        mus[0] = EPS1 / beta[0]
        mumax = abs(mus[0])

    else:

#   else
#       mu(1) = alpha(1)*nu(1)-alpha(j)*mu(1)
#       d = eps1*(dlapy2(alpha(j), beta(j)) + alpha(1)) + eps1*anorm
#       mu(1) = (mu(1) + dsign(d,mu(1))) / beta(j)
#       mumax = abs(mu(1))
        #pdb.set_trace()
        mus[0] = alpha[0] * nus[0] - alpha[j] * mus[0]
        d = EPS1 * (util.dlapy2(alpha[j], beta[j]) + alpha[0]) + \
            EPS1 * anorm
        mus[0] = (mus[0] + math.copysign(d, mus[0])) / beta[j]
        mumax = abs(mus[0])

#       do k=2,j-1
#          mu(k) = alpha(k)*nu(k) +beta(k-1)*nu(k-1)-alpha(j)*mu(k)
#          d = eps1*(dlapy2(alpha(j), beta(j)) + 
#    c         dlapy2(alpha(k), beta(k-1))) + eps1*anorm
#          mu(k) = (mu(k) + dsign(d,mu(k))) / beta(j)
#          mumax = max(mumax,abs(mu(k)))
#       enddo
        # FIXME dlapy2(alpha[j], beta[j]) is a loop invariant 
        #the_range = fortran_range(1, j)
        the_range = fortran_range(1, j - 1)
        # if not the_range:
        #     the_range = (1, )
        for k in the_range:
            mus[k] = (alpha[k] * nus[k]) + (beta[k - 1] * nus[k - 1]) - \
                     (alpha[j] * mus[k])
            d = EPS1 * (util.dlapy2(alpha[j], beta[j])        + \
                        util.dlapy2(alpha[k], beta[k - 1]))   + \
                EPS1 * anorm
            mus[k] = (mus[k] + math.copysign(d, mus[k])) / beta[j]
            mumax = max(mumax, abs(mus[k]))


#       mu(j) = beta(j-1)*nu(j-1)
#       d = eps1*(dlapy2(alpha(j), beta(j)) + 
#  c              dlapy2(alpha(j), beta(j-1))) + eps1*anorm
#       mu(j) = (mu(j) + sign(d,mu(j))) / beta(j)
#       mumax = max(mumax,abs(mu(j)))
        mus[j] = beta[j - 1] * nus[j - 1]
        d = EPS1 * (util.dlapy2(alpha[j], beta[j]) + util.dlapy2(alpha[j], beta[j - 1])) + \
            EPS1 * anorm
        mus[j] = (mus[j] + math.copysign(d, mus[j])) / beta[j]
        mumax = max(mumax, abs(mus[j]))


#   endif
#   mu(j+1) = one

    mus[j + 1] = 1.0


    # else case is formula = 0
#       else
#          if (j.eq.1) then
#             mu(1) = eps1/beta(1)
#             mumax = abs(mu(1))
#          else
#             mu(1) = alpha(1)*nu(1)-alpha(j)*mu(1)
#             mu(1) = (mu(1) + dsign(eps1,mu(1))) / beta(j)
#             mumax = abs(mu(1))
#             do k=2,j-1
#                mu(k) = alpha(k)*nu(k) + beta(k-1)*nu(k-1)-alpha(j)*mu(k)
#                mu(k) = (mu(k) + dsign(eps1,mu(k))) / beta(j)
#                mumax = max(mumax,abs(mu(k)))
#             enddo
#             mu(j) = beta(j-1)*nu(j-1)
#             mu(j) = (mu(j) + sign(eps1,mu(j))) / beta(j)
#             mumax = max(mumax,abs(mu(j)))
#          endif
#          mu(j+1) = one
#       endif
#       call second(t2)
#       tupdmu = tupdmu + (t2-t1)
# c      write (*,*) 'Exit update_mu'
#       end

    return mumax




# c
# c**********************************************************************
# c
#       subroutine dupdate_nu(numax,mu,nu,j,alpha,beta,anorm)
def dupdate_nu(mus, nus, j, alpha, beta, anorm):
    """Alters nus, returns numax"""

#   eps = dlamch('e')
#   eps1 = 100*eps

    # In the Fortran, forumla is hardcoded to 1 
#   if (formula.eq.1) then

#      if (j.gt.1) then
#         numax = zero
#         do k=1,j-1
#            nu(k) = beta(k)*mu(k+1) + alpha(k)*mu(k) -beta(j-1)*nu(k)
#            d = eps1*(dlapy2(alpha(k),beta(k)) +
#  c              dlapy2(alpha(j),beta(j-1))) + eps1*anorm
#            nu(k) = (nu(k) + dsign(d,nu(k))) / alpha(j)
#            numax = max(numax,abs(nu(k)))
#         enddo
#         nu(j) = one
#      endif
    numax = 0.0
    if j > 0:
        #pdb.set_trace()
        for k in fortran_range(0, j - 1):
            nus[k] = (beta[k] * mus[k + 1]) + (alpha[k] * mus[k]) - \
                     (beta[j - 1] * nus[k])
            d = EPS1 * (util.dlapy2(alpha[k], beta[k]) + \
                        util.dlapy2(alpha[j], beta[j - 1])
                       )
            d += EPS1 * anorm
            nus[k] = (nus[k] + math.copysign(d, nus[k])) / alpha[j]
            numax = max(numax, abs(nus[k]))

        nus[j] = 1.0
#   else
#      if (j.gt.1) then
#         numax = zero
#         do k=1,j-1
#            nu(k) = beta(k)*mu(k+1) + alpha(k)*mu(k) -beta(j-1)*nu(k)
#            nu(k) = (nu(k) + dsign(eps1,nu(k))) / alpha(j)
#            numax = max(numax,abs(nu(k)))
#         enddo
#         nu(j) = one
#      endif
#   endif

    return numax
