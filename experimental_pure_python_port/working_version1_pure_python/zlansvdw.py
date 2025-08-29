# Python modules
from __future__ import division
import ctypes
import math
import pprint
import pdb
pp = pprint.pprint

# 3rd party modules
import numpy as np
import numpy.linalg
import scipy
import scipy.linalg.blas
import scipy.linalg.lapack

# Our modules
import zget0w
import dlartg
import dbdsqr
import util
import zgemmina
import zlanbprow

MAX_SINGULAR_VALUES = 50


EPS = np.finfo(np.float64).eps
EPS34 = EPS ** (3.0 / 4.0)

# PS - These control flow in a few parts of the code. They're hardcoded in
# the Fortran we have, although they're probably passed as params in the 
# original PROPACK code.
JOBU = 'y'
JOBV = 'n'

# Variable name translations 
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# k          nsv_sought         # of singular values requested by caller
# kfit       nsv_found          # of singular values found
# lsinval    singular_values    Array containing the singular values (floats)
# ndp        n_data_points      
# kmax       MAX_SINGULAR_VALUES  Hardcoded to 50.
#           
#            options dict key
#            ----------------
# ioption(1) classical_gs       True ==> classical Gram-Schmidt reorth
#                               (see zreorth.py)
# ioption(2) extended_reorth_iterations  # of iterations in one case of 
#                               extended local reorth
# doption(1) delta              See zlanbprow() docstring
# doption(2) eta                See zlanbprow() docstring
# doption(3) anorm              See zlanbprow() docstring




# http://www.scipy.org/doc/api_docs/SciPy.lib.lapack.info.html

# Float BLAS functions
# drotg, = scipy.linalg.blas.get_blas_funcs( ['rotg'], 
#                                             np.array([0.0]) 
#                                           )

# Float LAPACK functions
dlamch, = scipy.linalg.lapack.get_lapack_funcs( ['lamch'], 
                                                np.array([np.float64(0.0)])
                                              )

# def zlansvdw(n_data_points, m, n, k, lambda_, trlambda, uuu):

#     libhlsvd = util.load_the_library()

#     f = libhlsvd.zlansvdw_python_

#     # if uuu is none, initialize it. Otherwise, split it into real/complex
#     # pairs in preparation for passing to zlansvdw_python().
#     UuuArrayType = ctypes.c_double * (m * (KMAX + 1))
#     # zlansvdw() is the call where U/uuu gets populated, so it's first 
#     # created here.
#     uuu_r = UuuArrayType()
#     uuu_i = UuuArrayType()
#     if uuu is not None:
#         for i, z in enumerate(uuu):
#             uuu_r[i] = z.real   
#             uuu_i[i] = z.imag


#     InputDoubleArrayType = ctypes.c_double * n_data_points
#     OutputDoubleArrayType = ctypes.c_double * MAX_SINGULAR_VALUES

#     # Input params
#     lambda_r   = InputDoubleArrayType()
#     lambda_i   = InputDoubleArrayType()
#     trlambda_r = InputDoubleArrayType()
#     trlambda_i = InputDoubleArrayType()


#     for i, z in enumerate(lambda_):
#         lambda_r[i] = z.real
#         lambda_i[i] = z.imag
#     for i, z in enumerate(trlambda):
#         trlambda_r[i] = z.real
#         trlambda_i[i] = z.imag


#     # Set up ctypes variables to call zlansvdw().
#     # In-line comments refer to the names of the corresponding variables
#     # in hlsvdpro.f.


#     nsv_sought = k

#     kmax = KMAX
#     ndpmax = NDPMAX

#     lsinval = OutputDoubleArrayType()               # sigma

#     for i, z in enumerate(lambda_):
#         lambda_r[i] = z.real   
#         lambda_i[i] = z.imag

#     #trlambda = trlambda.tolist()
#     for i, z in enumerate(trlambda):
#         trlambda_r[i] = z.real   
#         trlambda_i[i] = z.imag

#     n_data_points = ctypes.c_long(n_data_points)        # ndp
#     m = ctypes.c_long(m)                          # Lrow/m
#     n = ctypes.c_long(n)                          # mcoL/n
#     kmax = ctypes.c_long(kmax)                          # kmax
#     ndpmax = ctypes.c_long(ndpmax)                      # ndpmax

#     info = 0

#     nsv_sought = ctypes.c_long(nsv_sought)              # kuser/k/kfit
#     info = ctypes.c_long(info)                          # info (0, < 0, > 0)
#                                                         # 0 is good!

#     # Note - nsv_sought == k (alias kuser) going in. zlansvdw() changes this 
#     to the number of singular values it found (kfit).
#     libhlsvd.zlansvdw_python_(ctypes.pointer(n_data_points),
#                               ctypes.pointer(m), 
#                               ctypes.pointer(n), 
#                               ctypes.pointer(nsv_sought),
#                               ctypes.pointer(lsinval),
#                               ctypes.pointer(info),
#                               ctypes.pointer(lambda_r),
#                               ctypes.pointer(lambda_i),
#                               ctypes.pointer(trlambda_r),
#                               ctypes.pointer(trlambda_i),
#                               ctypes.pointer(uuu_r),
#                               ctypes.pointer(uuu_i),
#                              )

#     nsv_found = nsv_sought.value
#     info = info.value

#     uuu =     [complex(r, i) for r, i in zip(uuu_r, uuu_i)]
#     lambda_ = [complex(r, i) for r, i in zip(lambda_r, lambda_i)]
#     trlambda = [complex(r, i) for r, i in zip(trlambda_r, trlambda_i)]

#     # Convert this from a ctypes thingy to regular Python, and trim.
#     lsinval = [x for x in lsinval][:nsv_found]


#     return uuu, lsinval, info, nsv_found


def zlansvdw(n_data_points, m, n, nsv_sought, lambda_, trlambda, uuu, 
             tol=16e-16):

    # print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ zlansvdw() "
# c     %---------------------------------%
# c     | Set machine dependent constants |
# c     %---------------------------------%
#       eps = dlamch('e')
#       eps34 = eps**(3.0/4.0)
#       epsn = dble(max(m,n))*eps/2.0
#       epsn2 = sqrt(dble(max(m,n)))*eps/2.0
#       sfmin = dlamch('s')
    # PS - Fortran populates variable sfmin but otherwise doesn't use it.


    # options replaces the Fortran arrays ioption and doption.
    options = {}
    # Use classical Gram-Schmidt reorth instead of the modified version.
    options["classical_gs"] = True
    options["extended_reorth_iterations"] = 1
    options["delta"] = 1e-12
    options["eta"] = 1e-14
    options["anorm"] = 0.0

    #print m * (MAX_SINGULAR_VALUES + 1)
    
    EPSN = max(m,n) * EPS / 2.0

    if uuu is None:
        uuu = np.zeros( (m, MAX_SINGULAR_VALUES + 1), np.complex128)

    vvv = np.zeros( (n, MAX_SINGULAR_VALUES), np.complex128)

# c     %--------------------------------%
# c     | Guard against absurd arguments |
# c     %--------------------------------%
#       lanmax = min(n+1,m+1,kmax)
#       tol = min(one,max(16.0*eps,tolin))
#       anorm = zero
    lanmax = min(n + 1, m + 1, MAX_SINGULAR_VALUES)
    tol = min(1.0, max(16.0 * EPS, tol))
    anorm = 0.0

# c     %------------------------------%
# c     | Set pointers into work array |
# c     %------------------------------%
#   ibnd = 1
#   ib = ibnd + lanmax+1
#   ib1 = ib + 2*lanmax
#   ip = ib1 + 2*lanmax
#   iq = ip + (lanmax+1)**2
#   iwrk = iq + lanmax**2

    # The Fortran code uses multiple pointers into the arrays work and zwork
    # (and bnd). In Python I replace them with numpy arrays. Since I don't 
    # know what these arrays do, they have crummy names but at least they
    # harken back to the Fortran code.
    # Here's a map. The indices are a guide only; some of them will vary 
    # based on n_data_points.

    # Fortran        Fortran indices       Python
    # -------        ---------------       ------
    # work(ibnd)  =  work(   1 -   51)  =  work1
    # work(ib)    =  work(  52 -  151)  =  bbb_a, bbb_b
    # work(ib1)   =  work( 153 -  251)  =  bbb1_a, bbb1_b
    # work(ip)    =  work( 252 - 2852)  =  workp
    # work(iq)    =  work(2853 - 5352)  =  workq

    # Some of these lengths are adjusted by 1 to reflect the difference between
    # Fortran's 1-based arrays and Python's 0-based arrays. Since Python 
    # implements variable-length lists so nicely, it might be prefereable 
    # to convert some of these lists to dynamic sizing rather than 
    # preallocating them.
    len_work1  = lanmax + 1
    # In the Fortran code, work(ib) starts an array of length lanmax * 2. It's 
    # sometimes declared as matrix of dimension (lanmax, 2) but it's never used
    # as a matrix. It's really just two arrays (of length lanmax) bundled 
    # together under one variable name. It's confusing, and in Python I separate
    # them into bbb_a and bbb_b.
    # The chunk of work starting at work(ib1) is just used for a work copy of
    # work(ib), and has the same properties. I call the corresponding Python
    # arrays bbb1_a and bbb1_b.
    len_bbb = lanmax
    len_workp  = (lanmax + 1) ** 2
    # Note that before zlansdvw.f's second call to dbdsqr(), it writes a bunch
    # of 1s into work starting at index iq. It goes waaaay beyond lanmax**2
    # and into the work(iwrk) section. Intentional? I don't know. The size I
    # give to workq here should be about right, erring on the side of safety.
    len_workq  = n * (lanmax + 1)

# C     PS - changed from lwork to lzwrk
#       lwrk = lzwrk-iwrk+1
#       call dzero(7*lanmax + 2 + 2*lanmax**2,work,1)
#       call zzero(7*lanmax + 2 + 2*lanmax**2,zwork,1)


    work1  = np.zeros( (len_work1, ), dtype=np.float64)    # work(ibnd)
    workp  = np.zeros( (len_workp, ), dtype=np.float64)    # work(ip)
    workq  = np.zeros( (len_workq, ), dtype=np.float64)    # work(iq)

    bbb_a  = np.zeros( (len_bbb, ), dtype=np.float64)      # work(ib)
    bbb_b  = np.zeros( (len_bbb, ), dtype=np.float64)      # work(ib  + lanmax)
    bbb1_a = np.zeros( (len_bbb, ), dtype=np.float64)      # work(ib1)
    bbb1_b = np.zeros( (len_bbb, ), dtype=np.float64)      # work(ib1 + lanmax)


# c     %---------------------------------------------------------------%
# c     | Set up random starting vector if none is provided by the user |
# c     %---------------------------------------------------------------%
# c      rnorm = dnrm2(m,U(1,1),1)
#       rnorm =0.0d0
#       if (rnorm.eq.zero) then

    # FIXME this follows the convention in the Fortran code of hardcoding
    # the value for rnorm which makes the if statement below superfluous. 
    rnorm = 0.0
    if rnorm == 0:

#       print *, "calling zgetu0w()"
#           call zgetu0w('n',ndp, m,n,0,1,U,rnorm,U,ldu,
#      c        ierr,ioption(1),anorm,zwork(iwrk),lambda,
#      c        trlambda,planF,planB)

        rnorm, anorm, ierr = zget0w.zgetu0w('n', n_data_points, m, n, 0, 1, 
                                            uuu[:,0], 
                                            uuu, options["classical_gs"],
                                            lambda_, trlambda)
        # print "done zgetu0w(), rnorm = {}".format(rnorm)

    # info = 0
    # neig = 0
    # jold = 0
    # j = min(k+max(8,k)+1,lanmax)
    info = 0
    neig = 0
    jold = 0
    j = min(nsv_sought + max(8, nsv_sought) + 1, lanmax)

# c     %------------------------------%
# c     | Iterate until convergence... |
# c     %------------------------------%
#   do while (neig.lt.k)
    while neig < nsv_sought:
# c     %---------------------------------------------------%
# c     | Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j} |
# c     %---------------------------------------------------%
        # print "looptop, neig = {}, nsv_sought = {}".format(neig, nsv_sought)

        # for i in range(m):
        #     for j in range(KMAX + 1):
        #         print "uuu({}, {}) = {:.17G}".format(i + 1, j + 1, uuu[i, j])

        # for i in range(n):
        #     for j in range(KMAX):
        #         print "vvv({}, {}) = {:.17G}".format(i + 1, j + 1, vvv[i, j])

        # zlanbprow() modifies bbb_a and bbb_b in place.
        ierr, rnorm = zlanbprow.zlanbprow(n_data_points, m, n, jold, j, 
                                          uuu, vvv, bbb_a, bbb_b,
                                          rnorm, options, 
                                          lambda_, trlambda)


        # print "Python returning from zlanbprow.zlanbprow()"
#       jold = j
        jold = j

# c     %---------------------------------------------%
# c     | Compute and analyze SVD(B) and error bounds |
# c     %---------------------------------------------%
        
        # call dcopy(2*lanmax, work(ib),1,work(ib1),1)
        bbb1_a = bbb_a.copy()
        bbb1_b = bbb_b.copy()

        # call dzero(j+1,work(ibnd),1)
        work1[:j + 1] *= 0

        #      call dbdqr('N',j,work(ib1),work(ib1+lanmax),work(ibnd+j-1),
        # c         work(ibnd+j),work(ip),lanmax+1)
        # print "j =  {}".format(j)
        
        # for i, value in enumerate(workp):
        #     print "workp[{}] = {:.17E}".format(i + 1, workp[i])

        # print "^^^"

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {:.17E}".format(i + 1, work1[i])

        # print "^^^"

        # print "before dbdqr1(), j = {}".format(j)
        c1, c2 = dbdqr('N', j, bbb1_a, bbb1_b, workp, lanmax + 1)

        work1[j - 1] = c1
        work1[j] = c2

        # print "after dbdqr1(), c1 = {}, c2 = {}".format(c1, c2)

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {:.17E}".format(i + 1, work1[i])

        # print "^^^"

        # for i, value in enumerate(workp):
        #     print "workp[{}] = {:.17E}".format(i + 1, workp[i])

        # print "^^^"

        # print "before dbdsqr, j = {}" .format(j)

        # if work1[3]:
        #     for i in range(42):
        #         work1[i] = 0.0

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {}".format(i + 1, value)

        # for i, value in enumerate(bbb1_a):
        #     print "bbb1_a[{}] = {}".format(i + 1, value)

        # for i, value in enumerate(bbb1_b):
        #     print "bbb1_b[{}] = {}".format(i + 1, value)

        # Params are       uplo, n, ncvt, nru, ncc, ...

        #      call dbdsqr('u',j,0,1,0,work(ib1),work(ib1+lanmax),work,1,
        # c         work(ibnd),1,work,1,work(iwrk),info)

        xj      = j
        xbbb1_a = bbb1_a.copy()
        xbbb1_b = bbb1_b.copy()
        xwork1  = work1.copy()

#         r = dbdsqr.dbdsqr('u', j, 0, 1, 0, bbb1_a, bbb1_b, None, 1, work1, 1, None, 1)
#         info, bbb1_a, bbb1_b, _, work1, _ = r

        r = dbdsqr.dbdsqr_cython('u', xj, 0, 1, 0, bbb1_a, bbb1_b, None, 1, work1, 1, None, 1)
        info, bbb1_a, bbb1_b, _, work1, _ = r


        # print "after dbdsqr" 
        
        # for i, value in enumerate(work1):
        #     print "work1[{}] = {}".format(i + 1, value)

        # PS - dbdsqr() returns bbb1_a, bbb1_b, and work1 as lists that are 
        # shorter than they were on the way in. Subsequent code (refinebounds() 
        # and probably other code too) doesn't deal with this size change,
        # so I have to pad them back out to their original size.
        bbb1_a +=  [0.0] * (lanmax - len(bbb1_a) + 1)
        bbb1_b +=  [0.0] * (lanmax - len(bbb1_b) + 1)
        # print "len(bbb1_b) = {}".format(len(bbb1_b))
        bbb1_a = np.array(bbb1_a)
        bbb1_b = np.array(bbb1_b)
        work1 += [0.0] * (len_work1 - len(work1))
        work1 = np.array(work1)

        # if (j.gt.5) then
        #    anorm = work(ib1)
        # else
        #    anorm = max(anorm,work(ib1))
        # endif
        if j > 5:
            anorm = bbb1_a[0]
        else:
            anorm = max(anorm, bbb1_a[0])     

        # print "anorm = {}".format(anorm)

        # do i=1,j
        #    work(ibnd+i-1) = abs(rnorm*work(ibnd+i-1))
        #    work(ib1+lanmax+i-1) = work(ib1+lanmax+i-1)**2
        # enddo

        work1 *= rnorm
        np.absolute(work1, work1)
        bbb1_b **= 2

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {}".format(i + 1, value)

        # for i, value in enumerate(bbb1_b):
        #     print "bbb1_b[{}] = {}".format(i + 1, value)

# c     %---------------------------------------------%
# c     | Refine error bounds using the "Gap theorem" |
# c     %---------------------------------------------%

#          call refinebounds(j,work(ib1+lanmax),work(ibnd),
#      c        epsn*anorm,eps34)

        # print "before refinebounds"

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {}".format(i + 1, value)

        refinebounds(j, bbb1_b, work1, EPSN * anorm)

        # print "after refinebounds"

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {}".format(i + 1, value)

# c     %----------------------------------------------------%
# c     | Determine the number of converged singular values  |
# c     %----------------------------------------------------%
#          do i=1,min(j,k)
#             bnd(i) = work(ibnd+i-1)
#          enddo

        i = min(j, nsv_sought)
        bounds = work1[:i]

        # for i, value in enumerate(bounds):
        #     print "bounds[{}] = {}".format(i + 1, value)


        # i = 0
        # neig = 0
        # do while(i.lt.min(j,k))
        #    if (work(ibnd+i).le.tol*work(ib1+i)) then
        #       neig = neig + 1
        #       sigma(neig) = work(ib1+i)
        #       i = i+1
        #    else
        #       i = k
        #    endif
        # enddo
        sigma = [ ]
        i = 0
        # neig = 0
        # print "tol = {}".format(tol)
        while (i < min(j, nsv_sought)):
            # print "work[{}] = {}".format(i + 1, work1[i])
            # print "tol * bbb1_a[{}] = {}".format(i + 1, tol * bbb1_a[i])
            #print "sigma delta = {}".format((tol * bbb1_a[i]) - work1[i])
            if work1[i] <= (tol * bbb1_a[i]):
                # neig += 1
                sigma.append(bbb1_a[i])
                i += 1
            else:
                # Force loop exit
                i = nsv_sought

        # for i, value in enumerate(sigma):
        #     print "sigma[{}] = {}".format(i + 1, value)

        neig = len(sigma)


# c     %--------------------------------------------------%
# c     | Test if an invariant subspace have been found or |
# c     | the workspace has been exhausted.                |
# c     %--------------------------------------------------%
#          if (ierr.lt.0) then
#             goto 50               
#          endif
#          if (j.ge.lanmax) then
#             if (neig.lt.k) then
#                info = -1
#             endif
#             goto 50
#          endif
        if ierr < 0:
            # Bail out of loop
            # print "bailing, ierr = {}".format(ierr)
            break

        if j >= lanmax:
            if neig < nsv_sought:
                # print "bailing"
                # print "j = {}, lanmax = {}, neig = {}, nsv_sought = {}".format(j, lanmax, neig, nsv_sought)
                info = -1
            # Bail out of loop            
            break

# c     %----------------------------------------------------%
# c     | Increase the dimension of the Krylov subspace.     |
# c     | If any Ritz values have converged then try to      | 
# c     | estimate the average number of iterations per      |
# c     | converged Ritz value.                              |
# c     | Else increase the dimension with 50%.              |
# c     %----------------------------------------------------%
#          if (neig.gt.1) then
#             dj = min(j/2,((k-neig)*(j-6))/(2*neig+1))
#             dj = min(100,max(2,dj))
#          else
#             dj = j/2
#             dj = min(100,max(10,dj))
#         endif
#          j = min(j + dj,lanmax)

        # print "monkey, j = {}".format(j)
        # print "monkey, k = {}".format(nsv_sought)
        # print "monkey, neig = {}".format(neig)
        # print "monkey, form = {}".format(((nsv_sought - neig) * (j - 6)) // ( 2 * neig + 1))

        if neig > 1:
            dj = min(j // 2, ((nsv_sought - neig) * (j - 6)) // ( 2 * neig + 1))
            dj = min(100, max(2, dj))
        else:
            dj = j // 2
            dj = min(100, max(10, dj))

        j = min(j + dj, lanmax)

        # print "loop end, j = {}, dj = {}".format(j, dj)
        # print "loop end, neig = {}, k = {}".format(neig, nsv_sought)


##############             while loop has ended      ##################

    # print "Pre-calculate singular vectors"
    # print "neig = {}, nsv_sought = {}".format(neig, nsv_sought)

#  50   if (neig.ge.k) then
# c     %-----------------------------------------%
# c     | Calculate singular vectors if requested %
# c     %-----------------------------------------%
    if neig >= nsv_sought:
        # print "inside Calculate singular vectors"
        # j = jold
        j = jold

        # In Python, I don't attempt to modify 'k' as the Fortran code does.
        # I just return neig directly.
        # k = neig


        # if (lsame(jobu,'y') .or. lsame(jobv,'y')) then
        #    print *, "about to dcopy"
        #    call dcopy(2*lanmax, work(ib),1,work(ib1),1)
        #    print *, "done dcopy"

        if (JOBU == 'y') or (JOBV == 'y'):
            bbb1_a = bbb_a.copy()
            bbb1_b = bbb_b.copy()

            # %--------------------------------------------------------------
            # | The bidiagonal SVD is computed in a two-stage procedure:             
            # |                                                   
            # | 1. Compute a QR-factorization M^T*B = [R; 0] of the (k+1)-by-k  
            # |    lower bidiagonal matrix B.
            # | 2. Compute the SVD of the k-by-k upper bidiagonal matrix 
            # |    R = P*S*Q^T. The SVD of B is then (M*P)*S*Q^T.
            # %--------------------------------------------------------------

            # call dbdqr(jobu,j,work(ib1),work(ib1+lanmax),
            #            work(ibnd+j-1), work(ibnd+j),work(ip),lanmax+1)

            # print "before dbdqr2()"

            c1, c2 = dbdqr(JOBU, j, bbb1_a, bbb1_b, workp, lanmax + 1)
            work1[j - 1] = c1
            work1[j] = c2
            # print "after dbdqr2(), c1 = {}, c2 = {}".format(c1, c2)

            workp = workp.reshape(lanmax + 1, -1)

            # transpose to account for Fortran row-major order
            workp = workp.transpose().flatten()

            # for i, value in enumerate(workp.flatten()):
            #     print "workp[{}] = {}".format(252 + i, value)


            # if (lsame(jobu,'y')) then
            #    ncc = j+1
            # else
            #    ncc = 0
            # endif
            # if (lsame(jobv,'y')) then
            #    ncvt = j
            # else
            #    ncvt = 0
            # endif

            ncc =  (j + 1 if (JOBU == 'y') else 0)
            ncvt = (j     if (JOBV == 'y') else 0)

            # do i=1,n
            #    work(iq+(i-1)*(lanmax+1)) = one
            # enddo
            # print workq.shape
            for i in range(n - 1):
                # print i, i * (lanmax + 1)
                workq[i * (lanmax + 1)] = 1.0
            # pdb.set_trace()

            # THIS SHOULD POSSIBLY BE REPLACED BY A CALL TO THE FAST
            # DIVIDE-AND-CONQUER BIDIAGONAL SVD IN LAPACK 3 OR BY TRANSFORMING B
            # TO TRIDIAGONAL GOLUB-KAHAN FORM AND USING DHILLONS "HOLY GRAIL"
            # CODE.

            # %-----------------------------------------%
            # | R = P * S * Q^T, M^T <- P^T * M^T
            # %-----------------------------------------%
            # call dbdsqr('u',j,ncvt,0,ncc,work(ib1),work(ib1+lanmax),
            #      work(iq),lanmax, work,1, work(ip),lanmax+1,
            #      work(iwrk),info)

            # print "before dbdsqr2"
            # print "ncvt = {}, ncc = {}".format(ncvt, ncc)

            # for i, value in enumerate(workp):
            #     print "workp[{}] = {}".format(252 + i, value)

            # rbt = Really Big Tuple
            rbt = dbdsqr.dbdsqr('u', j, ncvt, 0, ncc, bbb1_a, bbb1_b, 
                                workq, lanmax, work1, 1, workp, lanmax + 1)

            info, bbb1_a, bbb1_b, _, work1, workp = rbt

            # print "after dbdsqr2"

            workp = np.array(workp)

            # for i, value in enumerate(workp):
            #     print "workp[{}] = {}".format(252 + i, value)

            # if (lsame(jobu,'y')) then
            #     %-----------------------------------------%
            #     | Form left Ritz-vectors
            #     | U = U * M * P
            #     %-----------------------------------------%            
            #   call zgemmina(m,j+1,j+1,U,ldu,
            #     work(ip),lanmax+1,zwork(iwrk),lwrk) 
            if JOBU == 'y':
                workp = workp.reshape(lanmax + 1, -1)

                # The Fortran code describes zgemmina() as A=AxB', which in 
                # this case is uuu * workp'. However, I found that I do not
                # need to transpose workp to get the correct result. I suspect
                # this is one bug cancelling out another. I *should* 
                # transpose workp, but since workp was generated in the call
                # to dbdsqr2() above, Fortran has already transposed it.
                # If I replace dbdsqr2() with a pure Python version of it,
                # I'll need to transpose it here as expected.

                # I have three different ways of doing what zgemmina() does.
                # The first is the fastest and most straightforward and 
                # produces good results.
                # The second exactly matches what the Fortran code does but
                # it's probably slower than the first method and definitely 
                # more confusing.
                # The third just calls the Fortran code
                zgemmina_option = 1

                if zgemmina_option == 1:
                    # Use ordinary matrix multiplication.
                    uuu = np.array(np.matrix(uuu) * np.matrix(workp))
                elif zgemmina_option == 2:
                    # Use partial matrix multiplication. This matches what
                    # the Fortran code does.
                    uuu_rows = uuu.shape[1]
                    uuu_copy = np.matrix(uuu[:, :j + 1].copy())
                    workp_copy = workp.copy()
                    workp_copy = np.matrix(workp_copy[:j + 1, :j + 1])
                    uuu_copy = uuu_copy * workp_copy
                    # uuu_copy is a truncated version of uuu. Here I append the
                    # columns of uuu that are not present in uuu_copy.
                    delta = uuu_rows - uuu_copy.shape[1]

                    uuu = np.hstack( (uuu_copy, uuu[:, -delta:]) )
                    uuu = np.array(uuu)
                elif zgemmina_option == 3:
                    # Call the Fortran code.
                    uuu, _ = zgemmina.zgemmina(m, j + 1, j + 1, uuu, workp, 
                                               lanmax + 1)


                # print "uuu[:1] --"
                # print uuu[:1]
                # print "uuu_copy[:1] --"
                # print uuu_copy[:1]
                # print "uuu[1:2] --"
                # print uuu[1:2]
                # print "uuu_copy[1:2] --"
                # print uuu_copy[1:2]

                # assert(np.allclose(uuu_copy, uuu[:m, :j + 1]))


                # for i in range(20):
                #     for j in range(20):
                #         print "uuu[{}, {}] = {}".format(i+1, j+1, uuu[i,j])

            # FIXME - didn't port case where jobv ='Y' since it is always
            # hardcoded to 'n' in my case


    return uuu, sigma, info, neig


#       subroutine dbdqr(jobq, n, D, E, c1, c2, Qt, ldq)
#       implicit none
def dbdqr(jobq, n, d, e, qt, ldq):
    jobq = jobq.upper()

    qt = qt.reshape(ldq, -1)
    #print "qt.shape = {}".format(qt.shape)

# c Compute QR factorization B = Q*R of (n+1) x n lower bidiagonal matrix 
# c with diagonal elements d(1)...d(n) and first subdiagonal elements
# c e(1)...e(n). On return [0 ... 0 c1 c2]' = Q'*[0 ... 0 1]'.
# c
# c If jobq=='Y' then on return Qt contains Q^T.

#       character*1 jobq
#       integer n,ldq
#       double precision D(*),E(*),c1,c2,Qt(ldq,*)
      
#       integer i,j
#       double precision cs,sn,r
#       logical lsame
#       external lsame

#       if (n.lt.1) return
    # import pdb
    # pdb.set_trace()

    # PS - tested for when jobq = 'N', not tested for jobq = 'Y'

    # print "inside dbdqr, n = {}".format(n)
    # for i, value in enumerate(d):
    #     print "d[{}] = {}".format(i + 1, d[i])
    # for i, value in enumerate(d):
    #     print "e[{}] = {}".format(i + 1, e[i])

    msg = "{} qt[{}, {}] = {}"
    pqt = lambda letter, i1, i2, value: msg.format(letter, i1 + 1, i2 + 1, 
                                                   value)

    if n >= 1:
#       if (lsame(jobq,'Y')) then
        if jobq == 'Y':
#           do j=1,n+1
#             do i=1,n+1
#                Qt(i,j) = 0.0
#             enddo
#           enddo
            qt[:n + 1, : n + 1] *= 0.0

#           do j=1,n+1
#             Qt(j,j) = 1.0
#           enddo
            for j in range(n + 1):
                # print "setting qt[{}, {}] = 1.0".format(j, j)
                qt[j, j] = 1.0
                # print pqt('O', j, j, 1.0)

#       endif

#       do i=1,n-1
        for i in range(n - 1):
#           call dlartg(d(i),e(i),cs,sn,r)
            # print "d[{}] = {}".format(i + 1, d[i])
            # print "e[{}] = {}".format(i + 1, e[i])
            cs, sn, r = dlartg.dlartg(d[i], e[i])
            #print "dlartg i = {}, cs = {}, sn = {}".format(i + 1, cs, sn)
            # print "cs = {}, sn = {}, r = {}".format(cs, sn, r)

            # cs = 0.0
            # sn = 0.0
            # r = 0.0
            # cs, sn, r = drotg(d[i], e[i])

#           d(i) = r
#           e(i) = sn*d(i+1)
#           d(i+1) = cs*d(i+1)
            d[i] = r
            e[i] = sn * d.item(i + 1)
            d[i + 1] = cs * d.item(i + 1)

#           if (lsame(jobq,'Y')) then
            if jobq == 'Y':
#               do j=1,i
#                   Qt(i+1,j) = -sn*Qt(i,j)
#                   Qt(i,j) = cs*Qt(i,j)
#               enddo
                # PS - This is the naive Python equivalent of the Fortran --
                # for j in range(i + 1):
                #     qt[i + 1, j] = -sn * qt.item(i, j)
                #     qt[i    , j] =  cs * qt.item(i, j)

                # PS - Here's the fast version --
                qt[i + 1, :i + 1] = -sn * qt[i, :i + 1]
                qt[i    , :i + 1] =  cs * qt[i, :i + 1]


#               Qt(i,i+1) = sn
#               Qt(i+1,i+1) = cs
#                print pqt('B', i, i + 1, sn)
                qt[i    , i + 1] = sn
                qt[i + 1, i + 1] = cs
#          endif
#       enddo

        # In fortran, a loop ends when the counter is one past the limit. For
        # example, at the end of the loop  `do i = 1,50`, i == 51. In Python,
        # this is not the case so I have to increment i here.
        i += 1

#       call dlartg(d(n),e(n),cs,sn,r)
        cs, sn, r = dlartg.dlartg(d[n - 1], e[n - 1])
#        print "i = {}, cs = {}, sn = {}".format(i + 1, cs, sn)

#       d(n) = r
#       e(n) = 0.0
#       c1 = sn
#       c2 = cs

        d[n - 1] = r
        e[n - 1] = 0.0
        c1 = sn
        c2 = cs

#       if (lsame(jobq,'Y')) then
        if jobq == 'Y':
#           do j=1,i
#               Qt(i+1,j) = -sn*Qt(i,j)
#               Qt(i,j) = cs*Qt(i,j)
#           enddo
            # PS - This is the naive Python equivalent of the Fortran --
            # for j in range(i + 1):
            #     qt[i + 1, j] = -sn * qt.item(i, j)
            #     qt[i    , j] =  cs * qt.item(i, j)

            # PS - Here's the fast version --
            qt[i + 1, :i + 1] = -sn * qt[i, :i + 1]
            qt[i    , :i + 1] =  cs * qt[i, :i + 1]

#           Qt(i,i+1) = sn
#           Qt(i+1,i+1) = cs
            qt[i    , i + 1] = sn
#            print pqt('D', i, i + 1, sn)
            qt[i + 1, i + 1] = cs

    # print "exiting dbdqr"
    # for i, value in enumerate(d):
    #     print "d[{}] = {}".format(i + 1, d[i])
    # for i, value in enumerate(d):
    #     print "e[{}] = {}".format(i + 1, e[i])

#   endif

    return c1, c2

#   end

def refinebounds(n, theta, bound, tol):
#       subroutine refinebounds(n,theta,bound,tol,eps34)
#       implicit none
#       integer n
#       double precision theta(*), bound(*), tol,eps34,gap
#       integer i,l
#       double precision dlapy2
#       external dlapy2
# c
# c     Refine Lanczos error bounds using the gap theorem.
# c     
#       if (n.le.1) return

    # print "len(theta) = {}".format(len(theta))
    # print "len(bound) = {}".format(len(bound))

    # print "inside refinebounds, n = {}".format(n)

    # for i, value in enumerate(theta):
    #     print "theta[{}] = {}".format(i + 1, theta[i])


    if n > 1:
#       do i=1,n
        for i in range(n):
#          do l=-1,1,2
            for j in (-1, 1):
#             if ((l.eq.1.and.i.lt.n) .or. (l.eq.-1.and.i.gt.1)) then
                if ((j == 1) and (i < n)) or ((j == -1) and (i > 1)):
#                if (abs(theta(i)-theta(i+l)) .lt. eps34*(theta(i))) then
#                    print "i = {}, j = {}, i + j = {}".format(i,j,i+j)

                    if abs(theta[i] - theta[i + j]) < (EPS34 * theta[i]):
#                   if (bound(i).gt.tol .and. bound(i+l).gt.tol) then
                        if (bound[i] > tol) and (bound[i + j] > tol):
#                      bound(i+l) = dlapy2(bound(i),bound(i+l))
                            bound[i + j] = util.dlapy2(bound[i], bound[i + j])
#                      bound(i) = 0.0
                            bound[i] = 0.0
#                   endif
#                endif
#             endif
#          enddo
#          gap = theta(i+1)-bound(i+1)-(theta(i)+bound(i))
            # print "theta[{}] = {}".format(i + 1 + 1, theta[i + 1])
            # print "bound[{}] = {}".format(i + 1 + 1, bound[i + 1])
            # print "theta[{}] = {}".format(i + 1, theta[i])
            # print "bound[{}] = {}".format(i + 1, bound[i])

            gap = theta[i + 1] - bound[i + 1] - (theta[i] + bound[i])
            # print "i = {}, gap = {}".format(i + 1, gap)
#          if (gap.gt.bound(i)) then
#             bound(i) = bound(i) * (bound(i)/gap)
#          endif
            if gap > bound[i]:
                # print "dong!"
                bound[i] *= (bound[i] / gap)

#       enddo
#       end


