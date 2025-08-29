# Python modules
from __future__ import division
import ctypes
import math
import os
import pprint
import re
import pdb
pp = pprint.pprint

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
import zgemmina

KMAX = 50
MAX_SINGULAR_VALUES = KMAX
NDPMAX = 8192

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

#     # Note - nsv_sought == kuser going in. zlansvdw() changes this to
#     # the number of singular values it found (kfit).
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


def zlansvdw_pure(n_data_points, lrow, mcol, kuser, lambda_,    
                  trlambda, uuu, tol=16e-16):

    print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ zlansvdw() "
# c     %---------------------------------%
# c     | Set machine dependent constants |
# c     %---------------------------------%
#       eps = dlamch('e')
#       eps34 = eps**(3.0/4.0)
#       epsn = dble(max(m,n))*eps/2.0
#       epsn2 = sqrt(dble(max(m,n)))*eps/2.0
#       sfmin = dlamch('s')

    jobu = 'y'
    jobv = 'n'

    ndp = n_data_points
    m = lrow
    n = mcol
    print m * (KMAX + 1)

    if uuu is None:
        uuu = np.zeros( (m, KMAX + 1), np.complex128)

    vvv = np.zeros( (n, KMAX), np.complex128)

    eps = np.finfo(np.float64).eps
    eps34 = eps**(3.0/4.0)
    epsn = max(m,n)*eps/2.0
    epsn2 = math.sqrt(epsn)
    # FIXME - would be nice to replace this with np.finfo()
    sfmin = dlamch('s')

    # print "eps = {}".format(eps)
    # print "eps34 = {}".format(eps34)
    # print "epsn = {}".format(epsn)
    # print "epsn2 = {}".format(epsn2)
    # print "sfmin = {}".format(sfmin)

# c     %--------------------------------%
# c     | Guard against absurd arguments |
# c     %--------------------------------%
#       lanmax = min(n+1,m+1,kmax)
#       tol = min(one,max(16.0*eps,tolin))
#       anorm = zero
    lanmax = min(n+1,m+1,KMAX)
    tol = min(1.0, max(16.0*eps, tol))
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

    ibnd = 1
    ib = ibnd + lanmax + 1
    ib1 = ib + (2 * lanmax)
    ip = ib1 + (2 * lanmax)
    iq = ip + ((lanmax + 1) ** 2)
    iwrk = iq + (lanmax ** 2)

    # These lengths are adjusted by 1 to reflect the difference between
    # Fortran's 1-based arrays and Python's 0-based arrays.
    len_work1  = lanmax + 1
    len_bbb = 2 * lanmax
    len_bbb1 = len_bbb
    len_workp  = (lanmax + 1) ** 2
    # Note that before zlansdvw.f's second call to dbdsqr(), it writes a bunch
    # of 1s into work starting at index iq. It goes waaaay beyond lanmax**2
    # and into the work(iwrk) section. Intentional? I don't know. The size I
    # give to workq here should be about right, erring on the side of safety.
    len_workq  = n * (lanmax + 1)

    msg = "lanmax = {}, ib = {}, ib1 = {}, ip = {}, iq = {}, iwrk = {}"
    print msg.format(lanmax, ib, ib1, ip, iq, iwrk)

 # lanmax =           50
 # ib =           52
 # ib1 =          152
 # ip =          252
 # iq =         2853
 # iwrk =         5353
# C     PS - changed from lwork to lzwrk
#       lwrk = lzwrk-iwrk+1
#       call dzero(7*lanmax + 2 + 2*lanmax**2,work,1)
#       call zzero(7*lanmax + 2 + 2*lanmax**2,zwork,1)

    # Fortran      Fortran indices     Python
    # -------      ---------------     ------
    # work(ibnd) = work(   1 -   51) = work1
    # work(ib)   = work(  52 -  151) = workb
    # work(ib1)  = work( 153 -  251) = workb1
    # work(ip)   = work( 252 - 2852) = workp
    # work(iq)   = work(2853 - 5352) = workq

    s = """
Fortran      Fortran indices     Python   Length
-------      ---------------     ------   ------
work(ibnd) = work(   1 -   51) = work1        50
work(ib)   = work(  52 -  151) = workb        99
work(ib1)  = work( 153 -  251) = workb1       98
work(ip)   = work( 252 - 2852) = workp      2600
work(iq)   = work(2853 - 5352) = workq      2499
"""
    print s

    work1  = np.zeros( (len_work1, ), dtype=np.float64)         # work(ibnd)
    bbb  = np.zeros( (len_bbb, ), dtype=np.float64)         # work(ib)
    bbb1 = np.zeros( (len_bbb1, ), dtype=np.float64)         # work(ib1)
    workp  = np.zeros( (len_workp, ), dtype=np.float64)        # work(ip)
    workq  = np.zeros( (len_workq, ), dtype=np.float64)        # work(iq)

    bbb.shape = (lanmax, 2)
    bbb1.shape = (lanmax, 2)


    # For some reason, zwork is declared to be of a certain size, but then
    # when functions are called only zwork(iwrk) is passed. The lower chunk 
    # of the array never gets used, so what's the point?
    # Here I make it big, probably much bigger than it needs to be.
    # FIXME trim this down to size later
    zwork = np.zeros( (m+n+32*m+7*KMAX+2+2*KMAX**2+5*ndp, ), np.complex128)



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

#           print *, "calling zgetu0w()"
#           call zgetu0w('n',ndp, m,n,0,1,U,rnorm,U,ldu,
#      c        ierr,ioption(1),anorm,zwork(iwrk),lambda,
#      c        trlambda,planF,planB)

        ioption = (1, 1)
        doption = (1e-12, 1e-14, 0.0)
        # uuu_copy = uuu.copy()
        rnorm, anorm, ierr = zget0w.zgetu0w('n', ndp, m, n, 0, 1, uuu[:,0], 
                                            uuu,
                                            ioption[0], lambda_, trlambda)
        print "done zgetu0w(), rnorm = {}".format(rnorm)

        # delta = uuu_copy - uuu
        # import pdb
        # pdb.set_trace()
        # for i in range(m):
        #     for j in range(KMAX + 1):
        #         if delta[i, j]:
        #             print "delta[{}, {}] = {}".format(i, j, delta[i, j])
#           print *, "done zgetu0w()"
#       endif

    # info = 0
    # neig = 0
    # jold = 0
    # j = min(k+max(8,k)+1,lanmax)
    info = 0
    neig = 0
    jold = 0
    j = min(kuser + max(8, kuser) + 1, lanmax)

# c     %------------------------------%
# c     | Iterate until convergence... |
# c     %------------------------------%
#   do while (neig.lt.k)
    libhlsvd = util.load_the_library()

    f = libhlsvd.zlanbprow_python_

    UuuArrayType = ctypes.c_double * (m * (KMAX + 1))
    uuu_r = UuuArrayType()
    uuu_i = UuuArrayType()

    VvvArrayType = ctypes.c_double * (n * KMAX)
    vvv_r = VvvArrayType()
    vvv_i = VvvArrayType()

    BbbArrayType = ctypes.c_double * (len_bbb)
    bbb_fortran = BbbArrayType()
    # bbb1_fortran = BbbArrayType()

    # Work1ArrayType = ctypes.c_double * (len_work1)
    # work1_fortran = Work1ArrayType()

    # WorkPArrayType = ctypes.c_double * (len_workp)
    # workp_fortran = WorkPArrayType()

    # WorkQArrayType = ctypes.c_double * (len_workq)
    # workq_fortran = WorkQArrayType()

    DoptionArrayType = ctypes.c_double * 3
    doption = DoptionArrayType()

    IoptionArrayType = ctypes.c_long * 2
    ioption = IoptionArrayType()

    ioption[0] = 1
    ioption[1] = 1
    doption[0] = 1e-12
    doption[1] = 1e-14
    doption[2] = 0.0


    LambdaDoubleArrayType = ctypes.c_double * n_data_points

    # Input params
    lambda_r   = LambdaDoubleArrayType()
    lambda_i   = LambdaDoubleArrayType()
    trlambda_r = LambdaDoubleArrayType()
    trlambda_i = LambdaDoubleArrayType()




    while neig < kuser:
# c     %---------------------------------------------------%
# c     | Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j} |
# c     %---------------------------------------------------%
        print "looptop, neig = {}, kuser = {}".format(neig, kuser)

        for i, z in enumerate(lambda_):
            lambda_r[i] = z.real
            lambda_i[i] = z.imag
        for i, z in enumerate(trlambda):
            trlambda_r[i] = z.real
            trlambda_i[i] = z.imag

        ierr = 0

        # for i in range(m):
        #     for j in range(KMAX + 1):
        #         print "uuu({}, {}) = {:.17G}".format(i + 1, j + 1, uuu[i, j])

        # for i in range(n):
        #     for j in range(KMAX):
        #         print "vvv({}, {}) = {:.17G}".format(i + 1, j + 1, vvv[i, j])

        # I call np.transpose() here to convert the 2D arrays into row major
        # order to make Fortran happy.

        uuu_fortran = np.transpose(uuu).flatten()

        for i, z in enumerate(uuu_fortran):
            # print "flat_uuu[{}] = {:.17E}".format(i + 1, z)
            uuu_r[i] = z.real
            uuu_i[i] = z.imag

        vvv_fortran = np.transpose(vvv).flatten()

        for i, z in enumerate(vvv_fortran):
            #print "flat_vvv[{}] = {:.17E}".format(i + 1, z)
            vvv_r[i] = z.real
            vvv_i[i] = z.imag


        # for i, value in enumerate(work1):
        #     work1_fortran[i] = value

        for i, value in enumerate(np.transpose(bbb).flatten()):
            # The call to float() here simply converts numpy floats to regular
            # Python floats.
            bbb_fortran[i] = float(value)
            #print "bbb_fortran[{}] = {}".format(i + 1, bbb_fortran[i])

        # for i, value in enumerate(np.transpose(bbb1).flatten()):
        #     # The call to float() here simply converts numpy floats to regular
        #     # Python floats.
        #     bbb1_fortran[i] = float(value)

        # for i, value in enumerate(workp):
        #     workp_fortran[i] = value

        # for i, value in enumerate(workq):
        #     workq_fortran[i] = value

        ldu = m
        ldv = n
        ldb = lanmax

        n_data_points = ctypes.c_long(ndp)        # ndp
        m = ctypes.c_long(m)
        n = ctypes.c_long(n)
        ierr = ctypes.c_long(ierr)
        jold = ctypes.c_long(jold)
        j = ctypes.c_long(j)
        lanmax = ctypes.c_long(lanmax)
        ldu = ctypes.c_long(ldu)
        ldv = ctypes.c_long(ldv)
        #ldb = ctypes.c_long(ldb)
        rnorm = ctypes.c_double(rnorm)
        #lwork = ctypes.c_long(len(work))

        # len_work1 = ctypes.c_long(len_work1)
        len_bbb = ctypes.c_long(len_bbb)
        # len_bbb1 = ctypes.c_long(len_bbb1)
        # len_workp = ctypes.c_long(len_workp)
        # len_workq = ctypes.c_long(len_workq)

        # for ips, value in enumerate(work1):
        #     kps = ips
        #     work[kps] = value
        # #     print "work[{}] = {}".format(kps, value)

        # # print "$$$"

        # for ips, value in enumerate(workb):
        #     kps = ips + lanmax + 1
        #     work[kps] = value
        # #     print "work[{}] = {}".format(kps, value)

        # # print "$$$"

        # for ips, value in enumerate(workb1):
        #     kps = ips + (2 * lanmax) + lanmax + 1
        #     work[kps] = value
        # #     print "work[{}] = {}".format(kps, value)

        # # print "$$$"

        # for ips, value in enumerate(workp):
        #     kps = ips + ((lanmax + 1) ** 2) + (2 * lanmax) + lanmax + 1
        #     work[kps] = value
        # #     print "work[{}] = {}".format(kps, value)

        # # print "$$$"

        # for ips, value in enumerate(workq):
        #     kps = ips + (lanmax ** 2) + ((lanmax + 1) ** 2) + (2 * lanmax) + lanmax + 1
        #     work[kps] = value
        # #     print "work[{}] = {}".format(kps, value)

        # # print "$$$"

        libhlsvd.zlanbprow_python_(ctypes.pointer(n_data_points),
                                  ctypes.pointer(m), 
                                  ctypes.pointer(n), 
                                  ctypes.pointer(jold),
                                  ctypes.pointer(j),
                                  ctypes.pointer(ierr),
                                  ctypes.pointer(uuu_r),
                                  ctypes.pointer(uuu_i),
                                  ctypes.pointer(ldu), 
                                  ctypes.pointer(vvv_r),
                                  ctypes.pointer(vvv_i),
                                  ctypes.pointer(ldv), 
                                  # ctypes.pointer(work1_fortran),
                                  # ctypes.pointer(len_work1),
                                  # bbb and bbb1 are declared flat but are 
                                  # 2-D in Fortran with the 2nd dim hardcoded
                                  # at 2. I only need to pass the length of the
                                  # leading dimension which is lanmax.
                                  ctypes.pointer(bbb_fortran),
                                  ctypes.pointer(lanmax),
                                  # ctypes.pointer(bbb1_fortran),
                                  # ctypes.pointer(workp_fortran),
                                  # ctypes.pointer(len_workp),
                                  # ctypes.pointer(workq_fortran),
                                  # ctypes.pointer(len_workq),
                                  ctypes.pointer(rnorm),
                                  ctypes.pointer(doption),
                                  ctypes.pointer(ioption),
                                  ctypes.pointer(lambda_r),
                                  ctypes.pointer(lambda_i),
                                  ctypes.pointer(trlambda_r),
                                  ctypes.pointer(trlambda_i),
                                 )

        print "Python returning from zlanbprow_python_()"

        # len_work1 = len_work1.value
        len_bbb = len_bbb.value
        # len_bbb1 = len_bbb1.value
        # len_workp = len_workp.value
        # len_workq = len_workq.value

        jold = jold.value
        j = j.value

        # pdb.set_trace()


        rnorm = rnorm.value
        m = m.value
        n = n.value
        ierr = ierr.value
        lanmax = lanmax.value

        uuu = [complex(r, i) for r, i in zip(uuu_r, uuu_i)]

        # for i, z in enumerate(uuu):
        #     print "uuu({}) = {:.17E}".format(i + 1, z)

        uuu = util.flat_fortran_to_2d_python(uuu, m, KMAX + 1)

        #print uuu

        # uuu_copy = [None] * len(uuu)

        # # print "copying..."
        # for i, z in enumerate(uuu):
        #     k = (i // m) + ((i % m) * (KMAX + 1))
        #     uuu_copy[k] = z

        # # for i, z in enumerate(uuu_copy):
        # #     print "uuu_copy({}) = {:.17E}".format(i + 1, z)

        # uuu = np.array(uuu_copy)
        # uuu.shape = (m, KMAX + 1)

        # for i in range(m):
        #     for j in range(KMAX + 1):
        #         print "uuu({}, {}) = {:.17E}".format(i + 1, j + 1, uuu[i, j])

        vvv = [complex(r, i) for r, i in zip(vvv_r, vvv_i)]

        vvv = util.flat_fortran_to_2d_python(vvv, n, KMAX)

        # for i in range(n):
        #     for j in range(KMAX):
        #         print "vvv({}, {}) = {:.17E}".format(i + 1, j + 1, vvv[i, j])

        bbb = util.flat_fortran_to_2d_python(bbb_fortran, lanmax, 2)

        # bbb1 = util.flat_fortran_to_2d_python(bbb1_fortran, lanmax, 2)

        # for i in range(lanmax):
        #     for j in range(2):
        #         print "bbb({}, {}) = {}".format(i + 1, j + 1, bbb[i, j])

        # work1 = np.array([value for value in work1])
        # workp = np.array([value for value in workp])
        # workq = np.array([value for value in workq])


        # work = [value for value in work]

        # start = 0
        # end = lanmax + 1
        # work1  = np.array(work[start:end])
        # start = end + 1
        # end += 2 * lanmax
        # workb  = np.array(work[start:end])
        # start = end + 1
        # end += 2 * lanmax
        # workb1 = np.array(work[start:end])
        # start = end + 1
        # end += (lanmax + 1) ** 2
        # workp  = np.array(work[start:end])
        # start = end + 1
        # end += lanmax ** 2
        # workq  = np.array(work[start:end])

        # import pdb
        # pdb.set_trace()

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {:.17E}".format(i + 1, work1[i])

        # print "^^^"

        # for i, value in enumerate(bbb):
        #     print "bbb[{}] = {:.17E}".format(i + 1, bbb.item(i))

        # print "^^^"
        
        # for i, value in enumerate(bbb1):
        #     print "workb1[{}] = {:.17E}".format(i + 1, bbb1.item(i))

        # print "^^^"
        
        # for i, value in enumerate(workp):
        #     print "workp[{}] = {:.17E}".format(i + 1, workp[i])

        # print "^^^"

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {:.17E}".format(i + 1, work1[i])

        # print "^^^"



     #     call zlanbprow(ndp, m, n, jold, j,  U, ldu, V, ldv,
     # c        work(ib),lanmax,rnorm,doption(1),ioption(1),
     # c        work(iwrk), iwork,ierr,  
     # c        zwork(iwrk), lambda,trlambda,planF,planB)
     #     jold = j

        jold = j

# c     %---------------------------------------------%
# c     | Compute and analyze SVD(B) and error bounds |
# c     %---------------------------------------------%

        
        # call dcopy(2*lanmax, work(ib),1,work(ib1),1)
        # for i in range(2 * lanmax):
        #     work[ib1 + i] = work[ib + i]
        bbb1 = bbb.copy()

        # call dzero(j+1,work(ibnd),1)

        for i in range(j + 1):
            work1[i] = 0.0

        #      call second(t2)
        #      call dbdqr('N',j,work(ib1),work(ib1+lanmax),work(ibnd+j-1),
        # c         work(ibnd+j),work(ip),lanmax+1)
        print "j =  {}".format(j)

        # for i, value in enumerate(workb):
        #     print "workb[{}] = {:.17E}".format(i + 1, workb[i])

        # print "^^^"
        
        # for i, value in enumerate(workb1):
        #     print "workb1[{}] = {:.17E}".format(i + 1, workb1[i])

        # print "^^^"
        
        # for i, value in enumerate(workp):
        #     print "workp[{}] = {:.17E}".format(i + 1, workp[i])

        # print "^^^"

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {:.17E}".format(i + 1, work1[i])

        # print "^^^"

        # import pdb
        # pdb.set_trace()

        # bbb1 (which starts life as an ordinary copy of bbb) is always 
        # manipulated as two distinct parts. 
        bbb1_a = bbb1.flatten()[::2]
        bbb1_b = bbb1.flatten()[1::2]


        # for value in (bbb1.flatten()[::2] + bbb1.flatten()[1::2]):
        #     print value

        print "before dbdqr1(), j = {}".format(j)
        c1, c2 = dbdqr('N', j, bbb1_a, bbb1_b, workp, lanmax+1)

        work1[j - 1] = c1
        work1[j] = c2

        #bbb1 = np.array(zip(aaa, ccc))

        print "after dbdqr1(), c1 = {}, c2 = {}".format(c1, c2)

        # for value in (bbb1.flatten()[::2] + bbb1.flatten()[1::2]):
        #     print value

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {:.17E}".format(i + 1, work1[i])

        # print "^^^"

        # for i in range(lanmax):
        #     for k in (0, 1):
        #         print "bbb1[{}, {}] = {:.17E}".format(i + 1, k + 1, bbb1[i, k])

        # print "^^^"

        # for i, value in enumerate(workp):
        #     print "workp[{}] = {:.17E}".format(i + 1, workp[i])

        # print "^^^"

        print "before dbdsqr, j = {}" .format(j)

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


        rbt = dbdsqr.dbdsqr('u', j,    0,   1,   0, bbb1_a, bbb1_b, None, 1, 
                            work1, 1, None, 1)

        info, bbb1_a, bbb1_b, _, work1, _ = rbt

        print "after dbdsqr" 
        
        # for i, value in enumerate(work1):
        #     print "work1[{}] = {}".format(i + 1, value)


        # FIXME this is a crappy hack to address the fact that code below 
        # (refinebounds) expects aaa and ccc to be larger than they are.
        # import pdb
        # pdb.set_trace()
        bbb1_a +=  [0.0] * (lanmax - len(bbb1_a) + 1)
        bbb1_b +=  [0.0] * (lanmax - len(bbb1_b) + 1)
        print "len(bbb1_b) = {}".format(len(bbb1_b))
        bbb1_a = np.array(bbb1_a)
        bbb1_b = np.array(bbb1_b)
        work1 += [0.0] * (len_work1 - len(work1))
        work1 = np.array(work1)

        #      call dbdsqr('u',j,0,1,0,work(ib1),work(ib1+lanmax),work,1,
        # c         work(ibnd),1,work,1,work(iwrk),info)

        # call  second(t3)
        # tbsvd = tbsvd + (t3-t2)
        # nbsvd = nbsvd + 1

        # if (j.gt.5) then
        #    anorm = work(ib1)
        # else
        #    anorm = max(anorm,work(ib1))
        # endif
        if j > 5:
            anorm = bbb1_a[0]
        else:
            anorm = max(anorm, bbb1_a[0])     

        print "anorm = {}".format(anorm)

         # do i=1,j
         #    work(ibnd+i-1) = abs(rnorm*work(ibnd+i-1))
         #    work(ib1+lanmax+i-1) = work(ib1+lanmax+i-1)**2
         # enddo

        work1 = work1 * rnorm
        work1 = np.absolute(work1)
        bbb1_b = bbb1_b**2

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

        refinebounds(j, bbb1_b, work1, epsn*anorm, eps34)

        print "after refinebounds"

        # for i, value in enumerate(work1):
        #     print "work1[{}] = {}".format(i + 1, value)

# c     %----------------------------------------------------%
# c     | Determine the number of converged singular values  |
# c     %----------------------------------------------------%
#          do i=1,min(j,k)
#             bnd(i) = work(ibnd+i-1)
#          enddo

        i = min(j, kuser)
        bounds = work1[:i]

        for i, value in enumerate(bounds):
            print "bounds[{}] = {}".format(i + 1, value)


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

    
        print "tol = {}".format(tol)
        while (i < min(j, kuser)):
            # print "work[{}] = {}".format(i + 1, work1[i])
            # print "tol * bbb1_a[{}] = {}".format(i + 1, tol * bbb1_a[i])
            #print "sigma delta = {}".format((tol * bbb1_a[i]) - work1[i])
            if work1[i] <= (tol * bbb1_a[i]):
                # neig += 1
                sigma.append(bbb1_a[i])
                i += 1
            else:
                # Force loop exit
                i = kuser

        for i, value in enumerate(sigma):
            print "sigma[{}] = {}".format(i + 1, value)

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
            print "bailing, ierr = {}".format(ierr)
            break

        if j >= lanmax:
            if neig < kuser:
                print "bailing"
                print "j = {}, lanmax = {}, neig = {}, kuser = {}".format(j, lanmax, neig, kuser)
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
        # print "monkey, k = {}".format(kuser)
        # print "monkey, neig = {}".format(neig)
        # print "monkey, form = {}".format(((kuser - neig) * (j - 6)) // ( 2 * neig + 1))

        if neig > 1:
            dj = min(j // 2, ((kuser - neig) * (j - 6)) // ( 2 * neig + 1))
            dj = min(100, max(2, dj))
        else:
            dj = j / 2
            dj = min(100, max(10, dj))

        j = min(j + dj, lanmax)

        #bbb1 = np.array(zip(aaa, ccc))
        print "loop end, j = {}, dj = {}".format(j, dj)
        print "loop end, neig = {}, k = {}".format(neig, kuser)


##############             while loop has ended      ##################

    print "Pre-calculate singular vectors"
    print "neig = {}, kuser = {}".format(neig, kuser)

#  50   if (neig.ge.k) then
# c     %-----------------------------------------%
# c     | Calculate singular vectors if requested %
# c     %-----------------------------------------%
    if neig >= kuser:
        print "inside Calculate singular vectors"
        # j = jold
        # k = neig

        j = jold
        kuser = neig

        # if (lsame(jobu,'y') .or. lsame(jobv,'y')) then
        #    print *, "about to dcopy"
        #    call dcopy(2*lanmax, work(ib),1,work(ib1),1)
        #    print *, "done dcopy"

        if (jobu == 'y') or (jobv == 'y'):
            bbb1 = bbb.copy()
            bbb1_a = bbb1.flatten()[::2]
            bbb1_b = bbb1.flatten()[1::2]

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

            print "before dbdqr2()"

            c1, c2 = dbdqr(jobu, j, bbb1_a, bbb1_b, workp, lanmax + 1)
            work1[j - 1] = c1
            work1[j] = c2
            print "after dbdqr2(), c1 = {}, c2 = {}".format(c1, c2)

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

            ncc =  (j + 1 if (jobu == 'y') else 0)
            ncvt = (j     if (jobv == 'y') else 0)

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

            print "before dbdsqr2"
            print "ncvt = {}, ncc = {}".format(ncvt, ncc)

            # for i, value in enumerate(workp):
            #     print "workp[{}] = {}".format(252 + i, value)

            rbt = dbdsqr.dbdsqr('u', j, ncvt, 0, ncc, bbb1_a, bbb1_b, 
                                workq, lanmax, work1, 1, workp, lanmax + 1)

            info, bbb1_a, bbb1_b, _, work1, workp = rbt

            print "after dbdsqr2"

            #pdb.set_trace()

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
            if jobu == 'y':
                #pdb.set_trace()

                workp = workp.reshape(lanmax + 1, -1)
                # trworkp = workp.transpose()
                # ppp = uuu.copy()
                # for i in range(lanmax + 1):
                #     for j in range(lanmax + 1):
                #         ppp[i, j] *= workp[i, j]

                # workp = workp.flatten()

                # The Fortran code describes zgemmina() as A=AxB', which in 
                # this case is uuu * workp'. However, I found that I do not
                # need to transpose workp to get the correct result. I suspect
                # this is one bug cancelling out another. I *should* 
                # transpose workp, but since workp was generated in the call
                # to dbdsqr2() above, Fortran has already transposed it.
                # If I replace dbdsqr2() with a pure Python version of it,
                # I'll need to transpose it here as expected.

                uuu_shape = uuu.shape
                uuu_rows = uuu.shape[1]
                uuu_copy = np.matrix(uuu[:, :j + 1].copy())
                #uuu_copy = np.matrix(uuu.copy())
                #workp_copy = workp.transpose()
                workp_copy = workp.copy()
                workp_copy = np.matrix(workp_copy[:j + 1, :j + 1])
                uuu_copy = uuu_copy * workp_copy
                #uuu_copy = uuu_copy.transpose()

                #uuu_copy = uuu_copy.transpose().reshape(m,j+1)

                #uuu = uuu_copy


                #uuu, _ = zgemmina.zgemmina(m, j + 1, j + 1, uuu, workp, lanmax + 1)


                # uuu = np.array(uuu)
                # uuu.shape = uuu_shape
                #uuu_copy = np.array(uuu_copy)

                # uuu_copy is a truncated version of uuu. Here I append the
                # columns of uuu that are not present in uuu_copy.
                delta = uuu_rows - uuu_copy.shape[1]

                uuu = np.hstack( (uuu_copy, uuu[:, -delta:]) )
                uuu = np.array(uuu)
                #uuu = np.append(uuu_copy, uuu[-delta:, :])

                # print "uuu[:1] --"
                # print uuu[:1]
                # print "uuu_copy[:1] --"
                # print uuu_copy[:1]
                # print "uuu[1:2] --"
                # print uuu[1:2]
                # print "uuu_copy[1:2] --"
                # print uuu_copy[1:2]

                # z = np.allclose(uuu_copy, uuu[:m, :j + 1])
                # assert(z)


                # for i in range(20):
                #     for j in range(20):
                #         print "uuu[{}, {}] = {}".format(i+1, j+1, uuu[i,j])

                # z = np.allclose(uuu,ppp)

                z = 42
            # FIXME - didn't port case where jobv ='Y' since it is always
            # hardcoded to 'n' in my case


    kuser = neig

    #pdb.set_trace()

    return uuu, sigma, info, kuser


#       subroutine dbdqr(jobq, n, D, E, c1, c2, Qt, ldq)
#       implicit none
def dbdqr(jobq, n, d, e, qt, ldq):
    jobq = jobq.upper()

    #pdb.set_trace()
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

    print "inside dbdqr, n = {}".format(n)
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
#          do j=1,n+1
#             do i=1,n+1
#                Qt(i,j) = 0.0
#             enddo
#          enddo
            for j in range(n + 1):
                for i in range(n + 1):
                    #print i, j
                    qt[i, j] = 0.0
                    # print pqt('Z', i, j, 0.0)

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
                # pdb.set_trace()
                #the_range = (range(i) if i else [0])
                for j in range(i + 1):
#                   Qt(i+1,j) = -sn*Qt(i,j)
#                   Qt(i,j) = cs*Qt(i,j)
#                    print pqt('A', i, j, qt.item(i, j))
                    qt[i + 1, j] = -sn * qt.item(i, j)
#                    print pqt('A', i + 1, j, qt[i + 1, j])
                    # if (i == 44) and (j == 0):
                    #     pdb.set_trace()
                    qt[i    , j] =  cs * qt.item(i, j)

#               enddo
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
            #the_range = (range(i) if i else [0])
            for j in range(i + 1):
            # for j in range(i):
#               Qt(i+1,j) = -sn*Qt(i,j)
#               Qt(i,j) = cs*Qt(i,j)
#                print pqt('C', i, j, qt.item(i, j))
                qt[i + 1, j] = -sn * qt.item(i, j)
                qt[i    , j] =  cs * qt.item(i, j)
#           enddo
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

#       endif

    return c1, c2

#       end

def refinebounds(n, theta, bound, tol, eps34):
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

    print "inside refinebounds, n = {}".format(n)

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

                    if abs(theta[i] - theta[i + j]) < (eps34 * theta[i]):
#                   if (bound(i).gt.tol .and. bound(i+l).gt.tol) then
                        if (bound[i] > tol) and (bound[i + j] > tol):
#                      bound(i+l) = dlapy2(bound(i),bound(i+l))
                            print "ding!"
                            bound[i + j] = math.sqrt(bound[i]**2 + bound[i + j] ** 2)
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
                print "dong!"
                bound[i] = bound[i] * (bound[i] / gap)

#       enddo
#       end


