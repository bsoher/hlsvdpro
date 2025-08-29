# Python modules
from __future__ import division
import ctypes
import math
import pdb

# 3rd party modules
import numpy as np

# Our modules
import zlansvdw

# Variable name translations 
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# kuser      nsv_sought         # of singular values requested by caller
# kfit       nsv_found          # of singular values found
# lsinval    singular_values    Array containing the singular values (floats)
# ndp        n_data_points      
# kmax       MAX_SINGULAR_VALUES  Hardcoded to 50.
# lrow       n_rows             # of rows in the very important uuu matrix.
#                               Once inside zlansvdw(), this is called 'm'.
# mcol       n_columns          # of rows? in the very important vvv matrix.
#                               Once inside zlansvdw(), this is called 'n'.
#                               Since it counts rows, I'm not sure why it was
#                               called 'mcol'.


def lanczopw(signals, n_rows, n_columns, nsv_sought):
# Call -- 
    # lanczopw(signals, len(signals), lrow, mcol, kuser, kmax, lsinval, u, v, 
    #          work, lwork, zwork, lzwork, zwork[ilambda], zwork[itrlambda]
    #          zwork[ifvect], ndiv)

    # Output: lsinval, uuu, ndiv (kfit == nsv_found)

# Signature --
    #   subroutine lanczopw(signal,ndp,m,n,k,kmax,sigma,U,V,
    #  c     work,lwrk,zwork,lzwrk,lambda,trlambda,fvect,info)
    #   implicit none

    #   double precision zero, one
    #   double complex zeroc
    #   parameter(zero = 0.0d0, one = 1.0d0, zeroc=(0.0d0,0.0d0))
    #   integer i,m,n,k,kmax,ioption(10),iwork(2*kmax+1),info,lwrk,lzwrk
    #   integer ind1,ndp
    #   integer*8 planF,planB
    #   double precision work(lwrk)
    #   double precision doption(10)
    #   double precision tol,pi,bnd(kmax)
    #   double precision dnrm2,sigma(kmax)
    #   complex*16 lambda(ndp),trlambda(ndp),
    #  c          U(m,kmax+1),V(n,kmax),zwork(lzwrk),
    #  c           fvect(ndp),signal(*)

    #libhlsvd = _load_the_library()

    # At this point libhlsvd should be a valid reference to the hlsvd library.
    # I test that assertion here as best as I can by generating a throwaway
    # reference to the function I'm going to call. 
    #f = libhlsvd.zlansvdw_python_


    # OK, all is well. I create and populate the variables I need as parameters
    # to the function.
    n_data_points = len(signals)

    # c
    # c  computation of the first column of the circulant matrix
    # c      

    # PS - There's a couple of odd loops here that set up fvect. I think the 'f'
    # might stand for 'fased' (phased)? These loops basically shove all the
    # points to the right by (ndp/2) + 1.
#   do i=1, ndp
#      fvect(i)=zeroc
#   end do
#   call zcopy(m,signal(n),1,fvect(1),1)

    fvect = np.roll(signals.copy(), (n_data_points // 2) + 1)

    # for i, z in enumerate(fvect):
    #     print "fvect[{}] = {:.17E}".format(i + 1, z)

    lambda_ = np.fft.fft(fvect)

#   call dfftw_execute_dft(planF, fvect, lambda)

    # for i, z in enumerate(lambda_):
    #     print "lambda[{}] = {}".format(i, z)

#   do i =1,ndp
#     fvect(i)=dconjg(fvect(i))
#   end do     
    fvect = np.conjugate(fvect)

    # PS - Reverse the elements of fvect, except for the first element which 
    # isn't changed.

#   ind1=ndp
#   do i=2,ndp
#      zwork(i)=fvect(ind1)
#      ind1=ind1-1
#   end do

#   do i=2,ndp
#      fvect(i)=zwork(i)
#  end do

    fvect = np.concatenate( (np.array([fvect[0]]), 
                             fvect[1:].copy()[::-1])
                          )
    # for i, z in enumerate(fvect):
    #     print "fvect[{}] = {:.17E}".format(i + 1, z)

    trlambda = np.fft.fft(fvect)

    # PS - lambda and trlambda never change from this point on.

    # call dfftw_execute_dft(planF, fvect, trlambda)

#   do i=1, ndp
#      lambda(i)=lambda(i)/dble(ndp)
#      trlambda(i)=trlambda(i)/dble(ndp) 
#   end do 

    lambda_  = lambda_  / n_data_points
    trlambda = trlambda / n_data_points

    # for i, z in enumerate(trlambda):
    #     print "trlambda[{}] = {}".format(i, z)

    uuu = None
    done = False
    while not done:
        uuu, singular_values, info, nsv_found = \
                zlansvdw.zlansvdw(n_data_points, n_rows, n_columns, nsv_sought, 
                                  lambda_, trlambda, uuu)
            
        if info == -1:
            if nsv_found > 0:
                # Need to call zlansvdw() again.
                nsv_sought = nsv_found
                # print "Python: calling zlansvdw again, nsv_sought = {}".format(nsv_sought)
                pass
            else:
                # FIXME - raise an error; PROPACK didn't converge.
                done = True
        else:
            done = True

    # print "lanczopw: done with zlansvdw, nsv_found = {}".format(nsv_found)
    # print "lanczopw: singular values"
    # pp(singular_values)


    # # Set up ctypes variables to call zlansvdw().
    # # In-line comments refer to the names of the corresponding variables
    # # in hlsvdpro.f.

    # InputDoubleArrayType = ctypes.c_double * n_data_points
    # OutputDoubleArrayType = ctypes.c_double * MAX_SINGULAR_VALUES
    # UuuArrayType = ctypes.c_double * (lrow * (KMAX + 1))

    # kmax = KMAX
    # ndpmax = NDPMAX

    # singular_values = OutputDoubleArrayType()               # sigma

    # # lwrk = lrow + mcol + (13*KMAX) + (8*KMAX**2) + (32 * lrow) + (n_data_points*KMAX)
    # # lwrk = ctypes.c_long(lwrk)

    # # lzwrk = lrow + mcol + (32*lrow) + (7*KMAX) + 2 + (2*KMAX**2) + (5 * n_data_points)
    # # lzwrk = ctypes.c_long(lzwrk)

    # # Input params
    # lambda_r   = InputDoubleArrayType()
    # lambda_i   = InputDoubleArrayType()
    # trlambda_r = InputDoubleArrayType()
    # trlambda_i = InputDoubleArrayType()
   
    # # import pdb
    # # pdb.set_trace()

    # #lambda_ = lambda_.tolist()
    # for i, z in enumerate(lambda_):
    #     lambda_r[i] = z.real   
    #     lambda_i[i] = z.imag

    # #trlambda = trlambda.tolist()
    # for i, z in enumerate(trlambda):
    #     trlambda_r[i] = z.real   
    #     trlambda_i[i] = z.imag

    # # zlansvdw() is the call where U/uuu gets populated, so it's first created
    # # here.
    # uuu_r = UuuArrayType()
    # uuu_i = UuuArrayType()

    # # for foo in uuu_r:
    # #     print foo

    # # kasdhfkljhsd

    # n_data_points = ctypes.c_long(n_data_points)        # ndp
    # lrow = ctypes.c_long(lrow)                          # Lrow/m
    # mcol = ctypes.c_long(mcol)                          # mcoL/n
    # kmax = ctypes.c_long(kmax)                          # kmax
    # ndpmax = ctypes.c_long(ndpmax)                      # ndpmax

    # info = 0

    # done = False
    # while not done:
    #     nsv_sought = ctypes.c_long(nsv_sought)              # kuser/k/kfit
    #     info = ctypes.c_long(info)                          # info (0, < 0, > 0)
    #                                                         # 0 is good!

    #     # Note - nsv_sought == kuser going in. zlansvdw() changes this to
    #     # the number of singular values it found (kfit).
    #     libhlsvd.zlansvdw_python_(ctypes.pointer(n_data_points),
    #                               ctypes.pointer(lrow), 
    #                               ctypes.pointer(mcol), 
    #                               ctypes.pointer(nsv_sought),
    #                               ctypes.pointer(singular_values),
    #                               ctypes.pointer(info),
    #                               ctypes.pointer(lambda_r),
    #                               ctypes.pointer(lambda_i),
    #                               ctypes.pointer(trlambda_r),
    #                               ctypes.pointer(trlambda_i),
    #                               ctypes.pointer(uuu_r),
    #                               ctypes.pointer(uuu_i),
    #                              )

    #     nsv_sought = nsv_sought.value
    #     info = info.value

    #     if (info == -1):
    #         if nsv_sought > 0:
    #             # Need to call zlansvdw() again.
    #             print "Python: calling zlansvdw again with kuser = {}".format(nsv_sought)
    #             pass
    #         else:
    #             # FIXME - raise an error; PROPACK didn't converge.
    #             done = True
    #     else:
    #         done = True

    # nsv_found = nsv_sought

    # #print "Python land: nsv_sought = {}".format(nsv_sought)

    # mcol = mcol.value
    # lrow = lrow.value

    # # Turn these into ordinary Python lists
    # uuu_r = [value for value in uuu_r]
    # uuu_i = [value for value in uuu_i]

    # # import vespa.common.util.math_ as util_math
    # # for j in range(lrow):
    # #     for i in range(KMAX + 1):
    # #         k = (i * lrow) + j
    # #         #if uuu_r[k] or uuu_i[k]:

    # #         real, imag = d[ (j + 1, i + 1) ]
    # #         try:
    # #             assert(util_math.eq(real, uuu_r[k]))
    # #             assert(util_math.eq(imag, uuu_i[k]))
    # #         except AssertionError:
    # #             print "AssertionError, k = {}, i = {}, j = {}".format(k, i ,j )
    # #             print real, imag, uuu_r[k], uuu_i[k]
    # #         except IndexError:
    # #             print "IndexError, k = {}, i = {}, j = {}".format(k, i ,j )
    # #         # print "U({},{}): ({:.17E},{:.17E})".format(j+1, i+1, uuu_r[k], uuu_i[k])    

    # uuu = [complex(uuu_r[i], uuu_i[i]) for i in range(lrow * (KMAX + 1))]

    # # for i, z in enumerate(uuu):
    # #     print "{}: {:.17E}".format(i, z)

    # # convert uuu into a numpy array and reorganize from Fortran's column-major
    # # order.
    # reorg = np.zeros( (lrow, KMAX + 1), np.complex128 )
    # for j in range(lrow):
    #     for i in range(KMAX + 1):
    #         k = (i * lrow) + j

    #         reorg[j][i] = uuu[k]

    # uuu = reorg

    # # for i in range(lrow):
    # #     for j in range(KMAX + 1):
    # #         print "U({}, {}) = {:.17E}".format(i, j, uuu[i][j])



    # # for i, z in enumerate(uuu):
    # #     # U(           1           1 ): ( 9.10583222219728122E-002,-6.27880386519005795E-002)
    # #     print "U({}): ({:.17E},{:.17E})".format(i, z.real, z.imag)
                               
    # # # I tease the returned variables into Python types before passing them
    # # # back to the caller. (Slicing a ctypes array returns a list.) The Fortran
    # # # code has already sorted them the way we like (largest singular value 
    # # # first).
    # # singular_values = singular_values[:nsv_found]
    # # frequencies     = frequencies[:nsv_found]
    # # damping_factors = damping_factors[:nsv_found]
    # # amplitudes      = amplitudes[:nsv_found]
    # # phases          = phases[:nsv_found]

    # # damping_factors = [1 / df for df in damping_factors]
    # # damping_factors = [df * dwell_time for df in damping_factors]

    # # frequencies = [frequency / dwell_time for frequency in frequencies]

    # # phases = [phase * constants.RADIANS_TO_DEGREES for phase in phases]

    # # Convert this from a ctypes thingy to regular Python, and trim.
    # singular_values = [x for x in singular_values][:nsv_found]

    return uuu, singular_values, nsv_found

