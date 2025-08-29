# Python modules
from __future__ import division
import math
import pdb
import pprint
pp=pprint.pprint

# 3rd party modules
import numpy as np
import scipy.linalg

# Our modules
import lanczopw

MAX_SINGULAR_VALUES = 50
# Max # of data points is set to 8192 in the Fortran. That's only necessary
# because of Fortran's lack (at the time) of dynamic array allocation. In 
# Python there's no such restriction.
#NDPMAX = 8192

# Complex LAPACK functions
zgelss, = scipy.linalg.lapack.get_lapack_funcs( ['gelss'], 
                                                np.array([1j])
                                              )

# Variable name translations 
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# kuser      nsv_sought         # of singular values requested by caller
# kfit       nsv_found          # of singular values found
# lsinval    singular_values    Array containing the singular values (floats)
# ndp        n_data_points      
# kmax       MAX_SINGULAR_VALUES  Hardcoded to 50.
# lrow       n_rows             # of rows in the critical work matrix 
# mcol       n_columns          # of columns in the critical work matrix 


#  subroutine hlsvdpw_python(signal_r, signal_i, ndp, Lrow, mcoL,
# +                          kuser, kfit, Lsinval, ampl, fase, damp, freq)
def hlsvdpro(signals, nsv_sought, step_size):

    # integer          ndpmax
    # parameter        (ndpmax=8192)
    # parameter        (kmax=50)
    # real*8           ampl(50),damp(50),freq(50),fase(50)
    # real*8           signal_r(ndpmax), signal_i(ndpmax), Lsinval(50)
    # complex*16       signal(ndpmax)
    # integer          kuser,kfit,Lrow,mcoL,ndp,i,kmax,lrwork,lzwork
    # parameter        (lrwork=20*ndpmax+13*kmax+8*kmax**2+ndpmax*kmax)
    # parameter        (lzwork=90*ndpmax+100*kmax+5*kmax**2+4*kmax*ndpmax+20)
    # real*8           rwork(lrwork)
    # complex*16       zwork(lzwork)

    # lsinval = output: Lanczos singular values of the Hankel matrix
    # kuser = nsv_sought
    # kfit = nsv_found
    # rwork = real (float) work array, lrwork = length of rwork
    # zwork = complex work array, lzwork = length of zwork
    # lrwork = 594090 and lzwork = 2393200
    # indU is the index of the first independent work array in zwork. The first
    #    indU elements of zwork are not dedicated to any array other than 
    #    zwork, as far as I can tell.
    # indU = 22*ndp + 7*kmax + 2*kmax**2 + 35 

    n_data_points = len(signals)
    n_columns = n_data_points // 2
    n_rows = n_columns + 1


    # PS Making a copy of the signal array is necessary in Fortran because 
    # zgelss() would overwrite signal otherwise. In Python's version of zgelss
    # we have the option not to overwrite the input array so copying signals
    # is not necessary.
    # do i = 1,ndp
    #    x(i) = signal(i)       ! x corrupted in amplitude computation
    # end do

# c     ******************************************************************
# c
# c     Lanczos SVD of hankel datamatrix.
# c     ================================
# c
# c     Lcwork = 2*kmax+ndp
#       Lcwork = 2*kmax+64*(ndp+kmax)
#       Lrwork = 5*kmax
#       ilambda = lzwrk+1
#       itrlambda = ilambda + ndp
#       ifvect = itrlambda + ndp

# 12      call lanczopw(signal,ndp,Lrow,mcoL,
#    c         kuser,kmax,Lsinval,U,V,
#    c     work,lwrk,zwork,lzwrk,zwork(ilambda),zwork(itrlambda),
#    c     zwork(ifvect),ndiv)


    # PS - lanczopw() is called in a loop until it returns a satisfactory 
    # answer. I moved that loop from this code into the call to lanczopw().
    uuu, singular_values, nsv_found = lanczopw.lanczopw(signals, n_rows, 
                                                        n_columns, nsv_sought)

    # for i in range(20):
    #     for j in range(20):
    #         print "uuu[{}, {}] = {}".format(i+1, j+1, uuu[i,j])

#    print "nsv_found = {}".format(nsv_found)

    # uuu = 0
    # ndiv = 0
    # zprime = 0
    # us = 0
    # unit = 0

    # lanczopw(signals, len(signals), lrow, mcol, kuser, kmax, lsinval, u, v, 
    #          work, lwork, zwork, lzwork, zwork[ilambda], zwork[itrlambda]
    #          zwork[ifvect], ndiv)


#       if (ndiv==-1) then 
# c     kuser = max(kuser - 5,0)
#         if (kuser>0) goto 12
#         print *,'PROPACK did not converge'
#         goto 999
#       endif


    # kfit = kuser                     ! kuser was provided by the user      
    # if (ndiv.gt.0) kfit = ndiv 

    # kfit = kuser
    # if ndiv > 0:
    #     kfit = ndiv

# c
# c     ******************************************************************
# c
# c     Calculation of the Z' matrix, with ...
# c     ======================================
# c      print *,'Computing matrix Z-prime with ...'
# c

#      call Zcalc(kfit,kmax,Lrow,U,zprime,us,unit)

    zprime = _zcalc(nsv_found, n_rows, uuu)

# c
# c     Diagonalization of Z' (=hx), yielding the 'signal'poles',
# c     ===========================================================
# c     also called roots.
# c     =================
# c
# c     zgeev - procedure, which calculates eigvals and eigvecs
# c     of a general matrix - only eigvals are used and computed
# c     into variable root
# c 
  
#       call zgeev('N','N',kfit,zprime,kmax,root,U,Lrow,U,Lrow,
#      +           cwork,Lcwork,rwork,info)

    # This is a LAPACK routine
    # zgeev('N', 'N', kfit, zprime, kmax, root, U, lrow, U, lrow, cwork, lcwork,
    #       rwork, info)

    roots, _ = np.linalg.eig(zprime)

    # Eigenvalues are returned unordered. I sort them here to make it easier
    # to compare them with output from the Fortran code.
    roots = np.array(sorted(roots))

    # for i, root in enumerate(roots):
    #     print "roots[{}] = {:.17E}".format(i + 1, root)



# c 
# c     ******************************************************************
# c
# c     Calculation of dampings (damp) and frequencies (freq) from roots
# c 
 #      pi         = 4.0d0*datan2(1.0d0,1.0d0)    
 #      do 77    i = 1,kfit
 #         damp(i) = dlog(cdabs(root(i)))
 #         freq(i) = datan2(dimag(root(i)),dble(root(i))) / (2.0d0*pi)   
 # 77   continue                  !10-8-99 dreal -> dble

    dampings = np.log(np.abs(roots))

    # dampings = sorted(dampings)
    # for i, damping in enumerate(dampings):
    #     print "damping[{}] = {:.17E}".format(i + 1, damping)

    frequencies = np.arctan2(roots.imag, roots.real) / (math.pi * 2)


    # frequencies = sorted(frequencies)
    # for i, frequencies in enumerate(frequencies):
    #     print "frequency[{}] = {:.17E}".format(i + 1, frequencies)


# c
# c     ******************************************************************
# c
# c     Calculation of complex-valued amplitudes , using the 
# c     pseudoinverse of the Lrow*kfit Vandermonde matrix zeta.
# c
# c     First calculation of zeta:
# c
#   call vanmon(ndp,kmax,kfit,root,zeta)

    zeta = _vanmon(len(signals), roots)

# c     zgells writes solution in space of vector x,
# c     but the vector 'signal' is preserved. 
# c     5-8-99: rcond was -1.0; g77 makes rank = 0!! Not with SUN.
# c             So I set rcond = +1.0d-10!!!!
# c
#       call zgelss(ndp,kfit,1,zeta,ndp,x,ndp,
#      +            sinval,1.0d-10,rank,cwork,Lcwork,rwork,info)

    #v,x,s,rank,info = zgelss(a,b,cond=-1.0,lwork=2*minmn+MAX(maxmn,nrhs),overwrite_a=0,overwrite_b=0)
    lwork = 2*MAX_SINGULAR_VALUES+64*(len(signals)+MAX_SINGULAR_VALUES)
    v, x, s, rank, _, info = zgelss(zeta, signals, 
                                     #cond=1e-10, 
                                     cond=-1.0,
                                     lwork=lwork,
                                     overwrite_a=False, overwrite_b=False)

    # FIXME this code, like the Fortran, ignores possible errors reported 
    # in the "info" return value.

    # print "zgelss return"
    # print x.shape

    # for i in range(4096):
    #     print "x[{}] = {:.17E}".format(i + 1, x[i])


    # print "------- V  ----------"
    # print v.shape, v
    # print "------- X  ----------"
    # print x.shape, x
    # print "------- S  ----------"
    # print s
    # print "------- Rank  ----------"
    # print rank
    # print "------- info  ----------"
    # print info

#       do  88  i = 1,kfit
#       ampl(i) = cdabs (x(i))
#    88 fase(i) = datan2(dimag(x(i)),dble(x(i))) ! 10-8-99 dreal -> dble

    # Discard the uneeded values of x.
    x = x[:nsv_found]

    amplitudes = np.abs(x)
    phases = np.arctan2(x.imag, x.real)

    # for i, value in enumerate(amplitudes):
    #     print "amplitudes[{}] = {}".format(i, value)

    # for i, value in enumerate(phases):
    #     print "phases[{}] = {}".format(i, value)

    return nsv_found, singular_values, frequencies, dampings, amplitudes, phases




#       subroutine vanmon(ndp,kmax,kfit,root,zeta)
# c     vanmon(.) calculates the ndp*kfit Vandermonde matrix zeta.
def _vanmon(n_data_points, roots):
      # integer    i,j,kfit,kmax,ndp
      # complex*16 one,root(*),rootj,temp,zeta(ndp,*)

    nsv_found = len(roots)

    zeta = np.zeros( (n_data_points, nsv_found), np.complex128)

#       one = dcmplx(1.0d0,0.0d0)
# c
# c     First row of zeta:
# c
#       do 10   j = 1,kfit
#    10 zeta(1,j) = one

    zeta[1, :nsv_found] = (1+0j)

# c
# c     Rest of zeta:
# c
#       do 20   j = 1,kfit
#           rootj = root(j)
#            temp = one
#       do 20   i = 2,ndp
#            temp = temp * rootj
#    20 zeta(i,j) = temp


    for j in range(nsv_found):
        root = roots[j]
        #print "vanmon, roots[{}] = {:.17E}".format(j + 1, root)
        temp = (1+0j)
        for i in range(1, n_data_points):
            temp *= root

            zeta[i, j] = temp
            
            #print "zeta[{},{}] = {:.17E}".format(i + 1, j + 1, temp)

    return zeta


def _zcalc(nsv_found, n_rows, uuu):
    # c=======================================================================
    #       subroutine Zcalc(kfit,kmax,Lrow,U,zprime,us,unit)

    # PS - zprime is the output of this function
    # PS - unit is only used in this function
    # PS - us is only used in this function

    uuu_sum = np.zeros( (nsv_found, nsv_found), np.complex128)

 #       integer    i,j,k,kfit,kmax,Lrow,m
# c
#       real*8     uot,zero
# c
#       complex*16 cero,sum,temp
#       complex*16 U(Lrow,*),us(kmax,*),unit(kmax,*),zprime(kmax,*)
# c
#       zero = 0.0d0
#       cero = dcmplx(zero,zero)
#       m    = Lrow

    m = n_rows
# c
#       do 60       i = 1,kfit
#          do 60    j = 1,kfit
#                 sum = cero
#             do 50 k = 1,m-1
#                 sum = sum+ conjg(U(k,i))*U(k+1,j)
# c!              sum = sum+dconjg(U(k,i))*U(k+1,j)
#    50 continue
#              us(i,j) = sum
#    60 continue

    for i in range(nsv_found):
        for j in range(nsv_found):
            # PS This is the Fortran equivalent, but slow --
            # sum_ = (0+0j)
            # for k in range(m - 1):
            #     sum_ += uuu[k][i].conjugate() * uuu[k + 1][j]

            # PS - This is the fast way to do the same --
            # FIXME - is there an off by one error here? Double check sums
            # against Fortran.
            sum_ = (uuu[:m - 1, i].conjugate() * uuu[1:m, j]).sum()

            #print "sum({}, {}) = {:.17E}".format(i + 1, j + 1, sum_)
            uuu_sum[i][j] = sum_



#           sum = cero
#       do 70 i = 1,kfit
# c         sum = sum+ conjg(U(m,i))*U(m,i)
#           sum = sum+dconjg(U(m,i))*U(m,i)
#    70 continue
#           uot = 1.0d0-dble (sum)

    # PS this is the simple Python equivalent of the fortran loop --
    # sum_ = (0+0j)
    # for i in range(kfit):
    #     sum_ += uuu[m - 1, i].conjugate() * uuu[m - 1, i]

    # Here's the fast way to do it -- 
    sum_ = (uuu[m - 1, :nsv_found].conjugate() * uuu[m - 1, :nsv_found]).sum()

    uot = 1.0 - sum_.real

    # print "sum = {:.17E}, uot = {:.17E}".format(sum_, uot)

# c
#       do 80         i = 1,kfit
#          do 80      j = 1,kfit
#             temp      = dcmplx(0.d0,0.d0)
#                               if (j.eq.i) 
#      +      temp      = dcmplx(1.d0,0.d0)
# c           unit(i,j) = temp +  conjg(U(m,i))*U(m,j)/uot
#             unit(i,j) = temp + dconjg(U(m,i))*U(m,j)/uot
#    80 continue
# c

    unit = np.zeros( (nsv_found, nsv_found), np.complex128)

    for i in range(nsv_found):
        for j in range(nsv_found):
            temp = ((1+0j) if j == i else (0+0j))
            unit[i, j] = temp + (uuu[m - 1][i].conjugate() * uuu[m - 1][j] / uot)
            #print "unit({}, {}) = {:.17E}".format(i + 1, j + 1, unit[i, j])

  #     do 100      i = 1,kfit
  #        do 100   j = 1,kfit
  #               sum = cero
  #           do 90 k = 1,kfit
  #               sum = sum + unit(i,k)*us(k,j)
  #  90       continue
  #       zprime(i,j) = sum
  # 100 continue

    zprime = np.zeros( (nsv_found, nsv_found), np.complex128)

    for i in range(nsv_found):
        for j in range(nsv_found):
            zprime[i, j] = (unit[i, :nsv_found] * uuu_sum[:nsv_found, j]).sum()
            #print "zprime({}, {}) = {:.17E}".format(i + 1, j + 1, zprime[i, j])

    return zprime
