# Python modules
from __future__ import division
import pdb

# 3rd party modules
import numpy as np


# Our modules




#  subroutine aprodw(transa,ndp,m,n,z,y,z1,y1,lambda,
# c                  trlambda,planF,planB)

def aprodw(transa, ndp, m, n, z, y, z1, y1, lambda_, trlambda):
# c     Computes y=TH*T*z where lambda contains the eigenvalues of the circulant
# c     matrix used for T and trlambda contains the eigenvalues of the circulant
# c     matrix used for TH.


#       parameter (zeroc=(0.0d0,0.0d0))
#       integer m,n,ndp
#       complex*16 z(*),y(*),lambda(*),trlambda(*),z1(ndp),y1(ndp)

#       character*1 transa
#       integer i,planF,planB
#       logical lsame
#       external lsame 


#       if (lsame(transa,'t')) then
    transa = transa.lower()
    if transa == 't':
#       do i=1,m
#          z1(i)=z(i)
#       end do
        z1[:m] = z[:m]

#       do i=m+1,ndp
#          z1(i)=zeroc
#       end do
        z1[m:ndp] *= 0j

#       call ztmultz(trlambda(1),ndp,z1(1),y1(1),planF,planB)
        z1_copy, y1_copy = ztmultz(trlambda, ndp, z1, y1)

        z1[:ndp] = z1_copy[:ndp]
        y1[:ndp] = y1_copy[:ndp]

#       do i=1,n
#         y(i)=y1(i)
#       end do
        y[:n] = y1[:n]
#   else
    else:
#       do i=1,n
#         z1(i)=z(i)
#       end do
        z1[:n] = z[:n]

#       do i=n+1,ndp
#         z1(i)=zeroc
#       end do
        z1[n:ndp] *= 0j

#       call ztmultz(lambda(1),ndp,z1(1),y1(1),planF,planB)
        # for ips, value in enumerate(z1[:n]):
        #     print "~ {:.17G}\t{:.17G}".format(value.real, value.imag)
        z1_copy, y1_copy = ztmultz(lambda_, ndp, z1, y1)
        z1[:ndp] = z1_copy[:ndp]
        y1[:ndp] = y1_copy[:ndp]

        # print "len(z, y, z1, y1) = {}, {}, {}, {}".format(len(z), len(y), 
        #                                                   len(z1), len(y1))

        # for i in range(len(z)):
        #     print "z[{}] = {:.17G}".format(i + 1, z[i])
        
        # for i in range(len(y)):
        #     print "y[{}] = {:.17G}".format(i + 1, y[i])

        # for i in range(len(z1)):
        #     print "z1[{}] = {:.17G}".format(i + 1, z1[i])
        
        # for i in range(len(y1)):
        #     print "y1[{}] = {:.17G}".format(i + 1, y1[i])

#       do i=1,m
#          y(i)=y1(i)
#       end do
        y[:m] = y1[:m]

#   end if

    return z, y, z1, y1

#   end



#       subroutine ztmultz(lambda,ndp,z,y,dummyF,dummyB)
# C                ztmultz(lambda,ndp,z,y,planF,planB)
# c
# c Multiplicerar en toepliz-matris med en vektor med hjalp av fft
# c
def ztmultz(lambda_, ndp, z, y):
    """Multiply a Toeplitz matrix by a vector with the help of FFT"""

#       integer m,n, ii
#       integer i,ndp,dummyF,dummyB
#       double complex lambda(*), z(*), y(*)
# C      double complex fz(m+n),temp

#   call dfftw_plan_dft_1d(planFF, ndp, z, y, 
#  +                        FFTW_FORWARD, FFTW_ESTIMATE)
#   call dfftw_execute_dft(planFF, z, y)   
#   call dfftw_destroy_plan(planFF)
    
    y = np.fft.fft(z)

#   do ii = 1, ndp
#      z(ii)=y(ii)*lambda(ii)
#   end do
    z = y * lambda_

#   call dfftw_plan_dft_1d(planBB, ndp, z, y, 
#  +                        FFTW_BACKWARD, FFTW_ESTIMATE)
#   call dfftw_execute_dft(planBB, z, y)   
#   call dfftw_destroy_plan(planBB)
    y = np.fft.ifft(z)

    # The FFTW doc says, "FFTW computes an unnormalized transform, in that 
    # there is no coefficient in front of the summation in the DFT. In other 
    # words, applying the forward and then the backward transform will multiply 
    # the input by n."
    # http://www.fftw.org/faq/section3.html#whyscaled
    # Apparently numpy's FFT is normalized, so in order to get the same 
    # results as FFTW computes, we have to multiply by n a.k.a. ndp.
    y *= ndp

    # for i in range(len(z)):
    #     print "z[{}] = {:.17G}".format(i + 1, z[i])
    # for i in range(len(y)):
    #     print "y[{}] = {:.17G}".format(i + 1, y[i])

# c      call dcfftb(m+n,fz(1),wsave(1))
# c      call zcopy(m,fz(1),1,y(1),1)
# cC      temp=cmplx(1.0d0/dble(m+n),0.0d0,kind(1.0d0))
# c      temp=cmplx(1.0d0/dble(m+n),0.0d0)
# c      call zscal(m,temp,y(1),1)
#       end

    return z, y
