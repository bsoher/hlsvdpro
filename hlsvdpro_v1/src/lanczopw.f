      subroutine lanczopw(signal,ndp,m,n,k,kmax,sigma,U,V,
     c     work,lwrk,zwork,lzwrk,lambda,trlambda,fvect,info)
      implicit none

      double precision zero, one
      double complex zeroc
      parameter(zero = 0.0d0, one = 1.0d0, zeroc=(0.0d0,0.0d0))
      integer i,m,n,k,kmax,ioption(10),iwork(2*kmax+1),info,lwrk,lzwrk
      integer ind1,ndp
      integer*8 planF,planB
      double precision work(lwrk)
      double precision doption(10)
      double precision tol,pi,bnd(kmax)
      double precision dnrm2,sigma(kmax)
      complex*16 lambda(ndp),trlambda(ndp),
     c          U(m,kmax+1),V(n,kmax),zwork(lzwrk),
     c           fvect(ndp),signal(*)

      real tt0,tt1
      logical DEBUG

c     This defines FFTW_XXX constants
      include "fftw3.f"

      external dnrm2

c
c     m      :  number of rows of the  Toeplitz matrix
c     n      :  number of columns of the Toeplitz matrix
c
 
      DEBUG = .false.

      tol = 16d-16

      ioption(1) = 1
      ioption(2) = 1
      doption(1) = 1d-12
      doption(2) = 1d-14
      doption(3) = zero

      pi = 2*acos(0d0)



c
c  Initialization of FFTW
c

c
c  computation of the first column of the circulant matrix
c      
 
      do i=1, ndp
         fvect(i)=zeroc
      end do

      call zcopy(m,signal(n),1,fvect(1),1)
       ind1=ndp
       do i=n-1,1,-1 
         fvect(ind1)=signal(i)
         ind1 =ind1-1
       end do
   
c
c computation of lambda and trlambda
c
        call dfftw_plan_dft_1d(planF, ndp, fvect, lambda, 
     +                        FFTW_FORWARD, FFTW_ESTIMATE)

       call dfftw_plan_dft_1d(planB, ndp, fvect, lambda, 
     +                        FFTW_BACKWARD, FFTW_ESTIMATE)

       call dfftw_execute_dft(planF, fvect, lambda)

       do i =1,ndp
         fvect(i)=dconjg(fvect(i))
       end do     

       ind1=ndp
       do i=2,ndp
          zwork(i)=fvect(ind1)
          ind1=ind1-1
       end do

       do i=2,ndp
          fvect(i)=zwork(i)
       end do

       call dfftw_execute_dft(planF, fvect, trlambda)
       do i=1, ndp
          lambda(i)=lambda(i)/dble(ndp)
          trlambda(i)=trlambda(i)/dble(ndp) 
       end do 


      call clearstat

       call second(tt0)
      call zlansvdw('y','n', ndp, m,n,k,kmax,U,m,sigma,bnd,V,n,
     c      tol,work,lwrk,iwork,doption,ioption,info,
     c      zwork,lambda,trlambda, planF,planB)
       call second(tt1)


       call dfftw_destroy_plan(planF)
       call dfftw_destroy_plan(planB)

      end

c
c
c
c
