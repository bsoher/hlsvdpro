cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c........1.........2.........3.........4.........5.........6.........7..
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine lanczopw_python(signal_r,signal_i,ndp,m,n,k,kmax,sigma,
     c                           uuu_r, uuu_i, lwrk,lzwrk,info)

c m == lrow
c n == mcol
c sigma = lsinval
c info = ndiv = nsv_found = kfit

      real*8     signal_r(ndp), signal_i(ndp)
      integer    ndp, m, n, k, kmax, i
      real*8     sigma(kmax)
      real*8     uuu_r(m, kmax + 1), uuu_i(m, kmax + 1)

c     Local variables.
      real*8     work(lwrk)
      complex*16 signal(ndp)
      complex*16 zwork(lzwrk)
      complex*16 lambda(ndp)
      complex*16 trlambda(ndp)
      complex*16 uuu(m, kmax + 1)
      complex*16 vvv(n, kmax)
      complex*16 fvect(ndp)

C       print *, "inside lanczopw_python()!"
       do i = 1,ndp
          signal(i) = dcmplx(signal_r(i), signal_i(i))
C            PRINT *, "signal(", i, "):", signal(i)
       end do

C        print *, "ndp = ", ndp
C        print *, "m (lrow) = ", m
C        print *, "n (mcol) = ", n
C        print *, "k (kuser, nsv_sought) = ", k
C        print *, "kmax = ", kmax
C        print *, "lwrk = ", lwrk
C        print *, "lzwrk = ", lzwrk


       do i = 1,m
          do j = 1, kmax + 1
C              PRINT *, i, ":", uuu_r(i), uuu_i(i)
              uuu(i,j) = dcmplx(uuu_r(i,j), uuu_i(i,j))
         end do
       end do

      call lanczopw(signal, ndp, m, n, k, kmax, sigma, uuu, vvv,
     c              work, lwrk, zwork, lzwrk, lambda, trlambda,
     c              fvect, info)

C        print *, "lanczopw_python, info = ", info

C      FIXME - maybe this is not necessary? signals might be input-only
C      In other words, lanczopw() might not modify it, so there's no 
C      need to recopy it here.
       do i = 1,ndp
C           PRINT *, i, ":", signal_r(i), signal_i(i)
          signal_r(i) = REALPART(signal(i))
          signal_i(i) = IMAGPART(signal(i))
       end do

       do i = 1, m
          do j = 1, kmax + 1
              uuu_r(i,j) = REALPART(uuu(i,j))
              uuu_i(i,j) = IMAGPART(uuu(i,j))
C               PRINT *, "uuu_i(", i, j, "):", uuu_i(i,j)
         end do
       end do

      end



      subroutine lanczopw(signal,ndp,m,n,k,kmax,sigma,U,V,
     c     work,lwrk,zwork,lzwrk,lambda,trlambda,fvect,info)
      implicit none

      double precision zero, one
      double complex zeroc
      parameter(zero = 0.0d0, one = 1.0d0, zeroc=(0.0d0,0.0d0))
      integer i,m,n,k,kmax,ioption(10),iwork(2*kmax+1),info,lwrk,lzwrk
      integer j
      integer ind1,ndp
      integer*8 planF,planB
      real*8 work(lwrk)
      real*8 doption(10)
      real*8 tol,pi,bnd(kmax)
      real*8 dnrm2,sigma(kmax)
      complex*16 lambda(ndp), trlambda(ndp)
      complex*16 U(m,kmax+1), V(n,kmax)
      complex*16 zwork(lzwrk)
      complex*16 fvect(ndp)
      complex*16 signal(*)
      integer ppoints

      real tt0,tt1
      logical DEBUG

c     This defines FFTW_XXX constants
      include "fftw3.f"

      external dnrm2

c
c     m      :  number of rows of the  Toeplitz matrix
c     n      :  number of columns of the Toeplitz matrix
c
 
C        do i = 1,ndp
C             PRINT *, "signal(", i, "):", signal(i)
C        end do

C        print *, "inside lanczopw()"

C        print *, "ndp = ", ndp
C        print *, "m (lrow) = ", m
C        print *, "n (mcol) = ", n
C        print *, "k (kuser, nsv_sought) = ", k
C        print *, "kmax = ", kmax
C        print *, "lwrk = ", lwrk
C        print *, "lzwrk = ", lzwrk

C        print *, "printing non-zero Us..."
C        do i = 1,m
C            do j = 1,kmax + 1
C                if (U(i,j).ne.0.0) then
C                   PRINT *, "U(", i,j, "):", U(i,j)
C                 end if
C            end do
C        end do


C        print *, "printing non-zero Vs..."
C        do i = 1,n
C            do j = 1,kmax
C                if (V(i,j).ne.0.0) then
C                   PRINT *, "V(", i,j, "):", V(i,j)
C                 end if
C            end do
C        end do

C        print *, "printing non-zero works..."
C        do i = 1,lwrk
C            if (work(i).ne.0.0) then
C                PRINT *, "work(", i, "):", work(i)
C            end if
C        end do

C        print *, "printing non-zero zworks..."
C        do i = 1,lzwrk
C            if (zwork(i).ne.0.0) then
C                PRINT *, "zwork(", i, "):", zwork(i)
C            end if
C        end do


C        print *, "printing non-zero lambdas..."
C        do i = 1,ndp
C             if (lambda(i).ne.0.0) then
C                PRINT *, "lambda(", i, "):", lambda(i)
C             end if
C        end do

C        print *, "printing non-zero trlambda..."
C        do i = 1,ndp
C             if (trlambda(i).ne.0.0) then
C                PRINT *, "trlambda(", i, "):", trlambda(i)
C             end if
C        end do

C        print *, "printing non-zero fvects..."
C        do i = 1,ndp
C             if (fvect(i).ne.0.0) then
C                PRINT *, "fvect(", i, "):", fvect(i)
C             end if
C        end do




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

C        print *, "about to call zcopy, n= ", n, "signal(n) = ", signal(n)
C     zcopy() copies m elements starting at signal(n) to fvect      
      call zcopy(m,signal(n),1,fvect(1),1)

C       do i =1, ndp 
C          	print *, "fvect(", i, ")=", fvect(i)
C       end do

       ind1=ndp

C        ppoints = 0
       do i=n-1,1,-1 
C        	 if (ppoints.eq.0) then
C        	 	print *, "first i = ", i
C        	 end if
C        	 ppoints = ppoints + 1
C        	 print *, "ind1 = ", ind1, "i = ", i
         fvect(ind1)=signal(i)
         ind1 =ind1-1
       end do

C        print *, "exited loop, ppoints = ", ppoints, "i = ", i

C       do i =1, ndp 
C          	print *, "fvect(", i, ")=", fvect(i)
C       end do



C       do i =1, ndp
C C          if (signal(i).ne.fvect(i)) then
C          	print *, i, ": signal = ", signal(i), ", fvect = ", fvect(i)
C C          end if
C       end do
   
c
c computation of lambda and trlambda
c
C        	print *, "about to call dfftw_plan_dft_1d"
        call dfftw_plan_dft_1d(planF, ndp, fvect, lambda, 
     +                        FFTW_FORWARD, FFTW_ESTIMATE)
C        	print *, "done dfftw_plan_dft_1d"

C        	print *, "about to call dfftw_plan_dft_1d"
       call dfftw_plan_dft_1d(planB, ndp, fvect, lambda, 
     +                        FFTW_BACKWARD, FFTW_ESTIMATE)
C        	print *, "done dfftw_plan_dft_1d"

C        	print *, "about to call dfftw_execute_dft"

C         print *, "printing non-zero fvects..."
C         do i = 1,ndp
C             if (fvect(i).ne.0.0) then
C                PRINT *, "fvect(", i, "):", fvect(i)
C             end if
C         end do

C         print *, "printing non-zero trlambda..."
C         do i = 1,ndp
C             if (trlambda(i).ne.0.0) then
C                PRINT *, "trlambda(", i, "):", trlambda(i)
C             end if
C         end do

       call dfftw_execute_dft(planF, fvect, lambda)
C         	print *, "done dfftw_execute_dft"

C         print *, "printing non-zero fvects..."
C         do i = 1,ndp
C             if (fvect(i).ne.0.0) then
C                PRINT *, "fvect(", i, "):", fvect(i)
C             end if
C         end do

C         print *, "printing non-zero lambda..."
C         do i = 1,ndp
C             if (lambda(i).ne.0.0) then
C                PRINT *, "lambda(", i, "):", lambda(i)
C             end if
C         end do

C        	print *, "about to call dconjg"
       do i =1,ndp
         fvect(i)=dconjg(fvect(i))
       end do     

C       do i =1, ndp 
C          	print *, "fvect(", i, ")=", fvect(i)
C       end do


C        	print *, "about to copy fvect"
       ind1=ndp
       do i=2,ndp
          zwork(i)=fvect(ind1)
          ind1=ind1-1
       end do

C        	print *, "about to copy zwork"
       do i=2,ndp
          fvect(i)=zwork(i)
       end do

C       do i =1, ndp 
C          	print *, "fvect(", i, ")=", fvect(i)
C       end do


C        	print *, "about to call dfftw_execute_dft"
       call dfftw_execute_dft(planF, fvect, trlambda)
C        	print *, "done dfftw_execute_dft"
       do i=1, ndp
          lambda(i)=lambda(i)/dble(ndp)
          trlambda(i)=trlambda(i)/dble(ndp) 
       end do 

C       do i =1, ndp 
C          	print *, "trlambda(", i, ")=", trlambda(i)
C       end do


C        	print *, "about to call clearstat"
      call clearstat

C        	print *, "done clearstat"

       call second(tt0)
C        	print *, "about to call zlansvdw"
C        	print *, "lwrk = ", lwrk, "lzwrk = ", lzwrk

      call zlansvdw('y','n', ndp, m,n,k,kmax,U,m,sigma,bnd,V,n,
     c      tol,work,lwrk,iwork,doption,ioption,info,
     c      zwork,lzwrk,lambda,trlambda, 0,0)
C        	print *, "done zlansvdw"

       call second(tt1)


       call dfftw_destroy_plan(planF)
       call dfftw_destroy_plan(planB)

C        print *, "exiting lanczopw, k = ", k, "info = ", info

C        print *, "At end of lanczopw:"
C        do i = 1,m
C            do j = 1,kmax+1
C                PRINT *, "U(", i,j, "):", U(i,j)
C C                if (U(i,j).ne.0.0) then
C C                   PRINT *, "U(", i,j, "):", U(i,j)
C C                 end if
C            end do
C        end do

      end

c
c
c
c
