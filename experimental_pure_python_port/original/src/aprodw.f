      subroutine aprodw(transa,ndp,m,n,z,y,z1,y1,lambda,
     c                  trlambda,planF,planB)

C     lambda and trlambda are input and are not altered
C     y, y1, z, and z1 are output
c
c     Computes y=TH*T*z where lambda contains the eigenvalues of the circulant
c     matrix used for T and trlambda contains the eigenvalues of the circulant
c     matrix used for TH.
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y == A * X.
c           If TRANSA.EQ.'T' then the function should compute the matrix-vector
c           product Y = A^T * X.
c           The arrays IPARM and DPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.

      parameter (zeroc=(0.0d0,0.0d0))
      integer m,n,ndp
      complex*16 z(*),y(*),lambda(*),trlambda(*),z1(ndp),y1(ndp)

C     z(m), y(n)

      character*1 transa
      integer i,planF,planB
      logical lsame
      external lsame 
      integer ips

      if (lsame(transa,'t')) then
        do i=1,m
           z1(i)=z(i)
        end do
        do i=m+1,ndp
           z1(i)=zeroc
        end do

C        print *, "calling ztmultz, passing ", trlambda(1)
        call ztmultz(trlambda(1),ndp,z1(1),y1(1),planF,planB)
        
        do i=1,n
           y(i)=y1(i)
        end do
      else
        do i=1,n
           z1(i)=z(i)
        end do
        do i=n+1,ndp
           z1(i)=zeroc
        end do

C         if (planF.eq.6) then
C             do ips=1, 10
C               print *, "lambda(", ips, ") = ", lambda(ips)
C             end do
C             do ips=1, n + 1
C               print *, "~ ", ips, " ~ = ", z1(ips)
C             end do
C             do ips=1, 10
C               print *, "y1(", ips, ") = ", y1(ips)
C             end do
C         end if

        call ztmultz(lambda(1),ndp,z1(1),y1(1),planF,planB)

C         if (planF.eq.6) then
C             print *, "yertz!, zeroc = ", zeroc

C             do ips=1, 10
C               print *, "z(", ips, ") = ", z(ips)
C             end do
C             print *, "-----"
C             do ips=1, 10
C               print *, "y(", ips, ") = ", y(ips)
C             end do
C             print *, "-----"
C             do ips=1, 10
C               print *, "z1(", ips, ") = ", z1(ips)
C             end do
C             print *, "-----"
C             do ips=1, 10
C               print *, "y1(", ips, ") = ", y1(ips)
C             end do
C             print *, "-----"

C             print *, "yertz!"
C         end if 

        do i=1,m
           y(i)=y1(i)
        end do

      end if

      end



C     PS - I changed this routine to ignore the FFT plan variables
C     The plan is computed locally in this routine.
C     lambda is input, y and z are output 

      subroutine ztmultz(lambda,ndp,z,y,dummyF,dummyB)
C                ztmultz(lambda,ndp,z,y,planF,planB)
c
c Multiplicerar en toepliz-matris med en vektor med hjalp av fft
c
      integer m,n, ii

      integer i,ndp,dummyF,dummyB
      integer ips, jps
      double complex lambda(*), z(*), y(*)
C      double complex fz(m+n),temp

      integer*8 planFF,planBB

      include "fftw3.f"
 


C     print *, "inside ztmultz, lambda  = ", lambda


C     PS - I assign values to these just to silence compiler 
C     warnings about unused variables.
C       dummyF = 42
C       dummyB = 42
      
c       do ii = 1, m + n
c          fz(ii)=(0.0d0,0.0d0)
c       end do
c      call zcopy(n,z(1),1,fz(1),1)
c      call dcfftf(m+n,fz(1),wsave(1))

c      y=FFT(z)

C     PS The plans are now computed locally. This is really inefficient,
C     but it allows me to use numpy's FFT in lancopw and FFTW here.
      call dfftw_plan_dft_1d(planFF, ndp, z, y, 
     +                        FFTW_FORWARD, FFTW_ESTIMATE)
C       if (dummyF.eq.6) then
C         do ips = 1, 10
C             print *, "zmultz: y(", ips, ") = ", y(ips)
C         end do
C       end if
      call dfftw_execute_dft(planFF, z, y)   
      call dfftw_destroy_plan(planFF)
c

C       if (dummyF.eq.6) then
C         do ips = 1, 10
C             print *, "zmultz: y(", ips, ") = ", y(ips)
C         end do
C       end if


      do ii = 1, ndp
         z(ii)=y(ii)*lambda(ii)
      end do
C       do ips = 1, ndp
C           print *, "zmultz: z(", ips, ") = ", z(ips)
C       end do

       call dfftw_plan_dft_1d(planBB, ndp, z, y, 
     +                        FFTW_BACKWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(planBB, z, y)   
      call dfftw_destroy_plan(planBB)


C         do i=1, ndp
C           print *, "z(", i, ") = ", z(i)
C         end do
C         print *, "-----"
C         do i=1, ndp
C           print *, "y(", i, ") = ", y(i)
C         end do
C         print *, "-----"

C         print *, "yertz!"


c      call dcfftb(m+n,fz(1),wsave(1))
c      call zcopy(m,fz(1),1,y(1),1)
cC      temp=cmplx(1.0d0/dble(m+n),0.0d0,kind(1.0d0))
c      temp=cmplx(1.0d0/dble(m+n),0.0d0)
c      call zscal(m,temp,y(1),1)
      end



