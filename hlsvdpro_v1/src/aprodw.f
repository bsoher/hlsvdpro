
      subroutine aprodw(transa,ndp,m,n,z,y,z1,y1,lambda,
     c                  trlambda,planF,planB)


c
c     Computes y=TH*T*z where lambda contains the eigenvalues of the circulant
c     matrix used for T and trlambda contains the eigenvalues of the circulant
c     matrix used for TH.
c

      parameter (zeroc=(0.0d0,0.0d0))
      integer m,n,ndp
      complex*16 z(*),y(*),lambda(*),trlambda(*),z1(ndp),y1(ndp)

      character*1 transa
      integer i,planF,planB
      logical lsame
      external lsame 

      if (lsame(transa,'t')) then
        do i=1,m
           z1(i)=z(i)
        end do
        do i=m+1,ndp
           z1(i)=zeroc
        end do

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

        call ztmultz(lambda(1),ndp,z1(1),y1(1),planF,planB)
        do i=1,m
           y(i)=y1(i)
        end do

      end if

      end



      subroutine ztmultz(lambda,ndp,z,y,planF,planB)
c
c Multiplicerar en toepliz-matris med en vektor med hjalp av fft
c
      integer m,n, ii

      integer i,ndp,planF,planB
      double complex lambda(*), z(*), y(*)
C      double complex fz(m+n),temp


      
c       do ii = 1, m + n
c          fz(ii)=(0.0d0,0.0d0)
c       end do
c      call zcopy(n,z(1),1,fz(1),1)
c      call dcfftf(m+n,fz(1),wsave(1))

c      y=FFT(z)

      call dfftw_execute_dft(planF, z, y)   


c
      do ii = 1, ndp
         z(ii)=y(ii)*lambda(ii)
      end do
      call dfftw_execute_dft(planB, z, y)   

c      call dcfftb(m+n,fz(1),wsave(1))
c      call zcopy(m,fz(1),1,y(1),1)
cC      temp=cmplx(1.0d0/dble(m+n),0.0d0,kind(1.0d0))
c      temp=cmplx(1.0d0/dble(m+n),0.0d0)
c      call zscal(m,temp,y(1),1)
      end

