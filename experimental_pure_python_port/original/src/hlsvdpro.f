cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c........1.........2.........3.........4.........5.........6.........7..
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine hlsvdpw_python(signal_r,signal_i,ndp,Lrow,mcoL,
     +                          kuser,kfit,Lsinval,ampl,fase,damp,freq)

C    This subroutine was added by Philip Semanchuk from the Vespa team. 
C    It's just a simplifying wrapper around hlsvdpw(). The main 
C    differences between this subroutine and hlsvdpw() are --
C      - The complex signal array is passed here as two correlated real 
C        and imaginary arrays. This makes it possible to call via 
C        Python's ctypes which doesn't know anything about passing 
C        complex variables.
C      - This subroutine provides the work arrays rwork, lrwork, zwork, 
C        and lzwork instead of asking the caller to pass them.
C      - This subroutine contains a number of debugging statements 
C        (commented out) to help troubleshoot problems related to 
C        calling from Python.


      implicit none  ! All variables have to be declared.

      integer          ndpmax
      parameter        (ndpmax=8192)
      real*8           ampl(50),damp(50),freq(50),fase(50)
      real*8           signal_r(ndpmax), signal_i(ndpmax), Lsinval(50)
      complex*16       signal(ndpmax)
      integer          kuser,kfit,Lrow,mcoL,ndp,i,kmax,lrwork,lzwork
      parameter        (kmax=50)
      parameter        (lrwork=20*ndpmax+13*kmax+8*kmax**2+ndpmax*kmax)
      parameter        (lzwork=90*ndpmax+100*kmax+5*kmax**2+4*kmax*
     +                  ndpmax+20)
      real*8           rwork(lrwork)
      complex*16       zwork(lzwork)

C       print *, "kuser,kfit,Lrow,mcoL,ndp=",kuser,kfit,Lrow,mcoL,ndp
      
       do i = 1,ndp
C           PRINT *, i, ":", signal_r(i), signal_i(i)
          signal(i) = dcmplx(signal_r(i), signal_i(i))
       end do
       

       call hlsvdpw(signal,ndp,Lrow,mcoL,kuser,kfit,kmax,Lsinval,
     +              ampl,fase,damp,freq,rwork,lrwork,zwork,lzwork)

C       do i = 1,kfit
C           print *, "hlsvdpw_wrapper: damp(i)", i, damp(i)
C       enddo
      
C       do i = 1,kfit
C           print *, "hlsvdpw_wrapper: freq(i)", i, freq(i)
C       enddo
      
C       do i = 1,kfit
C           print *, "hlsvdpw_wrapper: ampl(i)", i, ampl(i)
C       enddo
      
C       do i = 1,kfit
C           print *, "hlsvdpw_wrapper: fase(i)", i, fase(i)
C       enddo
      
       return
       end



      subroutine hlsvdpw(signal,ndp,Lrow,mcoL,kuser,kfit,kmax,Lsinval,
     +                   ampl,fase,damp,freq,rwork,lrwork,zwork,lzwork)
c
c    Intermediate function to deal with the workspace allocation.
c    
c      for the     LIST OF IMPORTANT VARIABLES   see hlsvdpro below
c                  ===========================
c
c     Input workspace arrays
c     ----------------------
c     
c     RWORK  = real*8 workspace of size lrwork
c     ZWORK  = complex*16 workspace of size lzwork
c     
c     LRWORK = integer, at least:  20*ndp+13*kmax+8*kmax**2
c     LZWORK = integer, at least:  90*ndp+100*kmax+5*kmax**2+4*kmax*ndp+20
c     
c     D. Sima, KUL, 2007-2010.    

c
c     Declarations of all variables
c     =============================
c
      implicit none  ! All variables have to be declared.
c
      integer          kmax
      integer          kuser,kfit,Lrow,mcoL,ndp,lrwork,lzwork,i
      integer          lwrk,lzwrk
      integer          indU,indV,indus,indunit,indz,indzeta,indc
      real*8           Lsinval(kmax)
      real*8           ampl(*),damp(*),freq(*),fase(*)

      real*8           rwork(lrwork)                 
c
      complex*16       signal(*)

      complex*16       zwork(lzwork)
c

      lwrk    = Lrow+mcoL+13*kmax+8*kmax**2+32*Lrow+ndp*kmax
      lzwrk   = Lrow+mcoL+32*Lrow+7*kmax+2+2*kmax**2+5*ndp
      
      indU    = lzwrk+3*ndp+1
      indV    = indU + 2*Lrow*kmax
      indus   = indV + mcoL*kmax
      indunit = indus + kmax**2
      indz    = indunit + kmax**2
      indzeta = indz + kmax**2
      indc    = indzeta + ndp*kmax

      if (lwrk.gt.lrwork) then
C        print *,'lwrk',lwrk
C        print *,'lrwork',lrwork
         print *,'Increase the dimension of the real workspace'
         goto 999
      endif
      if (indc+2*kmax+64*(ndp+kmax).gt.lzwork) then
         print *,'Increase the dimension of the complex workspace'
         goto 999
      endif

      call      hlsvdpro(signal,ndp,Lrow,mcoL,kuser,kfit,kmax,lsinval,
     +                ampl,fase,damp,freq,rwork,lwrk,zwork(1),lzwrk,
     +                zwork(indU),zwork(indV),zwork(indus),
     +                zwork(indunit),zwork(indz),zwork(indzeta), 
     +                zwork(indc)) 
c

 999  return 
      end
c 
c     &&&&&&&&&&&&&&&  end of subroutine hlsvdw  &&&&&&&&&&&&&&&&&&
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine hlsvdpro(signal,ndp,Lrow,mcoL,kuser,kfit,kmax,Lsinval,
     +                   ampl,fase,damp,freq,work,lwrk,zwork,lzwrk,
     +                   U,V,us,unit,zprime,zeta,cwork)  
c
c                  List of subroutine packages called
c     lanczos.f
c     hlsubrs.f
c     Lapack/blas 
c
c                  LIST OF IMPORTANT VARIABLES
c                  ===========================
c
c     Input from java:
c     ===============
c     signal  = complex*16 nmr signal (data)
c     ndp     = number of datapoints
c     Lrow    = number of rows of hankel data matrix hx, Lrow < ndp/2
c     mcoL    = number of columns of hx, mcoL = ndp+1-Lrow > ndp/2
c     kuser   = number of sinusoids initially wished by user 
c     kmax    = max number of iterations in PROPACK (recommended >= 100)
c
c     Output to java:
c     ==============
c     sinvals = Lrow singular values of hankel data matrix 
c              (Lrow < mcoL)
c     ampl    = kfitted amplitudes, estimated by htls
c     fase    = kfitted phases, estimated by htls
c     damp    = kfitted dampings, estimated by htls
c     freq    = kfitted frequencies, estimated by htls
c     kfit    = number of sinusoids actually fitted: identical to
c               number of sinvals found by Lanczos = ndiv
c
c     Other variables: 
c     =============== 
c     x       = complex*16 work-vector, initially for the signal
c     yr,yi   = real*4 real and imag parts of signal   
c     hx      = Lmax*(2*kmax) work-matrix for ztls and zgeev
c     lsinval = Lanczos sinvals of data-matrix
c     sinval  = work-vector for Lapack subroutines
c     sinvec  = Lmax*kmax left Lanczos sinvec matrix of data-matrix,
c               and, later on, right Lapack sinvec matrix of hx
c     vvt     = (2*kfit)*(2*kfit) right sinvecs of hx in ztls
c     root    = kfit complex-valued eigvals of Z' 
c     zeta    = Lrow*kfit Vandermonde matrix 
c     kmax    = max number of sinusoids
c     Lmax    = max number of rows  < ndp/2
c     Lrow    = number of rows
c     coef    = coefficient for optimal size of working arrays
c               for Lapack
c    
c     this routine differs computes the largest singular triplets
c     of the Hankel matrix by means of the Lanczos algorithm
c     with partial reorthogonalization
c     In particular it is  a modified version of the Larsen code
c     (PROPACK)  
c


c     Declarations of all variables
c     =============================
c
      implicit none  ! All variables have to be declared.
c
      integer          kmax
      integer          i,kuser,kfit,Lrow,mcoL,ndp,lwrk,lzwrk
      integer          ilambda,itrlambda,ifvect
      integer          info,rank,ndiv                       ! Lapack
      integer          Lcwork,Lrwork                        ! Lapack
c             (Lcwork>2*kmax+64*(ndp+kmax),Lrwork>5*kmax)   ! Lapack  zgelss

c
      real*4           tt0, tt1
      complex*16       us(kmax,kmax),unit(kmax,kmax)        ! Lanczos
      complex*16       zprime(kmax,kmax)                    !lanczos and Lapack
      
      real*8           sinval(2*kmax),Lsinval(kmax)
      real*8           ampl(*),damp(*),freq(*),fase(*)
      real*8           pi
      real*8           rwork(5*kmax),work(lwrk)             ! Lapack
c
      complex*16       signal(*),x(ndp)
      complex*16       U(Lrow,2*kmax)
      complex*16       cwork(2*kmax+64*(ndp+kmax))          ! Lapack
      complex*16       root(kmax)                           !eigvals of Z'
      complex*16       zeta(ndp,kmax)                       !vandermonde
      complex*16       V(mcoL,kmax),zwork(lzwrk+3*ndp)
      integer          ppi, ppj

c
c     ******************************************************************
c
c     Calculation of sinvals & sinvecs of the hankel data matrix.
c     ==========================================================
c
c     Copy complex*16 signal to real*4 yr and yi, and to complex*16 x.
c     ---------------------------------------------------------------
c
c
      if (kuser.GT.kmax) kuser = kmax

      do i = 1,ndp
         x(i) = signal(i)       ! x corrupted in amplitude computation
      end do

      
c
c computation of new ndp for the circulant matrix
c
c        call caldim(lrow,mcol,ndp)  

c
c     ******************************************************************
c
c     Lanczos SVD of hankel datamatrix.
c     ================================
c
c     Lcwork = 2*kmax+ndp
      Lcwork = 2*kmax+64*(ndp+kmax)
      Lrwork = 5*kmax
      ilambda = lzwrk+1
      itrlambda = ilambda + ndp
      ifvect = itrlambda + ndp
      call second(tt0) 
          
C       print *, "about to call lanczopw"
  12      call lanczopw(signal,ndp,Lrow,mcoL,
     c         kuser,kmax,Lsinval,U,V,
     c     work,lwrk,zwork,lzwrk,zwork(ilambda),zwork(itrlambda),
     c     zwork(ifvect),ndiv)
C       print *, "done calling lanczopw"

      if (ndiv==-1) then 
c     kuser = max(kuser - 5,0)
        if (kuser>0) goto 12
        print *,'PROPACK did not converge'
        goto 999
      endif

c 
       call second(tt1)

      kfit = kuser                     ! kuser was provided by the user      
      if (ndiv.gt.0) kfit = ndiv 


c
c     ******************************************************************
c
c     Calculation of the Z' matrix, with ...
c     ======================================
c      print *,'Computing matrix Z-prime with ...'
c
C       do ppi = 1, lrow
C       	do ppj = 1, 2*kmax
C       		print *, "U(", ppi, ",", ppj, ") = ", U(ppi, ppj)
C       	end do
C       end do      

      call Zcalc(kfit,kmax,Lrow,U,zprime,us,unit)
c
c     Diagonalization of Z' (=hx), yielding the 'signal'poles',
c     ===========================================================
c     also called roots.
c     =================
c
c     zgeev - procedure, which calculates eigvals and eigvecs
c     of a general matrix - only eigvals are used and computed
c     into variable root
c 
  
      call zgeev('N','N',kfit,zprime,kmax,root,U,Lrow,U,Lrow,
     +           cwork,Lcwork,rwork,info)

      call sortroot(kfit, root)
C       do i=1, kfit
C       	print *, "root(", i, ") = ", root(i)
C       end do
c 
c     ******************************************************************
c
c     Calculation of dampings (damp) and frequencies (freq) from roots
c 
      pi         = 4.0d0*datan2(1.0d0,1.0d0)    
      do 77    i = 1,kfit
         damp(i) = dlog(cdabs(root(i)))
         freq(i) = datan2(dimag(root(i)),dble(root(i))) / (2.0d0*pi)   
 77   continue                  !10-8-99 dreal -> dble

C       call sortreals(kfit, damp)
C       do i=1, kfit
C       	print *, "damp(", i, ") = ", damp(i)
C       end do


C       call sortreals(kfit, freq)
C       do i=1, kfit
C       	print *, "freq(", i, ") = ", freq(i)
C       end do

c
c     ******************************************************************
c
c     Calculation of complex-valued amplitudes , using the 
c     pseudoinverse of the Lrow*kfit Vandermonde matrix zeta.
c
c     First calculation of zeta:
c
      call vanmon(ndp,kmax,kfit,root,zeta)
c
c     zgells writes solution in space of vector x,
c     but the vector 'signal' is preserved. 
c     5-8-99: rcond was -1.0; g77 makes rank = 0!! Not with SUN.
c             So I set rcond = +1.0d-10!!!!
c
      call zgelss(ndp,kfit,1,zeta,ndp,x,ndp,
     +            sinval,1.0d-10,rank,cwork,Lcwork,rwork,info)
c
c
C       do i =1, 4096
C       	print *, "x(", i, ") = ", x(i) 
C       end do

      do  88  i = 1,kfit
      ampl(i) = cdabs (x(i))
   88 fase(i) = datan2(dimag(x(i)),dble(x(i))) ! 10-8-99 dreal -> dble
c
c  
C       do  i = 1,kfit
C         print *, "ampl(", i, ") = ", ampl(i)
C       end do

C       do  i = 1,kfit
C         print *, "fase(", i, ") = ", fase(i)
C       end do

c     ******************************************************************
c
999   return 
      end
c 
c     &&&&&&&&&&&&&&&  end of subroutine hltls/hlsvd  &&&&&&&&&&&&&&&&&&
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine vanmon(ndp,kmax,kfit,root,zeta)
c     vanmon(.) calculates the ndp*kfit Vandermonde matrix zeta.
c
      implicit   none
      integer    i,j,kfit,kmax,ndp
      complex*16 one,root(*),rootj,temp,zeta(ndp,*)
c
      one = dcmplx(1.0d0,0.0d0)
c
c     First row of zeta:
c
      do 10   j = 1,kfit
   10 zeta(1,j) = one
c
c     Rest of zeta:
c
      do 20   j = 1,kfit
          rootj = root(j)
C           print *, "root(", j, ") = ", root(j)
           temp = one
      do 20   i = 2,ndp
           temp = temp * rootj
   20 zeta(i,j) = temp


C       do j = 1, kfit
C       	do i = 2, ndp
C       	  print *, "zeta(", i, ",", j, ") = ", zeta(i, j)
C       	end do
C       end do

      return
      end


c
c=======================================================================
      subroutine Zcalc(kfit,kmax,Lrow,U,zprime,us,unit)
c
c     31-08-1999: Calculation of Z-matrix moved to separate subroutine.
c
      implicit none
c
      integer    i,j,k,kfit,kmax,Lrow,m
c
      real*8     uot,zero
c
      complex*16 cero,sum,temp
      complex*16 U(Lrow,*),us(kmax,*),unit(kmax,*),zprime(kmax,*)
c
      zero = 0.0d0
      cero = dcmplx(zero,zero)
      m    = Lrow
c
      do 60       i = 1,kfit
         do 60    j = 1,kfit
                sum = cero
            do 50 k = 1,m-1
                sum = sum+ conjg(U(k,i))*U(k+1,j)
c!              sum = sum+dconjg(U(k,i))*U(k+1,j)
   50 continue
C             print *, "sum(", i, ",", j, ") = ", sum
             us(i,j) = sum
   60 continue
c
          sum = cero
      do 70 i = 1,kfit
c         sum = sum+ conjg(U(m,i))*U(m,i)
          sum = sum+dconjg(U(m,i))*U(m,i)
   70 continue
          uot = 1.0d0-dble (sum)
c
C      print *, "sum = ", sum, ", uot = ", uot

      do 80         i = 1,kfit
         do 80      j = 1,kfit
            temp      = dcmplx(0.d0,0.d0)
                              if (j.eq.i) 
     +      temp      = dcmplx(1.d0,0.d0)
c           unit(i,j) = temp +  conjg(U(m,i))*U(m,j)/uot
            unit(i,j) = temp + dconjg(U(m,i))*U(m,j)/uot
C             print *, "unit(", i, ",", j, ") = ", unit(i,j)     
   80 continue
c
      do 100      i = 1,kfit
         do 100   j = 1,kfit
                sum = cero
            do 90 k = 1,kfit
                sum = sum + unit(i,k)*us(k,j)
   90       continue
        zprime(i,j) = sum
C         print *, "zprime(", i, ",", j, ") = ", zprime(i,j)     
  100 continue
c

      return
      end     

c=======================================================================
c
      subroutine sortroot(kfit, root)
c
c     Sorts roots using a crappy insertion sort. Written by Philip.

      integer     i,j
      complex*16  root(kfit), current

      do j=2, kfit
      	current = root(j)
      	do i = j-1, 1, -1
      		if (REALPART(root(i)) <= REALPART(current)) goto 10
      		root(i + 1) = root(i)
        end do
        i = 0
   10   root(i+1) = current
      end do

      end

c=======================================================================
c
      subroutine sortreals(laaa, aaa)
c
c     Sorts an array of real numbers using a crappy insertion sort. 
C     Written by Philip.

      integer     i,j
      real*8      aaa(laaa), current

      do j=2, laaa
      	current = aaa(j)
      	do i = j-1, 1, -1
      		if (aaa(i) <= current) goto 10
      		aaa(i + 1) = aaa(i)
        end do
        i = 0
   10   aaa(i+1) = current
      end do

      end

c=======================================================================
c
      subroutine sortfreq(kfit,ampl,fase,damp,freq)
c
c     Sort sinusoid parameters according to increasing frequency
c
c     Last change: 01-09-1999     
c
      implicit none
c
      integer     k,kfit,kstop,ktop
      real*8      damp(*),freq(*),ampl(*),fase(*)
      real*8      damp1,damp2,freq1,freq2,ampl1,ampl2,fase1,fase2  
c
           ktop = kfit
   10     kstop = 0
      do 20   k = 1, ktop-1
          damp1 = damp(k)
          damp2 = damp(k+1)
          freq1 = freq(k)
          freq2 = freq(k+1)
          ampl1 = ampl(k)
          ampl2 = ampl(k+1)
          fase1 = fase(k)
          fase2 = fase(k+1)
      if (freq1.gt.freq2) then
        damp(k) = damp2
      damp(k+1) = damp1
        freq(k) = freq2
      freq(k+1) = freq1
        ampl(k) = ampl2
      ampl(k+1) = ampl1
        fase(k) = fase2
      fase(k+1) = fase1
          kstop = k
      end if
   20 continue
c
          ktop = kstop
      if (kstop.gt.1) go to 10    
c
c
      return
      end




c
c
c=======================================================================
c      
      subroutine savev(v1,n1,k1,hh,nry,kk)
c
c-----------------------------------------------------------------------
c modifications by a. van den boogaart
c last update: 18-02-1991
c-----------------------------------------------------------------------
c
      integer*4 n1,k1,nry,kk
      real*8    v1(*)
      real*8    hh(nry,*)

      integer*4 i

      do 10 i=1,n1-k1+1
         hh(i,kk)=v1(k1+i)
         hh(n1-k1+1+i,kk)=v1(n1+k1+1+i)
   10 continue

      return
      end

c
c=======================================================================
c      
      subroutine evscan(nik,n,ma,alpha,beta,nev,low,up,div,ndiv,
     +                  tk,tk1,ierr,eps)
c
c-----------------------------------------------------------------------
c this subprogram is portable to other types of computer
c
c subroutines required:
c   imtql1,mondis,scanmp,comint
c
c modifications by a. van den boogaart
c last update: 04-03-1991
c-----------------------------------------------------------------------
c
      integer*4 nik,n,ma,nev,ndiv,ierr
      real*8    eps
      real*8    alpha(*),beta(*)
      real*8    tk(ma,2),tk1(ma,2)
      logical*4 low,up,div

      integer*4 i,j,k,k1,m1
      real*8    epsold

      nev=0
      low = .false.
      up = .false.
      ndiv = 0
      do 10 i=2,nik
         j = nik - i + 1
         tk(j+1,2) = beta(j)
   10 continue
      do 20 i=1,nik
         tk(i,1) = alpha(i)
         tk1(i,1) = alpha(i)
         tk1(i,2) = tk(i,2)
   20 continue
      m1 = nik - 1

      call imtql1(m1,tk,tk(1,2),ierr)

      if (ierr.ne.0) return
      tk(nik,1) = 0.0d0

      call imtql1(nik,tk1,tk1(1,2),ierr)

      if (ierr.ne.0) return

      call mondis(tk,tk1,m1,eps)

      eps = dble(n)*dmax1(eps,2.0d0**(-52))
c
c if one wants a larger value of eps, this value should be inserted
c here
c
      do 30 i=1,m1
         tk(i,2) = tk(i,1)
         tk1(i,2) = tk1(i,1)
   30 continue
      tk1(nik,2) = tk1(nik,1)
      k = m1
      k1 = nik

      call scanmp(ma,k,k1,tk,tk1,eps)

   40 continue
      epsold = eps


      call comint(ma,tk,tk1,k,k1,eps,low,up,nev,div,ndiv)

c     write(*,*) tk(2,2)
c     write(*,*) !26-08-1999

      if (eps.gt.epsold) go to 40
      if (.not.low) div = .false.

      return
      end
c
c=======================================================================
c      
      subroutine imtql1(n,d,e,ierr)
c
c-----------------------------------------------------------------------
c this subroutine is a translation of the algol procedure imtql1,
c num. math. 12, 377-383(1968) by martin and wilkinson,
c as modified in num. math. 15, 450(1970) by dubrulle.
c handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).
c
c this subroutine finds the eigenvalues of a symmetric
c tridiagonal matrix by the implicit ql method.
c
c on input-
c   n     is the order of the matrix,
c   d     contains the diagonal elements of the input matrix,
c   e     contains the subdiagonal elements of the input matrix
c         in its last n-1 positions.  e(1) is arbitrary.
c
c on output-
c   d     contains the eigenvalues in ascending order.  if an
c         error exit is made, the eigenvalues are correct and
c         ordered for indices 1,2,...ierr-1, but may not be
c         the smallest eigenvalues,
c   e     has been destroyed,
c   ierr  is set to-
c                  zero   for normal return,
c                  j      if the j-th eigenvalue has not been
c                         determined after 30 iterations.
c
c questions and comments should be directed to b. s. garbow,
c applied mathematics division, argonne national laboratory
c
c modifications by a. van den boogaart
c last update: 11-12-1990
c-----------------------------------------------------------------------
c
      integer*4 n,ierr
      real*8    d(n),e(n)

      integer*4 i,ii,j,l,m,mml
      real*8    b,c,f,g,machep,p,r,s
c
c machep is a machine dependent parameter specifying
c the relative precision of floating point arithmetic.
c
      machep = 2.0d0**(-50) 
 
      ierr = 0
      if (n.eq.1) go to 140
 
      do 10 i=2,n
         e(i-1) = e(i)
   10 continue
 
      e(n) = 0.0d0
 
      do 120 l=1,n
         j = 0
c
c look for small sub-diagonal element
c
   20    do 30 m=l,n
            if (m.eq.n) go to 40
            if (dabs(e(m)).le.machep*(dabs(d(m))+dabs(d(m+1)))) goto 40
   30    continue
 
   40    p = d(l)
         if (m.eq.l) go to 80
         if (j.eq.30) go to 130
         j = j + 1
c
c form shift
c
         g = (d(l+1)-p)/(2.0d0*e(l))
         r = dsqrt(g*g+1.0d0)
         g = d(m) - p + e(l)/(g+dsign(r,g))
         s = 1.0d0
         c = 1.0d0
         p = 0.0d0
         mml = m - l
c
c for i=m-1 step -1 until l do --
c
         do 70 ii=1,mml
            i = m - ii
            f = s*e(i)
            b = c*e(i)
            if (dabs(f).lt.dabs(g)) go to 50
            c = g/f
            r = dsqrt(c*c+1.0d0)
            e(i+1) = f*r
            s = 1.0d0/r
            c = c*s
            go to 60
   50       s = f/g
            r =dsqrt(s*s+1.0d0)
            e(i+1) = g*r
            c = 1.0d0/r
            s = s*c
   60       g = d(i+1) - p
            r = (d(i)-g)*s + 2.0*c*b
            p = s*r
            d(i+1) = g + p
            g = c*r - b
   70    continue
 
         d(l) = d(l) - p
         e(l) = g
         e(m) = 0.0d0
         go to 20
c
c order eigenvalues
c
   80    if (l.eq.1) go to 100
c
c for i=l step -1 until 2 do --
c
         do 90 ii=2,l
            i = l + 2 - ii
            if (p.ge.d(i-1)) go to 110
            d(i) = d(i-1)
   90    continue

  100    i = 1
  110    d(i) = p
  120 continue

      go to 140
c
c set error -- no convergence to an eigenvalue after 30 iterations
c
  130 ierr = l
  140 return
c
c last card of imtql1
c date      01/10/76
c end of deck
c
      end
c
c=======================================================================
c
      subroutine mondis(tk,tk1,k,eps)
c
c-----------------------------------------------------------------------
c modifications by a. van den boogaart
c last update: 11-12-1990
c-----------------------------------------------------------------------
c
      integer*4 k
      real*8    eps
      real*8    tk(k),tk1(*)

      integer*4 i
      real*8    ep,low,up,val
c
c disturbance of monotony of two successive rows of eigenvalues
c given in and tk1.
c tk1 has one eigenvalue more than tk (k+1 and k resp.).
c
      up = tk(1)
      val = tk1(1)
      eps = 0.0d0
      if (val.le.up) go to 10
      eps = dabs(val-up)/(abs(val)+1.0d0)
   10 continue
      do 30 i=2,k
         low = up
         val = tk1(i)
         up = tk(i)
         if (val.ge.low) go to 20
         ep =dabs(val-low)/(dabs(val)+1.0d0)
         if (ep.gt.eps) eps = ep
   20    if (val.le.up) go to 30
         ep =dabs(up-val)/(dabs(val)+1.0d0)
         if (ep.gt.eps) eps = ep
   30 continue
      val = tk1(k+1)
      if (val.ge.up) go to 40
      ep =dabs(val-up)/(dabs(val)+1.0d0)
      if (ep.gt.eps) eps = ep
   40 continue
c
c eps is the maximum of the relative distance in case of disturbance
c of the monotony criterium.
c
      return
      end
c
c=======================================================================
c
      subroutine scanmp(m,k,k1,tk,tk1,eps)
c
c-----------------------------------------------------------------------
c subroutines required:
c   mltplt
c
c modifications by a. van den boogaart
c last update: 11-12-1990
c-----------------------------------------------------------------------
c
      integer*4 m,k,k1
      real*8    eps
      real*8    tk(m,2),tk1(m,2)

      integer*4 kold,k1old
      real*8    epsold
c
c with the value of eps, two rows (each apart) are scanned for
c multiplet this may deliver a new value of eps and the process is 
c repeated as long as eps changes.
c this results in two multiplet-free rows.
c
c m is the number of rows in tk and tk1.
c k is the number of intervals in tk.
c k1 is the number of intervals in tk1.
c
   10 continue
      kold = k
      k1old = k1
      epsold = eps

      call mltplt(m,k,epsold,tk)

      call mltplt(m,k1,eps,tk1)

      eps = dmax1(epsold,eps)
      if ((k.ne.kold) .or. (k1.ne.k1old)) go to 10

      return
      end
c
c=======================================================================
c      
      subroutine mltplt(m,nint,eps,a)
c
c-----------------------------------------------------------------------
c this subprogram is portable to other types of computer
c
c modifications by a. van den boogaart
c last update: 11-12-1990
c-----------------------------------------------------------------------
c
      integer*4 m,nint
      real*8    eps
      real*8    a(m,2)

      integer*4 i,j
      real*8    epsm,epsmul
c
c the part a(nint,2) contains the multiplet-intervals.
c
      epsm = eps
      j = 1

      do 20 i=2,nint
         if (dabs(a(i,1)-a(j,2))/(1.0d0+dabs(a(i,1))).gt.eps) go to 10
         a(j,2) = a(i,2)
c
c the upperbound of the j-th multiplet is updated.
c
         epsmul = dabs(a(j,2)-a(j,1))/(1.0d0+dabs(a(j,1)))
         if (epsmul.gt.epsm) epsm = epsmul
         go to 20
   10    j = j + 1
         a(j,1) = a(i,1)
         a(j,2) = a(i,2)
   20 continue
      nint = j
      eps = epsm

      return
      end
c
c=======================================================================
c
      subroutine comint(m,tk,tk1,k,k1,eps,lower,upper,jdef,div,jdiv)
c
c-----------------------------------------------------------------------
c subroutines required:
c   epsspn
c
c modifications by a. van den boogaart
c last update: 11-12-1990
c-----------------------------------------------------------------------
c
      integer*4 m,k,k1,jdef,jdiv
      real*8    eps
      real*8    tk(m,2),tk1(m,2)
      logical*4 lower,upper,div

      integer*4 idis,j,j1
      real*8    e,eold,f,fold
      logical*4 disjct

c
c comparision of two rows of intervals given in tk and tk1.
c if an interval in one row is close to an interval in the other row,
c than the span of both intervals is recorded in tk.
c so the values of tk have been overwritten.
c
c lower is .true. : convergence at the lower end of the spectrum.
c up    is .true. : convergence at the upper end of the spectrum.
c div   is .true. : a hole in the spectrum is detected.
c jdef is the number of the computed eigenvalue intervals.
c jdiv gives the number of eigenvalue intervals at the lower end.
c
      div = .false.
      idis = 0
      lower = .false.
      upper = .false.
      fold = dmin1(tk1(1,1),tk(1,1))
      eold = (1.0d0-dsign(0.1d0,fold))*fold
      fold = eold
      j = 1
      j1 = 1
      jdef = 0
   10 continue

      call epsspn(tk(j,1),tk(j,2),tk1(j1,1),tk1(j1,2),eps,e,f,disjct)

      if (disjct) go to 50
      idis = 0
      if (j.eq.1 .or. j1.eq.1) lower = .true.
      if (j.eq.k .or. j1.eq.k1) upper = .true.
      if (f.le.fold) go to 20
      if (e.gt.eold) then
         jdef=jdef+1
      endif
      tk(jdef,1) = e
      eold = e
      tk(jdef,2) = f
      fold = f
   20 if (tk(j,2).le.tk1(j1,2)) go to 40
   30 continue
      j = j + 1
      if (j.le.k) go to 10
      j = j - 1
      j1 = j1 + 1
      if (j1.gt.k1) go to 80
   40 continue
      j1 = j1 + 1
      if (j1.le.k1) go to 10
      j1 = j1 - 1
      j = j + 1
      if (j.gt.k) go to 80
   50 continue
c
c if successively 6 disjunct intervals are met, a hole in the
c spectrum is supposed.
c
      idis = idis + 1
      if (.not. div .and. (idis.eq.6)) go to 60
      go to 70
   60 jdiv = jdef
      div = .true.
   70 if (tk1(j1,1).gt.tk(j,1)) go to 30
      go to 40
   80 continue

      if (.not.(lower .and. (.not.div))) return
      div = .true.
      jdiv = jdef

      return
      end
c
c=======================================================================
c      
      subroutine epsspn(a,b,c,d,eps,e,f,disjct)
c
c-----------------------------------------------------------------------
c modifications by a. van den boogaart
c last update: 11-12-1990
c-----------------------------------------------------------------------
c
      real*8    a,b,c,d,eps,e,f
      logical*4 disjct
 
      real*8    epss
c
c if two intervals [a,b] and [c,d] are relatively close with respect
c to eps, [e,f] gives the span of both, and disjct is set to .false.
c if they are not close, disjct is set to .true..
c
      disjct = .true.
      epss = 0.0d0
      if (a.le.c) go to 10
      if (d.ge.a) go to 20
      if (dabs(d-a)/(1.0d0+dabs(d)).le.eps) go to 20
      return
   10 if (c.le.b) go to 20
      if (dabs(c-b)/(1.0d0+dabs(c)).le.eps) go to 20
      return
   20 continue
      disjct = .false.
      e = dmin1(a,c)
      f = dmax1(b,d)
      epss =( f - e)/(1.0d0+dabs(e))
      if (epss.gt.eps) eps = epss

      return
      end
c
c=======================================================================
c      
      subroutine dgtsl(n,b,c,d,e,info)
c
c-----------------------------------------------------------------------
c modifications by a. van den boogaart
c last update: 11-12-1990
c-----------------------------------------------------------------------
c
      integer*4 n,info
      real*8    b(*),c(*),d(*),e(*)

      integer*4 k,kb,kp1,nm1,nm2
      real*8    t

      info =0
      c(1)=d(1)
      nm1=n-1
      if (nm1 .lt. 1) goto 40
      d(1)=e(1)
      e(1)=0.0d0
      e(n)=0.0d0
      do 30 k=1,nm1
         kp1=k+1
         if (dabs(c(kp1)) .lt. dabs(c(k))) goto 10
         t=c(kp1)
         c(kp1)=c(k)
         c(k)=t
         t=d(kp1)
         d(kp1)=d(k)
         d(k)=t
         t=e(kp1)
         e(kp1)=e(k)
         e(k)=t
         t=b(kp1)
         b(kp1)=b(k)
         b(k)=t
   10    continue
         if (c(k) .ne. 0.0d0) goto 20
         info=k
         goto 100
   20    continue
         t=-c(kp1)/c(k)
         c(kp1)=d(kp1)+t*d(k)
         d(kp1)=e(kp1)+t*e(k)
         e(kp1)=0.0d0
         b(kp1)=b(kp1)+t*b(k)
   30 continue
   40 continue
      if (c(n) .ne. 0.0d0) goto 50
      info=n
      goto 90
   50 continue
      nm2=n-2
      b(n)=b(n)/c(n)
      if (n .eq. 1) goto 80
      b(nm1)=(b(nm1)-d(nm1)*b(n))/c(nm1)
      if (nm2 .lt. 1 ) goto 70
      do 60 kb=1,nm2
         k=nm2 -kb +1
         b(k)=(b(k)-d(k)*b(k+1)-e(k)*b(k+2))/c(k)
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue

      return
      end
c
c=======================================================================



       subroutine caldim(m,n,ndp)
       integer m,n,ndp,ind, mn1


       mn1 = m + n - 1
       ind = 0
       ndp = 1
       do while (ndp .lt. mn1)
          ind = ind+1
          ndp = ndp * 2
       end do

 
       end 


c=======================================================================



      Program hlsvdmain
c
c     fort77 hlsvdmain.f lanczos.f hlsubrs.f lapack... -o hlsvd.prg
c     g77 hlsvdmain.f lanczos.f hlsubrs.f -llapack -lblas -o hlsvd.prg
c
c     31-08-1999: Lanczos svd. Rest with Lapack.
c     Last change: 31-08-1999, ~12 h
c     01-09-1999: Frequency sort added in Main.
c     19-09-2000: Comment adapted for Sabine
c
c
c     List of variables in MAIN
c     =========================
c     signal  = complex-valued signal
c     ndp     = number of data points
c     Lrow    = number of rows in Hankel matrix, < ndp/2
c     mcoL    = number of columns in Hankel matrix, > ndp/2 + 1
c     kuser   = number of sinusoids initially wished by user
c     kfit    = number of sinusoids fitted (sinvals found by Lanczos)
c     lsinval = lanczos sinvals of the Hankel data-matrix
c     ampli   = estimated amplitudes
c     fase    = estimated phases
c     damp    = estimated dampings
c     freq    = estimated frequencies
c
      implicit none
c
      integer    kmax,ndpmax
      parameter  (kmax=100,ndpmax=12512)
      integer    i                                   ! running integer
      integer    kfit,kuser,Lrow,mcoL,ndp              ! parameters
      real*4     expon,test                            ! fftcheck
      real*8     Lsinval(kmax)                           ! singular values
      real*8     ampl(kmax),fase(kmax),damp(kmax), freq(kmax)
      complex*16 signal(ndpmax)                          ! signal
      real*8 datiin(2),stdev      
      character*40 nomefile

      integer    lwrk,lzwrk
      real*8     rwork(20*ndpmax+13*kmax+8*kmax**2+kmax*ndpmax)
      complex*16 zwork(20*ndpmax+7*kmax+2*kmax**2+2*kmax*ndpmax)
      complex*16 us(kmax,kmax),unit(kmax,kmax),zprime(kmax,kmax),
     +     zeta(ndpmax*kmax),cwork(2*kmax+64*(ndpmax+kmax)) 

      real*8        r1work(20*ndpmax+13*kmax+8*kmax**2+kmax*ndpmax)
      complex*16    z1work(84*ndpmax+100*kmax+5*kmax**2+4*kmax*ndpmax)
c
c     Some input & output:
c     Read: no. of components kfit, 
c           no. of data points ndp,
c           complex-valued signal x(ndp). 
c
c     Open files:
C       write(*,*)'Please, introduce the name of the file'
C       read(*,*)nomefile
C       write(*,*)'Model order'
C       read(*,*)kuser
      nomefile = 'press_cp0.txt'
      kuser = 20
      kfit=kuser  

      open     ( 8,file='hlsvdppro.out'  ,status='unknown')
c      open     ( 9,file='signal.dat',status='old')
      open     ( 9,file=nomefile,status='old')
c       open     ( 9,file='data1024.txt',status='old')
c
c     ******************************************************************
C      read (9,'(1F26.51X1F26.5)') datiin(1),datiin(2)
      read (9,*) datiin(1),datiin(2)
      ndp=int(real(datiin(1)))
C       print *, "ndp=", ndp
      stdev=datiin(2)

      do i = 1,ndp
C   10 read (9,'(d26.161Xd26.16)') signal(i)
            read (9,*) datiin(1),datiin(2)
            signal(i) = dcmplx(datiin(1),datiin(2))
C            print *, "signal(", i, ")=", signal(i)
      enddo
      close (9)
      
c
c     Adapt ndp to the mixed-radix fft used in the Lanczos-algorithm.
c     --------------------------------------------------------------
          expon = log(real(ndp)) / log(2.)
          test  = mod(100. * expon,100.)
c          if (test .ne. 0.) then
c              call check(ndp) ! In subrd.f; fftcheck.fil used too.
c          endif 
c

      mcol=int(ndp/2)
      lrow =ndp-mcol+1
c      Lrow = 512             ! 7-8-99 
c      mcol = ndp + 1 - Lrow !  ndp-kfit-9 > mcol > ndp/2
c
      lwrk  = Lrow+mcoL+13*kmax+8*kmax**2+32*Lrow+kmax*ndp
      lzwrk = Lrow+mcoL+32*Lrow+7*kmax+2+2*kmax**2
c
      call      hlsvdpro(signal,ndp,Lrow,mcoL,kuser,kfit,kmax,lsinval,
     +               ampl,fase,damp,freq,rwork(1),lwrk,zwork(1),lzwrk,
     +               zwork(1+lzwrk),zwork(1+lzwrk+lrow*2*kmax),
     +               us,unit,zprime,zeta,cwork) 
c
c
c     Write results to screen & disc
c     ------------------------------
c
      call sortfreq(kfit,ampl,fase,damp,freq)
c   
      print *,'freq'
      do 20 i = 1,kfit
         print *,freq(i)
   20 continue
      print *,'ampl'
      do 30 i = 1,kfit
         print *,ampl(i)
   30 continue
      print *,'damp'
      do 40 i = 1,kfit
         print *,damp(i)
   40 continue

c     TEMP DEBUG CHECK: hlsvdro call above is the same as :

c      call  hlsvdpw(signal,ndp,Lrow,mcoL,kuser,kfit,kmax,Lsinval,
c     +     ampl,fase,damp,freq,r1work,20*ndp+13*kmax+8*kmax**2+kmax*ndp,
c     +     z1work,84*ndp+100*kmax+5*kmax**2+4*kmax*ndp+20)

c      call sortfreq(kfit,ampl,fase,damp,freq)
c      print *,'freq'
c      do 50 i = 1,kfit
c         print *,freq(i)
c   50 continue
c      print *,'ampl'
c      do 60 i = 1,kfit
c         print *,ampl(i)
c   60 continue
c      print *,'damp'
c      do 70 i = 1,kfit
c         print *,damp(i)
c   70 continue
      close (8)
      stop 
      end
c
c     &&&&&&&&&&&&&&&&&&&&&&&  END OF MAIN  &&&&&&&&&&&&&&&&&&&&&&&&
c
