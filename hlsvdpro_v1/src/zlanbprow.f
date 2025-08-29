      subroutine zlanbprow( ndp, m, n, k0, k,  U, ldu, V, ldv, B, ldb,
     c     rnorm, doption, ioption, work, iwork,
     c      ierr,zwork,lambda,trlambda,planF,planB)
c
c
c
c     doption(1) = delta
c     doption(2) = eta
c     doption(3) = anorm
c
c     ioption(1) = cgs
c     ioption(2) = elr
c
c     work  :  double precision array
c     iwork :  integer array
c     (C) Rasmus Munk Larsen, Stanford, 1999
c


c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer m, n, k0, k, ldb, ldu, ldv, ierr
      integer ioption(*), iwork(*),ndp,PLANF,PLANB
      double precision rnorm,B(ldb,*), doption(*), work(*)
      complex*16  U(ldu,*),V(ldv,*), zwork(*), czero,s,cone,
     c            lambda(*),trlambda(*)
      external aprodw

c     %------------%
c     | Parameters |
c     %------------%
      logical DEBUG
      integer MGS
      double precision one, zero, FUDGE, kappa
      parameter(one = 1.0d0, zero = 0.0d0, FUDGE = 1.01)
c      parameter (DEBUG=.TRUE., MGS=0, kappa = 0.717)
      parameter (DEBUG=.FALSE., MGS=1, kappa = 0.717)
      parameter (czero=(0.0d0,0.0d0),cone=(1.0d0,0.0d0))
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,iv,iu,inu,imu,is,iidx,j0,iy1,iz1
      double precision eps,eps34,epsn2,epsn,sfmin,delta,eta,anorm
      double precision mumax,numax,alpha,beta,a1,b1,amax,anormest
      logical force_reorth,full_reorth,debito
      real t1,t2,t3


c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external zgetu0w,zreorth,zsafescal,zzero,izero,zcopy,zaxpy
      external dcompute_int,dupdate_nu,dupdate_mu,dzero

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch,dznrm2,dlapy2
      double complex zdotc
      external dznrm2,zdotc
      external dlamch,dlapy2

c-------------------- Here begins executable code ---------------------
      call second(t1)

c     %---------------------------------%
c     | Set machine dependent constants |
c     %---------------------------------%
      eps = dlamch('e')
      eps34 = eps**(3d0/4d0)
      epsn = dble(max(m,n))*eps
      epsn2 = sqrt(dble(max(m,n)))*eps
      sfmin = dlamch('s')

c      kappa =1.0d0/dsqrt(2.0d0)
      debito = .false.

c     %------------------------%
c     | Set default parameters |
c     %------------------------%
      if (doption(1).lt.zero) then
         delta = sqrt(eps/k)
      else
         delta = doption(1)
      endif
      if (doption(2).lt.zero) then
         eta = eps34/sqrt(dble(k))
      else
         eta = doption(2)
      endif
      if (delta.le.eta .or. delta.eq.zero) then
         full_reorth = .true.
      else
         full_reorth = .false.
      endif
      if (doption(3).gt.zero) then
         anorm = doption(3)

      else if (k0.gt.0) then
         anorm = dlapy2(B(1,1),B(1,2))
         if (anorm.le.zero) then
            ierr = -1
            goto 9999
         endif
      else
         anorm = zero
      endif
c               delta =3.107106895105165e-09
c              eta = 3.792855096563922e-13

c     %---------------------%
c     | Get starting vector |
c     %---------------------%


      if (rnorm .eq. zero) then
c         call zgetu0('n',m, n, k0, 3, U(1,k0+1) , rnorm, U,ldu, aprod,
c     c         ierr, ioption(1), anormest,zwork,lambda,trlambda)
c         anorm = max(anorm,anormest)
         call zgetu0w('n',ndp,m, n, k0, 3, U(1,k0+1) , rnorm, U,ldu,
     c         ierr, ioption(1), anormest,zwork,lambda,
     c           trlambda,planF,planB) ! it was zwork instead of zwork(is) ! diana

ccc         call zgetu0('n',m,n,0,1,U,rnorm,U,ldu,
ccc     c        ierr,ioption(1),anorm,zwork(iwrk),lambda,
ccc     c        trlambda,wsave)

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         anorm = max(anorm,anormest)
      endif

c     %------------------------------%
c     | Set pointers into work array |
c     %------------------------------%
      iu = 1
      iv = iu+m
      imu = iv+n
c      iv = iu+ndp
c      imu = iv+ndp
      inu = imu+k+1
      is = inu+k+1
      iz1 = is+3*ndp
      iy1 = iz1+ndp
      iidx = 1
      call dzero(m+n+max(m,n)+2*k+2,work,1)
      call zzero(m+n+max(m,n)+2*k+2,zwork,1)
      call izero(k,iwork,1)

c     %---------------------------%
c     | Prepare Lanczos iteration |
c     %---------------------------%
c      write(*,*)k0,'******************************** k0'

      if (k0.eq.0) then
         amax = zero
         alpha = zero
         beta = rnorm
         force_reorth = .false.
c     %-----------------------------------------------%
c     | Compute ||A x|| / ||x|| for a random vector x |
c     | to make it less likely that ||A|| is grossly  |
c     | underestimated.                               |
c     %-----------------------------------------------%
         if (n.gt.m) then
c            call zgetu0('n',m,n,0,1,zwork(iu),s,U,ldu,aprod,
c     c           ierr,ioption(1),anormest,zwork(is),lambda,
c     c           trlambda,wsave)
            call zgetu0w('n',ndp,m,n,0,1,zwork(iu),s,U,ldu,
     c           ierr,ioption(1),anormest,zwork(is),lambda,
     c           trlambda,planF,planB)
         else

c            call zgetu0('y',m,n,0,1,zwork(iv),s,V,ldv,aprod,
c     c           ierr,ioption(1),anormest,zwork(is),lambda,
c     c           trlambda,wsave)
            call zgetu0w('t',ndp,m,n,0,1,zwork(iv),s,V,ldv,
     c           ierr,ioption(1),anormest,zwork(is),lambda,
     c           trlambda,planF,planB)
         endif

         anorm = max(anorm,FUDGE*anormest)
         j0 = 1
         if (beta.ne.zero) then
            call donothing()
            call zsafescal(m,beta,U(1,k0+1))
         endif
         call zcopy(m,U(1,k0+1),1,zwork(iu),1)

         work(imu) = one
         work(inu) = one
         zwork(imu) = cone
         zwork(inu) = cone


      else
         alpha = B(k0,1)
         call zcopy(n,V(1,k0),1,zwork(iv),1)
         beta = rnorm
         force_reorth = .true.


         if (k0.lt.k .and. beta*delta.lt.anorm*eps) then
            full_reorth = .true.
            ierr = k0
         endif
         iwork(iidx) = 1
         iwork(iidx+1) = k0
         iwork(iidx+2) = k0+1

         call second(t2)


c         if (MGS.eq.0) then
c            call zreorth(m,k0,U,ldu,U(1,k0+1),s,iwork(iidx),kappa,
c     c           zwork(is),ioption(1))
c         else
c            call zreorth2(m,k0,U,ldu,U(1,k0+1),s,iwork(iidx))
            call zreorth2(m,k0,U,ldu,U(1,k0+1),rnorm,iwork(iidx))
c         endif
         call second(t3)
         treorthu = treorthu+(t3-t2)

         call dset_mu(k0,work(imu),iwork(iidx),epsn2)


c     %--------------------------------------%
c     | Estimate ||B||_2^2 as ||B^T * B||_1  |
c     %--------------------------------------%

c
c         beta = s*beta
c

c*******************************
c         beta = dble(s)*beta
c*******************************
         B(k0,2) = beta
         amax = zero
         do j=1,k0
            amax = max(amax,B(j,1),B(j,2))
            if (j.eq.1) then
               anorm = max(anorm,FUDGE*alpha)
            else if (j.eq.2) then
               a1 = B(1,2)/amax
               a1 = FUDGE*amax*sqrt((B(1,1)/amax)**2 + a1**2 +
     c              B(2,1)/amax*a1)
               anorm = max(anorm,a1)
            else
               a1 = B(j-1,1)/amax
               b1 = B(j-1,2)/amax
               a1 = FUDGE*amax*sqrt( a1**2 + b1**2 +
     c              a1*B(j-2,2)/amax + B(j,1)/amax*b1)
               anorm = max(anorm,a1)
            endif
            work(imu+j-1) = epsn2
            work(inu+j-1) = epsn2
         enddo
         j0 = k0+1
         call zcopy(m,U(1,k0+1),1,zwork(iu),1)
      endif
      numax = zero
      mumax = zero



c     %-------------------------------------------%
c     | Start Lanczos bidiagonalization iteration |
c     %-------------------------------------------%
      if (DEBUG) then
         call flush(6)
      endif

      do j=j0,k


         if (DEBUG .and. mod(j,10).eq.0) then
c            write (*,*) 'alpha, beta = ',alpha,beta
            call flush(6)
         endif

         if (j.eq.1) then
            call second(t2)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            call aprodw('t',ndp,m,n,zwork(iu),zwork(iv),zwork(iz1),
     c                   zwork(iy1),lambda,trlambda,
     c                   planF,planB )


c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            call second(t3)
            tmvopx = tmvopx + (t3-t2)
            nopx = nopx+1
            alpha = dznrm2(n,zwork(iv),1)
            anorm = max(anorm,FUDGE*alpha)

         else

c     %---------------------------------------------%
c     | alpha_{j} v_{j} = A'*u_{j} - beta_{j} v_{j} |
c     %---------------------------------------------%
            call second(t2)
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            call aprodw('t',ndp,m,n,zwork(iu),zwork(is),zwork(iz1),
     c                   zwork(iy1),lambda,trlambda,
     c                     planF,planB )

c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            call second(t3)
            tmvopx = tmvopx + (t3-t2)
            nopx = nopx+1

            do i=0,n-1
               zwork(iv+i) = zwork(is+i) - beta*zwork(iv+i)
            enddo

c     %------------------------------------%
c     | Extended local reorthogonalization |
c     %------------------------------------%
            alpha = dznrm2(n,zwork(iv),1)
            call second(t2)
            if (j.gt.1 .and. ioption(2).gt.0 .and.
     c           alpha.lt.kappa*beta) then
               do i=1,2
                  s = zdotc(n,V(1,j-1),1,zwork(iv),1)

                  call zaxpy(n,-s,V(1,j-1),1,zwork(iv),1)
                  if (beta .ne. zero) then
                     beta = beta + dble(s)
                     B(j-1,2) = beta
                  endif
                  s = dznrm2(n,zwork(iv),1)
c                  if (dsqrt(dble(conjg(s)*s)) .ge. kappa*alpha) goto 10
                   if (dble(s) .ge. kappa*alpha) goto 10
               enddo
 10            work(inu+j-2) = eps
c               alpha =dsqrt( dble(conjg(s)*s))
               alpha =dble(s)
            endif
            call second(t3)
            telrv = telrv + (t3-t2)

            B(j,1) = alpha
            amax = max(amax,alpha)
c     %---------------------------%
c     | Update estimate of ||A||_2 |
c     %---------------------------%
            if (j.eq.2) then
               a1 = B(1,2)/amax
               a1 = FUDGE*amax*sqrt((B(1,1)/amax)**2 + a1**2 +
     c              B(2,1)/amax*a1)
            else
               a1 = B(j-1,1)/amax
               b1 = B(j-1,2)/amax
               a1 = FUDGE*amax*sqrt( a1**2 + b1**2 +
     c              a1*B(j-2,2)/amax + B(j,1)/amax*b1)
            endif
            anorm = max(anorm,a1)
         endif

c     %--------------------------%
c     | Update the nu recurrence |
c     %--------------------------%
         if (.not.full_reorth .and. alpha.ne.zero) then
            call dupdate_nu(numax,work(imu),work(inu),j,
     c           B(1,1),B(1,2),anorm)
         endif

c     %------------------------------%
c     | Reorthogonalize if necessary |
c     %------------------------------%
         if ( .true. .and.
     c        (full_reorth .or. numax.gt.delta .or. force_reorth)
     c        .and. alpha.ne.zero) then
            if (full_reorth .or. eta.eq.zero) then
               iwork(iidx) = 1
               iwork(iidx+1) = j-1
               iwork(iidx+2) = j
            else if (.not. force_reorth) then
               call dcompute_int(work(inu),j-1,delta,eta,iwork(iidx))
            endif
            call second(t2)

            if (MGS.eq.0) then
               call zreorth(n,j-1,V,ldv,zwork(iv),alpha,iwork(iidx),
     c              kappa,zwork(is),ioption(1))
            else
               call zreorth2(n,j-1,V,ldv,zwork(iv),alpha,iwork(iidx))
            endif
            call second(t3)
            treorthv = treorthv+(t3-t2)

            call dset_mu(j-1,work(inu),iwork(iidx),epsn2)
            numax = eta
            if (force_reorth) then
               force_reorth = .false.
            else
               force_reorth = .true.
            endif
         endif




c     %-----------------------------------------------%
c     | Check whether an invariant subspace was found |
c     %-----------------------------------------------%
         if (alpha .lt. anorm*epsn .and. j.lt.k) then
            rnorm = alpha
            alpha = zero
c     %------------------------------------------------%
c     | Try to build and orthogonal subspace, starting |
c     | with a random vector.                          |
c     %------------------------------------------------%
          call zgetu0w('t',ndp, m, n, j-1, 3, zwork(iv), alpha, V, ldv,
     c             ierr,ioption(1),anormest,
     c           zwork(is),lambda,trlambda,planF,planB)
            if (alpha .eq. zero) then
               k = j-1
               ierr = -j
               goto 9999
            else
               nrstrt = nrstrt + 1
               call donothing()
               call zsafescal(n,alpha,zwork(iv))
               alpha = zero
               force_reorth = .true.
               if (delta.gt.zero) then
                  full_reorth = .false.
               endif
            endif
         else if (j.gt.1 .and. .not. full_reorth .and. j.lt.k .and.
     c           (delta*alpha .lt. anorm*eps)) then
c            full_reorth = .true.
c            full_reorth = .true.
            ierr = j
         endif
         B(j,1) = alpha

         if (alpha.ne.zero) then
            call donothing()
            call zsafescal(n,alpha,zwork(iv))
         endif
         call zcopy(n,zwork(iv),1,V(1,j),1)

c     %------------------------------------------------%
c     | beta_{j+1} u_{j+1} = A*v_{j} - alpha_{j} u_{j} |
c     %------------------------------------------------%
         call second(t2)


         call aprodw('n',ndp,m,n,zwork(iv),zwork(is),zwork(iz1),
     c                   zwork(iy1),lambda,trlambda,
     c                       planF,planB )

         call second(t3)
         tmvopx = tmvopx + (t3-t2)
         nopx = nopx+1

         do i=0,m-1
            zwork(iu+i) = zwork(is+i) - alpha*zwork(iu+i)

         enddo



c     %------------------------------------%
c     | Extended local reorthogonalization |
c     %------------------------------------%
         call second(t2)
         do i=1,ioption(2)
            s = zdotc(m,U(1,j),1,zwork(iu),1)


            call zaxpy(m,-s,U(1,j),1,zwork(iu),1)
            if (alpha .ne. zero) then
               alpha = alpha + dble(s)

               B(j,1) = alpha
            endif
            work(imu+j-1) = eps
         enddo
         call second(t3)
         telru = telru + (t3-t2)

         beta= dznrm2(m,zwork(iu),1)

         B(j,2) = beta

         amax = max(amax,beta)

c     %---------------------------%
c     | Update estimate of ||A||_2 |
c     %---------------------------%
         if (j.eq.1) then
            a1 = dlapy2(B(1,1), B(1,2))
         else
            a1 = B(j,1)/amax
            a1 = amax*sqrt(a1**2 + (B(j,2)/amax)**2 +
     c           a1*B(j-1,2)/amax)
         endif
         anorm=max(anorm,a1)

c     %--------------------------%
c     | Update the mu recurrence |
c     %--------------------------%
         if (.not.full_reorth .and. beta.ne.zero) then
            call dupdate_mu(mumax,work(imu),work(inu),j,B(1,1),
     c           B(1,2),anorm)
         endif

c     %--------------------------------------%
c     | Reorthogonalize u_{j+1} if necessary |
c     %--------------------------------------%
         if ( .true. .and.
     c        (full_reorth .or. mumax.gt.delta .or. force_reorth)
     c        .and. beta.ne.zero) then
            if (full_reorth .or. eta.eq.zero) then
               iwork(iidx) = 1
               iwork(iidx+1) = j
               iwork(iidx+2) = j+1
            else if (.not. force_reorth) then
               call dcompute_int(work(imu),j,delta,eta,iwork(iidx))
            else
               do i=1,j+1
                  if (iwork(iidx+i-1).eq.j) then
                     iwork(iidx+i-1) = j+1
                     goto 25
                  endif
               enddo
            endif

 25         call second(t2)

c 25         if (MGS.eq.0) then
c 25       call zreorth(m,j,U,ldu,zwork(iu),beta,iwork(iidx),
c     c              kappa, zwork(is),ioption(1))
c            else
               call zreorth2(m,j,U,ldu,zwork(iu),beta,iwork(iidx))
c            endif
            call second(t3)
            treorthu = treorthu+(t3-t2)

            call dset_mu(j,work(imu),iwork(iidx),epsn2)
            mumax = eta
            if (force_reorth) then
               force_reorth = .false.
            else
               force_reorth = .true.
            endif
         endif




c     %-----------------------------------------------%
c     | Check whether an invariant subspace was found |
c     %-----------------------------------------------%
         if (beta .lt. anorm*epsn .and. j.lt.k) then
            rnorm = beta
            beta = zero
c     %------------------------------------------------%
c     | Try to build an orthogonal subspace, starting |
c     | with a random vector.                          |
c     %------------------------------------------------%
            call zgetu0w('n',ndp, m, n, j, 3, zwork(iu), beta, U, ldu,
     c             ierr,ioption(1),anormest,zwork(is),
     c             lambda,trlambda,planF,planB)
            if (beta .eq. zero) then
               k = j
               ierr = -j
               goto 9999
            else
               nrstrt = nrstrt + 1
               call donothing()
               call zsafescal(m,beta,zwork(iu))
               beta = zero
               force_reorth = .true.
               if (delta .gt. zero) then
                  full_reorth = .false.
               endif
            endif
         else if (.not.full_reorth .and. j.lt.k .and.
     c           (delta*beta .lt. anorm*eps)) then
c            full_reorth = .true.

            debito=.true.
            ierr = j
         endif

         B(j,2) = beta
c         write (*,*) 'beta1 = ',beta
         if (beta.ne.zero .and. beta.ne.one) then
            call donothing()
            call zsafescal(m,beta,zwork(iu))
         endif
         call zcopy(m,zwork(iu),1,U(1,j+1),1)


         rnorm = beta
         call second(t2)
      enddo
 9999 doption(3) = anorm
      call second(t2)
      tlanbpro = tlanbpro + (t2-t1)
 100  format(1i6,6e12.3)
 101  format(1i6,4e12.3)
      return
      end

c
c**********************************************************************
c

      subroutine dset_mu(k,mu,index,val)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer k,index(*)
      double precision mu(*), val

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,p,q

      i=1
      do while(index(i).le.k .and. index(i).gt.0)
         p = index(i)
         q = index(i+1)
         do j=p,q
            mu(j) = val
         enddo

         i = i+2
      enddo
      end
c
c**********************************************************************
c
      subroutine dcompute_int(mu,j,delta,eta,index)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer j,index(*)
      double precision mu(*)
      double precision delta,eta

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,k,s,ip
      real t1,t2



      call second(t1)

      if (delta.lt.eta) then
         return
      endif

c      write (*,*) 'Enter compute_int, j=',j
      ip = 0
      index(1) = 0
      i=0
      do while(i.lt.j)
c     find the next mu(k), k>i where abs(mu(k)) > delta
         do k=i+1,j
c            write(*,*) i,k,abs(mu(k))
            if (abs(mu(k)).gt.delta) goto 10
         enddo
         goto 40

c     find smallest i<k such that for all j=i,..,k, m(j) >= eta
 10      do s=k,max(i,1),-1
            if (abs(mu(s)).lt.eta) goto 20
         enddo
 20      ip= ip+1
         index(ip) = s+1
         do i=s+1,j
            if (abs(mu(i)).lt.eta) goto 30
         enddo
 30      ip= ip+1
         index(ip) = i-1
      enddo
 40   ip = ip+1
      index(ip) = j+1
c      write (*,*) 'i  index(i)'
c      do i=1,j
c         write(*,*) i,index(i)
c         if (index(i).gt.j) goto 50
c      enddo
c 50   write (*,*) 'Exit compute_int, ip=',ip
      call second(t2)
      tintv = tintv + (t2-t1)
      end
c
c**********************************************************************
c
      subroutine dupdate_mu(mumax,mu,nu,j,alpha,beta,anorm)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer j,formula
      parameter(formula=1)
      double precision mumax,eps1,eps,anorm
      double precision mu(*),nu(*),alpha(*),beta(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero, FUDGE
      parameter(one = 1.0d0, zero = 0.0d0, FUDGE = 1.01)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision d
      integer k
      real t1,t2

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch,dlapy2
      external dlamch,dlapy2

c      write (*,*) 'Enter update_mu'
      call second(t1)
      eps = dlamch('e')
      eps1 = 100*eps
      if (formula.eq.1) then
         if (j.eq.1) then
            d = eps1*(dlapy2(alpha(j), beta(j)) + alpha(1)) + eps1*anorm
            mu(1) = eps1/beta(1)
            mumax = abs(mu(1))
         else
            mu(1) = alpha(1)*nu(1)-alpha(j)*mu(1)
            d = eps1*(dlapy2(alpha(j), beta(j)) + alpha(1)) + eps1*anorm
            mu(1) = (mu(1) + dsign(d,mu(1))) / beta(j)
            mumax = abs(mu(1))
            do k=2,j-1
               mu(k) = alpha(k)*nu(k) +beta(k-1)*nu(k-1)-alpha(j)*mu(k)
               d = eps1*(dlapy2(alpha(j), beta(j)) +
     c              dlapy2(alpha(k), beta(k-1))) + eps1*anorm
               mu(k) = (mu(k) + dsign(d,mu(k))) / beta(j)
               mumax = max(mumax,abs(mu(k)))
            enddo
            mu(j) = beta(j-1)*nu(j-1)
            d = eps1*(dlapy2(alpha(j), beta(j)) +
     c           dlapy2(alpha(j), beta(j-1))) + eps1*anorm
            mu(j) = (mu(j) + sign(d,mu(j))) / beta(j)
            mumax = max(mumax,abs(mu(j)))
         endif
         mu(j+1) = one
      else
         if (j.eq.1) then
            mu(1) = eps1/beta(1)
            mumax = abs(mu(1))
         else
            mu(1) = alpha(1)*nu(1)-alpha(j)*mu(1)
            mu(1) = (mu(1) + dsign(eps1,mu(1))) / beta(j)
            mumax = abs(mu(1))
            do k=2,j-1
               mu(k) = alpha(k)*nu(k) + beta(k-1)*nu(k-1)-alpha(j)*mu(k)
               mu(k) = (mu(k) + dsign(eps1,mu(k))) / beta(j)
               mumax = max(mumax,abs(mu(k)))
            enddo
            mu(j) = beta(j-1)*nu(j-1)
            mu(j) = (mu(j) + sign(eps1,mu(j))) / beta(j)
            mumax = max(mumax,abs(mu(j)))
         endif
         mu(j+1) = one
      endif
      call second(t2)
      tupdmu = tupdmu + (t2-t1)
c      write (*,*) 'Exit update_mu'
      end
c
c**********************************************************************
c

      subroutine dupdate_nu(numax,mu,nu,j,alpha,beta,anorm)
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      integer j,formula
      parameter(formula=1)
      double precision numax,eps1,anorm
      double precision mu(*),nu(*),alpha(*),beta(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero, FUDGE
      parameter(one = 1.0d0, zero = 0.0d0, FUDGE = 1.01)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      double precision d,eps
      integer k
      real t1,t2

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch,dlapy2
      external dlamch,dlapy2

c      write (*,*) 'Enter update_nu'
      call second(t1)
      eps = dlamch('e')
      eps1 = 100*eps
      if (formula.eq.1) then
         if (j.gt.1) then
            numax = zero
            do k=1,j-1
               nu(k) = beta(k)*mu(k+1) + alpha(k)*mu(k) -beta(j-1)*nu(k)
               d = eps1*(dlapy2(alpha(k),beta(k)) +
     c              dlapy2(alpha(j),beta(j-1))) + eps1*anorm
               nu(k) = (nu(k) + dsign(d,nu(k))) / alpha(j)
               numax = max(numax,abs(nu(k)))
            enddo
            nu(j) = one
         endif
      else
         if (j.gt.1) then
            numax = zero
            do k=1,j-1
               nu(k) = beta(k)*mu(k+1) + alpha(k)*mu(k) -beta(j-1)*nu(k)
               nu(k) = (nu(k) + dsign(eps1,nu(k))) / alpha(j)
               numax = max(numax,abs(nu(k)))
            enddo
            nu(j) = one
         endif
      endif
      call second(t2)
      tupdnu = tupdnu + (t2-t1)
c      write (*,*) 'Exit update_nu'
      end

      subroutine donothing()

c     Added by the Vespa team to address ticket #49:
c     http://scion.duhs.duke.edu/vespa/analysis/ticket/49
c
c     It writes "" to /dev/null which should do nothing but apparently
c     it has an important side effect. (Read the ticket.)

      character*255 envname, envvalue, bitbucket
      integer envlength, i

c     Determine the location of the bit bucket based on OS. The
c     OS environment variable exists on Windows but not elsewhere.
      envname = "OS"
      call get_environment_variable(envname, envvalue, envlength)

c     Convert OS name (if it exists) to lower case
      do i=1, envlength
         if (LGE(envvalue(i:i),'A') .AND. LLE(envvalue(i:i),'Z')) THEN
            envvalue(i:i) = ACHAR(IACHAR(envvalue(i:i)) + 32)
         endif
      enddo

      i = index(envvalue, "windows")
      if (i.ne.0) then
c        We're on Windows
         bitbucket = "NUL"
      else
c        We're on *nix
         bitbucket = "/dev/null"
      endif

      open(10, file=bitbucket, status="old")
      write(10, "()")
      close(10)

      end

