
      subroutine zlansvdw(jobu,jobv,ndp,m,n,k,kmax,U,ldu,Sigma,bnd,V,
     c     ldv,tolin,work,lwork,iwork,doption,ioption,info,
     c       zwork,lambda,trlambda, planF,planB)


c     DLANSVD: Compute the leading singular triplets of a large and 
c     sparse matrix A by Lanczos bidiagonalization with partial
c     reorthogonalization.
c
c     Parameters:
c
c     JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
c           Otherwise the array U is not touched.
c     JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular 
c           vectors. Otherwise the array V is not touched.
c     M: INTEGER. Number of rows of A.
c     N: INTEGER. Number of columns of A.
c     K: INTEGER. Number of desired singular triplets. K <= MIN(KMAX,M,N)
c     KMAX: INTEGER. maximal number of iterations / maximal dimension of
c           Krylov subspace.
c     APROD: Subroutine defining the linear operator A. 
c            APROD should be of the form:
c
c           SUBROUTINE DAPROD(TRANSA,M,N,X,Y,DPARM,IPARM)
c           CHARACTER*1 TRANSA
c           INTEGER M,N,IPARM(*)
c           DOUBLE PRECISION X(*),Y(*),DPARM(*)
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y == A * X.
c           If TRANSA.EQ.'T' then the function should compute the matrix-vector
c           product Y = A^T * X.
c           The arrays IPARM and DPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.
c     U(LDU,K): DOUBLE PRECISION array. On return the first K columns of U
c               will contain approximations to the left singular vectors 
c               corresponding to the K largest singular values of A.
c               On entry the first column of U contains the starting vector
c               for the Lanczos bidiagonalization. A random starting vector
c               is used if U is zero.
c     LDU: INTEGER. Leading dimension of the array U. LDU >= M.
c     SIGMA(K): DOUBLE PRECISION array. On return Sigma contains approximation
c               to the K largest singular values of A.
c     BND(K)  : DOUBLE PRECISION array. Error estimates on the computed 
c               singular values. The computed SIGMA(I) is within BND(I)
c               of a singular value of A.
c     V(LDV,K): DOUBLE PRECISION array. On return the first K columns of V
c               will contain approximations to the right singular vectors 
c               corresponding to the K largest singular values of A.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     TOLIN: DOUBLE PRECISION. Desired relative accuracy of computed singular 
c            values. The error of SIGMA(I) is approximately 
c            MAX( 16*EPS*SIGMA( 1), TOLIN*SIGMA(I) )
c     WORK(LWORK): DOUBLE PRECISION array. Workspace. Dimension must be
c                at least M+7*KMAX+2 + 2*KMAX**2 + MAX(M,N,4*KMAX+4))
c     LWORK: INTEGER. Dimension of WORK. Should be at least
c                M+7*KMAX+2 + 2*KMAX**2 + MAX(M,N,4*KMAX+4))
c     IWORK: INTEGER array. Integer workspace. Dimension must be 
c            at least 2*K+1.
c     DOPTION: DOUBLE PRECISION array. Parameters for LANBPRO.
c        doption(1) = delta. Level of orthogonality to maintain among
c          Lanczos vectors.
c        doption(2) = eta. During reorthogonalization, all vectors with
c          with components larger than eta along the latest Lanczos vector
c          will be purged.
c        doption(3) = anorm. Estimate of || A ||.
c
c     IOPTION: INTEGER array. Parameters for LANBPRO.
c        ioption(1) = CGS.  If CGS.EQ.1 then reorthogonalization is done
c          using iterated classical GRAM-SCHMIDT. IF CGS.EQ.0 then 
c          reorthogonalization is done using iterated modified Gram-Schmidt.
c        ioption(2) = ELR. If ELR.EQ.1 then extended local orthogonality is
c          enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
c
c     INFO: INTEGER. 
c         INFO = 0  : The K largest singular triplets were computed succesfully
c         INFO = J>0, J<K: An invariant subspace of dimension J was found.
c         INFO = -1 : K singular triplets did not converge within KMAX
c                     iterations.   
c
c
c     (C) Rasmus Munk Larsen, Stanford, 1999
c


c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      character*1 jobu,jobv
      integer m,n,k,kmax,lanmax,ldu,ldv,iwork(*),lwork,info,ioption(*)
      double precision Sigma(*),bnd(*),work(*)
      double precision tolin,doption(*)
      double complex U(ldu,*),V(ldv,*),zwork(*),
     c               lambda(*),trlambda(*)


c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero, FUDGE,DEBUG     
      parameter(one = 1.0d0, zero = 0.0d0, FUDGE = 1.01d0, DEBUG=0.0d0)
            
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,dj,jold,ibnd,ib,ib1,iwrk,ierr,ip,iq,neig,lwrk
      integer ncvt,ncc,ndp, planF,planB
      double precision eps,eps34,epsn2,epsn,sfmin,anorm,rnorm,tol

      real t0,t1,t2,t3

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external dzero,izero,dcopy,daxpy,dbdsqr,zgemmina,zzero
      
c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      double precision dlamch,dnrm2,ddot,dlapy2
      external dnrm2,ddot,lsame
      external dlamch,dlapy2,zgetu0w

c-------------------- Here begins executable code ---------------------
      
c     %-------------%
c     | Start timer |
c     %-------------%
      call second(t0)

 

c     %---------------------------------%
c     | Set machine dependent constants |
c     %---------------------------------%
      eps = dlamch('e')
      eps34 = eps**(3.0/4.0)
      epsn = dble(max(m,n))*eps/2.0
      epsn2 = sqrt(dble(max(m,n)))*eps/2.0
      sfmin = dlamch('s')
      

c     %--------------------------------%
c     | Guard against absurd arguments |
c     %--------------------------------%
      lanmax = min(n+1,m+1,kmax)
      tol = min(one,max(16.0*eps,tolin))
      anorm = zero

c     %------------------------------%
c     | Set pointers into work array |
c     %------------------------------%
      ibnd = 1
      ib = ibnd + lanmax+1
      ib1 = ib + 2*lanmax
      ip = ib1 + 2*lanmax
      iq = ip + (lanmax+1)**2
      iwrk = iq + lanmax**2

      lwrk = lwork-iwrk+1
      call dzero(7*lanmax + 2 + 2*lanmax**2,work,1)
      call zzero(7*lanmax + 2 + 2*lanmax**2,zwork,1)

c     %---------------------------------------------------------------%
c     | Set up random starting vector if none is provided by the user |
c     %---------------------------------------------------------------%
c      rnorm = dnrm2(m,U(1,1),1)
      rnorm =0.0d0
      if (rnorm.eq.zero) then
          call zgetu0w('n',ndp, m,n,0,1,U,rnorm,U,ldu,
     c        ierr,ioption(1),anorm,zwork(iwrk),lambda,
     c        trlambda,planF,planB)

c         call zgetu0('n',m,n,0,1,U,rnorm,U,ldu,
c     c        ierr,ioption(1),anorm,zwork(iwrk),lambda,
c     c        trlambda,wsave)
  
      endif

      info = 0
      neig = 0
      jold = 0
      j = min(k+max(8,k)+1,lanmax)
       
c     %------------------------------%
c     | Iterate until convergence... |
c     %------------------------------%
      do while (neig.lt.k)
c     %---------------------------------------------------%
c     | Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j} |
c     %---------------------------------------------------%
         call zlanbprow(ndp, m, n, jold, j,  U, ldu, V, ldv,
     c        work(ib),lanmax,rnorm,doption(1),ioption(1),
     c        work(iwrk), iwork,ierr,  
     c        zwork(iwrk), lambda,trlambda,planF,planB)
         jold = j
 

c     %---------------------------------------------%
c     | Compute and analyze SVD(B) and error bounds |
c     %---------------------------------------------%

         call dcopy(2*lanmax, work(ib),1,work(ib1),1)
         call dzero(j+1,work(ibnd),1)

         call second(t2)
         call dbdqr('N',j,work(ib1),work(ib1+lanmax),work(ibnd+j-1),
     c        work(ibnd+j),work(ip),lanmax+1)
         call dbdsqr('u',j,0,1,0,work(ib1),work(ib1+lanmax),work,1,
     c        work(ibnd),1,work,1,work(iwrk),info)

         call  second(t3)
         tbsvd = tbsvd + (t3-t2)
         nbsvd = nbsvd + 1

         if (j.gt.5) then
            anorm = work(ib1)
         else
            anorm = max(anorm,work(ib1))
         endif
         do i=1,j
            work(ibnd+i-1) = abs(rnorm*work(ibnd+i-1))
            work(ib1+lanmax+i-1) = work(ib1+lanmax+i-1)**2
         enddo

c     %---------------------------------------------%
c     | Refine error bounds using the "Gap theorem" |
c     %---------------------------------------------%

         call refinebounds(j,work(ib1+lanmax),work(ibnd),
     c        epsn*anorm,eps34)
         

c     %----------------------------------------------------%
c     | Determine the number of converged singular values  |
c     %----------------------------------------------------%
         do i=1,min(j,k)
            bnd(i) = work(ibnd+i-1)
         enddo
         i = 0
         neig = 0
         do while(i.lt.min(j,k))
            if (work(ibnd+i).le.tol*work(ib1+i)) then
               neig = neig + 1
               sigma(neig) = work(ib1+i)
               i = i+1
            else
               i = k
            endif
         enddo
        
c     %--------------------------------------------------%
c     | Test if an invariant subspace have been found or |
c     | the workspace has been exhausted.                |
c     %--------------------------------------------------%
         if (ierr.lt.0) then
            goto 50               
         endif
         if (j.ge.lanmax) then
            if (neig.lt.k) then
               info = -1
            endif
            goto 50
         endif

c     %----------------------------------------------------%
c     | Increase the dimension of the Krylov subspace.     |
c     | If any Ritz values have converged then try to      | 
c     | estimate the average number of iterations per      |
c     | converged Ritz value.                              |
c     | Else increase the dimension with 50%.              |
c     %----------------------------------------------------%
         if (neig.gt.1) then
            dj = min(j/2,((k-neig)*(j-6))/(2*neig+1))
            dj = min(100,max(2,dj))
         else
            dj = j/2
            dj = min(100,max(10,dj))
        endif
         j = min(j + dj,lanmax)
      enddo

 50   if (neig.ge.k) then
c     %-----------------------------------------%
c     | Calculate singular vectors if requested %
c     %-----------------------------------------%
         j = jold
         k = neig

         if (lsame(jobu,'y') .or. lsame(jobv,'y')) then
            call dcopy(2*lanmax, work(ib),1,work(ib1),1)

c     %----------------------------------------------------------------------%
c     | The bidiagonal SVD is computed in a two-stage procedure:             
c     |                                                   
c     | 1. Compute a QR-factorization M^T*B = [R; 0] of the (k+1)-by-k lower 
c     |    bidiagonal matrix B.
c     | 2. Compute the SVD of the k-by-k upper bidiagonal matrix 
c     |    R = P*S*Q^T. The SVD of B is then (M*P)*S*Q^T.
c     %----------------------------------------------------------------------%


c     %-----------------------------------------%
c     | M^T*B = [R; 0]
c     %-----------------------------------------%
C            do ii =1,j+1
C              do i=1,m
C              
C                write(22,'(2d24.16)')U(i,ii)
C             end do
C            end do   
C            write(*,*)' Matrice B',lanmax 
C            do i=1,j
C               write(*,'(2d24.16)')work(ib1+i-1),work(ib1+lanmax+i-1)
C            end do
C            pause 

            call dbdqr(jobu,j,work(ib1),work(ib1+lanmax),
     c           work(ibnd+j-1), work(ibnd+j),work(ip),lanmax+1)
 
Cc            do i=1,j+1
Cc               do ii= 1,j+1
Cc                  write(16,'(d24.16)')work(ip +(i-1)*(lanmax+1) +ii-1)
Cc               end do
Cc            end do
               

            if (lsame(jobu,'y')) then
               ncc = j+1
            else
               ncc = 0
            endif
            if (lsame(jobv,'y')) then
               ncvt = j
            else
               ncvt = 0
            endif
            
            do i=1,n
               work(iq+(i-1)*(lanmax+1)) = one
            enddo

c     THIS SHOULD POSSIBLY BE REPLACED BY A CALL TO THE FAST
c     DIVIDE-AND-CONQUER BIDIAGONAL SVD IN LAPACK 3 OR BY TRANSFORMING B
c     TO TRIDIAGONAL GOLUB-KAHAN FORM AND USING DHILLONS "HOLY GRAIL"
c     CODE.

c     %-----------------------------------------%
c     | R = P * S * Q^T, M^T <- P^T * M^T
c     %-----------------------------------------%
            call dbdsqr('u',j,ncvt,0,ncc,work(ib1),work(ib1+lanmax),
     c           work(iq),lanmax, work,1, work(ip),lanmax+1,
     c           work(iwrk),info)
            call  second(t3)
            tbsvd = tbsvd + (t3-t2)
            nbsvd = nbsvd + 1
C              do i=1,j+1
C               do ii= 1,j+1
C                  write(16,'(d24.16)')work(ip +(i-1)*(lanmax+1) +ii-1)
C               end do
C            end do
                         
            if (lsame(jobu,'y')) then
c     %-----------------------------------------%
c     | Form left Ritz-vectors
c     | U = U * M * P
c     %-----------------------------------------%
              
               call zgemmina(m,j+1,j+1,U,ldu,
     c           work(ip),lanmax+1,zwork(iwrk),lwrk) 
c              call zgemmina(m,j+1,j+1,U,ldu,
c     c           work(ip),lanmax+1,zwork(iwrk),lwrk)
            endif
            
 
C            do ii =1,j+1
C              do i=1,m
C              
C                write(23,'(2d24.16)')U(i,ii)
C             end do
C            end do  

            if (lsame(jobv,'y')) then

c     %-----------------------------------------%
c     | Form right Ritz-vectors
c     | V = V * Q
c     %-----------------------------------------%
               call zgemmina(n,j,j,V,ldv,
     c           work(iq),lanmax,zwork(iwrk),lwrk) 
c              call zgemmina(n,j,j,V,ldv,
c     c           work(iq),lanmax,zwork(iwrk),lwrk)
c               call zgemm_ovwr_left('t',n,j,j,(1.0d0,0.0d0),
c     c           V,ldv,(0.0d0,0.0d0),
c     c           zwork(iq),lanmax,zwork(iwrk),lwrk)
            endif               
         endif
      endif
c        do i=1,n
c                write(14,'(10d24.16)')( V(i,j),j=1,5)
c            end do 

 100  k = neig
      nlandim = j
      call second(t1)
      tlansvd = t1-t0
      end
  



      subroutine refinebounds(n,theta,bound,tol,eps34)
      implicit none
      integer n
      double precision theta(*), bound(*), tol,eps34,gap
      integer i,l
      double precision dlapy2
      external dlapy2
c
c     Refine Lanczos error bounds using the gap theorem.
c     
      if (n.le.1) return
      do i=1,n
         do l=-1,1,2
            if ((l.eq.1.and.i.lt.n) .or. (l.eq.-1.and.i.gt.1)) then
               if (abs(theta(i)-theta(i+l)) .lt. eps34*(theta(i))) then
                  if (bound(i).gt.tol .and. bound(i+l).gt.tol) then
                     bound(i+l) = dlapy2(bound(i),bound(i+l))
                     bound(i) = 0.0
                  endif
               endif
            endif
         enddo
         gap = theta(i+1)-bound(i+1)-(theta(i)+bound(i))
         if (gap.gt.bound(i)) then
            bound(i) = bound(i) * (bound(i)/gap)
         endif
      enddo
      end



      subroutine dbdqr(jobq, n, D, E, c1, c2, Qt, ldq)
      implicit none

c Compute QR factorization B = Q*R of (n+1) x n lower bidiagonal matrix 
c with diagonal elements d(1)...d(n) and first subdiagonal elements
c e(1)...e(n). On return [0 ... 0 c1 c2]' = Q'*[0 ... 0 1]'.
c
c If jobq=='Y' then on return Qt contains Q^T.

      character*1 jobq
      integer n,ldq
      double precision D(*),E(*),c1,c2,Qt(ldq,*)
      
      integer i,j
      double precision cs,sn,r
      logical lsame
      external lsame

      if (n.lt.1) return
      if (lsame(jobq,'Y')) then
         do j=1,n+1
            do i=1,n+1
               Qt(i,j) = 0.0
            enddo
         enddo
         do j=1,n+1
            Qt(j,j) = 1.0
         enddo
      endif
      do i=1,n-1
         call dlartg(d(i),e(i),cs,sn,r)
         d(i) = r
         e(i) = sn*d(i+1)
         d(i+1) = cs*d(i+1)
         if (lsame(jobq,'Y')) then
            do j=1,i
               Qt(i+1,j) = -sn*Qt(i,j)
               Qt(i,j) = cs*Qt(i,j)
            enddo
            Qt(i,i+1) = sn
            Qt(i+1,i+1) = cs
         endif
      enddo
      call dlartg(d(n),e(n),cs,sn,r)
      d(n) = r
      e(n) = 0.0
      c1 = sn
      c2 = cs
      if (lsame(jobq,'Y')) then
         do j=1,i
            Qt(i+1,j) = -sn*Qt(i,j)
            Qt(i,j) = cs*Qt(i,j)
         enddo
         Qt(i,i+1) = sn
         Qt(i+1,i+1) = cs
      endif
      end

