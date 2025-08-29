cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c........1.........2.........3.........4.........5.........6.........7..
c23456789012345678901234567890123456789012345678901234567890123456789012
c
      subroutine zlansvdw_python(ndp, m, n, k, sigma, 
     c                           info, 
     c                           lambda_r,  lambda_i, 
     c                           trlambda_r,  trlambda_i,
     c                           uuu_r, uuu_i)
      implicit none

C output = k, info, uuu, sigma

c m == lrow
c n == mcol
c sigma = lsinval (output)
c info = ndiv = nsv_found = kfit
c uuu = output

C     params
      integer ndp, m, n, k, kmax
      parameter (kmax=50)
      real*8     sigma(kmax)
      integer info
      real*8 lambda_r(ndp), lambda_i(ndp)
      real*8 trlambda_r(ndp), trlambda_i(ndp)
      real*8 uuu_r(m, kmax + 1), uuu_i(m, kmax + 1)


c     Local variables.
      integer i, j
      real*8     tol, pi, bnd(kmax)
      complex*16 uuu(m, kmax + 1)
      complex*16 lambda(ndp)
      complex*16 trlambda(ndp)
      integer lwrk,lzwrk
C       integer ndpmax
C       parameter (ndpmax=8192)

      integer ioption(10), iwork(2*kmax+1)
      real*8 doption(10)
      complex*16 vvv(n, kmax)

      real*8, allocatable :: work(:)
      complex*16, allocatable :: zwork(:)


      tol = 16d-16
      pi = 2*acos(0d0)

      ioption(1) = 1
      ioption(2) = 1
      doption(1) = 1d-12
      doption(2) = 1d-14
      doption(3) = 0.0d0

C     PS - I allocate work & zwork dynamically,  otherwise I'll have
C     to delcare them much earlier (like in hlsvdpro()).
      lwrk=m+n+13*kmax+8*kmax**2+32*m+ndp*kmax
      lzwrk=m+n+32*m+7*kmax+2+2*kmax**2+5*ndp

      allocate(work(lwrk))
      allocate(zwork(lzwrk))


C     PS - copying this isn't really necessary. ctypes and Fortran both
C     initialize them to 0, so I'm copying a bunch of zeroes here.      
      do i = 1, m
        do j = 1, kmax + 1
          uuu(i,j) = dcmplx(uuu_r(i,j), uuu_i(i,j))
C            PRINT *, "uuu_i(", i, j, "):", uuu_i(i,j)
       end do
      end do

      do i = 1, ndp
         lambda(i) = dcmplx(lambda_r(i), lambda_i(i))
C          print *, "lambda(", i, ") = ", lambda(i)
      end do

      do i = 1, ndp
         trlambda(i) = dcmplx(trlambda_r(i), trlambda_i(i))
C          print *, "trlambda(", i, ") = ", trlambda(i)
      end do



      call zlansvdw('y','n', ndp, m, n, k, kmax, uuu, m, sigma, 
     c              bnd, vvv, n,
     c              tol, work, lwrk, iwork, doption, ioption, info,
     c              zwork, lzwrk, lambda, trlambda, 0, 0)
C        print *, "done zlansvdw()"
       do i = 1, m
          do j = 1, kmax + 1
C               PRINT *, "uuu(", i, j, "):", uuu(i,j)
              uuu_r(i,j) = REALPART(uuu(i,j))
              uuu_i(i,j) = IMAGPART(uuu(i,j))
         end do
       end do


      deallocate(work)
      deallocate(zwork)

      end


      subroutine zlansvdw(jobu,jobv,ndp,m,n,k,kmax,U,ldu,Sigma,bnd,V,
     c     ldv,tolin,work,lwork,iwork,doption,ioption,info,
     c       zwork,lzwrk,lambda,trlambda, planF,planB)


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
      integer lzwrk
      integer m,n,k,kmax,lanmax,ldu,ldv,iwork(*),lwork,info,ioption(*)
C     len iwork = 2*kmax+1      
      double precision Sigma(*),bnd(*),work(*)
      double precision tolin,doption(*)
      double complex U(ldu,*),V(ldv,*),zwork(lzwrk),
     c               lambda(*),trlambda(*)
      complex*16, DIMENSION(:, :), allocatable :: uuu_copy
      complex*16 delta


c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero, FUDGE,DEBUG     
      parameter(one = 1.0d0, zero = 0.0d0, FUDGE = 1.01d0, DEBUG=0.0d0)
            
c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,j,dj,jold,ibnd,ib,ib1,iwrk,ierr,ip,iq,neig,lwrk
      integer ips, jps
      integer ncvt,ncc,ndp, planF,planB
      double precision eps,eps34,epsn2,epsn,sfmin,anorm,rnorm,tol
      integer ppi, ppj

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

C       print *, "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ zlansvdw() "
C       call FLUSH()

      allocate(uuu_copy(ldu, k))

C     print scalars for debugging
C         print *, "jobu, jobv = ", jobu, jobv
C         print *, "ndp, m, n = ", ndp, m, n
C         print *, "k, kmax = ", k, kmax
C         print *, "ldu, ldv = ", ldu, ldv
C         print *, "tolin, lwork = ", tolin, lwork
C         print *, "info = ", info
C         print *, "lzwrk = ", lzwrk
C         print *, "planF,planB = ", planF,planB


C       print *, "printing non-zero UUUs..."
C       do i=1,ldu
C           do j = 1,kmax+1
C               if (U(i,j).ne.(0.0d0,0.0d0)) then
C                   print *, "U(", i, ", ", j, ") = ", U(i,j)
C               end if
C           end do
C       end do 

C       do i=1,ndp
C         print *, "lambda(", i, ") = ", lambda(i)
C       end do 
 
C       do i=1,ndp
C         print *, "trlambda(", i, ") = ", trlambda(i)
C       end do 
 

c     %---------------------------------%
c     | Set machine dependent constants |
c     %---------------------------------%
      eps = dlamch('e')
      eps34 = eps**(3.0/4.0)
      epsn = dble(max(m,n))*eps/2.0
      epsn2 = sqrt(dble(max(m,n)))*eps/2.0
      sfmin = dlamch('s')
      


C       print *, "eps = ", eps
C       print *, "eps34 = ", eps34
C       print *, "epsn = ", epsn
C       print *, "epsn2 = ", epsn2
C       print *, "sfmin = ", sfmin


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

C       print *, "lanmax = ", lanmax
C       print *, "ib = ", ib
C       print *, "ib1 = ", ib1
C       print *, "ip = ", ip
C       print *, "iq = ", iq
C       print *, "iwrk = ", iwrk
C       print *, "lzwrk = ", lzwrk
C       print *, "lzwrk - iwrk = ", lzwrk - iwrk




C     PS - changed from lwork to lzwrk
C       print *, "lwork - iwrk + 1 = ", lwork - iwrk + 1
C       print *, "lzwrk - iwrk + 1 = ", lzwrk - iwrk + 1
      lwrk = lzwrk-iwrk+1
      call dzero(7*lanmax + 2 + 2*lanmax**2,work,1)
      call zzero(7*lanmax + 2 + 2*lanmax**2,zwork,1)

C       print *, "lwrk (calc) = ", lwrk
C       print *, "lwork (param) = ", lwork
C       print *, "lzwrk = ", lzwrk
C       print *, "iwrk = ", iwrk
C       print *, "lanmax = ", lanmax
C       print *, "z/dzero param = ", 7*lanmax + 2 + 2*lanmax**2
C       print *, "PlanF/B = ", planF, planB




c     %---------------------------------------------------------------%
c     | Set up random starting vector if none is provided by the user |
c     %---------------------------------------------------------------%
c      rnorm = dnrm2(m,U(1,1),1)
      rnorm =0.0d0
      if (rnorm.eq.zero) then
C          print *, "calling zgetu0w()"
C           do i=1, ldu
C             do j = 1, k
C                 uuu_copy(i, j) = U(i, j)
C             end do
C           end do

          call zgetu0w('n',ndp, m,n,0,1,U,rnorm,U,ldu,
     c        ierr,ioption(1),anorm,zwork(iwrk),lambda,
     c        trlambda,planF,planB)
C           print *, "done zgetu0w()"
C           print *, "rnorm = ", rnorm, "anorm = ", anorm


C           print *, "printing uuu_copy deltas"          
C           do i=1, ldu
C             do j = 1, k
C                 delta = uuu_copy(i, j) - U(i, j)
C                 if (delta.ne.(0.0d0,0.0d0)) then
C                     print *, "delta(", i, ",", j, ") = ", delta
C                 end if
C             end do
C           end do

C           print *, "stifle!"

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
C          print *, "looptop, neig = ", neig, "k = ", k
C          print *, "before zlanbprow, jold = ", jold, ", j = ", j,
C      c            ", rnorm = ", rnorm
C         call FLUSH()
         call zlanbprow(ndp, m, n, jold, j,  U, ldu, V, ldv,
     c        work(ib),lanmax,rnorm,doption(1),ioption(1),
     c        work(iwrk), iwork,ierr,  
     c        zwork(iwrk), lambda,trlambda,planF,planB)
C          print *, "zlansvdw: done zlanbprow, k = ", k
C          print *, "after zlanbprow, jold = ", jold, ", j = ", j,
C      c            ", rnorm = ", rnorm
C          call FLUSH()
         jold = j

C          print *, "printing non-zero work"
C          do ips = 1, lwork
C           if (work(ips).ne.0.0d0) then
C                print *, "work(", ips, ") = ", work(ips)    
C             end if
C          end do 

C          print *, "####"

C       do i=1,ldu
C           do j = 1,kmax+1
C               print *, "U(", i, ", ", j, ") = ", U(i,j)
C           end do
C       end do 


c     %---------------------------------------------%
c     | Compute and analyze SVD(B) and error bounds |
c     %---------------------------------------------%

         call dcopy(2*lanmax, work(ib),1,work(ib1),1)
         call dzero(j+1,work(ibnd),1)

C          print *, "printing non-zero work"
C          do ips = 1, lwork
C           if (work(ips).ne.0.0d0) then
C                print *, "work(", ips, ") = ", work(ips)    
C             end if
C          end do 

C          print *, "####"

         call second(t2)
C          print *, "before dbdqr1(), j = ", j

C         print *, "work(1b1) --"
C         print *, "ib1 = ", ib1, ", ib1+ lanmax = ", ib1+ lanmax
C         do ips = 1, (lanmax * 2)
C           print *, "work(", ips - 1 + ib1, ") = ", work(ips - 1 + ib1)
C         end do
C          print *, "####"

C          print *, "ibnd+j-1 = ", ibnd+j-1
C          print *, "ibnd+j = ", ibnd+j

         call dbdqr('N',j,work(ib1),work(ib1+lanmax),work(ibnd+j-1),
     c        work(ibnd+j),work(ip),lanmax+1)

C          print *, "after dbdqr1(), c1 = ", work(ibnd+j-1),
C      c            "c2 = ", work(ibnd+j)

C         do ips = 1, (lanmax * 2)
C           print *, "work(", ips - 1 + ib1, ") = ", work(ips - 1 + ib1)
C         end do
C          print *, "####"


C         do ips = 1, 52
C           print *, "work(", ips, ") = ", work(ips)
C         end do
C         print *, "####"

C         do ips = 1, lanmax **2
C           print *, "work(", ip - 1 + ips, ") = ", work(ip - 1 + ips)
C         end do
C         print *, "####"



C          print *, "printing non-zero work"
C          do ips = 1, lwork
C           if (work(ips).ne.0.0d0) then
C                print *, "work(", ips, ") = ", work(ips)    
C             end if
C          end do 

C          print *, "####"

C          do ips = 1, lwork
C             print *, "work(", ips, ") = ", work(ips)    
C          end do 


C          print *, "ib1 = ", ib1
C          print *, "ib1+lanmax = ", ib1+lanmax
C          print *, "ibnd+j-1 = ", ibnd+j-1
C          print *, "ibnd+j = ", ibnd+j
C          print *, "ip = ", ip

C        UPLO = 'u'
C        N = j
C        NCVT = 0, NRU = 1, NCC = 0
C        D = work(ib1),            (input/output)
C        E = work(ib1+lanmax)      (input/output)
C        VT = work (not referenced because NCVT == 0)
C        LDVT = 1  (meaningless because NCVT == 0)
C        U = work(ibnd) (input/output)
C        LDU = 1
C        C = work (not referenced because NCC == 0)
C        LDC = 1  (meaningless because NCC == 0)
C        WORK = work(iwrk)  (len  = 4 * N)
C        INFO = info (output)

C        dbdsqr "destroys" the contents of E (work(ib1+lanmax)).
C        Really, only D (work(ib1) and U (work(ibnd)) are significant
C        when this call completes, right?
C          print *, "before dbdsqr, j = ", j

C          do ips = 1, lanmax + 1
C             jps = ibnd + ips - 1
C             print *, "work(", jps, ") = ", work(jps)
C          end do


C          print *, "work(ib1+lanmax) --"
C          do ips = 1, lanmax
C             jps = ips - 1 + ib1 + lanmax
C            print *, "work(", jps, ") = ", work(jps)
C          end do

C          print *, "####"

         call dbdsqr('u',j,0,1,0,work(ib1),work(ib1+lanmax),work,1,
     c        work(ibnd),1,work,1,work(iwrk),info)

C           print *, "after dbdsqr"
C          do ips = 1, lanmax + 1
C           jps = ibnd + ips - 1
C           print *, "work(", jps, ") = ", work(jps)
C          end do


C          do ips = 1, lanmax + 1
C             jps = ibnd + ips - 1
C             print *, "work(", jps, ") = ", work(jps)
C          end do

C          print *, "work(ib1) --"
C          do ips = 1, lanmax
C            print *, "work(", ips - 1 + ib1, ") = ", work(ips - 1 + ib1)
C          end do
C          print *, "####"

C          print *, "work(ibnd == 1) --"
C          do ips = 1, lanmax + 1
C            print *, "work(", ips, ") = ", work(ips)
C          end do
C          print *, "####"

C          print *, "work(ib1+lanmax) --"
C          do ips = 1, lanmax + 1
C             jps = ips + ib1+lanmax - 1
C            print *, "work(", jps, ") = ", work(jps)
C          end do
C          print *, "####"

         call  second(t3)
         tbsvd = tbsvd + (t3-t2)
         nbsvd = nbsvd + 1

         if (j.gt.5) then
            anorm = work(ib1)
         else
            anorm = max(anorm,work(ib1))
         endif

C          print *, "anorm = ", anorm

         do i=1,j
            work(ibnd+i-1) = abs(rnorm*work(ibnd+i-1))
            work(ib1+lanmax+i-1) = work(ib1+lanmax+i-1)**2
         enddo
C          do i=1,j
C             print *, "work(", ibnd+i-1, ") = ", work(ibnd+i-1)
C          enddo
C          do i=1,j
C             print *, "work(",ib1+lanmax+i-1,") = ",work(ib1+lanmax+i-1)
C          enddo

c     %---------------------------------------------%
c     | Refine error bounds using the "Gap theorem" |
c     %---------------------------------------------%

C          print *, "before refinebounds"

C          do ips = 1, lanmax + 1
C           jps = ibnd + ips - 1
C           print *, "work(", jps, ") = ", work(jps)
C          end do

         call refinebounds(j,work(ib1+lanmax),work(ibnd),
     c        epsn*anorm,eps34)


C         print *, "after refinebounds"

C          do ips = 1, lanmax + 1
C           jps = ibnd + ips - 1
C           print *, "work(", jps, ") = ", work(jps)
C          end do
         

c     %----------------------------------------------------%
c     | Determine the number of converged singular values  |
c     %----------------------------------------------------%
C          print *, "converge time, j = ", j, "k = ", k

         do i=1,min(j,k)
            bnd(i) = work(ibnd+i-1)
C             print *, "bnd(", i, ") = ", bnd(i)
         enddo


         i = 0
         neig = 0
C         print *, "tol = ", tol
         do while(i.lt.min(j,k))
            ips = ibnd + i
C             print *, "work(", ips, ") = ", work(ips)
            ips = ib1 + i
C             print *, "tol * bbb1_a(",ips,") = ", tol * work(ips)
C            print *,"sigma delta = ", (tol * work(ib1+i)) - work(ibnd+i)
            if (work(ibnd+i).le.tol*work(ib1+i)) then
               neig = neig + 1
               sigma(neig) = work(ib1+i)
C              print *, "sigma(", neig, ") = ", sigma(neig)
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
C             print *, "ierr = ", ierr, ", going to 50"
            goto 50               
         endif
         if (j.ge.lanmax) then
            if (neig.lt.k) then
               info = -1
            endif
C             print *, "going to 50"
C             print *, "j=",j,",lanmax=",lanmax,",neig=",neig,"k=",k
            goto 50
         endif

c     %----------------------------------------------------%
c     | Increase the dimension of the Krylov subspace.     |
c     | If any Ritz values have converged then try to      | 
c     | estimate the average number of iterations per      |
c     | converged Ritz value.                              |
c     | Else increase the dimension with 50%.              |
c     %----------------------------------------------------%
C          print *, "monkey, j = ", j
C          print *, "monkey, j/2 = ", j/2
C          print *, "monkey, k = ", k
C          print *, "monkey, neig = ", neig
C          print *, "monkey, form = ", ((k-neig)*(j-6))/(2*neig+1)

         if (neig.gt.1) then
            dj = min(j/2,((k-neig)*(j-6))/(2*neig+1))
            dj = min(100,max(2,dj))
         else
            dj = j/2
            dj = min(100,max(10,dj))
        endif
         j = min(j + dj,lanmax)
C          print *, "loop end, j = ", j, "dj = ", dj
C          print *, "loop end, neig = ", neig, "k = ", k

      enddo

 50   continue
C       print *, "Pre-calculate singular vectors"
C       print *, "neig = ", neig, ", k = ", k

C  50   if (neig.ge.k) then
      if (neig.ge.k) then
c     %-----------------------------------------%
c     | Calculate singular vectors if requested %
c     %-----------------------------------------%
C          print *, "inside Calculate singular vectors"
         j = jold
         k = neig

         if (lsame(jobu,'y') .or. lsame(jobv,'y')) then
C             print *, "about to dcopy"
            call dcopy(2*lanmax, work(ib),1,work(ib1),1)
C             print *, "done dcopy"

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

C             print *, "before dbdqr2"


C             print *, "work(1b1) --"
C             print *, "ib1 = ", ib1, ", ib1+ lanmax = ", ib1+ lanmax
C             do ips = 1, (lanmax * 2)
C               print *, "work(", ips + ib1, ") = ", work(ips + ib1)
C             end do

            call dbdqr(jobu,j,work(ib1),work(ib1+lanmax),
     c           work(ibnd+j-1), work(ibnd+j),work(ip),lanmax+1)
 
C             print *, "after dbdqr2(), c1 = ", work(ibnd+j-1),
C      c               "c2 = ", work(ibnd+j)
C             do ips=1, (lanmax + 1)**2
C                 jps = ips + ip - 1
C                 print *, "work(", jps, ") = ", work(jps)
C             end do


C             print *, "work(1b1) --"
C             do ips = 1, (lanmax * 2)
C               print *, "work(", ips + ib1, ") = ", work(ips + ib1)
C             end do

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
            
C             print *, "about to set work"
            do i=1,n
                jps = iq+(i-1)*(lanmax+1)
C               print *, "Setting work(", jps, ") = 1.0"
                work(jps) = one
            enddo
C             print *, "done setting work"

c     THIS SHOULD POSSIBLY BE REPLACED BY A CALL TO THE FAST
c     DIVIDE-AND-CONQUER BIDIAGONAL SVD IN LAPACK 3 OR BY TRANSFORMING B
c     TO TRIDIAGONAL GOLUB-KAHAN FORM AND USING DHILLONS "HOLY GRAIL"
c     CODE.

c     %-----------------------------------------%
c     | R = P * S * Q^T, M^T <- P^T * M^T
c     %-----------------------------------------%
C             print *, "about to call dbdsqr2"

C             do ips=1, n * (lanmax + 1)
C               jps = iq + ips - 1
C               print *, "work(", jps, ") = ", work(jps)
C             enddo

C             do ips=1, (lanmax + 1)**2
C                 jps = ips + ip - 1
C                 print *, "work(", jps, ") = ", work(jps)
C             end do

C           Note - dbdsqr() changes the content of work(ip)
C           work(ip) is array C inside dbdsqr


            call dbdsqr('u',j,ncvt,0,ncc,work(ib1),work(ib1+lanmax),
     c           work(iq),lanmax, work,1, work(ip),lanmax+1,
     c           work(iwrk),info)
C             print *, "done dbdsqr2"

C             do ips=1, (lanmax + 1)**2
C                 jps = ips + ip - 1
C                 print *, "work(", jps, ") = ", work(jps)
C             end do

            call  second(t3)
            tbsvd = tbsvd + (t3-t2)
            nbsvd = nbsvd + 1
C              do i=1,j+1
C               do ii= 1,j+1
C                  write(16,'(d24.16)')work(ip +(i-1)*(lanmax+1) +ii-1)
C               end do
C            end do

C             do ips = 1, lanmax + 1
C               jps = ib1 + ips
C               print *, "work(", jps, "( = ", work(jps)
C             end do
                         
C             print *, "########"

C             do ips = 1, lanmax + 1
C               jps = ib1 + lanmax + ips
C               print *, "work(", jps, "( = ", work(jps)
C             end do

C             print *, "########"

            if (lsame(jobu,'y')) then
c     %-----------------------------------------%
c     | Form left Ritz-vectors
c     | U = U * M * P
c     %-----------------------------------------%
              
C                 print *, "about to call zgemmina"
C                 print *, "m = ", m
C                 print *, "j = ", j
C                 print *, "ldu = ", ldu
C                 print *, "ip = ", ip
C                 print *, "work(ip) = ", work(ip)
C                 print *, "lanmax = ", lanmax
C                 print *, "iwrk = ", iwrk
C                 print *, "lwrk = ", lwrk
C                print *, "U(1, 1) = ", U(1, 1)
C                 do ips=1, (lanmax + 1)**2
C                     jps = ips + ip - 1
C                     print *, "work(", jps, ") = ", work(jps)
C                 end do

C                do ips = 1,5
C                  do jps = 1, 5
C                     print *, "U(", ips, ",", jps, ") = ", U(ips, jps)
C                  end do
C                end do

               call zgemmina(m,j+1,j+1,U,ldu,
     c           work(ip),lanmax+1,zwork(iwrk),lwrk) 
C                 print *, "done zgemmina"
c              call zgemmina(m,j+1,j+1,U,ldu,
c     c           work(ip),lanmax+1,zwork(iwrk),lwrk)
C                do ips = 1,20
C                  do jps = 1, 20
C                     print *, "U(", ips, ",", jps, ") = ", U(ips, jps)
C                  end do
C                end do

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

      deallocate(uuu_copy)

C        print *, "At end of zlansvdw()"
C       do i=1,ldu
C         do j = 1,kmax+1
C             print *, "U(", i, ", ", j, ") = ", U(i,j)
C         end do
C       end do 



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

C       print *, "inside refinebounds"

C       do i=1,n
C         print *, "theta(", i, ") = ", theta(i)
C       end do

      if (n.le.1) return
      do i=1,n
         do l=-1,1,2
C               print *, "i = ", i, "l = ", l

            if ((l.eq.1.and.i.lt.n) .or. (l.eq.-1.and.i.gt.1)) then
               if (abs(theta(i)-theta(i+l)) .lt. eps34*(theta(i))) then
                  if (bound(i).gt.tol .and. bound(i+l).gt.tol) then
                     bound(i+l) = dlapy2(bound(i),bound(i+l))
                     bound(i) = 0.0
C                      print *, "1: bound(", i + l, ") = ", bound(i + l)
                  endif
               endif
            endif
         enddo
C          print *, "theta(", i + 1, ") = ", theta(i+1)
C          print *, "bound(", i + 1, ") = ", bound(i+1)
C          print *, "theta(", i, ") = ", theta(i)
C          print *, "bound(", i, ") = ", bound(i)
         gap = theta(i+1)-bound(i+1)-(theta(i)+bound(i))
C          print *, "i = ", i, "gap = ", gap
         if (gap.gt.bound(i)) then
            bound(i) = bound(i) * (bound(i)/gap)
C             print *, "2: bound(", i , ") = ", bound(i)
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
      integer ips, jps
      double precision D(*),E(*),c1,c2,Qt(ldq,*)
      
      integer i,j
      double precision cs,sn,r
      logical lsame
      external lsame

C        print *, "inside dbdqr, n = ", n
C       do ips = 1, 50
C         print *, "d(", ips, ") = ", d(ips)
C       end do
C       do ips = 1, 50
C         print *, "e(", ips, ") = ", e(ips)
C       end do

C       print *, "####"
          
      if (n.lt.1) return
      if (lsame(jobq,'Y')) then
         do j=1,n+1
            do i=1,n+1
               Qt(i,j) = 0.0
C                print *, "Z Qt(", i, ",", j, ") = ", Qt(i,j)
            enddo
         enddo
         do j=1,n+1
            Qt(j,j) = 1.0
C              print *, "O Qt(", j, ",", j, ") = ", Qt(j,j)
         enddo
      endif
      do i=1,n-1
C          print *, "calling dlartg, i = ", i
C          print *, "d(", i, ") = ", d(i)
C          print *, "e(", i, ") = ", e(i)
         call dlartg(d(i),e(i),cs,sn,r)
C         print *, "i = ", i, ", cs = ", cs, "sn = ", sn, "r = ", r
C         print *, "dlartg i = ", i, ", cs = ", cs, "sn = ", sn
         d(i) = r
         e(i) = sn*d(i+1)
C          print *, "setting d(", i + 1, ") = ", cs*d(i+1)
         d(i+1) = cs*d(i+1)
         if (lsame(jobq,'Y')) then
            do j=1,i
C                 print *, "A Qt(", i, ",", j, ") = ", Qt(i,j)
               Qt(i+1,j) = -sn*Qt(i,j)
               Qt(i,j) = cs*Qt(i,j)
            enddo
            Qt(i,i+1) = sn
C             print *, "B Qt(", i, ",", i +1, ") = ", Qt(i,i+1)
            Qt(i+1,i+1) = cs
         endif
      enddo
      call dlartg(d(n),e(n),cs,sn,r)
C       print *, "i = ", i, ", cs = ", cs, "sn = ", sn
      d(n) = r
      e(n) = 0.0
      c1 = sn
      c2 = cs
      if (lsame(jobq,'Y')) then
         do j=1,i
C             print *, "C Qt(", i, ",", j, ") = ", Qt(i,j)
            Qt(i+1,j) = -sn*Qt(i,j)
            Qt(i,j) = cs*Qt(i,j)
         enddo
         Qt(i,i+1) = sn
C           print *, "D Qt(", i, ",", i +1, ") = ", Qt(i,i+1)
         Qt(i+1,i+1) = cs
      endif

C       print *, "exiting dbdqr"
C       do ips = 1, 50
C         print *, "d(", ips, ") = ", d(ips)
C       end do
C       do ips = 1, 50
C         print *, "e(", ips, ") = ", e(ips)
C       end do

C       print *, "####"


      end

