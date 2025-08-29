      subroutine zlansvd_irl_aprod(A,which,jobu,jobv,m,n,dim,p,neig,
     c     maxiter,U,ldu,Sigma,bnd,V,ldv,tolin,work,lwork,zwork,
     c     lzwrk,iwork,liwork,doption,ioption,info,zparm,iparm)

c     ZLANSVD_IRL_APROD: Compute the leading singular triplets of a large
c     and sparse matrix A by implicitly restarted Lanczos bidiagonalization 
c     with partial reorthogonalization.
c
c     Parameters:
c
c     A:    DOUBLE COMPLEX array.
c     WHICH: CHARACTER*1. Decides which singular triplets to compute. 
c            If WHICH.EQ.'L' then compute triplets corresponding to the K
c            largest singular values. 
c            If WHICH.EQ.'S' then compute triplets corresponding to the K
c            smallest singular values. 
c     JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
c           Otherwise the array U is not touched.
c     JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular 
c           vectors. Otherwise the array V is not touched.
c     M:    INTEGER. Number of rows of A.
c     N:    INTEGER. Number of columns of A.
c     DIM:  INTEGER. Dimension of the Krylov subspace.
c     P:    INTEGER. Number of shift per restart.
c     NEIG: INTEGER. Number of desired singular triplets.  
c           NEIG <= MIN(DIM-P,M,N)
c     MAXITER: INTEGER. Maximum number of restarts.
c     APROD: Subroutine defining the linear operator A. 
c            APROD should be of the form:
c
c           SUBROUTINE ZAPROD(TRANSA,M,N,X,Y,DPARM,IPARM)
c           CHARACTER*1 TRANSA
c           INTEGER M,N,IPARM(*)
c           DOUBLE COMPLEX X(*),Y(*),DPARM(*)
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y = A * X.
c           If TRANSA.EQ.'T' then the function should compute the matrix-vector
c           product Y = A^T * X.
c           The arrays IPARM and DPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.
c     U(LDU,KMAX+1): DOUBLE COMPLEX array. On return the first K columns of U
c               will contain approximations to the left singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH)  singular values of A.
c               On entry the first column of U contains the starting vector
c               for the Lanczos bidiagonalization. A random starting vector
c               is used if U is zero.
c     LDU: INTEGER. Leading dimension of the array U. LDV >= M.
c     SIGMA(K): DOUBLE PRECISION array. On return Sigma contains approximation
c               to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     BND(K)  : DOUBLE PRECISION array. Error estimates on the computed 
c               singular values. The computed SIGMA(I) is within BND(I)
c               of a singular value of A.
c     V(LDV,KMAX): DOUBLE COMPLEX array. On return the first K columns of V
c               will contain approximations to the right singular vectors 
c               corresponding to the K largest or smallest (depending on the 
c               value of WHICH) singular values of A.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     TOLIN: DOUBLE PRECISION. Desired relative accuracy of computed singular 
c            values. The error of SIGMA(I) is approximately 
c            MAX( 16*EPS*SIGMA(1), TOLIN*SIGMA(I) )
c     WORK(LWORK): DOUBLE PRECISION array. Workspace of dimension LWORK.
c     LWORK: INTEGER. Dimension of WORK. 
c            If JOBU.EQ.'N' and JOBV.EQ.'N' then  LWORK should be at least
c            M + N + 10*KMAX + 2*KMAX**2 + 5 + MAX(M,N,4*KMAX+4).
c            If JOBU.EQ.'Y' or JOBV.EQ.'Y' then LWORK should be at least
c            M + N + 10*KMAX + 5*KMAX**2 + 4 + 
c            MAX(3*KMAX**2+4*KMAX+4, NB*MAX(M,N)), where NB>0 is a block 
c            size, which determines how large a fraction of the work in
c            setting up the singular vectors is done using fast BLAS-3 
c            operation. 
c     ZWORK: DOUBLE COMPLEX array of dimension ????.
c     IWORK: INTEGER array. Integer workspace of dimension LIWORK.
c     LIWORK: INTEGER. Dimension of IWORK. Should be at least 8*KMAX if
c             JOBU.EQ.'Y' or JOBV.EQ.'Y' and at least 2*KMAX+1 otherwise.
c     DOPTION: DOUBLE PRECISION array. Parameters for LANBPRO.
c        doption(1) = delta. Level of orthogonality to maintain among
c          Lanczos vectors.
c        doption(2) = eta. During reorthogonalization, all vectors with
c          with components larger than eta along the latest Lanczos vector
c          will be purged.
c        doption(3) = anorm. Estimate of || A ||.
c        doption(4) = min relgap. Smallest relgap allowed between any shift
c                     the smallest requested Ritz value.
c
c     IOPTION: INTEGER array. Parameters for LANBPRO.
c        ioption(1) = CGS.  If CGS.EQ.1 then reorthogonalization is done
c          using iterated classical Gram-Schmidt. IF CGS.EQ.0 then 
c          reorthogonalization is done using iterated modified Gram-Schmidt.
c        ioption(2) = ELR. If ELR.EQ.1 then extended local orthogonality is
c          enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
c
c     INFO: INTEGER. 
c         INFO = 0  : The K largest or smallest (depending on the 
c                     value of WHICH) singular triplets were computed 
c                     succesfully.
c         INFO = J>0, J<K: An invariant subspace of dimension J was found.
c         INFO = -1 : K singular triplets did not converge within KMAX
c                     iterations.   
c     ZPARM: DOUBLE COMPLEX array. Array used for passing data to the APROD
c         function.   
c     IPARM: INTEGER array. Array used for passing data to the APROD
c         function.   
c
c
c     (C) Rasmus Munk Larsen, Stanford University, 2000,2004
c


c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      
      complex*16, intent(in),dimension(m,n),target :: A
      character*1, intent(in) :: which
      character*1, intent(in) :: jobu
      character*1, intent(in) :: jobv
      integer, intent(in)  :: m
      integer, intent(in)  :: n
      integer, intent(in)  :: dim
      integer, intent(in)  :: p
      integer, intent(in)  :: neig
      integer, intent(in)  :: maxiter
      complex*16, intent(inout),dimension(ldu,dim+1) :: U
      integer ldu
      double precision, intent(out),dimension(neig) :: Sigma
      double precision, intent(out),dimension(neig) :: bnd
      complex*16, intent(inout),dimension(ldv,dim) :: V
      integer ldv
      double precision, intent(in)  :: tolin
      double precision, intent(in),dimension(lwork) :: work(lwork)
      integer lwork
      complex*16, intent(in),dimension(lzwrk) :: zwork(lzwrk)
      integer lzwrk
      integer, intent(in),dimension(liwork) :: iwork(liwork)
      integer liwork
      double precision, intent(in),dimension(4) :: doption
      integer, intent(in),dimension(2) :: ioption
      integer, intent(out) :: info
      complex*16, intent(in),dimension(*) :: zparm(*)
      integer, intent(in),dimension(*) :: iparm(*)

      integer i,j

c     %------------------------%
c     | Common block variables |
c     %------------------------%

      complex*16, pointer :: pA(:,:)  
      common /csvdp/ pA     

      external zaprod, zlansvd_irl


c-------------------- Here begins executable code ---------------------
      
c     We have to point the common block var pA at A then move on ...

      pA=>A

c     Here we call zlansvd_irl

      call zlansvd_irl(which,jobu,jobv,m,n,dim,p,neig,maxiter,
     c     zaprod,U,ldu,Sigma,bnd,V,ldv,tolin,work,lwork,zwork,lzwrk,
     c     iwork,liwork,doption,ioption,info,zparm,iparm)

      end
