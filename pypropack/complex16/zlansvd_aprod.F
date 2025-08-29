      subroutine zlansvd_aprod(A,jobu,jobv,m,n,k,kmax,U,ldu,
     c     Sigma,bnd,V,ldv,tolin,work,lwork,zwork,lzwrk,iwork,liwork,
     c     doption,ioption,info,zparm,iparm)

c     ZLANSVD_APROD: Compute the leading singular triplets of a large and
c     sparse matrix by Lanczos bidiagonalization with partial
c     reorthogonalization.
c
c     Parameters:
c
c     A(M,N): DOUBLE COMPLEX array.
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
c           SUBROUTINE DAPROD(TRANSA,M,N,X,Y,ZPARM,IPARM)
c           CHARACTER*1 TRANSA
c           INTEGER M,N,IPARM(*)
c           DOUBLE PRECISION X(*),Y(*),ZPARM(*)
c
c           If TRANSA.EQ.'N' then the function should compute the matrix-vector
c           product Y = A * X.
c           If TRANSA.EQ.'T' then the function should compute the matrix-vector
c           product Y = A^T * X.
c           The arrays IPARM and ZPARM are a means to pass user supplied
c           data to APROD without the use of common blocks.
c     U(LDU,KMAX+1): DOUBLE COMPLEX array. On return the first K columns of U
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
c     V(LDV,KMAX): DOUBLE COMPLEX array. On return the first K columns of V
c               will contain approximations to the right singular vectors 
c               corresponding to the K largest singular values of A.
c     LDV: INTEGER. Leading dimension of the array V. LDV >= N.
c     TOLIN: DOUBLE PRECISION. Desired relative accuracy of computed singular 
c            values. The error of SIGMA(I) is approximately 
c            MAX( 16*EPS*SIGMA(1), TOLIN*SIGMA(I) )
c     WORK(LWORK): DOUBLE PRECISION array. Workspace of dimension LWORK.
c     LWORK: INTEGER. Dimension of WORK.
c            If JOBU.EQ.'N' and JOBV.EQ.'N' then  LWORK should be at least
c            M + N + 9*KMAX + 2*KMAX**2 + 4 + MAX(M,N,4*KMAX+4).
c            If JOBU.EQ.'Y' or JOBV.EQ.'Y' then LWORK should be at least
c            M + N + 9*KMAX + 5*KMAX**2 + 4 + 
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
c     IOPTION: INTEGER array. Parameters for LANBPRO.
c        ioption(1) = CGS.  If CGS.EQ.1 then reorthogonalization is done
c          using iterated classical GRAM-SCHMIDT. IF CGS.EQ.0 then 
c          reorthogonalization is done using iterated modified Gram-Schmidt.
c        ioption(2) = ELR. If ELR.EQ.1 then extended local orthogonality is
c          enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
c     INFO: INTEGER. 
c         INFO = 0  : The K largest singular triplets were computed succesfully
c         INFO = J>0, J<K: An invariant subspace of dimension J was found.
c         INFO = -1 : K singular triplets did not converge within KMAX
c                     iterations.   
c     ZPARM: DOUBLE COMPLEX array. Array used for passing data to the APROD
c         function.   
c     IPARM: INTEGER array. Array used for passing data to the APROD
c         function.   
c
c     (C) Brian J Soher, 2020 
c


c     %-----------%
c     | Arguments |
c     %-----------%
c
c     All variables have to be declared.
      
      implicit none
      
      include 'stat.h'
      
      complex*16, intent(in),dimension(m,n),target :: A
      character*1, intent(in) :: jobu
      character*1, intent(in) :: jobv
      integer, intent(in)  :: m
      integer, intent(in)  :: n
      integer, intent(in)  :: k
      integer, intent(in)  :: kmax
      complex*16, intent(inout),dimension(ldu,kmax+1) :: U
      integer ldu
      double precision, intent(out),dimension(k) :: Sigma
      double precision, intent(out),dimension(k) :: bnd
      complex*16, intent(inout),dimension(ldv,kmax) :: V
      integer ldv
      double precision, intent(in)  :: tolin
      double precision, intent(in),dimension(lwork) :: work(lwork)
      integer lwork
      complex*16, intent(in),dimension(lzwrk) :: zwork(lzwrk)
      integer lzwrk
      integer, intent(in),dimension(liwork) :: iwork(liwork)
      integer liwork
      double precision, intent(in),dimension(3) :: doption
      integer, intent(in),dimension(2) :: ioption
      integer, intent(out) :: info
      complex*16, intent(in),dimension(*) :: zparm(*)
      integer, intent(in),dimension(*) :: iparm(*)

      integer i,j, xn, idist, iseed(4)
      complex*16 bob(16)

c     %------------------------%
c     | Common block variables |
c     %------------------------%

      complex*16, pointer :: pA(:,:)  
      common /csvdp/ pA     

      external zaprod, zlansvd, zlarnv_h, zlarnvbs_h


c-------------------- Here begins executable code ---------------------

c     We have to point the common block var pA at A then move on ...

      pA=>A

c     Here we call zlansvd

      call zlansvd(jobu,jobv,m,n,k,kmax,zaprod,U,ldu,Sigma,bnd,
     c     V,ldv,tolin,work,lwork,zwork,lzwrk,iwork,liwork,doption,
     c     ioption,info,zparm,iparm)

      end
