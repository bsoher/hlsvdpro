
      subroutine zgetu0w(transa, ndp, m, n, j, ntry, u0, u0norm, U, ldu,
     c      ierr, icgs, anormest, work,lambda,trlambda,planF,planB)
c
c Modified by Nicola Mastronardi
c
c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      include 'stat.h'
      character*1 transa
      integer m, n, j, i, ntry, ldu, ierr,icgs,ndp
      integer ips, jps
      integer planF,planB
      double precision u0norm,anormest
      
ccc      double precision u0(*),U(*),work(*)
ccc      external aprod      
      double complex  u0(*),U(*),work(*),lambda(*),
     c                trlambda(*) 


c     %------------%
c     | Parameters |
c     %------------%
      integer MGS
      double precision kappa
      parameter(kappa = 0.717,MGS=1)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer itry,idist,iseed(4),rsize,usize,index(3)
      double precision nrm
      real t1,t2,t3

c     %--------------------%
c     | External Functions |
c     %--------------------%
      logical lsame
      double precision dznrm2
      external dznrm2,lsame
      external aprodw
c-------------------- Here begins executable code ---------------------
      call second(t1)
      iseed(1) = 1
      iseed(2) = 3
      iseed(3) = 5
      iseed(4) = 7
c      iseed(3) = 11
c      iseed(4) = 7
      if (lsame(transa,'n')) then
c     %-------------------------%
c     | u0 is to be an m-vector |
c     %-------------------------%
         rsize = n
         usize = m
      else
c     %-------------------------%
c     | u0 is to be an n-vector |
c     %-------------------------%
         rsize = m
         usize = n
      endif

      idist = 2


      ierr = 0
      do itry=1,ntry
C       	 print *, "rsize = ", rsize

         call zlarnv(idist, iseed, rsize, work)   ! work of size rsize

C          do ips=1,rsize
C             print *, "work(", ips, ") = ", work(ips)
C          end do

         nrm = dznrm2(rsize,work,1)
C          print *, "zgetu0w: nrm = ", nrm
         call second(t2)
C            do i=1, usize
C               print *, "u0(", i, ") = ", u0(i)
C            end do

         call aprodw(transa,ndp,m,n,work,u0,work(ndp+1),
     c        work(2*ndp+1),lambda,trlambda,planF,planB)
C          print *, "zgetu0w: after aprodw"
C          do ips = 1, 513
C             print *, "u0(", ips, ") = ", u0(ips)
C          end do

c   ! work of size 3*ndp

C            do i=1, rsize
C               print *, "work(", i, ") = ", work(i)
C            end do

C            print *, "--------"
C            do i=1, usize
C               print *, "u0(", i, ") = ", u0(i)
C            end do

C            print *, "--------"
C            do i=ndp + 1, 2*ndp + 1
C               print *, "work(", i, ") = ", work(i)
C            end do

C            print *, "--------"
C            do i=2*ndp + 1, 3*ndp + 1
C               print *, "work(", i, ") = ", work(i)
C            end do

C            print *, "--------"
C            print *, "stifle!"


         call second(t3)
         tmvopx = tmvopx + (t3-t2)
         nopx = nopx+1

         u0norm = dznrm2(usize,u0,1)
         anormest = u0norm/nrm

C          print *, "zgetu0w: u0norm = ",u0norm,"anormest = ",anormest


         index(1) = 1
         index(2) = j
         index(3) = j+1
c
c         if  (MGS.eq.0) then
            call zreorth(usize,j,U,ldu,u0,u0norm,index,kappa,work,
     c           icgs)

C          do i=1, n
C             print *, "u0(", i, ") = ", u0(i)
C          end do
C          print *, "--------"

C          print *, "u0norm = ", u0norm
C          print *, "--------"
C          print *, "stifle!"



c         else
c            call zreorth2(usize,j,U,ldu,u0,u0norm,index)
c         endif
         if (u0norm.gt.0) goto 9999
      enddo
      ierr = -1
 9999 call second(t2)
      tgetu0 = tgetu0 + (t2-t1)
      return
      end
      
