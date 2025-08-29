c
c     (C) Brian J Soher, April 3,2020
c


ccccccccccccccccccccc dense format cccccccccccccccccccccccccccc
c
c     Note: We assume the matrix is always stored in major 
c     column order.
c

      subroutine daprod(transa,m,n,x,y,dparm,iparm)

      implicit none

      character*1 transa
      integer m,n,iparm(*)
      double precision x(*),y(*),dparm(*)

      double precision, pointer :: pA(:,:)  
      common /csvdp/ pA     

      external dgemv_h
     
      call dgemv_h(transa,m,n,1d0,pA,m,x,1,0d0,y,1)

      end


