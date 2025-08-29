c
c     (C) Brian J Soher, April 3,2020
c


ccccccccccccccccccccc dense format cccccccccccccccccccccccccccc
c
c     Note: We assume the matrix is always stored in major 
c     column order.
c

      subroutine caprod(transa,m,n,x,y,cparm,iparm)

      implicit none

      character*1 transa
      integer m,n,iparm(*)
      complex x(*),y(*),cparm(*)

      complex, pointer :: pA(:,:)  
      common /csvdp/ pA     

      external cgemv_h
     
      call cgemv_h(transa,m,n,cmplx(1d0,0d0),pA,m,
     c     x,1,cmplx(0d0,0d0),y,1)

      end


