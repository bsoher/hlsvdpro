c
c     (C) Brian J Soher, April 3,2020
c


ccccccccccccccccccccc dense format cccccccccccccccccccccccccccc
c
c     Note: We assume the matrix is always stored in major 
c     column order.
c

      subroutine zaprod(transa,m,n,x,y,zparm,iparm)

      implicit none

      character*1 transa
      integer m,n,iparm(*)
      complex*16 x(*),y(*),zparm(*)

      complex*16, pointer :: pA(:,:)  
      common /csvdp/ pA     

      external zgemv_h
     
      call zgemv_h(transa,m,n,dcmplx(1d0,0d0),pA,m,
     c     x,1,dcmplx(0d0,0d0),y,1)

      end


