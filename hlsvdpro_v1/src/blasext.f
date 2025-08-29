      subroutine szero(n, x , incx)
      implicit none
      integer n, incx
      real x(*),zero
      parameter (zero = 0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end

                 
      subroutine dzero(n, x , incx)
      implicit none
      integer n, incx
      double precision x(*),zero
      parameter (zero = 0.0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end


      subroutine czero(n, x , incx)
      implicit none
      integer n, incx
      complex x(*),zero
      parameter (zero = 0.0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end

      subroutine zzero(n, x , incx)
      implicit none
      integer n, incx
      double complex x(*),zero
      parameter (zero = 0.0)
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end

      subroutine izero(n, x , incx)
      implicit none
      integer n, incx
      integer x(*),zero
      parameter (zero = 0)          
      integer i

      if ((n.gt.0).and.(incx.ne.0)) then         
         if (incx.eq.1) then
            do i=1,n
               x(i) = zero
            enddo
         else
            do i=1,n
               x(1+(i-1)*incx) = zero
            enddo
         endif
      endif
      return
      end

         
