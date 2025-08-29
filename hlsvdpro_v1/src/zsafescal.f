      subroutine zsafescal(n,alpha,x)

c     %-----------%
c     | Arguments |
c     %-----------%
      implicit none
      integer n
      double precision alpha
      complex*16  x(*)

c     %------------%
c     | Parameters |
c     %------------%
      double precision one, zero
      parameter(one = 1.0, zero = 0.0)

c     %-----------------%
c     | Local variables |
c     %-----------------%
      integer i,info
      double precision sfmin

c     %----------------------%
c     | External Subroutines |
c     %----------------------%
      external zscal,zlascl

c     %--------------------%
c     | External Functions |
c     %--------------------%
      double precision dlamch
      external dlamch


      sfmin = dlamch('s')

      if (abs(alpha).ge.sfmin) then
         call zscal(n,one/alpha, x, 1)
      else
         call zlascl('General',i,i,alpha,one,n,1,x,n,info)
      endif
      
      end
