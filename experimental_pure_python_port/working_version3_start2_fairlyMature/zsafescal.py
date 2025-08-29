# Python modules
from __future__ import division

# 3rd party modules

# Our modules


#       subroutine zsafescal(n,alpha,x)

def zsafescal(n, alpha, x):
    # PS - This code replicates the Fortran in the case where 
    # (abs(alpha) > SFMIN), and SFMIN = dlamch('s'). When that's not the case,
    # the Fortran code calls the LAPACK routine zlascl() which produces the 
    # same result but contains some special code to avoid under/overflow.
    # scipy doesn't provide access to zlascl(), so I just use the simple 
    # implementation for all cases.

    x /= alpha



# c     %-----------%
# c     | Arguments |
# c     %-----------%
#       implicit none
#       integer n
#       double precision alpha
#       complex*16  x(*)

# c     %------------%
# c     | Parameters |
# c     %------------%
#       double precision one, zero
#       parameter(one = 1.0, zero = 0.0)

# c     %-----------------%
# c     | Local variables |
# c     %-----------------%
#       integer i,info
#       double precision sfmin

# c     %----------------------%
# c     | External Subroutines |
# c     %----------------------%
#       external zscal,zlascl

# c     %--------------------%
# c     | External Functions |
# c     %--------------------%
#       double precision dlamch
#       external dlamch


#       sfmin = dlamch('s')

#       if (abs(alpha).ge.sfmin) then
#          call zscal(n,one/alpha, x, 1)
#       else
#          call zlascl('General',i,i,alpha,one,n,1,x,n,info)
#       endif
      
#       end
