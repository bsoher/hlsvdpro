#
#     (C) Brian J Soher, 2020
#



def zsafescal(n,alpha,x):
    """ Scale the vector x by 1/alpha avoiding unnecessary under- and overflow. """
    
    # PS - This code replicates the Fortran in the case where 
    # (abs(alpha) > SFMIN), and SFMIN = dlamch('s'). When that's not the case,
    # the Fortran code calls the LAPACK routine zlascl() which produces the 
    # same result but contains some special code to avoid under/overflow.
    # scipy doesn't provide access to zlascl(), so I just use the simple 
    # implementation for all cases.
    #
    # FIXME bjs cython_lapack bindings now provide c/z/d/slascl, could revisit

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    x[:n] /= alpha0


#       subroutine zsafescal(n,alpha,x)
# c
# c     Scale the vector x by 1/alpha avoiding unnecessary under- and overflow.
# c
# 
# c     %-----------%
# c     | Arguments |
# c     %-----------%
#       implicit none
#       integer n
#       double precision alpha
#       complex*16 x(*)
# 
# c     %------------%
# c     | Parameters |
# c     %------------%
#       double precision one, zero
#       parameter(one = 1.0, zero = 0.0)
# 
# c     %-----------------%
# c     | Local variables |
# c     %-----------------%
#       integer i,info
#       double precision sfmin
# 
# c     %----------------------%
# c     | External Subroutines |
# c     %----------------------%
#       external pzdscal,zlascl
# 
# c     %--------------------%
# c     | External Functions |
# c     %--------------------%
#       double precision dlamch
#       external dlamch
# 
# c     %-----------------%
# c     | Data statements |
# c     %-----------------%
#       save
#       data sfmin /-1d0/
#       
#       if (sfmin.eq.-1d0) then         
#          sfmin = dlamch('s')
#       endif
# 
#       if (abs(alpha).ge.sfmin) then
#          call pzdscal(n,one/alpha, x, 1)
#       else
#          call zlascl('General',i,i,alpha,one,n,1,x,n,info)
#       endif
#       
#       end
