#
#  (C) Brian J Soher, 2020
#

from __future__ import division, print_function

import scipy.linalg.blas as blas


#****************************************************************************

def pdnrm2(n, x, incx):
    
    r = blas.dnrm2(x[:n],incx=incx)
    
    return r

#       double precision function pdnrm2(n, x, incx)
#       implicit none
#       integer n, incx
#       double precision x(*), dnrm2
#       external dnrm2
# 
#       pdnrm2 = dnrm2(n, x, incx)
#       end


#****************************************************************************

def pdscal(n, alpha, x , incx):

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    x[:n] = blas.dscal(alpha0,x[:n],incx=incx)
    
    return 
    
#       subroutine pdscal(n, alpha, x , incx)
#       implicit none
#       integer n, incx
#       double precision alpha,x(*)
# 
#       call dscal(n, alpha, x , incx)
#       end


#****************************************************************************

def pdcopy(n, x , incx, y, incy):
    
    y[:n] = blas.dcopy(x[:n],y[:n],incx=incx,incy=incy)
    
    return 
    
    
#       subroutine pdcopy(n, x , incx, y, incy)
#       implicit none
#       integer n, incx, incy
#       double precision x(*),y(*)
# 
#       call dcopy(n, x , incx, y, incy)
#       end


#****************************************************************************

def pdaxpy(n, alpha, x , incx, y, incy):

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    y[:n] = blas.daxpy(x[:n],y[:n],a=alpha0,incx=incx,incy=incy)
    
    return 

#       subroutine pdaxpy(n, alpha, x , incx, y, incy)
#       implicit none
#       integer n, incx, incy
#       double precision alpha,x(*),y(*)
#      
#       call daxpy(n, alpha, x , incx, y, incy)
#       end


#****************************************************************************

def pddot(n, x , incx, y, incy):      
    
    r = blas.ddot(x[:n],y[:n],incx=incx,incy=incy)
    
    return r

#       double precision function pddot(n, x , incx, y, incy)
#       implicit none
#       integer n, incx, incy
#       double precision x(*),y(*),ddot
#       external ddot
#       
#       pddot = ddot(n, x , incx, y, incy)
#       end


#**************************************************************************** 

def pdzero(n, x , incx):    
    
    if n > 0 and incx != 0:
        if incx == 1:
            x[:n] *= 0.0
        else:
            for i in range(n):
                x[i*incx] = 0.0 

    return 
     
#       subroutine pdzero(n, x , incx)
#       implicit none
#       integer n, incx
#       double precision x(*)
#       
#       integer i
# 
#       if ((n.gt.0).and.(incx.ne.0)) then         
#          if (incx.eq.1) then
#             do i=1,n
#                x(i) = 0d0
#             enddo
#          else
#             do i=1,n
#                x(1+(i-1)*incx) = 0d0
#             enddo
#          endif
#       endif
#       return
#       end


#****************************************************************************

def pizero(n, x , incx):
    
    if n > 0 and incx != 0:
        if incx == 1:
            x[:n] *= 0
        else:
            for i in range(n):
                x[i*incx] = 0 

    return 
            
#       subroutine pizero(n, x , incx)
#       implicit none
#       integer n, incx
#       integer x(*)
#       
#       integer i
# 
#       if ((n.gt.0).and.(incx.ne.0)) then         
#          if (incx.eq.1) then
#             do i=1,n
#                x(i) = 0
#             enddo
#          else
#             do i=1,n
#                x(1+(i-1)*incx) = 0
#             enddo
#          endif
#       endif
#       return
#       end


#****************************************************************************

def pdset(n, alpha, x, incx):

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    for i in range(n):
        x[i*incx] = alpha0
    
    return 

#       subroutine pdset(n, alpha, x , incx)
#       implicit none
#       integer n, incx
#       double precision alpha,x(*)
#       
#       integer i
# 
#       if ((n.gt.0).and.(incx.ne.0)) then         
#          if (incx.eq.1) then
#             do i=1,n
#                x(i) = alpha
#             enddo
#          else
#             do i=1,n
#                x(1+(i-1)*incx) = alpha
#             enddo
#          endif
#       endif
#       return
#       end


#****************************************************************************
 
def pdaxpby(n,alpha,x,incx,beta,y,incy): 
    """ Y = alpha*X + beta*Y """

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    try:
        beta0 = beta[0]
    except:
        beta0 = beta

    if n<=0 or incy==0 or incx==0:    # bjs-check, should I return x?
        return y
    
    if alpha0 == 0 and beta0 == 0:
        if incy == 1:
            y[:n] *= 0.0
        else:
            for i in range(n):
                y[i*incy] = 0.0
    else:
        if incx==1 and incy==1:
            y[:n] = alpha0*x[:n] + beta0*y[:n]
        else:
            for i in range(n):
                y[i*incy] = alpha0*x[i*incx] + beta0*y[i*incy]

    return     

    
#       subroutine pdaxpby(n,alpha,x,incx,beta,y,incy)
# #
# #     Y = alpha*X + beta*Y
# #     
# 
#       implicit none
#       double precision one,zero
#       parameter(one = 1d0,zero = 0d0)
#       integer n,incx,incy,i
#       double precision alpha,beta,x(n),y(n)
# 
#       if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
#       if (alpha.eq.zero .and. beta.eq.zero) then
#          if (incy.eq.1) then
#             do i=1,n
#                y(i) = zero
#             enddo
#          else
#             do i=1,n
#                y(1+(i-1)*incy) = zero
#             enddo
#          endif
#          
#       else if (alpha.eq.zero .and. beta.ne.zero) then
#          
#          call pdscal(n,beta,y,incy)
# 
#       else if (alpha.ne.zero .and. beta.eq.zero) then
# 
#          if (alpha.eq.one) then
#             call pdcopy(n,x,incx,y,incy)
#          else
#             if (incx.eq.1 .and. incy.eq.1) then
#                do i=1,n
#                   y(i) = alpha*x(i)
#                enddo
#             else
#                do i=1,n
#                   y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)
#                enddo
#             endif
#          endif
# 
#       else
# 
#          if (beta.eq.one) then
# c DAXPY
#             call pdaxpy(n,alpha,x,incx,y,incy)
#          else
#             if (incx.eq.1 .and. incy.eq.1) then
#                do i=1,n
#                   y(i) = alpha*x(i) + beta*y(i)
#                enddo
#             else
#                do i=1,n
#                   y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx) + 
#      c                 beta*y(1+(i-1)*incy)
#                enddo
#             endif
#          endif
#       endif
#       return
#       end


#****************************************************************************

def pdaxty(n,alpha,x,incx,y,incy):
    """ Y = alpha*X*Y """

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    if n<=0 or incy==0 or incx==0:    # bjs-check, should I return x?
        return y
    
    if alpha0 == 0:
        if incy == 1:
            y[:n] *= 0.0
        else:
            for i in range(n):
                y[i*incy] = 0.0
    else:
        
        if alpha0 == 1:
            if incx == 1 and incy == 1:         # bjs-check, optimize, could shrink logic
                y[:n] = x[:n] * y[:n]
            else:
                for i in range(n):
                    y[i*incy] = x[i*incx]*y[i*incy]
        else:
            if incx == 1 and incy == 1:
                y[:n] = alpha0*x[:n]*y[:n]
            else:
                for i in range(n):
                    y[i*incy] = alpha0*x[i*incx]*y[i*incy]

    return y    
    

#       subroutine pdaxty(n,alpha,x,incx,y,incy)
# #
# #     Y = alpha*X*Y
# #     
# 
#       implicit none
#       double precision one,zero
#       parameter(one = 1d0,zero = 0d0)
#       integer n,incx,incy,i
#       double precision alpha,x(n),y(n)
# 
#       if (n.le.0 .or. incy.eq.0 .or. incx.eq.0) return
#       if (alpha.eq.zero) then
#          if (incy.eq.1) then
#            do i=1,n
#                y(i) = zero
#             enddo
#          else
#            do i=1,n
#                y(1+(i-1)*incy) = zero
#             enddo
#          endif
#          
#       else if (alpha.ne.zero) then
# 
#          if (alpha.eq.one) then
#             if (incx.eq.1 .and. incy.eq.1) then
#                do i=1,n
#                   y(i) = x(i)*y(i)
#                enddo
#             else
#                do i=1,n
#                   y(1+(i-1)*incy) = x(1+(i-1)*incx)*y(1+(i-1)*incy)
#                enddo
#             endif
# 
#          else
#             if (incx.eq.1 .and. incy.eq.1) then
#                do i=1,n
#                   y(i) = alpha*x(i)*y(i)
#                enddo
#             else
#                do i=1,n
#                   y(1+(i-1)*incy) = alpha*x(1+(i-1)*incx)*
#      c                 y(1+(i-1)*incy)
#                enddo
#             endif
#          endif
#       endif
#       return
#       end


#****************************************************************************

def szero(n, x , incx):
    
    x[:n] = dzero(n,x[:n],incx)
    
    return x
    
#       subroutine szero(n, x , incx)
#       implicit none
#       integer n, incx
#       real x(*),zero
#       parameter (zero = 0)          
#       integer i
# 
#       if ((n.gt.0).and.(incx.ne.0)) then         
#          if (incx.eq.1) then
#             do i=1,n
#                x(i) = zero
#             enddo
#          else
#             do i=1,n
#                x(1+(i-1)*incx) = zero
#             enddo
#          endif
#       endif
#       return
#       end


#****************************************************************************
                 
def dzero(n, x , incx):

    if n > 0 and incx != 0:         
        if incx == 1:
            x[:n] *= 0.0
        else:
            for i in range(n):
                x[i*incx] = 0.0
    return 

#       subroutine dzero(n, x , incx)
#       implicit none
#       integer n, incx
#       double precision x(*),zero
#       parameter (zero = 0.0)          
#       integer i
# 
#       if ((n.gt.0).and.(incx.ne.0)) then         
#          if (incx.eq.1) then
#             do i=1,n
#                x(i) = zero
#             enddo
#          else
#             do i=1,n
#                x(1+(i-1)*incx) = zero
#             enddo
#          endif
#       endif
#       return
#       end


#****************************************************************************

def izero(n, x , incx):

    if n > 0 and incx != 0:         
        if incx == 1:
            x[:n] *= 0
        else:
            for i in range(n):
                x[i*incx] = 0
    return x 

#       subroutine izero(n, x , incx)
#       implicit none
#       integer n, incx
#       integer x(*),zero
#       parameter (zero = 0)          
#       integer i 
#
#       if ((n.gt.0).and.(incx.ne.0)) then         
#          if (incx.eq.1) then
#             do i=1,n
#                x(i) = zero
#             enddo
#          else
#             do i=1,n
#                x(1+(i-1)*incx) = zero
#             enddo
#          endif
#       endif
#       return
#       end

         
