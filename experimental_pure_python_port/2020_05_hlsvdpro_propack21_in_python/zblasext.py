#
#  (C) Brian J Soher, 2020
#

from __future__ import division, print_function

import numpy as np

import scipy.linalg.blas as blas



def pdznrm2(n, x, incx):

    # bjs - these give wrong answer in 15th decimal place
    #       don't know why the 'copied' blas code is different, but
    #       I recreate it here to compensate
    # r1 = np.sqrt((x[:n].conjugate() * x[:n]).sum()).real
    r2 = blas.dznrm2(x[:n], incx=incx)

#    r3 = my_dznrm2(n, x, incx)
    
    return r2

def my_dznrm2(n,x,incx):

    scale = 0.0
    ssq = 1.0
    for ix in range(n):
        if x[ix].real != 0.0:
            temp = np.abs(x[ix].real)
            if scale < temp:
              ssq = 1.0 + ssq*(scale/temp)**2
              scale = temp
            else:
                ssq = ssq + (temp/scale)**2

        if x[ix].imag != 0.0:
            temp = np.abs(x[ix].imag)
            if scale < temp:
                ssq = 1.0 + ssq*(scale/temp)**2
                scale = temp
            else:
                ssq = ssq + (temp/scale)**2

    norm = scale*np.sqrt(ssq)

    return norm

#            DO 10 IX = 1,1 + (N-1)*INCX,INCX
#               IF (DBLE(X(IX)).NE.ZERO) THEN
#                   TEMP = ABS(DBLE(X(IX)))
#                   IF (SCALE.LT.TEMP) THEN
#                       SSQ = ONE + SSQ* (SCALE/TEMP)**2
#                       SCALE = TEMP
#                   ELSE
#                       SSQ = SSQ + (TEMP/SCALE)**2
#                   END IF
#               END IF
#               IF (DIMAG(X(IX)).NE.ZERO) THEN
#                   TEMP = ABS(DIMAG(X(IX)))
#                   IF (SCALE.LT.TEMP) THEN
#                       SSQ = ONE + SSQ* (SCALE/TEMP)**2
#                       SCALE = TEMP
#                   ELSE
#                       SSQ = SSQ + (TEMP/SCALE)**2
#                   END IF
#               END IF
#    10     CONTINUE
#           NORM = SCALE*SQRT(SSQ)
#       END IF
# *
#       DZNRM2 = NORM

#----------------------------------------------------------------------
 
def pzscal(n, alpha, x , incx):

    x[:] = blas.zscal(alpha, x[:n], incx=incx)
    
    return 


#----------------------------------------------------------------------

def pzdscal(n, alpha, x , incx):

    x[:] = blas.zdscal(alpha, x[:n], incx=incx)

    return 


#----------------------------------------------------------------------

def pzcopy(n, x , incx, y, incy):

    y[:n] = blas.zcopy(x[:n],y[:n],incx=incx,incy=incy)

    return 


#----------------------------------------------------------------------

def pzaxpy(n, alpha, x , incx, y, incy):

    y[:n] = blas.zaxpy(x[:n],y[:n],a=alpha,incx=incx,incy=incy)

    return 


#----------------------------------------------------------------------

def pzdaxpy(n, alpha, x , incx, y, incy):

    if n>0 and incx!=0 and incy!=0:
        if incx==1 and incy==1:
            y[:n] = alpha*x[:n] + y[:n]
        else:
            for i in range(n):
                y[i*incy] = alpha*x[i*incx] + y[i*incy]     # FIXME bjs, list comp
    return


#----------------------------------------------------------------------
    
def pzdotc(n, x , incx, y, incy):

#    sum = 0+0j
#    for i in range(n):
# bjs        print('   dconjg(x[i]), y(i) = ',np.conj(x[i]),y[i])
#        sum = sum + np.conjugate(x[i]) * y[i]
# bjs        print('sum = ',sum)
#    return sum

#    return np.vdot(x[::incx],y[::incy])

    return blas.zdotc(x[:n],y[:n],incx=incx,incy=incy)




#----------------------------------------------------------------------

def pzdotu(n, x , incx, y, incy):

    return blas.zdotu(x[:n],y[:n],incx=incx,incy=incy)



#----------------------------------------------------------------------
       
def pzzero(n, x, incx):

    if n>0 and incx!=0:
        if incx==1:
            x[:n] *= 0.0+0.0j
        else:
            for i in range(n):
                x[i*incx] = 0.0+0.0j
    return


#----------------------------------------------------------------------

def pzset(n, alpha, x, incx):

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    if n>0 and incx!=0:
        if incx==1:
            x[:n] = x[:n]*0 + alpha0
        else:
            for i in range(n):
                x[(i-1)*incx] = alpha0
    return


#----------------------------------------------------------------------
  
def pzdaxpby(n,alpha,x,incx,beta,y,incy):
    """ Y = alpha*X + beta*Y """    

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    try:
        beta0 = beta[0]
    except:
        beta0 = beta

    if n <= 0 or incy == 0 or incx == 0:
        return
    
    if alpha0 == 0 and beta0 == 0:
        if incy==1:
            y[:n] *= 0.0+0.0j
        else:
            for i in range(n):
                y[i*incy] = 0.0+0.0j
    else:
        for i in range(n):
            y[i*incy] = alpha0*x[i*incx] + beta0*y[i*incy]
    
    return


#----------------------------------------------------------------------

def pzaxpby(n,alpha,x,incx,beta,y,incy):
    """ Y = alpha*X + beta*Y """    

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    try:
        beta0 = beta[0]
    except:
        beta0 = beta

    if n <= 0 or incy == 0 or incx == 0:
        return
    
    if alpha0 == 0.0+0.0j and beta0 == 0.0+0.0j:
        if incy==1:
            y[:n] *= 0.0+0.0j
        else:
            for i in range(n):
                y[i*incy] = 0.0+0.0j
    else:
        for i in range(n):
            y[i*incy] = alpha0*x[i*incx] + beta0*y[i*incy]
    
    return


#----------------------------------------------------------------------

def pzaxty(n,alpha,x,incx,y,incy):
    """ Y = alpha*X*Y """    

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    if n <= 0 or incy == 0 or incx == 0:
        return
     
    if alpha0 == 0.0+0.0j:
        if incy==1:
            for i in range(n):
                y[:n] *= 0.0+0.0j
        else:
            for i in range(n):
                y[i*incy] = 0.0+0.0j
    else: 
        if incx==1 and incy==1:
            y[:n] = alpha0*x[:n]*y[:n]
        else:
            for i in range(n):
                y[i*incy] = alpha0*x[i*incx]* y[i*incy]
    
    return


#----------------------------------------------------------------------

def pzdaxty(n,alpha,x,incx,y,incy):
    """ Y = alpha*X*Y """

    try:
        alpha0 = alpha[0]
    except:
        alpha0 = alpha

    if n <= 0 or incy == 0 or incx == 0:
        return
    
    if alpha0 == 0:
        if incy==1:
            y[:n] *= 0.0+0.0j
        else:
            for i in range(n):
                y[i*incy] = 0.0+0.0j
    else:
        if incx==1 and incy==1:
            y[:n] = alpha0*x[:n]*y[:n]
        else:
            for i in range(n):
                y[i*incy] = alpha0*x[i*incx]* y[i*incy]
    return

#----------------------------------------------------------------------
             
def zzero(n, x , incx):
    
    if n>0 and incx!=0:
        if incx==1:
            x[:n] *= 0.0+0.0j
        else:
            for i in range(n):
                x[i*incx] = 0.0+0.0j
    return


