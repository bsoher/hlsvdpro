#
#     (C) Brian J Soher, 2020
#

from __future__ import division, print_function

import math
import numpy as np
import scipy.linalg.blas as blas
import scipy.linalg.lapack as lapack

#from dlartg import dlartg

EPS = np.finfo(np.float64).eps
EPS34 = EPS ** (3.0 / 4.0)



def dbsvdstep(jobu,jobv,m,n,k,sigma,D,E,U,ldu,V,ldv):
    """
    Perform one implicit LQ SVD sweep with shift SIGMA.

    Used in 'IRL' SVD code
    """

    # FIXME bjs, NOT PORTED FULLY, index checks, etc.

    if k<=1: return

    dou = jobu == 'y'
    dov = jobv == 'y'
    
    #     Compute the initial rotation based on B*B^T-sigma^2
    x = D[0]*D[0] - sigma*sigma                                     
    y = E[0]*D[0]
    
    #     Chase the "bulge" down the lower bidiagonal with Givens rotations.
    #     Below 'y' is the "bulge" and 'x' is the element used to eliminate it.
    for i in range(k-1):

        if i>0:
            c,s,r = lapack.dlartg(x,y)
            E[i-1] = r
        else:
            c,s,r = lapack.dlartg(x,y)
            
        x      =  c*D[i] + s*E[i]
        E[i]   = -s*D[i] + c*E[i]
        D[i]   = x
        y      = s*D[i+1]
        D[i+1] = c*D[i+1]
        
        if dou and m>0:
            #             drot(n,dx,  incx, dy,  incy,c,s)
            U[:m,i], U[:m,i+1] = blas.drot(U[:m,i],U[:m,i+1],c,s)

        # call dlartg(x,y,c,s,D[i])
        c,s,r =  lapack.dlartg(x,y)
        D[i] = r
        
        x      =  c*E[i] + s*D[i+1]
        D[i+1] = -s*E[i] + c*D[i+1]
        E[i]   = x
        y      = s*E[i+1]
        E[i+1] = c*E[i+1]
        
        if dov and n>0:
            # drot(n,V[1,i],1,V[1,i+1],1,c,s)
            V[:n,i], V[:n,i+1] = blas.drot(V[:n,i],V[:n,i+1],c,s)

    # call dlartg(x,y,c,s,E[k-1])
    c,s,r = lapack.dlartg(x,y)
    E[k-1] = r
    
    x    =  c*D[k] + s*E[k]
    E[k] = -s*D[k] + c*E[k]
    D[k] = x

    if (dou and m>0):
        # drot(m,U[1,k],1,U[1,k+1],1,c,s)
        U[:m,k], U[:m,k+1] = blas.drot(m,U[:,k],1,U[:,k+1],1,c,s)

    return 



def dbdqr(ignorelast, jobq, n, d, e, c1, c2, Qt0, ldq):
    """
    c Compute QR factorization B = Q*R of (n+1) x n lower bidiagonal matrix 
    c with diagonal elements d(1)...d(n) and first subdiagonal elements
    c e(1)...e(n). On return [0 ... 0 c1 c2]' = Q'*[0 ... 0 1]'.
    c If ignorelast.eq..true. then e(n) is assumed to be zero.
    c
    c If jobq=='Y' then on return Qt contains Q^T.

    NB. bjs, changed D to d and E to e in params to work in code.

    """
    if len(Qt0.shape) == 1:
        Qt = np.zeros((n+1,n+1), dtype=np.float64)
        write_back = True
    else:
        Qt = Qt0
        write_back = False

    if n >=1 :

        if jobq in ['y','Y']:
            Qt[:n+1, :n+1] *= 0.0
            for j in range(n+1):
                Qt[j, j] = 1.0

        for i in range(n-1):

            cs,sn,r = lapack.dlartg(d[i],e[i])

            d[i]   = r
            e[i]   = sn*d[i+1]      # e[i] = sn * d.item(i + 1)
            d[i+1] = cs*d[i+1]      # d[i + 1] = cs * d.item(i + 1)
            if jobq in ['y','Y']:
                # for j in range(i+1):
                #     Qt[i+1,j] = -sn*Qt[i,j]
                #     Qt[i,j]   =  cs*Qt[i,j]

                # PS - Here's the fast version --
                Qt[i+1, :i+1] = -sn*Qt[i, :i+1]
                Qt[i  , :i+1] =  cs*Qt[i, :i+1]

                Qt[i,  i+1] = sn
                Qt[i+1,i+1] = cs

        # In fortran, a loop ends when the counter is one past the limit. For
        # example, at the end of the loop  `do i = 1,50`, i == 51. In Python,
        # this is not the case so I have to increment i here.#
        i += 1

        if not ignorelast:
            # call dlartg(d[n],e[n],cs,sn,r)
            cs,sn,r = lapack.dlartg(d[n-1],e[n-1])
            d[n-1] = r
            e[n-1] = 0.0
            c1 = sn
            c2 = cs
            if jobq in ['y','Y']:
                # for j in range(i+1):
                #     Qt[i+1,j] = -sn*Qt[i,j]
                #     Qt[i,j]   =  cs*Qt[i,j]

                # PS - Here's the fast version --
                Qt[i+1, :i+1] = -sn*Qt[i, :i+1]
                Qt[i  , :i+1] =  cs*Qt[i, :i+1]

                Qt[i,  i+1] = sn
                Qt[i+1,i+1] = cs

    if write_back:
        for j in range(n + 1):
            for i in range(n + 1):
                Qt0[i+(n+1)*j] = Qt[i,j]

    return c1, c2



def dbdqr1(jobq, n, d, e, qt, ldq):
    """
    from v1


    # c Compute QR factorization B = Q*R of (n+1) x n lower bidiagonal matrix
    # c with diagonal elements d(1)...d(n) and first subdiagonal elements
    # c e(1)...e(n). On return [0 ... 0 c1 c2]' = Q'*[0 ... 0 1]'.
    # c
    # c If jobq=='Y' then on return Qt contains Q^T.

    """
    jobq = jobq.upper()
    qt = qt.reshape(ldq, -1)

    # PS - tested for when jobq = 'N', not tested for jobq = 'Y'


    msg = "{} qt[{}, {}] = {}"
    pqt = lambda letter, i1, i2, value: msg.format(letter, i1 + 1, i2 + 1,
                                                   value)

    if n >= 1:
        if jobq == 'Y':
            qt[:n + 1, : n + 1] *= 0.0
            for j in range(n + 1):
                qt[j, j] = 1.0

        for i in range(n - 1):
            cs, sn, r = lapack.dlartg(d[i], e[i])

            d[i] = r
            e[i] = sn * d.item(i + 1)
            d[i + 1] = cs * d.item(i + 1)

            if jobq == 'Y':
                # PS - Here's the fast version --
                qt[i+1, :i+1] = -sn * qt[i, :i+1]
                qt[i, :i+1] = cs * qt[i, :i+1]

                qt[i, i+1] = sn
                qt[i+1, i+1] = cs

        # In fortran, a loop ends when the counter is one past the limit. For
        # example, at the end of the loop  `do i = 1,50`, i == 51. In Python,
        # this is not the case so I have to increment i here.
        i += 1

        cs, sn, r = lapack.dlartg(d[n - 1], e[n - 1])

        d[n - 1] = r
        e[n - 1] = 0.0
        c1 = sn
        c2 = cs

        if jobq == 'Y':

            # PS - Here's the fast version --
            qt[i+1, :i+1] = -sn * qt[i, :i+1]
            qt[i, :i+1] = cs * qt[i, :i+1]

            qt[i, i+1] = sn
            qt[i+1, i+1] = cs

    return c1, c2


#------------------------------------------------------------------------------

def refinebounds(n, theta, bound, tol, eps34):
    """
         Refine Lanczos error bounds using the gap theorem.

    """
    if n > 1:
        for i in range(n):
            for j in (-1, 1):
                if ((j == 1) and (i < n)) or ((j == -1) and (i > 1)):

                    if abs(theta[i] - theta[i + j]) < (eps34 * theta[i]):
                        if (bound[i] > tol) and (bound[i + j] > tol):
                            bound[i + j] = np.sqrt(bound[i]*bound[i]+bound[i+j]*bound[i+j])
                            bound[i] = 0.0

            gap = theta[i + 1] - bound[i + 1] - (theta[i] + bound[i])

            if gap > bound[i]:
                bound[i] *= (bound[i] / gap)



def drefinebounds(n,k,theta,bound,tol,eps34):
    """ Refine Lanczos error bounds using the gap theorem.  
     
    c     Input arguments: 
    c              n:     smallest dimension of original matrix
    c              k:     number of Ritz values to refine
    c              theta: array of Ritz values
    c              bound: array of unrefined error bounds
    c              tol:   clustering tolerance
    c              eps34: machine epsilon to the power 3/4.

    """
    if k > 1:
    
        for i in range(k):
            for l in [-1,1]:
                if (l==1 and i<k) or (l==-1 and i>1):
                    if (np.abs(theta[i]-theta[i+l]) < eps34*(theta[i])):
                        if (bound[i]>tol and bound[i+l]>tol):
                            bound[i+l] = np.sqrt(bound[i]*bound[i] + bound[i+l]*bound[i+l])
                            bound[i]   = 0.0

        for i in range(k):
            if i<k or k==n:

                # We cannot compute a reliable value for the gap of the last
                # Ritz value unless we know it is an approximation to the
                # smallest singular value (k.eq.n). In this case we can take the
                # distance to the next bigger one as the gap, which can really
                # save us from getting stuck on matrices with a single isolated tiny
                # singular value.

                if (i==0):
                    gap = np.abs(theta[i]-theta[i+1]) - max(bound[i],bound[i+1])
                elif (i==(n-1)):
                    gap = np.abs(theta[i-1]-theta[i]) - max(bound[i-1],bound[i])
                else:
                    gap = np.abs(theta[i]-theta[i+1]) - max(bound[i],bound[i+1])
                    gap = min(gap,np.abs(theta[i-1] - theta[i]) - max(bound[i-1],bound[i]))

                if gap>bound[i]:
                    bound[i] *=  (bound[i]/gap)
      
    return



