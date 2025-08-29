# Python modules
from __future__ import division
import ctypes
import math
import pprint
import pdb
pp = pprint.pprint

# 3rd party modules
import numpy as np
import numpy.linalg
import scipy
import scipy.linalg.blas
import scipy.linalg.lapack

# Our modules
import zget0w
#import dlartg
import dbdsqr
import util
import zgemmina
import zlanbprow

MAX_SINGULAR_VALUES = 50
EPS = np.finfo(np.float64).eps
EPS34 = EPS ** (3.0 / 4.0)

# PS - These control flow in a few parts of the code. They're hardcoded in
# the Fortran we have, although they're probably passed as params in the 
# original PROPACK code.

JOBU = 'y'
JOBV = 'n'

# Variable name translations 
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# k          nsv_sought         # of singular values requested by caller
# kfit       nsv_found          # of singular values found
# lsinval    singular_values    Array containing the singular values (floats)
# ndp        n_data_points      
# kmax       MAX_SINGULAR_VALUES  Hardcoded to 50.
#           
#            options dict key
#            ----------------
# ioption(1) classical_gs       True ==> classical Gram-Schmidt reorth
#                               (see zreorth.py)
# ioption(2) extended_reorth_iterations  # of iterations in one case of 
#                               extended local reorth
# doption(1) delta              See zlanbprow() docstring
# doption(2) eta                See zlanbprow() docstring
# doption(3) anorm              See zlanbprow() docstring


# Float LAPACK functions
dlamch, dlartg = scipy.linalg.lapack.get_lapack_funcs( ['lamch', 'lartg'], np.array([0.0]))




def zlansvdw(n_data_points, m, n, nsv_sought, lambda_, trlambda, uuu, tol=16e-16):

    # PS - Fortran populates variable sfmin but otherwise doesn't use it.

    # options replaces the Fortran arrays ioption and doption.
    options = {}

    # Use classical Gram-Schmidt reorth instead of the modified version.
    options["classical_gs"] = True
    options["extended_reorth_iterations"] = 1
    options["delta"] = 1e-12
    options["eta"]   = 1e-14
    options["anorm"] = 0.0

    EPSN = max(m,n) * EPS / 2.0

    if uuu is None:
        uuu = np.zeros( (m, MAX_SINGULAR_VALUES + 1), np.complex128)

    vvv = np.zeros( (n, MAX_SINGULAR_VALUES), np.complex128)

    #  %--------------------------------%
    #  | Guard against absurd arguments |
    #  %--------------------------------%

    lanmax = min(n + 1, m + 1, MAX_SINGULAR_VALUES)
    tol = min(1.0, max(16.0 * EPS, tol))
    anorm = 0.0

    #  %------------------------------%
    #  | Set pointers into work array |
    #  %------------------------------%

    # The Fortran code uses multiple pointers into the arrays work and zwork
    # (and bnd). In Python I replace them with numpy arrays. Since I don't 
    # know what these arrays do, they have crummy names but at least they
    # harken back to the Fortran code.
    # Here's a map. The indices are a guide only; some of them will vary 
    # based on n_data_points.

    # Fortran        Fortran indices       Python
    # -------        ---------------       ------
    # work(ibnd)  =  work(   1 -   51)  =  work1
    # work(ib)    =  work(  52 -  151)  =  bbb_a, bbb_b
    # work(ib1)   =  work( 153 -  251)  =  bbb1_a, bbb1_b
    # work(ip)    =  work( 252 - 2852)  =  workp
    # work(iq)    =  work(2853 - 5352)  =  workq

    # Some of these lengths are adjusted by 1 to reflect the difference between
    # Fortran's 1-based arrays and Python's 0-based arrays. Since Python 
    # implements variable-length lists so nicely, it might be prefereable 
    # to convert some of these lists to dynamic sizing rather than 
    # preallocating them.
    
    len_work1  = lanmax + 1
    
    # In the Fortran code, work(ib) starts an array of length lanmax * 2. It's 
    # sometimes declared as matrix of dimension (lanmax, 2) but it's never used
    # as a matrix. It's really just two arrays (of length lanmax) bundled 
    # together under one variable name. It's confusing, and in Python I separate
    # them into bbb_a and bbb_b.
    # The chunk of work starting at work(ib1) is just used for a work copy of
    # work(ib), and has the same properties. I call the corresponding Python
    # arrays bbb1_a and bbb1_b.
    
    len_bbb = lanmax
    len_workp  = (lanmax + 1) ** 2
    
    # Note that before zlansdvw.f's second call to dbdsqr(), it writes a bunch
    # of 1s into work starting at index iq. It goes waaaay beyond lanmax**2
    # and into the work(iwrk) section. Intentional? I don't know. The size I
    # give to workq here should be about right, erring on the side of safety.
    
    len_workq  = n * (lanmax + 1)

    # PS - changed from lwork to lzwrk
    #  lwrk = lzwrk-iwrk+1
    #  call dzero(7*lanmax + 2 + 2*lanmax**2,work,1)
    #  call zzero(7*lanmax + 2 + 2*lanmax**2,zwork,1)

    work1  = np.zeros( (len_work1, ), dtype=np.float64)    # work(ibnd)
    workp  = np.zeros( (len_workp, ), dtype=np.float64)    # work(ip)
    workq  = np.zeros( (len_workq, ), dtype=np.float64)    # work(iq)

    bbb_a  = np.zeros( (len_bbb, ), dtype=np.float64)      # work(ib)
    bbb_b  = np.zeros( (len_bbb, ), dtype=np.float64)      # work(ib  + lanmax)
    bbb1_a = np.zeros( (len_bbb, ), dtype=np.float64)      # work(ib1)
    bbb1_b = np.zeros( (len_bbb, ), dtype=np.float64)      # work(ib1 + lanmax)


    #  %---------------------------------------------------------------%
    #  | Set up random starting vector if none is provided by the user |
    #  %---------------------------------------------------------------%

    # FIXME this follows the convention in the Fortran code of hardcoding
    # the value for rnorm which makes the if statement below superfluous. 

    rnorm = 0.0
    if rnorm == 0:

        rnorm, anorm, ierr = zget0w.zgetu0w('n', n_data_points, m, n, 0, 1, 
                                            uuu[:,0], 
                                            uuu, options["classical_gs"],
                                            lambda_, trlambda)
    info = 0
    neig = 0
    jold = 0
    j = min(nsv_sought + max(8, nsv_sought) + 1, lanmax)

    # %------------------------------%
    # | Iterate until convergence... |
    # %------------------------------%

    while neig < nsv_sought:
        # %---------------------------------------------------%
        # | Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j} |
        # %---------------------------------------------------%

        # zlanbprow() modifies bbb_a and bbb_b in place.

        ierr, rnorm = zlanbprow.zlanbprow(n_data_points, m, n, jold, j, 
                                          uuu, vvv, bbb_a, bbb_b,
                                          rnorm, options, 
                                          lambda_, trlambda)

        jold = j

        # %---------------------------------------------%
        # | Compute and analyze SVD(B) and error bounds |
        # %---------------------------------------------%
        
        bbb1_a = bbb_a.copy()
        bbb1_b = bbb_b.copy()

        work1[:j + 1] *= 0

        c1, c2 = dbdqr('N', j, bbb1_a, bbb1_b, workp, lanmax + 1)

        work1[j - 1] = c1
        work1[j] = c2

        xj      = j
        xbbb1_a = bbb1_a.copy()
        xbbb1_b = bbb1_b.copy()
        xwork1  = work1.copy()

        # THIS ONE CALLS ORIGINAL HLSVDPRO FORTRAN LIB TO ACCESS DBDSQR
        # r = dbdsqr.dbdsqr('u', j, 0, 1, 0, bbb1_a, bbb1_b, None, 1, work1, 1, None, 1)
        # info, bbb1_a, bbb1_b, _, work1, _ = r

        r = dbdsqr.dbdsqr_cython('u', xj, 0, 1, 0, bbb1_a, bbb1_b, None, 1, work1, 1, None, 1)
        info, bbb1_a, bbb1_b, _, work1, _ = r

        # PS - dbdsqr() returns bbb1_a, bbb1_b, and work1 as lists that are 
        # shorter than they were on the way in. Subsequent code (refinebounds() 
        # and probably other code too) doesn't deal with this size change,
        # so I have to pad them back out to their original size.
        
        bbb1_a +=  [0.0] * (lanmax - len(bbb1_a) + 1)
        bbb1_b +=  [0.0] * (lanmax - len(bbb1_b) + 1)

        bbb1_a = np.array(bbb1_a)
        bbb1_b = np.array(bbb1_b)
        work1 += [0.0] * (len_work1 - len(work1))
        work1 = np.array(work1)

        if j > 5:
            anorm = bbb1_a[0]
        else:
            anorm = max(anorm, bbb1_a[0])     

        work1 *= rnorm
        np.absolute(work1, work1)
        bbb1_b **= 2

        # %---------------------------------------------%
        # | Refine error bounds using the "Gap theorem" |
        # %---------------------------------------------%

        refinebounds(j, bbb1_b, work1, EPSN * anorm)

        # %----------------------------------------------------%
        # | Determine the number of converged singular values  |
        # %----------------------------------------------------%

        i = min(j, nsv_sought)
        bounds = work1[:i]

        sigma = [ ]
        i = 0

        while (i < min(j, nsv_sought)):

            if work1[i] <= (tol * bbb1_a[i]):
                sigma.append(bbb1_a[i])
                i += 1
            else:
                # Force loop exit
                i = nsv_sought

        neig = len(sigma)


        # %--------------------------------------------------%
        # | Test if an invariant subspace have been found or |
        # | the workspace has been exhausted.                |
        # %--------------------------------------------------%

        if ierr < 0:
            # Bail out of loop
            break

        if j >= lanmax:
            if neig < nsv_sought:
                # print "bailing"
                info = -1
            # Bail out of loop            
            break

        # %----------------------------------------------------%
        # | Increase the dimension of the Krylov subspace.     |
        # | If any Ritz values have converged then try to      | 
        # | estimate the average number of iterations per      |
        # | converged Ritz value.                              |
        # | Else increase the dimension with 50%.              |
        # %----------------------------------------------------%

        if neig > 1:
            dj = min(j // 2, ((nsv_sought - neig) * (j - 6)) // ( 2 * neig + 1))
            dj = min(100, max(2, dj))
        else:
            dj = j // 2
            dj = min(100, max(10, dj))

        j = min(j + dj, lanmax)



    ##############      while loop has ended      ##################

    # print "Pre-calculate singular vectors"

    # %-----------------------------------------%
    # | Calculate singular vectors if requested %
    # %-----------------------------------------%
    
    if neig >= nsv_sought:

        j = jold

        # In Python, I don't attempt to modify 'k' as the Fortran code does.
        # I just return neig directly.
        # k = neig

        if (JOBU == 'y') or (JOBV == 'y'):
            bbb1_a = bbb_a.copy()
            bbb1_b = bbb_b.copy()

            # %--------------------------------------------------------------
            # | The bidiagonal SVD is computed in a two-stage procedure:             
            # |                                                   
            # | 1. Compute a QR-factorization M^T*B = [R; 0] of the (k+1)-by-k  
            # |    lower bidiagonal matrix B.
            # | 2. Compute the SVD of the k-by-k upper bidiagonal matrix 
            # |    R = P*S*Q^T. The SVD of B is then (M*P)*S*Q^T.
            # %--------------------------------------------------------------

            c1, c2 = dbdqr(JOBU, j, bbb1_a, bbb1_b, workp, lanmax + 1)
            work1[j - 1] = c1
            work1[j] = c2

            workp = workp.reshape(lanmax + 1, -1)

            # transpose to account for Fortran row-major order
            workp = workp.transpose().flatten()

            ncc =  (j + 1 if (JOBU == 'y') else 0)
            ncvt = (j     if (JOBV == 'y') else 0)

            for i in range(n - 1):
                workq[i * (lanmax + 1)] = 1.0

            # THIS SHOULD POSSIBLY BE REPLACED BY A CALL TO THE FAST
            # DIVIDE-AND-CONQUER BIDIAGONAL SVD IN LAPACK 3 OR BY TRANSFORMING B
            # TO TRIDIAGONAL GOLUB-KAHAN FORM AND USING DHILLONS "HOLY GRAIL"
            # CODE.

            # %-----------------------------------------%
            # | R = P * S * Q^T, M^T <- P^T * M^T
            # %-----------------------------------------%

            # rbt = Really Big Tuple
            rbt = dbdsqr.dbdsqr('u', j, ncvt, 0, ncc, bbb1_a, bbb1_b, 
                                workq, lanmax, work1, 1, workp, lanmax + 1)

            info, bbb1_a, bbb1_b, _, work1, workp = rbt

            # print "after dbdsqr2"

            workp = np.array(workp)

            # %-----------------------------------------%
            # | Form left Ritz-vectors
            # | U = U * M * P
            # %-----------------------------------------%            
            #   call zgemmina(m,j+1,j+1,U,ldu, work(ip),lanmax+1,zwork(iwrk),lwrk) 
            
            if JOBU == 'y':
                workp = workp.reshape(lanmax + 1, -1)

                # The Fortran code describes zgemmina() as A=AxB', which in 
                # this case is uuu * workp'. However, I found that I do not
                # need to transpose workp to get the correct result. I suspect
                # this is one bug cancelling out another. I *should* 
                # transpose workp, but since workp was generated in the call
                # to dbdsqr2() above, Fortran has already transposed it.
                # If I replace dbdsqr2() with a pure Python version of it,
                # I'll need to transpose it here as expected.

                # I have three different ways of doing what zgemmina() does.
                # The first is the fastest and most straightforward and 
                # produces good results.
                # The second exactly matches what the Fortran code does but
                # it's probably slower than the first method and definitely 
                # more confusing.
                # The third just calls the Fortran code
                zgemmina_option = 1

                if zgemmina_option == 1:
                    # Use ordinary matrix multiplication.
                    uuu = np.array(np.matrix(uuu) * np.matrix(workp))
                elif zgemmina_option == 2:
                    # Use partial matrix multiplication. This matches what
                    # the Fortran code does.
                    uuu_rows = uuu.shape[1]
                    uuu_copy = np.matrix(uuu[:, :j + 1].copy())
                    workp_copy = workp.copy()
                    workp_copy = np.matrix(workp_copy[:j + 1, :j + 1])
                    uuu_copy = uuu_copy * workp_copy
                    # uuu_copy is a truncated version of uuu. Here I append the
                    # columns of uuu that are not present in uuu_copy.
                    delta = uuu_rows - uuu_copy.shape[1]
 
                    uuu = np.hstack( (uuu_copy, uuu[:, -delta:]) )
                    uuu = np.array(uuu)
                elif zgemmina_option == 3:
                    # Call the Fortran code.
                    uuu, _ = zgemmina.zgemmina(m, j + 1, j + 1, uuu, workp, lanmax + 1)

                # print "uuu[:1] --"
                # print uuu[:1]
                # print "uuu_copy[:1] --"
                # print uuu_copy[:1]
                # print "uuu[1:2] --"
                # print uuu[1:2]
                # print "uuu_copy[1:2] --"
                # print uuu_copy[1:2]

                # assert(np.allclose(uuu_copy, uuu[:m, :j + 1]))


            # FIXME - didn't port case where jobv ='Y' since it is always
            # hardcoded to 'n' in my case


    return uuu, sigma, info, neig


def dbdqr(jobq, n, d, e, qt, ldq):

    jobq = jobq.upper()

    qt = qt.reshape(ldq, -1)

    # Compute QR factorization B = Q*R of (n+1) x n lower bidiagonal matrix 
    # with diagonal elements d(1)...d(n) and first subdiagonal elements
    # e(1)...e(n). On return [0 ... 0 c1 c2]' = Q'*[0 ... 0 1]'.
    # 
    # If jobq=='Y' then on return Qt contains Q^T.

    # PS - tested for when jobq = 'N', not tested for jobq = 'Y'


    msg = "{} qt[{}, {}] = {}"
    pqt = lambda letter, i1, i2, value: msg.format(letter, i1 + 1, i2 + 1, value)

    if n >= 1:
        if jobq == 'Y':
            qt[:n + 1, : n + 1] *= 0.0

            for j in range(n + 1):
                qt[j, j] = 1.0

        for i in range(n - 1):

            #cs, sn, r = dlartg.dlartg(d[i], e[i])
            cs, sn, r = dlartg(d[i], e[i])

            d[i] = r
            e[i] = sn * d.item(i + 1)
            d[i + 1] = cs * d.item(i + 1)

            if jobq == 'Y':

                # PS - This is the naive Python equivalent of the Fortran --
                # for j in range(i + 1):
                #     qt[i + 1, j] = -sn * qt.item(i, j)
                #     qt[i    , j] =  cs * qt.item(i, j)

                # PS - Here's the fast version --
                qt[i + 1, :i + 1] = -sn * qt[i, :i + 1]
                qt[i    , :i + 1] =  cs * qt[i, :i + 1]


                qt[i    , i + 1] = sn
                qt[i + 1, i + 1] = cs

        # In fortran, a loop ends when the counter is one past the limit. For
        # example, at the end of the loop  `do i = 1,50`, i == 51. In Python,
        # this is not the case so I have to increment i here.
        i += 1

        #cs, sn, r = dlartg.dlartg(d[n - 1], e[n - 1])
        cs, sn, r = dlartg(d[n - 1], e[n - 1])


        d[n - 1] = r
        e[n - 1] = 0.0
        c1 = sn
        c2 = cs

        if jobq == 'Y':
            # PS - This is the naive Python equivalent of the Fortran --
            # for j in range(i + 1):
            #     qt[i + 1, j] = -sn * qt.item(i, j)
            #     qt[i    , j] =  cs * qt.item(i, j)
            # PS - Here's the fast version --

            qt[i + 1, :i + 1] = -sn * qt[i, :i + 1]
            qt[i    , :i + 1] =  cs * qt[i, :i + 1]

            qt[i    , i + 1] = sn
            qt[i + 1, i + 1] = cs

    return c1, c2


def refinebounds(n, theta, bound, tol):
    # 
    # Refine Lanczos error bounds using the gap theorem.
    #      

    if n > 1:
        for i in range(n):
            for j in (-1, 1):
                if ((j == 1) and (i < n)) or ((j == -1) and (i > 1)):

                    if abs(theta[i] - theta[i + j]) < (EPS34 * theta[i]):
                        if (bound[i] > tol) and (bound[i + j] > tol):
#                            bound[i + j] = util.dlapy2(bound[i], bound[i + j])

                            bound[i + j] = math.sqrt(bound[i]**2 + bound[i + j]**2)
                            
                            bound[i] = 0.0

            gap = theta[i + 1] - bound[i + 1] - (theta[i] + bound[i])

            if gap > bound[i]:
                # print "dong!"
                bound[i] *= (bound[i] / gap)



