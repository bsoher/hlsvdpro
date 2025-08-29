#
#     (C) Brian J Soher, 2020
#

from __future__ import division, print_function

import math
import numpy as np

from zgetu0 import zgetu0
from zreorth import zreorth
from zsafescal import zsafescal
from zblasext import pdznrm2, pzdotc, pzdaxpy, pzdscal, pzaxpy


MAX_SINGULAR_VALUES = 50
FUDGE = 1.01
KAPPA = 0.717

# Machine- and compiler-dependent constants
# On my Mac, EPS in Python is 2.22044604925e-16. In Fortran
# it is half as big (1.11022302462515654E-016).
EPS   = np.finfo(np.float64).eps
EPS34 = EPS ** (3.0 / 4.0)
EPS1  = 100 * EPS


# I create a shortcut for sqrt() since I use it a lot.
sqrt = math.sqrt


def fortran_range(a, b):
    the_range = list(range(a, b))
    if b >= a:
        the_range.append(b)
    return the_range


def zlanbpro( m, n, k0, k, aprod, U, ldu, V, ldv, bbb_a, bbb_b, rnorm, options, parm):
    """
    This docstring is incomplete.

        # FIXME complete the docstring

        bbb_a and bbb_b are modified in-place.

        This function accepts a dict of options that exert some control over how
        it behaves. The dict keys and their meaning are as follows --
           - classical_gs: True ==> classical Gram-Schmidt reorth (see zreorth)
           - extended_reorth_iterations: # of iterations in one case of
                 extended local reorth
           - delta: used if >= 0, otherwise a reasonable default is calculated.
           - eta:   used if >= 0, otherwise a reasonable default is calculated.
           - anorm: used if >  0, otherwise a reasonable default is calculated.

        None of these options are changed by this function EXCEPT FOR "anorm"
        which is rewritten every time (presumably so it can be passed back in
        on the next iteration).




    DLANBPRO: Computes K steps of the Lanczos bidiagonalization (LBD) 
    algorithm with partial reorthogonalization (BPRO) with M-by-1 starting 
    vector U(:,k0+1), producing a lower bidiagonal K+1-by-K matrix B_k, an 
    N-by-K matrix V_k, an M-by-K+1 matrix U_{k+1} such that
          A*V_k = U_{k+1}*B_k
    Partial reorthogonalization is used to keep the columns of V_K and U_k
    semiorthogonal to a level prescribed in DOPTION(1), i.e.
          MAX(DIAG((EYE(K) - V_K'*V_K))) <= DOPTION(1)
    and
          MAX(DIAG((EYE(K) - U_K'*U_K))) <= DOPTION(1).

    If K0>0 and K>K0 an existing K0-step LBD of stored in U, V and B is 
    extended to a K-step LBD.

    Parameters:

    M: INTEGER. Number of rows of A.
    N: INTEGER. Number of columns of A.
    K0: INTEGER. The dimension of the previously computed Lanczos
                 bidiagonalization stored in U, V, and B.
    K: INTEGER. On entry: The desired dimension of the Lanczos 
          bidiagonalization. On exit: the actual size of the LBD computed.
          This can be smaller than the input value if an invariant subspace
          is computed.
    APROD: Subroutine defining the linear operator A. 
           APROD should be of the form:

          SUBROUTINE DAPROD(TRANSA,M,N,X,Y,DPARM,IPARM)
          CHARACTER*1 TRANSA
          INTEGER M,N,IPARM(*)
          DOUBLE PRECISION X(*),Y(*),DPARM(*)

          If TRANSA.EQ.'N' then the function should compute the matrix-vector
          product Y = A * X.
          If TRANSA.EQ.'C' then the function should compute the matrix-vector
          product Y = A^H * X.
          The arrays IPARM and DPARM are a means to pass user supplied
          data to APROD without the use of common blocks.
    U(LDU,K+1): DOUBLE PRECISION array. On return the first K+1 columns of
              U will contain the left Lanczos vectors.
              On entry: 
                 If K0==0 the first column of U contains the starting
                 vector for the Lanczos bidiagonalization. A random 
                 starting vector is used if the first column is U is zero.
                 If K0>0 the first K0+1 columns of U are assumed to 
                 contain the first K0+1 left Lanczos vectors of an
                 existing LBD.

    LDU: INTEGER. Leading dimension of the array U. LDU >= M.
    V(LDV,K): DOUBLE PRECISION array. On return the first K columns of
              V will contain the right Lanczos vectors.
              On entry: 
                 If K0>0 the first K0 columns of V are assumed to 
                 contain the first K0 right Lanczos vectors of an
                 existing LBD.
    LDV: INTEGER. Leading dimension of the array V. LDV >= N.
    B(K,2): DOUBLE PRECISION array. On return the first columns of
              B will contain the K diagonal elements of B_k, and
              the second column of B will contain the K elements
              of the first sub-diagonal of B_k. 
    LDB: INTEGER. Leading dimension of the array B. LDB >= K.
    RNORM: DOUBLE PRECISION. On entry RNORM must contain the norm of
              the K0+1st column of U.
              On exit RNORM contains the value of the (K+1,K) element
              of B_k.
    DOPTION: DOUBLE PRECISION array. 
       doption(1 -> 0) = delta. Level of orthogonality to maintain among
         Lanczos vectors.
       doption(2 -> 1) = eta. During reorthogonalization, all vectors with
         with components larger than eta along the latest Lanczos vector
         will be purged.
       doption(3 -> 2) = anorm. Estimate of || A ||.
    IOPTION: INTEGER array. 
       ioption(1 -> 0) = CGS.  If CGS.EQ.1 then reorthogonalization is done
         using iterated classical GRAM-SCHMIDT. IF CGS.EQ.0 then 
         reorthogonalization is done using iterated modified Gram-Schmidt.
       ioption(2 -> 1) = ELR. If ELR.EQ.1 then extended local orthogonality is
         enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
    DWORK(LWORK): DOUBLE PRECISION array of dimension >= m+n+2*k+2+max(m,n).
    ZWORK(max(M,N)): DOUBLE COMPLEX array of dimension max(M,N): Workspace.
    IWORK(2*K+1): INTEGER ARRAY. Integer workspace.
    DPARM: DOUBLE PRECISION array. Array used for passing data to the APROD
        function.   
    IPARM: INTEGER array. Array used for passing data to the APROD
        function.   
    IERR: INTEGER. Error status code.
        IERR < 0  : An invariant subspace of dimension -J was found.
        IERR == 0 : The computation succeeded.
        IERR > 0  : The computation succeeded, but the algorithm 
                    came close to computing an invariant subspace after 
                    IERR steps. In a previous version this would have caused 
                    the algorithm to switch to full reorthogonalization 
                    after IERR steps, but that is no longer the case. 
                    It is probably safe to ignore.
    
    """

    ndp, _ = aprod.shape

    # %---------------------------------%
    # | Set machine dependent constants |
    # %---------------------------------%
    
    #     eps = dlamch('e')
    #     eps34 = eps**(3d0/4d0)
    #     epsn = dble(max(m,n))*eps
    #     epsn2 = sqrt(dble(max(m,n)))*eps
    #     sfmin = dlamch('s')

    EPSN = max(m, n) * EPS
    EPSN2 = sqrt(max(m, n)) * EPS

    ierr = 0

    # %------------------------%
    # | Set default parameters |
    # %------------------------%

    delta = sqrt(EPS / k) if (options["delta"] < 0) else options["delta"]
    
    eta = EPS34 / sqrt(k) if (options["eta"] < 0) else options["eta"]

    full_reorth = ((delta <= eta) or (delta == 0.0))

    if options["anorm"] > 0:
        anorm = options["anorm"]
    elif k0 > 0:
        #anorm = dlapy2.dlapy2(bbb_a[0], bbb_b[0])
        anorm = sqrt(bbb_a[0]*bbb_a[0]+bbb_b[0]*bbb_b[0])
        if anorm <= 0:
            ierr = -1
            options['anorm'] = anorm
            # FIXME raise an error instead?
            return ierr, 0.0, k
    else:
        anorm = 0.0


    # %---------------------%
    # | Get starting vector |
    # %---------------------%
    if not rnorm:
        rnorm, anormest, ierr = zgetu0('n', m, n, k0, 3, U[:,k0], U, ldu, aprod, parm, options['cgs'])
        anorm = max(anorm,anormest)


    # # %------------------------------%
    # # | Set pointers into work array |
    # # %------------------------------%
    # imu  = 1
    # inu  = imu+k[0]+1     # FIXME bjs, check off by one
    # iidx = 1
    # isx  = 1

    zwork   = np.zeros((m,), np.complex128)         # replaces zwork[iu] a.k.a. zwork
    zworkv  = np.zeros((n,), np.complex128)         # replaces zwork[iv]
    zworks  = np.zeros((3 * ndp,), np.complex128)   # replaces zwork[is]
    zworkz1 = np.zeros((ndp,), np.complex128)       # replaces zwork[iz1]
    # I'm not 100% sure about the size of zworky1 since it is the last pointer
    # in the Fortran zwork array. Since it doesn't have a next neighbor in
    # Fortran, I can't tell where iy1 stops. Based on how it's used, I think
    # a length of ndp is OK.
    zworky1 = np.zeros((ndp,), np.complex128)  # replaces zwork[iy1]

    # mus and nus track the mu and nu values. They replace the work array in
    # Fortran which appears to be much bigger than it needs to be.
    mus = np.zeros((k + 1), np.float64)
    nus = np.zeros((k + 1), np.float64)

    # indices replaces iwork from the Fortran code. It might not need to be this
    # big.
    indices = np.zeros(((2 * MAX_SINGULAR_VALUES) + 1,), np.int)

    indices -= 1




    # %---------------------------%
    # | Prepare Lanczos iteration |
    # %---------------------------%

    if k0 == 0:
        amax  = 0.0
        alpha = 0.0
        beta  = rnorm
        force_reorth = False 

        # %---------------------------------------------------%
        # | Compute ||A x|| / ||x|| for a random vector x     |
        # | to make it less likely that ||A|| is grossly      |
        # | underestimated at the beginning of the iteration. |
        # %---------------------------------------------------%
        if n > m:            # -1 index for j here is flag to not do certain code
            u0norm, anormest, ierr = zgetu0('n',m, n, -1, 1, zworks, U, ldu, aprod, parm, options['cgs'])
        else:
            u0norm, anormest, ierr = zgetu0('c', m, n, -1, 1, zworks, V, ldv, aprod, parm, options['cgs'])

        anorm = max(anorm, FUDGE*anormest)

        j0 = 0
        if beta != 0.0:
            zsafescal(m, beta, U[:,0])

        mus[0] = 1        # bjs - these are indices?
        nus[0] = 1

    else:
        force_reorth = True 
        alpha = bbb_a[k0 - 1]
        beta  = rnorm

        if (k0 < k) and (beta * delta < anorm * EPS):
            full_reorth = True
            ierr = k0

        indices[0] = 0
        indices[1] = k0 - 1
        indices[2] = k0

        pzdscal(m, rnorm, U[:,k0], 1)

        rnorm = zreorth(m, k0, U, ldu, U[:,k0], rnorm, indices, KAPPA, options['cgs'])

        zsafescal(m, rnorm, U[:,k0])

        dset_mu(k0, mus, indices, EPSN2)        # dset_mu(k0,dwork[imu:],iwork[iidx:],epsn2)
        dset_mu(k0, nus, indices, EPSN2)        # dset_mu(k0,dwork[inu:],iwork[iidx:],epsn2)

        beta = rnorm        # duplicate to line above, because rnorm not changed between?

        # %--------------------------------------%
        # | Estimate ||B||_2^2 as ||B^T * B||_1  |   
        # %--------------------------------------%

        bbb_b[k0 - 1] = beta
        amax = 0.0

        for j in range(k0):
            amax = max(amax, bbb_a[j], bbb_b[j])  # FIXME bjs, check off by one
            if j == 0:
                anorm = max(anorm, FUDGE * alpha)
            elif j == 1:
                a1 = bbb_b[0] / amax
                a1 = FUDGE * amax * sqrt((bbb_a[0] / amax)**2 + a1**2 + \
                                          bbb_a[1] / amax * a1)
                anorm = max(anorm, a1)
            else:
                a1 = bbb_a[j - 1] / amax        # FIXME bjs, check off by one
                b1 = bbb_b[j - 1] / amax        # FIXME bjs, check off by one
                a1 = FUDGE * amax * sqrt(a1**2 + b1**2 + a1*bbb_b[j-2] / amax + \
                                         bbb_a[j] / amax * b1
                                         )
                anorm = max(anorm, a1)

        j0 = k0

    numax = 0.0
    mumax = 0.0



    # %-------------------------------------------%
    # | Start Lanczos bidiagonalization iteration |
    # %-------------------------------------------%            

    for j in range(j0,k):       # FIXME bjs, check for off by one error

        # %---------------------------------------------%
        # | alpha_{j} v_{j} = A'*u_{j} - beta_{j} v_{j} |
        # %---------------------------------------------%

        aprod('c', m, n, U[:, j], V[:, j], parm)

        if j == 0:
            alpha = pdznrm2(n, V[:, j], 1)
            anorm = max(anorm, FUDGE * alpha)
        else:
            pzdaxpy(n, -beta, V[:, j-1], 1, V[:, j], 1)
            alpha = pdznrm2(n, V[:, j], 1)

#            # Fast Python equivalent --
#            zworkv[:n] = zworks[:n] - (beta * zworkv[:n])


            # %------------------------------------%
            # | Extended local reorthogonalization |
            # %------------------------------------%

            if (j > 0) and (options["exri"] > 0) and (alpha < (KAPPA * beta)):

                for i in range(1):
                    #  do i=1,ioption(2) which is 0 or 1 so "do i=1,1"

                    s = pzdotc(n,V[:,j-1],1,V[:,j],1)

                    pzaxpy(n,-s,V[:,j-1],1,V[:,j],1)

                    ## This is a fast version --        # from v1, not that much time to save
                    # zworkv[:n] += (-s * V[:n, j-1])   # vis cProfile

                    nrm = pdznrm2(n,V[:,j],1)

                    if nrm >= (KAPPA * alpha):
                        # Bail out of loop
                        break
                    alpha = nrm

                nus[j-1] = EPS
                alpha = nrm

            bbb_a[j] = alpha
            amax = max(amax, alpha)


            # %----------------------------%
            # | Update estimate of ||A||_2 |         
            # %----------------------------%

            if j == 1:
                a1 = bbb_b[0] / amax
                a1 = FUDGE * amax * sqrt((bbb_a[0] / amax)**2 + (a1**2) + \
                                         bbb_a[1] / amax * a1 )
            else:
                a1 = bbb_a[j - 1] / amax
                b1 = bbb_b[j - 1] / amax
                a1 = FUDGE * amax * sqrt(a1**2 + b1**2 + \
                                         a1 * bbb_b[j - 2] / amax + \
                                         bbb_a[j] / amax * b1 )
            anorm = max(anorm, a1)

        # %--------------------------%
        # | Update the nu recurrence |
        # %--------------------------%

        if (not full_reorth) and (alpha != 0.0):
            numax = dupdate_nu(mus, nus, j, bbb_a, bbb_b, anorm, EPSN2)

        # %------------------------------%
        # | Reorthogonalize if necessary |
        # %------------------------------%

        if (full_reorth or force_reorth or (numax > delta)) and (alpha != 0.0):

            if full_reorth or (eta == 0.0):
                indices[0] = 0
                indices[1] = j-1
                indices[2] = j
                indices[3] = -1     # not in v2.1 not in v1 fortran?  PS maybe
                indices[4] = -1     # not in v2.1 not in v1 fortran?  PS maybe
            elif not force_reorth:
                indices = dcompute_int(nus, j-1, delta, eta)

            alpha = zreorth(n, j-1, V, ldv, V[:,j], alpha, indices, KAPPA, options['cgs'])

            dset_mu(j-1, nus, indices, EPS)
            numax = eta
            force_reorth = not force_reorth


        # %-----------------------------------------------%
        # | Check whether an invariant subspace was found |
        # %-----------------------------------------------% 

        if (alpha < (anorm * EPSN)) and (j < k-1):  # FIXME BJS, check k-1?

            rnorm = alpha
            alpha = 0.0

            # %------------------------------------------------%
            # | Try to build an orthogonal subspace, starting  |
            # | with a random vector.                          |
            # %------------------------------------------------%

            alpha, anormest, ierr = zgetu0('c', m, n, j-1, 3, V[:,j], V, ldv, aprod, parm, options["cgs"])

            if alpha == 0.0:
                # %------------------------------------------------%
                # | We failed to generate a new random vector      |
                # | in span(A^T) orthogonal to span(V(:,1:j-1)).   |
                # | Most likely span(V(:,1:j-1)) is an invariant   |
                # | subspace.                                      |
                # %------------------------------------------------%

                k = j-1
                ierr = -j
                return ierr, rnorm, k

            else:
                # %-------------------------------------------------%
                # | We have managed to generate a random vector     |
                # | in span(A^T) orthogonal to V(:,1:j-1), so we    |
                # | can continue the LBD and "deflate" the subspace |
                # | by setting alpha_{j} = 0.                       |
                # %-------------------------------------------------%

                zsafescal(n, alpha, V[:,j])
                alpha = 0.0
                force_reorth = True
                if delta > 0.0:
                    full_reorth = False

        elif (j > 1) and (not full_reorth) and (j < k-1) and (delta * alpha < anorm * EPS):
             ierr = j

        bbb_a[j] = alpha

        if alpha != 0.0:
            zsafescal(n, alpha, V[:,j])

        # %------------------------------------------------%
        # | beta_{j+1} u_{j+1} = A*v_{j} - alpha_{j} u_{j} |
        # %------------------------------------------------%

        aprod('n',m,n,V[:,j],U[:,j+1],parm)

        pzdaxpy(m,-alpha,U[:,j],1,U[:,j+1],1)

        ## PS - This is the fast version of pzdaxpy()--
        ##
        #U[:m,j+1] = U[:m,j+1] - (alpha * U[:m,j])

        beta = pdznrm2(m, U[:, j+1], 1)


        # %------------------------------------%
        # | Extended local reorthogonalization |
        # %------------------------------------%

        if (options["exri"] > 0 and (beta < KAPPA*alpha)):
            for i in range(options["exri"]):

                s = pzdotc(m, U[:,j], 1, U[:,j+1], 1)

                pzaxpy(m, -s, U[:,j], 1, U[:,j+1], 1)

                nrm = pdznrm2(m, U[:,j+1], 1)

                if nrm >= KAPPA*beta:
                    break
                beta = nrm

            mus[j] = EPS        # line 774 .f
            beta = nrm

        bbb_b[j] = beta
        amax = max(amax, beta)

        # %----------------------------%
        # | Update estimate of ||A||_2 |         
        # %----------------------------%

        if j == 0:
            a1 = sqrt(bbb_a[0]*bbb_a[0]+bbb_b[0]*bbb_b[0])
        else:
            a1 = bbb_a[j] / amax
            a1 = amax * sqrt(a1**2 + \
                             (bbb_b[j] / amax)**2 +
                             a1 * bbb_b[j - 1] / amax )
        anorm = max(anorm, a1)


        # %--------------------------%
        # | Update the mu recurrence |
        # %--------------------------%
        if (not full_reorth) and (beta != 0.0):
            mumax = dupdate_mu(mus, nus, j, bbb_a, bbb_b, anorm, EPSN2)


        # %--------------------------------------%
        # | Reorthogonalize u_{j+1} if necessary |
        # %--------------------------------------%

        if (beta != 0.0) and (full_reorth or force_reorth or (mumax > delta)):
            if full_reorth or (eta == 0.0):
                indices[0] = 0
                indices[1] = j
                indices[2] = j+1
                indices[3] = -1
                indices[4] = -1
            elif not force_reorth:
                indices = dcompute_int(mus, j, delta, eta)
            else:
                for i in range(2*j+1):
                    if indices[i] == j:
                        indices[i] = j + 1
                        break

            beta = zreorth(m, j, U, ldu, U[:,j+1], beta, indices, KAPPA, options['cgs'])

            dset_mu(j, mus, indices, EPS)
            mumax = eta
            force_reorth = not force_reorth


        # %-----------------------------------------------%
        # | Check whether an invariant subspace was found |
        # %-----------------------------------------------%

        if (beta < (anorm * EPSN)) and (j < k-1):       # FIXME bjs, check if k-1? or k
            rnorm = beta
            beta = 0.0

            # %-----------------------------------------------%
            # | Try to build an orthogonal subspace, starting |
            # | with a random vector.                         |
            # %-----------------------------------------------%
            beta, anormest, ierr = zgetu0('n', m, n, j, 3, U[:,j+1], U, ldu, aprod, parm, options["cgs"])

            if beta == 0.0:
                # %-----------------------------------------------%
                # | We failed to generate a new random vector     |
                # | in span(A) orthogonal to span(U(:,1:j)).      |
                # | Most likely span(U(:,1:j)) is an invariant    |
                # | subspace.                                     |
                # %-----------------------------------------------%
                k = j
                ierr = -j
                return ierr, rnorm, k
            else:
                # %------------------------------------------------%
                # | We have managed to generate a random vector    |
                # | in span(A) orthogonal to U(:,1:j), so we can   |
                # | continue the LBD and "deflate" the subspace by |
                # | setting beta_{j+1} = 0.                        |
                # %------------------------------------------------%
                zsafescal(n, beta, U[:,j+1])
                beta = 0.0
                force_reorth = True
                if delta > 0:
                    full_reorth = False

        elif (not full_reorth) and (j < k-1) and ((delta * beta) < (anorm * EPS)):
            ierr = j

        bbb_b[j] = beta

        if (beta != 0) and (beta != 1):
            zsafescal(m, beta, U[:,j+1])

        rnorm = beta

    return ierr, rnorm, k






#------------------------------------------------------------------------------

def dset_mu(k, mus, indices, val):

    i = 0
    while (indices[i] <= k) and (indices[i] >= 0):
        # p = indices[i]
        # q = indices[i + 1]
        # # FIXME check for off-by-one error
        # for j in fortran_range(p, q):   # bjs fastest
        #     mus[j] = val
        for j in fortran_range(indices[i], indices[i+1]):   # bjs fastest
            mus[j] = val
        i += 2



#------------------------------------------------------------------------------

def dcompute_int(mus, j, delta, eta):
    """
    I think "int" is short for "intervals". returns indices.

    """
    indices = [ ]

    ip = 0
    i = -1

    if delta >= eta:
        while i < j:
            # find the next mu(k), k>i where abs(mu(k)) > delta.
            # If there is no such mu(k), set k = i + 1.
            k = i + 1

            while (k < j) and (abs(mus[k]) <= delta):
                k += 1

            # FIXME the 'if' block below seemed like a good idea but
            # doesn't actually work. Que pasa???
            # if k >= j:
            #     # If the loop above completes without finding
            #     # (abs(mus[k]) <= delta), it should skip the rest. Fortran
            #     # implements this with a goto.
            #     break

            # find smallest s<k such that for all j=s,..,k, m(j) >= eta
            # the_range = [value for value in reversed(range(max(i, 0), k))]
            # the_range will be something like [3, 2, 1, 0] and .pop() will
            # remove entries from the right-hand side.
            s = k
            done = False
            while not done:
                s -= 1
                if s < max(i, 0):
                    # Oooops, we're off the end of our search range.
                    done = True
                else:
                    if abs(mus[s]) < eta:
                        done = True
            # We actually want to record the index one to the right of
            # the position we found.
            s += 1
            indices.append(s)

            # Now search to the right
            found = False
            for i in range(s, j):
                if abs(mus[i]) < eta:
                    found = True
                    break

            # When a loop iterator is exhausted in fortran, the iterator
            # variable is set to the end of the range + 1. This code mimics
            # that behavior.
            if not found:
                i += 1

            indices.append(i)

        indices.append(j + 1)
        # Append end-of-list flags
        indices.append(-1)
        indices.append(-1)

    return indices



#------------------------------------------------------------------------------

def dupdate_mu(mus, nus, j, alpha, beta, anorm, epsi1):
    """
    Alters mus, returns mumax. I'm not sure if this is truly the max value
    in abs(mus) or just the max of the subset that this function touches.

    """
    if j == 0:
        # d = epsi1 * (sqrt(alpha[j]*alpha[j]+beta[j]*beta[j]) + alpha[0]) + \
        #     epsi1 * anorm
        mus[0] = epsi1 / beta[0]
        mumax = abs(mus[0])

    else:
        const1 = sqrt(alpha[j]*alpha[j]+beta[j]*beta[j])

        mus[0] = alpha[0] * nus[0] - alpha[j] * mus[0]
        d = epsi1 * (const1 + alpha[0]) + epsi1 * anorm
        mus[0] = (mus[0] + math.copysign(d, mus[0])) / beta[j]
        mumax  = abs(mus[0])

        the_range = fortran_range(1, j-1)
        for k in the_range:
            mus[k] = (alpha[k] * nus[k]) + (beta[k-1] * nus[k-1]) - (alpha[j] * mus[k])
            d = epsi1 * ((const1 + sqrt(alpha[k]*alpha[k]+beta[k-1]*beta[k-1])) + anorm)
            mus[k] = (mus[k] + math.copysign(d, mus[k])) / beta[j]
            mumax = max(mumax, abs(mus[k]))

        mus[j] = beta[j - 1] * nus[j - 1]
        d = epsi1 * ((const1 + sqrt(alpha[j]*alpha[j] + beta[j-1]*beta[j-1])) + anorm)
        mus[j] = (mus[j] + math.copysign(d, mus[j])) / beta[j]
        mumax = max(mumax, abs(mus[j]))

    mus[j + 1] = 1.0

    return mumax


#**********************************************************************

def dupdate_nu(mus, nus, j, alpha, beta, anorm, epsi1):
    """
    Alters nus, returns numax

    """
    numax = 0.0

    if j > 0:

        const1 = sqrt(alpha[j]*alpha[j]+beta[j-1]*beta[j-1])

        for k in fortran_range(0, j - 1):
            nus[k] = (beta[k] * mus[k+1]) + (alpha[k]*mus[k]) - (beta[j-1] * nus[k])
            d = epsi1 * ((sqrt(alpha[k]*alpha[k]+beta[k]*beta[k]) + const1) + epsi1)
            nus[k] = (nus[k] + math.copysign(d, nus[k])) / alpha[j]
            numax = max(numax, abs(nus[k]))

        # for k in fortran_range(0, j - 1):
        #     nus[k] = (beta[k] * mus[k + 1]) + (alpha[k] * mus[k]) - \
        #              (beta[j - 1] * nus[k])
        #     d = epsi1 * (sqrt(alpha[k]*alpha[k]+beta[k]*beta[k]) + \
        #                 sqrt(alpha[j]*alpha[j]+beta[j-1]*beta[j-1])
        #                )
        #     d += epsi1 * anorm
        #     nus[k] = (nus[k] + math.copysign(d, nus[k])) / alpha[j]
        #     numax = max(numax, abs(nus[k]))

        nus[j] = 1.0

    return numax



