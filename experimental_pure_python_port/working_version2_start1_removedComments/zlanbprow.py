# Python modules
from __future__ import division
import ctypes
import math
import pdb

# 3rd party modules
import numpy as np
import numpy.linalg
import scipy
import scipy.linalg.blas
import scipy.linalg.lapack

# Our modules
import vespa.common.util.misc as util_misc
import zget0w
import util
import aprodw
import zreorth
import zsafescal


MAX_SINGULAR_VALUES = 50
FUDGE = 1.01
KAPPA = 0.717

# Machine- and compiler-dependent constants
# On my Mac, EPS in Python is 2.22044604925e-16. In Fortran
# is is half as big (1.11022302462515654E-016).
EPS = np.finfo(np.float64).eps
EPS34 = EPS ** (3.0 / 4.0)
EPS1 = 100 * EPS


# In the Fortran, the constant MGS (==> modified Gram-Schmidt) is hardcoded 
# to 1. It's not present in the original PROPACK code. Since it's always 1,
# I've opted to ignore it.

# I create a shortcut for sqrt() since I use it a lot.
sqrt = math.sqrt

# http://www.scipy.org/doc/api_docs/SciPy.lib.lapack.info.html

# Float BLAS functions
dznrm2, = scipy.linalg.blas.get_blas_funcs( ['znrm2'], np.array([0.0]))

# Complex BLAS functions
zdotc, zaxpy = scipy.linalg.blas.get_blas_funcs( ['dotc', 'axpy'], np.array([0j]))


# Float LAPACK functions
dlamch, = scipy.linalg.lapack.get_lapack_funcs( ['lamch'], np.array([np.float64(0.0)]))


def fortran_range(a, b):
    the_range = range(a, b)
    if b >= a:
        the_range.append(b)
    return the_range


# cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
# c........1.........2.........3.........4.........5.........6.........7..
# c23456789012345678901234567890123456789012345678901234567890123456789012
# c
#       subroutine zlanbprow_python(ndp, m, n, k0, k, ierr,
#      c                            uuu_r, uuu_i, ldu,
#      c                            vvv_r, vvv_i, ldv,
#      c                            bbb, ldb,
#      c                            rnorm, doption, ioption, 
#      c                            lambda_r, lambda_i,
#      c                            trlambda_r, trlambda_i)


def zlanbprow(ndp, m, n, k0, k, uuu, vvv, bbb_a, bbb_b, rnorm,
              options, lambda_, trlambda):
    # FIXME complete the docstring
    """This docstring is incomplete.

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

    """

    # c     %---------------------------------%
    # c     | Set machine dependent constants |
    # c     %---------------------------------%
    # PS - Fortran populates variable sfmin but otherwise doesn't use it.
    
    EPSN = max(m, n) * EPS
    EPSN2 = sqrt(max(m, n)) * EPS

    ierr = 0

    # c     %------------------------%
    # c     | Set default parameters |
    # c     %------------------------%

    delta = sqrt(EPS / k) if (options["delta"] < 0) else options["delta"]

    eta = EPS34 / sqrt(k) if (options["eta"] < 0) else options["eta"]

    full_reorth = ((delta <= eta) or (delta == 0.0))

    if options["anorm"] > 0:
        anorm = options["anorm"]
    elif k0 > 0:
#        anorm = util.dlapy2(bbb_a[0], bbb_b[0])
        anorm = sqrt(bbb_a[0]**2 + bbb_b[0]**2)
        
        if anorm <= 0:
            ierr = -1
            # FIXME raise an error instead?
            return
    else:
        anorm = 0.0

    # c     %---------------------%
    # c     | Get starting vector |
    # c     %---------------------%

    if not rnorm:
        rnorm, anorm, ierr = zget0w.zgetu0w('n', ndp, m, n, k0, 3, uuu[:,k0:],
                                               uuu, options["classical_gs"],
                                               lambda_, trlambda)

    # c %------------------------------%
    # c | Set pointers into work array |
    # c %------------------------------%

    zwork  = np.zeros( (m, ), np.complex128)  # replaces zwork[iu] a.k.a. zwork
    zworkv = np.zeros( (n, ), np.complex128)  # replaces zwork[iv]
    zworks = np.zeros( (3 * ndp, ), np.complex128)  # replaces zwork[is]
    zworkz1 = np.zeros( (ndp, ), np.complex128)  # replaces zwork[iz1]

    # I'm not 100% sure about the size of zworky1 since it is the last pointer
    # in the Fortran zwork array. Since it doesn't have a next neighbor in 
    # Fortran, I can't tell where iy1 stops. Based on how it's used, I think 
    # a length of ndp is OK.
    
    zworky1 = np.zeros( (ndp, ), np.complex128)  # replaces zwork[iy1]

    # mus and nus track the mu and nu values. They replace the work array in
    # Fortran which appears to be much bigger than it needs to be.
    
    mus = np.zeros( (k + 1), np.float64)
    nus = np.zeros( (k + 1), np.float64)

    # indices replaces iwork from the Fortran code. It might not need to be this big.
    
    indices = np.zeros( ((2 * MAX_SINGULAR_VALUES) + 1, ), np.int)
    indices -= 1

    # The Fortran code does something fishy in the first calls to zgetu0w()
    # below. It passes a complex
    # variable 's' in the slot that zgetu0w() calls u0norm. u0norm
    # is declared as a real in zgetu0w() despite being passed a complex
    # by zlanbprow. Whether or not this is valid Fortran I don't know.
    # In any case, the practical result is that the real part of s
    # is altered while the imaginary part remains unchanged. 
    # In Fortran, s is initially populated with stack garbage. I "declare"
    # it here to simulate Fortran's behavior.
    
    s = 0j

    # c %---------------------------%
    # c | Prepare Lanczos iteration |
    # c %---------------------------%

    if not k0:
        # print "zlanbprow: k0 is zero"
        amax = 0.0
        alpha = 0.0
        beta = rnorm
        force_reorth = False

        # c     %-----------------------------------------------%
        # c     | Compute ||A x|| / ||x|| for a random vector x |
        # c     | to make it less likely that ||A|| is grossly  |
        # c     | underestimated.                               |
        # c     %-----------------------------------------------%

        if n > m:
            # print "zlanbprow: n > m"
            s_real, anormest, ierr = zget0w.zgetu0w('n', ndp, m, n, 0, 1, zwork,
                                               uuu, options["classical_gs"],
                                               lambda_, trlambda)
            s.real = s_real
        else:
            # print "zlanbprow: n <= m"
            s_real, anormest, ierr = zget0w.zgetu0w('t', ndp, m, n, 0, 1, 
                                                zworkv,
                                                vvv, options["classical_gs"],
                                                lambda_, trlambda)
            s = complex(s_real, s.imag)

        anorm = max(anorm, FUDGE * anormest)
        j0 = 0
        if beta:
            zsafescal.zsafescal(m, beta, uuu[:, 0])

        zwork = uuu[:, 0].copy()

        mus[0] = 1.0
        nus[0] = 1.0

        # I suspect that the two lines below from the Fortran code are a 
        # harmless bug. In every other case where they're used, imu and inu 
        # are indices into work, not zwork. I'm just going to ignore these
        # two lines.
        #       zwork(imu) = cone
        #       zwork(inu) = cone
        
    else:
        # k0 != 0

        # for ips, value in enumerate(bbb_a):
        #     print "bbb_a[{}] = {}".format(ips + 1, value)
        # print "boom, k0 = {}".format(k0)

        alpha = bbb_a[k0 - 1]
        zworkv = vvv[:, k0].copy()
        beta = rnorm
        force_reorth = True

        if (k0 < k) and (beta * delta < anorm * EPS):
            full_reorth = True
            ierr = k0

        indices[0] = 0
        indices[1] = k0 - 1 
        indices[2] = k0 

        # PS - In the Fortran, MGS is hardcoded to 1 and the MGS == 0 case is 
        # commented out. I didn't port the dead code.
        zreorth.zreorth2(m, k0, uuu, uuu[:, k0], rnorm, indices)

        dset_mu(k0, mus, indices, EPSN2)


        # c     %--------------------------------------%
        # c     | Estimate ||B||_2^2 as ||B^T * B||_1  |   
        # c     %--------------------------------------%

        bbb_b[k0 - 1] = beta
        amax = 0.0

        for j in range(k0):
            amax = max(amax, bbb_a[j], bbb_b[j])
            if j == 0:
                anorm = max(anorm, FUDGE * alpha)
            elif j == 1:
                a1 = bbb_b[0] / amax
                a1 = FUDGE * amax * sqrt((bbb_a[0] / amax) ** 2 + a1 ** 2 + bbb_a[1] / amax * a1)
            else:
                a1 = bbb_a[j - 1] / amax
                b1 = bbb_b[j - 1] / amax
                a1 = FUDGE * amax * sqrt(a1**2 + b1**2 + a1 * bbb_b[j - 2] / amax + bbb_a[j] / amax * b1)
                anorm = max(anorm, a1)
            mus[j - 1] = EPSN2
            nus[j - 1] = EPSN2

        j0 = k0
        zwork = uuu[:, k0].copy()

    numax = 0.0
    mumax = 0.0


    
    # c  %-------------------------------------------%
    # c  | Start Lanczos bidiagonalization iteration |
    # c  %-------------------------------------------%            

    # FIXME check for off-by-one error
    try:
        range(j0, k)
    except:
        pdb.set_trace()


    for j in range(j0, k):
        if j == 0:
            aprodw.aprodw('t', ndp, m, n, zwork, zworkv, zworkz1, zworky1, lambda_, trlambda)

            alpha = dznrm2(zworkv)
            # bjs alpha = scipy.linalg.norm(zworkv)
            anorm = max(anorm, FUDGE * alpha)

        else:
            
            # c        %---------------------------------------------%
            # c        | alpha_{j} v_{j} = A'*u_{j} - beta_{j} v_{j} |
            # c        %---------------------------------------------%

            aprodw.aprodw('t', ndp, m, n, zwork, zworks, zworkz1, zworky1, lambda_, trlambda)

            # PS - simple (but slow) Python equivalent -- 
            # for i in range(n):
            #     zworkv[i] = zworks[i] - (beta * zworkv[i])

            # Fast Python equivalent -- 
            zworkv[:n] = zworks[:n] - (beta * zworkv[:n])            


            # c         %------------------------------------%
            # c         | Extended local reorthogonalization |
            # c         %------------------------------------%
            alpha = dznrm2(zworkv)
            # bjs alpha = scipy.linalg.norm(zworkv)

            # FIXME I think 'j > 0' is always True -- look at the surrounding if statement.

            if (j > 0) and (options["extended_reorth_iterations"] > 0) and (alpha < (KAPPA * beta)):
                for i in (0, 1):

                    s = np.vdot(vvv[:, j - 1], zworkv)

#                   call zaxpy(n,-s,V(1,j-1),1,zwork(iv),1)

                    # This is a slow inline Python equivalent of the call 
                    # to zaxpy --
                    # for k in range(n):
                    #     zworkv[k] += (-s * vvv[k, j - 1])

                    # This is a fast version --
                    zworkv[:n] += (-s * vvv[:n, j - 1])
                    # bjs zworkv[:n] = zaxpy(-s,vvv(1,j-1),n)
                    # bjs zworkv[:n] = zaxpy(-s,vvv(1,j-1))


                    if beta != 0.0:
                        beta += s
                        bbb_b[j - 1] = beta.real
                    s = scipy.linalg.norm(zworkv)
                    # bjs s = dznrm2(zworkv)

                    if s >= KAPPA * alpha:
                        # Bail out of for loop
                        break

                nus[j - 1] = EPS
                alpha = s

            bbb_a[j] = alpha
            amax = max(amax, alpha)

            # c     %---------------------------%
            # c     | Update estimate of ||A||_2 |         
            # c     %---------------------------%

            if j == 1:
                a1 = bbb_b[0] / amax
                a1 = FUDGE * amax * sqrt( (bbb_a[0] / amax) ** 2 + (a1 ** 2) + bbb_a[1] / amax * a1)
            else:
                a1 = bbb_a[j - 1] / amax
                b1 = bbb_b[j - 1] / amax
                a1 = FUDGE * amax * sqrt(a1**2 + b1**2 + a1 * bbb_b[j - 2] / amax + bbb_a[j] / amax * b1)
            anorm = max(anorm, a1)

        # c     %--------------------------%
        # c     | Update the nu recurrence |
        # c     %--------------------------%

        if (not full_reorth) and (alpha != 0.0):

            numax = dupdate_nu(mus, nus, j, bbb_a, bbb_b, anorm)

        # c     %------------------------------%
        # c     | Reorthogonalize if necessary |
        # c     %------------------------------%

        if (full_reorth or force_reorth or (numax > delta)) and (alpha != 0.0):
            # print "Reorthogonalizing..."

            if full_reorth or (eta == 0.0):
                indices[0] = 0
                indices[1] = j - 1
                indices[2] = j 
                indices[3] = -1
                indices[4] = -1
            elif not force_reorth:
                indices = dcompute_int(nus, j - 1, delta, eta)

            # Note: MGS is hardcoded to 1
            zreorth.zreorth2(n, j - 1, vvv, zworkv, alpha, indices)

            dset_mu(j - 1, nus, indices, EPSN2)
            numax = eta

            force_reorth = not force_reorth

        # c     %-----------------------------------------------%
        # c     | Check whether an invariant subspace was found |
        # c     %-----------------------------------------------% 

        if (alpha < (anorm * EPSN)) and (j < k):
            rnorm = alpha
            alpha = 0.0

            # c     %------------------------------------------------%
            # c     | Try to build and orthogonal subspace, starting |
            # c     | with a random vector.                          |
            # c     %------------------------------------------------%
            
            alpha, anormest, ierr = zget0w.zgetu0w('t', ndp, m, n, j - 1, 3, zworkv, vvv, options["classical_gs"], lambda_, trlambda)

            if alpha == 0.0:
                k = j - 1
                ierr = -j
                # FIXME this should merely bail out of the loop, not raise an 
                # error.
                raise ValueError
            else:
                zsafescal.zsafescal(n, alpha, zworkv)
                alpha = 0.0
                force_reorth = True
                if delta > 0.0:
                    full_reorth = False

        elif (j > 1) and (not full_reorth) and (j < k) and (delta * alpha < anorm * EPS):
            ierr = j

        bbb_a[j] = alpha

        if alpha != 0.0:
            zsafescal.zsafescal(n, alpha, zworkv)

        # PS - This is the slow but simple Python equivalent of the fortran 
        # call to zcopy() --
        # for i, value in enumerate(zworkv):
        #     vvv[i, j] = value

        # This is the fast version --
        vvv[:, j] = zworkv.copy()

 
        # c     %------------------------------------------------%
        # c     | beta_{j+1} u_{j+1} = A*v_{j} - alpha_{j} u_{j} |
        # c     %------------------------------------------------%

             
        aprodw.aprodw('n', ndp, m, n, zworkv, zworks, zworkz1, zworky1, lambda_, trlambda)

        zwork[:m] = zworks[:m] - (alpha * zwork[:m])


        # c     %------------------------------------%
        # c     | Extended local reorthogonalization |
        # c     %------------------------------------%

        #       do i=1,ioption(2)
        # C        PS - call to zdotc() does not modify U
        #          s = zdotc(m,U(1,j),1,zwork(iu),1)
        
        # C        PS - I don't think this call to zaxpy() modifies U
        #          call zaxpy(m,-s,U(1,j),1,zwork(iu),1)
        #          if (alpha .ne. zero) then
        #             alpha = alpha + dble(s)
        
        #             B(j,1) = alpha
        #          endif
        #          work(imu+j-1) = eps
        #       enddo

        for i in range(options["extended_reorth_iterations"]):

            s = np.vdot(uuu[:, j], zwork)

            # PS - This is the naive (and slow) equivalent of the fortran 
            # call to zaxpy() --
            # for ips in range(m):
            #     zwork[ips] = zwork[ips] + (-s * uuu[ips, j])

            # PS - This is the fast version --
            zwork[:m] += (-s * uuu[:m, j])

            if alpha != 0.0:
                alpha += s
                bbb_a[j] = alpha.real

            # FIXME - why is this inside a for loop? It is invariant. The
            # Fortran code does the same.
            mus[j] = EPS

        beta = scipy.linalg.norm(zwork)
        # bjs beta = dznrm2(zwork)

        bbb_b[j] = beta
        amax = max(amax, beta)
      
        # c     %---------------------------%
        # c     | Update estimate of ||A||_2 |         
        # c     %---------------------------%

        if j == 0:
#            a1 = util.dlapy2(bbb_a[0], bbb_b[0])
            a1 = sqrt(bbb_a[0]**2 + bbb_b[0]**2)
        else:
            a1 = bbb_b[j] / amax
            a1 = amax * sqrt(a1 ** 2 + (bbb_b[j] / amax) ** 2 + a1 * bbb_b[j - 1] / amax)

        anorm = max(anorm, a1)

        # c     %--------------------------%
        # c     | Update the mu recurrence |
        # c     %--------------------------%

        if (not full_reorth) and (beta != 0.0):
            mumax = dupdate_mu(mus, nus, j, bbb_a, bbb_b, anorm)

        # c     %--------------------------------------%
        # c     | Reorthogonalize u_{j+1} if necessary |
        # c     %--------------------------------------%

        if (beta != 0.0) and (full_reorth or force_reorth or (mumax > delta)):
            if full_reorth or (eta == 0.0):
                indices[0] = 0
                indices[1] = j
                indices[2] = j + 1
                indices[3] = -1
                indices[4] = -1
            elif not force_reorth:
                indices = dcompute_int(mus, j, delta, eta)
            else:
                # ==> force_reorth = True
                for i in range(j + 1):
                    if indices[i] == j:
                        indices[i] = j + 1
                        break

            zreorth.zreorth2(m, j, uuu, zwork, beta, indices)
            
            dset_mu(j, mus, indices, EPSN2)       

            mumax = eta
            force_reorth = not force_reorth


        # c     %-----------------------------------------------%
        # c     | Check whether an invariant subspace was found |
        # c     %-----------------------------------------------%

        if (beta < (anorm * EPSN)) and (j < k):
            rnorm = beta
            beta = 0.0
            # c         %------------------------------------------------%
            # c         | Try to build an orthogonal subspace, starting |
            # c         | with a random vector.                          |
            # c         %------------------------------------------------%

            beta, anormest, ierr = zget0w.zgetu0w('n', ndp, m, n, j, 3, zwork, uuu, options["classical_gs"], lambda_, trlambda)
            if beta == 0.0:
                k = j
                ierr = -j
                # FIXME This should merely exit the loop, not raise an error
                raise ValueError
            else:
                zsafescal.zsafescal(n, beta, zwork)
                beta = 0.0
                force_reorth = True
                if delta > 0:
                    full_reorth = False

        elif not full_reorth and (j < k) and ((delta * beta) < (anorm * EPS)):
            ierr = j 

        bbb_b[j] = beta

        if (beta != 0) and (beta != 1):
            zsafescal.zsafescal(m, beta, zwork)

        # call zcopy(m,zwork(iu),1,U(1,j+1),1)
        #
        # PS - this is the naive equivalent of the call to zcopy() -- 
        # for ips in range(m):
        #     uuu[ips, j + 1] = zwork[ips]
        #
        # PS - this is the fast version
        uuu[:m, j + 1] = zwork[:m].copy()

        rnorm = beta

    options["anorm"] = anorm

    # print "zlanbprow done"

    return ierr, rnorm

   
#==============================================================================

def dset_mu(k, mus, indices, val):

    i = 0
    while (indices[i] <= k) and (indices[i] >= 0):
        p = indices[i]
        q = indices[i + 1]
        # FIXME check for off-by-one error
        for j in fortran_range(p, q):
            mus[j] = val

        i += 2



#==============================================================================

def dcompute_int(mus, j, delta, eta):
    """I think "int" is short for "intervals". returns indices."""

    indices = [ ]

    ip = 0
    i  = -1

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
    # else:
    #     delta < eta, so raise an error?

    
    return indices



#==============================================================================

def dupdate_mu(mus, nus, j, alpha, beta, anorm):
    """
    Alters mus, returns mumax. I'm not sure if this is truly the max value
    in abs(mus) or just the max of the subset that this function touches.
    
    """
    # In the Fortran, forumla is hardcoded to 1 
    if j == 0:
        # FIXME is this a bug in the fortran? d is computed but not used.
#        d = EPS1 * (util.dlapy2(alpha[j], beta[j]) + alpha[0]) + EPS1 * anorm
        
        d = EPS1 * (sqrt(alpha[j]**2 + beta[j]**2) + alpha[0]) + EPS1 * anorm
        
        mus[0] = EPS1 / beta[0]
        mumax = abs(mus[0])

    else:

        mus[0] = alpha[0] * nus[0] - alpha[j] * mus[0]
#        d = EPS1 * (util.dlapy2(alpha[j], beta[j]) + alpha[0]) + EPS1 * anorm
        
        d = EPS1 * (sqrt(alpha[j]**2 + beta[j]**2) + alpha[0]) + EPS1 * anorm
        
        mus[0] = (mus[0] + copysign(d, mus[0])) / beta[j]
        mumax = abs(mus[0])

        # FIXME dlapy2(alpha[j], beta[j]) is a loop invariant - maybe calc once and reuse?
        the_range = fortran_range(1, j - 1)
        for k in the_range:
            mus[k] = (alpha[k] * nus[k]) + (beta[k - 1] * nus[k - 1]) - (alpha[j] * mus[k])
#            d  = EPS1 * (util.dlapy2(alpha[j], beta[j]) + util.dlapy2(alpha[k], beta[k - 1])) 
            
            d  = EPS1 * (sqrt(alpha[j]**2 + beta[j]**2) + sqrt(alpha[k]**2 + beta[k - 1]**2))
            
            d += EPS1 * anorm
            mus[k] = (mus[k] + math.copysign(d, mus[k])) / beta[j]
            mumax = max(mumax, abs(mus[k]))

        mus[j] = beta[j - 1] * nus[j - 1]
#        d  = EPS1 * (util.dlapy2(alpha[j], beta[j]) + util.dlapy2(alpha[j], beta[j - 1])) 
        
        d  = EPS1 * (sqrt(alpha[j]**2 + beta[j]**2) + sqrt(alpha[j]**2 + beta[j - 1]**2))
        
        d += EPS1 * anorm
        mus[j] = (mus[j] + math.copysign(d, mus[j])) / beta[j]
        mumax = max(mumax, abs(mus[j]))

    mus[j + 1] = 1.0

    return mumax




#==============================================================================

def dupdate_nu(mus, nus, j, alpha, beta, anorm):
    """Alters nus, returns numax"""

    # In the Fortran, forumla is hardcoded to 1 
    numax = 0.0
    if j > 0:
        for k in fortran_range(0, j - 1):
            nus[k] = (beta[k] * mus[k + 1]) + (alpha[k] * mus[k]) - (beta[j - 1] * nus[k])
#            d = EPS1 * (util.dlapy2(alpha[k], beta[k]) + util.dlapy2(alpha[j], beta[j - 1]))
            
            d = EPS1 * (sqrt(alpha[k]**2 + beta[k]**2) + sqrt(alpha[j]**2 + beta[j - 1]**2))
            
            d += EPS1 * anorm
            nus[k] = (nus[k] + math.copysign(d, nus[k])) / alpha[j]
            numax = max(numax, abs(nus[k]))

        nus[j] = 1.0

    return numax
