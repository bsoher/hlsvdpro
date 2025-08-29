#
#     (C) Brian J Soher, 2020
#
from __future__ import division, print_function

import numpy as np

from dbsvd         import dbdqr, drefinebounds
from zblasext      import pdznrm2
from dbdsqr        import dbdsqr
from zritzvec      import zritzvec
from zgetu0        import zgetu0
from zlanbpro      import zlanbpro

# PS - Machine- and compiler-dependent constants
# On my Mac, EPS in Python is 2.22044604925e-16. In Fortran
# it is half as big (1.11022302462515654E-016).

#EPS   = 1.1102230246251565E-016  #np.finfo(np.float64).eps
#EPS34 = EPS ** (3.0 / 4.0)
#EPS1  = 100 * EPS

# Variable name translations
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# k          nsv_sought         # of singular values requested by caller
# kfit       nsv_found          # of singular values found
# lsinval    singular_values    Array containing the singular values (floats)
# ndp        n_data_points
# kmax       kmax               passed from wrapper function, default 5*k
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




#def zlansvdw(n_data_points, m, n, nsv_sought, lambda_, trlambda, uuu, tol=16e-16):


def zlansvd(jobu,jobv,m,n,k,kmax,aprod,U,V,tolin, doption,ioption,parm):
    """
    DLANSVD: Compute the leading singular triplets of a large and
    sparse matrix by Lanczos bidiagonalization with partial
    reorthogonalization.

    Parameters:

    JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
          Otherwise the array U is not touched.
    JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular
          vectors. Otherwise the array V is not touched.
    M: INTEGER. Number of rows of A.
    N: INTEGER. Number of columns of A.
    K: INTEGER. Number of desired singular triplets. K <= MIN(KMAX,M,N)
    KMAX: INTEGER. maximal number of iterations / maximal dimension of
          Krylov subspace.
    APROD: Subroutine defining the linear operator A.
           APROD should be of the form:

          SUBROUTINE DAPROD(TRANSA,M,N,X,Y,ZPARM,IPARM)
          CHARACTER*1 TRANSA
          INTEGER M,N,IPARM(*)
          DOUBLE PRECISION X(*),Y(*),ZPARM(*)

          If TRANSA.EQ.'N' then the function should compute the matrix-vector
          product Y = A * X.
          If TRANSA.EQ.'T' then the function should compute the matrix-vector
          product Y = A^T * X.
          The arrays IPARM and ZPARM are a means to pass user supplied
          data to APROD without the use of common blocks.
    U(LDU,KMAX+1): DOUBLE COMPLEX array. On return the first K columns of U
              will contain approximations to the left singular vectors
              corresponding to the K largest singular values of A.
              On entry the first column of U contains the starting vector
              for the Lanczos bidiagonalization. A random starting vector
              is used if U is zero.
    LDU: INTEGER. Leading dimension of the array U. LDU >= M.
    SIGMA(K): DOUBLE PRECISION array. On return Sigma contains approximation
              to the K largest singular values of A.
    BND(K)  : DOUBLE PRECISION array. Error estimates on the computed
              singular values. The computed SIGMA(I) is within BND(I)
              of a singular value of A.
    V(LDV,KMAX): DOUBLE COMPLEX array. On return the first K columns of V
              will contain approximations to the right singular vectors
              corresponding to the K largest singular values of A.
    LDV: INTEGER. Leading dimension of the array V. LDV >= N.
    TOLIN: DOUBLE PRECISION. Desired relative accuracy of computed singular
           values. The error of SIGMA(I) is approximately
           MAX( 16*EPS*SIGMA(1), TOLIN*SIGMA(I) )
    WORK(LWORK): DOUBLE PRECISION array. Workspace of dimension LWORK.
    LWORK: INTEGER. Dimension of WORK.
           If JOBU.EQ.'N' and JOBV.EQ.'N' then  LWORK should be at least
           M + N + 9*KMAX + 2*KMAX**2 + 4 + MAX(M,N,4*KMAX+4).
           If JOBU.EQ.'Y' or JOBV.EQ.'Y' then LWORK should be at least
           M + N + 9*KMAX + 5*KMAX**2 + 4 +
           MAX(3*KMAX**2+4*KMAX+4, NB*MAX(M,N)), where NB>0 is a block
           size, which determines how large a fraction of the work in
           setting up the singular vectors is done using fast BLAS-3
           operation.
    ZWORK: DOUBLE COMPLEX array of dimension ????.
    IWORK: INTEGER array. Integer workspace of dimension LIWORK.
    LIWORK: INTEGER. Dimension of IWORK. Should be at least 8*KMAX if
            JOBU.EQ.'Y' or JOBV.EQ.'Y' and at least 2*KMAX+1 otherwise.
    DOPTION: DOUBLE PRECISION array. Parameters for LANBPRO.
       doption(1 -> 0) = delta. Level of orthogonality to maintain among
         Lanczos vectors.
       doption(2 -> 1) = eta. During reorthogonalization, all vectors with
         with components larger than eta along the latest Lanczos vector
         will be purged.
       doption(3 -> 2) = anorm. Estimate of || A ||.
    IOPTION: INTEGER array. Parameters for LANBPRO.
       ioption(1 -> 0) = CGS.  If CGS.EQ.1 then reorthogonalization is done
         using iterated classical GRAM-SCHMIDT. IF CGS.EQ.0 then
         reorthogonalization is done using iterated modified Gram-Schmidt.
       ioption(2 -> 1) = ELR. If ELR.EQ.1 then extended local orthogonality is
         enforced among u_{k}, u_{k+1} and v_{k} and v_{k+1} respectively.
    INFO: INTEGER.
        INFO = 0  : The K largest singular triplets were computed succesfully
        INFO = J>0, J<K: An invariant subspace of dimension J was found.
        INFO = -1 : K singular triplets did not converge within KMAX
                    iterations.
    ZPARM: DOUBLE COMPLEX array. Array used for passing data to the APROD
        function.
    IPARM: INTEGER array. Array used for passing data to the APROD
        function.

    (C) Brian J Soher, 2020

    """

    ldu = 1
    ldv = 1

    m,n = aprod.shape
    
    # options replaces the Fortran arrays ioption and doption.
    options = {}
    # Use classical Gram-Schmidt reorth instead of the modified version.
    options["cgs"]   = ioption[0]   # def = 1
    options["exri"]  = ioption[1]   # extended reorth iterations def. 1
    options["delta"] = doption[0]   # 1e-12
    options["eta"]   = doption[1]   # 1e-14
    options["anorm"] = doption[2]   # 0.0

    EPS = 1.1102230246251565E-016  # np.finfo(np.float64).eps
    EPS34 = EPS ** (3.0 / 4.0)
    EPSN = max(m,n) * EPS / 2.0

    # print('eps, eps34, epsn = ',EPS, EPS34, EPSN, )

    if U is None:
        U = np.zeros( (m, kmax + 1), np.complex128)
    if V is None:
        V = np.zeros( (n, kmax), np.complex128)

    bnd = np.ndarray([k,], dtype='d', order='F')
    
    # %--------------------------------%
    # | Guard against absurd arguments |
    # %--------------------------------%

    lanmax = min(n + 1, m + 1, kmax)
    tol = min(1.0, max(16.0 * EPS, tolin))
    anorm = 0.0

    # # %------------------------------%
    # # | Set pointers into work array |
    # # %------------------------------%
    # ibnd = 1
    # ib   = ibnd + lanmax+1
    # ib1  = ib   + 2*lanmax
    # ip   = ib1  + 2*lanmax
    # iq   = ip   + (lanmax+1)**2
    # iwrk = iq   + lanmax**2
    # lwrk = lwork-iwrk+1

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
    len_work1 = lanmax + 1

    # In the Fortran code, work(ib) starts an array of length lanmax * 2. It's
    # sometimes declared as matrix of dimension (lanmax, 2) but it's never used
    # as a matrix. It's really just two arrays (of length lanmax) bundled
    # together under one variable name. It's confusing, and in Python I separate
    # them into bbb_a and bbb_b.
    # The chunk of work starting at work(ib1) is just used for a work copy of
    # work(ib), and has the same properties. I call the corresponding Python
    # arrays bbb1_a and bbb1_b.
    len_bbb = lanmax
    len_workp = (lanmax + 1) ** 2

    # Note that before zlansdvw.f's second call to dbdsqr(), it writes a bunch
    # of 1s into work starting at index iq. It goes waaaay beyond lanmax**2
    # and into the work(iwrk) section. Intentional? I don't know. The size I
    # give to workq here should be about right, erring on the side of safety.
    len_workq = n * (lanmax + 1)

    # C     PS - changed from lwork to lzwrk

    work1 = np.zeros((len_work1,), dtype=np.float64)  # work(ibnd)
    workp = np.zeros((len_workp,), dtype=np.float64)  # work(ip)
    workq = np.zeros((len_workq,), dtype=np.float64)  # work(iq)

    bbb_a  = np.zeros((len_bbb,), dtype=np.float64)  # work(ib)
    bbb_b  = np.zeros((len_bbb,), dtype=np.float64)  # work(ib  + lanmax)
    bbb1_a = np.zeros((len_bbb,), dtype=np.float64)  # work(ib1)
    bbb1_b = np.zeros((len_bbb,), dtype=np.float64)  # work(ib1 + lanmax)

    # %---------------------------------------------------------------%
    # | Set up random starting vector if none is provided by the user |
    # %---------------------------------------------------------------%

    rnorm = pdznrm2(m,U[:,0],1)
    if rnorm == 0:
        rnorm, anorm, ierr = zgetu0('n', m, n, 0, 1, U[:,0], U, ldu, aprod, parm, options["cgs"])

    info = 0
    neig = 0
    jold = 0
    j = min(k+max(8, k)+1, lanmax)

    # %------------------------------%
    # | Iterate until convergence... |
    # %------------------------------%
    while neig < k:

        # %---------------------------------------------------%
        # | Compute bidiagonalization A*V_{j} = U_{j+1}*B_{j} |
        # %---------------------------------------------------%

        ierr, rnorm, jnew = zlanbpro( m, n, jold, j, aprod, U, ldu, V, ldv, bbb_a, bbb_b, rnorm, options, parm)

        jold = j
        j = jnew

        # %---------------------------------------------%
        # | Compute and analyze SVD(B) and error bounds |
        # %---------------------------------------------%
        bbb1_a = bbb_a.copy()
        bbb1_b = bbb_b.copy()
        work1[:j+1] *= 0

        # v2.1 In dbsvd.py dbdqr() has a condition where c1 and c2 are not
        # calculated, so I need to recreate the values for work1[j] and
        # work1[j-1] to remain the same if this happens. I will just set
        # c1 and c2 to the current work1 values going in!

        c1 = work1[j-1]
        c2 = work1[j]
        c1, c2 = dbdqr((j==min(m,n)), 'N', j, bbb1_a, bbb1_b, c1, c2, workp, lanmax+1)

        work1[j-1] = c1
        work1[j] = c2

        r = dbdsqr('u', j, 0, 1, 0, bbb1_a, bbb1_b, None, 1, work1, 1, None, 1)
        info, bbb1_a, bbb1_b, _, work1, _ = r

        # PS - dbdsqr() returns bbb1_a, bbb1_b, and work1 as lists that are
        # shorter than they were on the way in. Subsequent code (refinebounds()
        # and probably other code too) doesn't deal with this size change,
        # so I have to pad them back out to their original size.
        bbb1_a +=  [0.0] * (lanmax - len(bbb1_a) + 1)
        bbb1_b +=  [0.0] * (lanmax - len(bbb1_b) + 1)
        # print("len(bbb1_b) = {}".format(len(bbb1_b)))
        bbb1_a = np.array(bbb1_a)
        bbb1_b = np.array(bbb1_b)
        work1 += [0.0] * (len_work1 - len(work1))
        work1 = np.array(work1)

        if j > 5:
            anorm = bbb1_a[0]
        else:
            anorm = max(anorm, bbb1_a[0])

        work1 *= rnorm
        np.absolute(work1, work1)   # copies results in place


        # %---------------------------------------------%
        # | Refine error bounds using the "Gap theorem" |
        # %---------------------------------------------%

        drefinebounds(min(m,n), j, bbb1_a, work1, EPSN*anorm, EPS34)


        # %----------------------------------------------------%
        # | Determine the number of converged singular values  |
        # %----------------------------------------------------%

        for i in range(min(j, k)):
            bnd[i] = work1[i]       # need for loop to copy in place bnd[] array

        sigma = []
        i = 0

        while i < min(j, k):
            if work1[i] <= (tol * bbb1_a[i]):
                sigma.append(bbb1_a[i])
                i += 1
            else:
                i = k       # Force loop exit

        neig = len(sigma)

        # %--------------------------------------------------%
        # | Test if an invariant subspace has been found or |
        # | if the workspace has been exhausted.             |
        # %--------------------------------------------------%

        if ierr < 0:
            if j < k:
                print( 'WARNING: Invariant subspace found.',' Dimension = ',j)
                info = j
            break       # Bail out of loop

        if j >= lanmax:
            if neig < k:
                print( 'WARNING (py): Maximum dimension of Krylov',
                       ' subspace exceeded prior to convergence.',
                       ' Try increasing KMAX.')
                print( 'neig = ',neig)
                info = -1
            break       # Bail out of loop

        # %----------------------------------------------------%
        # | Increase the dimension of the Krylov subspace.     |
        # | If any Ritz values have converged then try to      | 
        # | estimate the average number of iterations per      |
        # | converged Ritz value.                              |
        # | Else increase the dimension by 50%.                |
        # %----------------------------------------------------%

        if neig > 1:
            dj = min(j // 2, ((k - neig) * (j - 6)) // ( 2 * neig + 1))
            dj = min(100, max(2, dj))
        else:
            dj = j // 2
            dj = min(100, max(10, dj))

        j = min(j + dj, lanmax)


##############      while loop has ended      ##################


    # %-----------------------------------------%
    # | Calculate singular vectors if requested %
    # %-----------------------------------------%

    if ((neig >= k or info > 0) and (jobu=='y' or jobv=='y')):

        zritzvec('L', jobu, jobv, m, n, neig, jold, bbb_a, bbb_b, bbb1_a, U, ldu, V, ldv, kmax)

    return U, sigma, V, info, bnd



