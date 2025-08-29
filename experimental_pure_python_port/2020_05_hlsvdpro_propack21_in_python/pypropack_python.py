#
#     (C) Brian J Soher, 2020
#

from __future__ import division, print_function, absolute_import

__all__ = ['svdp']

import cProfile

# 3rd party modules
import numpy as np
from scipy.sparse.linalg import aslinearoperator

import propack_python.zlansvd as zlansvd


class _AProd(object):
    """
    Wrapper class for linear operator

    The call signature of the __call__ method matches the callback of
    the PROPACK routines.
    """

    def __init__(self, A):
        try:
            self.aa = A.copy()
            self.A = aslinearoperator(A)
        except TypeError:
            self.aa = np.asarray(A.copy())
            self.A = aslinearoperator(np.asarray(A))

    def __call__(self, transa, m, n, x, y, parm, iparm):
        if transa == 'n':
            y[:m] = self.A.matvec(x[:n])
            # res = my_zgemv(transa,m,n,1+0j,self.aa,x,1,0+0j,y,1)
        else:
            y[:n] = self.A.rmatvec(x[:m])

    @property
    def shape(self):
        return self.A.shape

    @property
    def dtype(self):
        try:
            return self.A.dtype
        except AttributeError:
            return self.A.matvec(np.zeros(self.A.shape[1])).dtype


def my_zgemv(trans, m, n, alpha, a, x, incx, beta, y, incy):
    """
    ZGEMV(TRANS,M,N,ALPHA,A,LDA,X,INCX,BETA,Y,INCY)

    """

    if m == 0 or n == 0 or ((alpha == 0 + 0j) and (beta == 1 + 0j)):
        return

    noconj = trans in ['t', 'T']

    #     Set  LENX  and  LENY, the lengths of the vectors x and y, and set
    #     up the start points in  X  and  Y.

    if trans in ['n', 'N']:
        lenx = n
        leny = m
    else:
        lenx = m
        leny = n

    if incx > 0:
        kx = 0
    else:
        kx = 0 - (lenx - 1) * incx

    if incy > 0:
        ky = 0
    else:
        ky = 0 - (leny - 1) * incy

    # Start the operations. In this version the elements of A are
    # accessed sequentially with one pass through A.
    #
    #  First form  y := beta*y.

    if beta != 1 + 0j:
        if incy == 1:
            if beta == 0 + 0j:
                for i in range(leny):
                    y[i] = 0 + 0j
            else:
                for i in range(leny):
                    y[i] = beta * y[i]
        else:
            iy = ky
            if beta == 0 + 0j:
                for i in range(leny):
                    y[iy] = 0 + 0j
                    iy = iy + incy
            else:
                for i in range(leny):
                    y[iy] = beta * y[iy]
                    iy = iy + incy

    if alpha == 0 + 0j:
        return

    if trans in ['n', 'N']:

        # Form  y := alpha*A*x + y.
        jx = kx
        if incy == 1:
            for j in range(n):
                temp = alpha * x[jx]
                for i in range(m):
                    y[i] = y[i] + temp * a[i, j]
                jx = jx + incx
        else:
            for j in range(n):
                temp = alpha * x[jx]
                iy = ky
                for i in range(m):
                    y[iy] = y[iy] + temp * a[i, j]
                    iy = iy + incy
                jx = jx + incx
    else:
        # Form  y := alpha*A**T*x + y  or  y := alpha*A**H*x + y.
        jy = ky
        if incx == 1:
            for j in range(n):
                temp = 0 + 0j
                if noconj:
                    for i in range(m):
                        temp = temp + a[i, j] * x[i]
                else:
                    for i in range(m):
                        temp = temp + np.conjugate(a[i, j]) * x[i]

                y[jy] = y[jy] + alpha * temp
                jy = jy + incy
        else:
            for j in range(n):
                temp = 0 + 0j
                ix = kx
                if noconj:
                    for i in range(m):
                        temp = temp + a[i, j] * x[ix]
                        ix = ix + incx
                else:
                    for i in range(m):
                        temp = temp + np.conjugate(a[i, j]) * x[ix]
                        ix = ix + incx

                y[jy] = y[jy] + alpha * temp
                jy = jy + incy
    return y


def svdp(A, k, kmax=None,
         compute_u=True, compute_v=True, v0=None, full_output=False, tol=0,
         delta=None, eta=None, anorm=0, elr=True, blocksize=1,
         min_relgap=0.002, shifts=100, maxiter=1000):
    """Compute the singular value decomposition of A using PROPACK

    This function assumes that the user will create an external APROD function
    that the Fortran lansvd() or lansvd_irl() calls will access.

    Use the svdp_aprod() call if you don't want to provide aprod()

    Parameters
    ----------
    A : array_like, sparse matrix, or LinearOperator
        Operator for which svd will be computed.  If A is a LinearOperator
        object, it must define both ``matvec`` and ``rmatvec`` methods.
    k : int
        number of singular values/vectors to compute
    which : NOTE. not implemented in this version
        string (optional)
        which singluar triplets to compute:
        - 'L': compute triplets corresponding to the k largest singular values
        - 'S': compute triplets corresponding to the k smallest singular values
        which='S' requires irl=True.
        Computes largest singular values by default.
    irl_mode : NOTE. not implemented in this version
        boolean (optional)
        If True, then compute SVD using iterative restarts.  Default is False.
    kmax : int (optional)
        maximal number of iterations / maximal dimension of Krylov subspace.
        default is 5 * k
    compute_u : bool (optional)
        if True (default) then compute left singular vectors u
    compute_v : bool (optional)
        if True (default) then compute right singular vectors v
    tol : float (optional)
        The desired relative accuracy for computed singular values.
        If not specified, it will be set based on machine precision.
    v0 : starting vector (optional)
        Starting vector for iterations: should be of length A.shape[0].
        If not specified, PROPACK will generate a starting vector.
    full_output : boolean (optional)
        If True, then return info and sigma_bound.  Default is False
    delta : float (optional)
        Level of orthogonality to maintain between Lanczos vectors.
        Default is set based on machine precision.
    eta : float (optional)
        Orthogonality cutoff.  During reorthogonalization, vectors with
        component larger than eta along the Lanczos vector will be purged.
        Default is set based on machine precision.
    anorm : float (optional)
        estimate of ||A||.  Default is zero.
    cgs : NOTE. Not implemented at this moment. Issue with Fortran indexing.
        boolean (optional)
        If True, reorthogonalization is done using classical Gram-Schmidt.
        If False (default), it is done using modified Gram-Schmidt.
    elr : boolean (optional)
        If True (default), then extended local orthogonality is enforced
        when obtaining singular vectors.
    blocksize : int (optional)
        If computing u or v, blocksize controls how large a fraction of the
        work is done via fast BLAS level 3 operations.  A larger blocksize
        may lead to faster computation, at the expense of greater memory
        consumption.  blocksize must be >= 1; default is 1.
    min_relgap : float (optional)
        The smallest relative gap allowed between any shift in irl mode.
        Default = 0.001.  Accessed only if irl_mode == True.
    shifts : int (optional)
        Number of shifts per restart in irl mode.  Default is 100
        Accessed only if irl_mode == True.
    maxiter : int (optional)
        maximum number of restarts in irl mode.  Default is 1000
        Accessed only if irl_mode == True.

    Returns
    -------
    u : ndarray
        The top k left singular vectors, shape = (A.shape[0], 3),
        returned only if compute_u is True.
    sigma : ndarray
        The top k singular values, shape = (k,)
    vt : ndarray
        The top k right singular vectors, shape = (3, A.shape[1]),
        returned only if compute_v is True.
    info : integer
        convergence info, returned only if full_output is True
        - INFO = 0  : The K largest singular triplets were computed succesfully
        - INFO = J>0, J<K: An invariant subspace of dimension J was found.
        - INFO = -1 : K singular triplets did not converge within KMAX
                     iterations.
    sigma_bound : ndarray
        the error bounds on the singular values sigma, returned only if
        full_output is True


    """
    irl_mode = False  # hard set, not implemented in this version
    which = 'L'  # hard set, requires irl_mode to be implemented
    cgs = False  # hard set, not implemented in this version

    which = which.upper()
    if which not in ['L', 'S']:
        raise ValueError("`which` must be either 'L' or 'S'")
    if not irl_mode and which == 'S':
        raise ValueError("`which`='S' requires irl_mode=True")

    aprod = _AProd(A)
    typ = aprod.dtype.char

    if typ != np.dtype(complex).char:
        raise ValueError("A array is not complex128 type")

    try:
        lansvd = zlansvd.zlansvd
    except:
        raise ValueError("Can't assign lansvd function")

    m, n = aprod.shape

    if (k < 1) or (k > min(m, n)):
        raise ValueError("k must be positive and not greater than m or n")

    if kmax is None:
        kmax = 5 * k

    # guard against unnecessarily large kmax
    kmax = min(m + 1, n + 1, kmax)

    if kmax < k:
        raise ValueError("kmax must be greater than or equal to k")

    if compute_u:
        jobu = 'y'
    else:
        jobu = 'n'

    if compute_v:
        jobv = 'y'
    else:
        jobv = 'n'

    # these will be the output arrays
    u = np.zeros((m, kmax + 1), order='F', dtype=typ)
    v = np.zeros((n, kmax), order='F', dtype=typ)

    # Specify the starting vector.  if v0 is all zero, PROPACK will generate
    # a random starting vector: the random seed cannot be controlled in that
    # case, so we'll instead use numpy to generate a random vector

    # b = np.array([0.39574246 + 8.64960398e-04j, -0.92272058 - 9.16567150e-01j,
    #             0.11759638 - 2.99626252e-01j, 0.90382696 - 2.50451048e-01j,
    #             0.33224741 - 2.90239220e-01j, 0.65236429 + 8.80067891e-01j,
    #             -0.71063568 - 8.31851118e-01j, -0.34365925 + 4.61225337e-01j,
    #             -0.58001341 - 7.88448776e-02j, 0.02353656 + 5.01997418e-01j,
    #             -0.8337296 + 4.95846704e-01j, -0.24152437 + 6.93621382e-01j,
    #             -0.88483768 + 3.52059160e-01j, 0.12335653 - 6.49370645e-01j,
    #             -0.57815739 + 9.14450396e-01j, -0.64556854 + 5.31735418e-01j,
    #             -0.8890779 + 2.50473465e-01j, -0.94626406 + 1.15594701e-01j,
    #             -0.0032114 - 1.93704580e-01j, 0.982657 - 5.06614451e-01j,
    #             0.2551178 + 9.00256757e-01j, 0.49246474 + 5.54539609e-01j,
    #             0.81632678 + 8.81443234e-01j, -0.53894361 + 5.25947376e-01j,
    #             -0.42540597 - 2.05983196e-01j, 0.51321294 + 3.45785800e-01j,
    #             -0.12606551 + 4.75812526e-01j, 0.75708181 - 3.97062401e-01j,
    #             0.93164061 - 6.82958942e-01j, -0.32918585 + 2.68986279e-01j,
    #             -0.99988349 - 9.08570727e-01j, -0.35626236 - 2.65835788e-01j,
    #             -0.51331877 - 3.47527898e-01j, -0.33542511 - 1.01409189e-01j,
    #             0.40823508 - 2.08655536e-01j, 0.52845909 - 2.42148510e-01j,
    #             -0.96029183 - 2.64879512e-01j, 0.6032914 + 1.74623095e-01j,
    #             -0.86712007 + 8.03955829e-01j, 0.67776572 - 9.12164531e-01j,
    #             0.4714103 - 5.51504957e-01j, -0.14340699 + 6.44959103e-01j,
    #             0.74969417 + 6.83657363e-01j, 0.14652763 - 3.56386870e-01j,
    #             0.90337919 + 8.76818220e-01j, 0.68816882 - 1.47016472e-01j,
    #             0.28900988 + 6.16928026e-01j, 0.0566198 - 5.04284020e-01j,
    #             0.84521551 + 2.90981578e-01j, 0.32737165 - 7.25811194e-01j])

    if v0 is None:
        u[:, 0] = np.random.random(m)
        if np.iscomplexobj(np.zeros(0, dtype=typ)):  # complex type
            u[:, 0] += 1j * np.random.random(m)
    # # bjs debug
    #     # orig rnd array u[:,0] = np.array([0.25662954773130986+0.57062648579991326*1j, 0.20366349962301655+0.17718181865061933*1j, 0.13613447425116998+0.34298286063336469*1j, 0.62030760743804181+0.60769243986903876*1j, 0.31760778194161077+0.70415718758479073*1j, 0.29932099778923460+0.66607309269575221*1j])
    #     u[:,0] = b[:m]
    else:
        try:
            u[:, 0] = v0
        except:
            raise ValueError("v0 must be of length %i" % m)

    # process options for the fit
    if delta is None:
        delta = np.finfo(typ).eps ** 0.5
    if eta is None:
        eta = np.finfo(typ).eps ** 0.75

    if irl_mode:
        doption = np.array([delta, eta, anorm, min_relgap], dtype=typ.lower())
    else:
        doption = np.array([delta, eta, anorm], dtype=typ.lower())

    ioption = np.array([int(bool(cgs)), int(bool(elr))], dtype='i')

    # Determine lwork & liwork:
    # the required lengths are specified in the PROPACK documentation
    if compute_u or compute_v:
        lwork = (m + n + 10 * kmax + 5 * kmax ** 2 + 4 +
                 max(3 * kmax ** 2 + 4 * kmax + 4,
                     max(1, int(blocksize)) * max(m, n)))
        liwork = 8 * kmax
    else:
        lwork = (m + n + 10 * kmax + 2 * kmax ** 2 + 5 +
                 max(m + n, 4 * kmax + 4))
        liwork = 2 * kmax + 1
    work = np.zeros(lwork, dtype=typ.lower())
    iwork = np.zeros(liwork, dtype='i')

    # dummy arguments: these are passed to aprod, and not used in this wrapper
    parm = np.zeros(1, dtype=typ.lower())

    u, sigma, v, info, bnd = lansvd(jobu, jobv, m, n, k, kmax, aprod, u, v, tol, doption, ioption, parm)

    neig = len(sigma)

    # construct return tuple
    ret = ()
    if compute_u:
        ret = ret + (u[:, :k],)
    ret = ret + (sigma,)
    if compute_v:
        ret = ret + (v[:, :k].conj().T,)
    if full_output:
        ret = ret + (info, bnd)

    return ret


#------------------------------------------------------------------------------

def test():

    import scipy

    k = 10
    M = 18
    X = [0.41896704+0.99980257j, 0.42153555+0.12091503j, 0.21960718+0.68795215j, 0.5484147 +0.1769976j ,
         0.26129899+0.56660992j, 0.56709373+0.55560421j, 0.89536965+0.19176404j, 0.65056708+0.69196265j,
         0.73759826+0.41295692j, 0.54597543+0.23774541j, 0.17429778+0.17341309j, 0.90913384+0.62080491j,
         0.47303659+0.01698477j, 0.19599911+0.53031064j, 0.27059218+0.49050355j, 0.4968034 +0.44224764j,
         0.87622486+0.61846116j, 0.75811665+0.59743874j, 0.36398074+0.51103615j, 0.29445307+0.44056691j,
         0.00661976+0.76269535j, 0.61738681+0.99455311j, 0.4132981 +0.87536523j, 0.41927287+0.42516633j,
         0.4473831 +0.14482903j, 0.46747749+0.24649163j, 0.25166896+0.21109008j, 0.73257403+0.94693629j,
         0.7184696 +0.3033616j , 0.64551375+0.58885844j, 0.5199645 +0.02910201j, 0.82232501+0.51692768j,
         0.73125811+0.7325245j , 0.96816824+0.14838221j, 0.34608904+0.96139791j, 0.77807039+0.19527549j,
         0.02245328+0.92376054j, 0.79903798+0.32057331j, 0.37806417+0.89394476j, 0.38891312+0.89303679j]

    sing0 = np.array([14.53067919,  3.13876274,  2.88112334,  2.48764528,  2.33695611,
                       2.14050356,  1.99244337,  1.79766561,  1.62749004,  1.46056729,
                       1.35906006,  1.22204062,  1.08330475,  0.93928118,  0.72308523,
                       0.6785563 ,  0.64182857,  0.54659049,  0.44160335 ])


    # A = np.random.random((m, n)) + 1j * np.random.random((m, n))
    X = np.array(X)
    X = X.astype(np.complex128)

    N = len(X)
    L = N - M - 1
    A = scipy.linalg.hankel(X[:L + 1], X[L:])

    u0, sigma0, v0 = np.linalg.svd(A, full_matrices=True)

    # compute SVD via propack and lapack
    r1 = svdp(A, k, kmax=k*8)

    u1, sigma1, v1, = r1[0:3]
    neig1 = len(sigma1)

    # print the results
    np.set_printoptions(suppress=True, precision=6)

    print('')
    print('Back in Python')

    print('SingVals svd linalg = ', sigma0[:k])
    print('SingVals svdp       = ', sigma1)
    print('SingVal Diffs (linalg vs svdp_aprod    ) = ', sigma0[:k] - sigma1)

    utmp1 = abs(np.dot(u1.conjugate().T,u1))
    utmp1[utmp1<1e-6] = 0.0
    vtmp1 = abs(np.dot(v1, v1.conjugate().T))
    vtmp1[vtmp1<1e-3] = 0.0
    print('  np.clip(abs(np.dot(u1.conjugate().T, u1))) = ')
    print(utmp1)
    print('  np.clip(abs(np.dot(v1, v1.conjugate().T))) = ')
    print(vtmp1)

    # print(u1)
    # print('')
    # print(v1)


    print('')


if __name__ == "__main__":

    test()