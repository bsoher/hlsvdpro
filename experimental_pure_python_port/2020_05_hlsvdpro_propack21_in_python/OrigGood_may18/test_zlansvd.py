#
#     (C) Brian J Soher, 2020
#

from __future__ import division, print_function, absolute_import

__all__ = ['svdp']
import sys
import xml.etree.cElementTree as ElementTree
import cProfile
import os
import pprint

pp = pprint.pprint

# 3rd party modules
import numpy as np
import scipy
from scipy.sparse.linalg import aslinearoperator

import zlansvd

# Vespa modules
import vespa.common.util.xml_ as util_xml



class _AProd(object):
    """
    Wrapper class for linear operator

    The call signature of the __call__ method matches the callback of
    the PROPACK routines.
    """
    def __init__(self, A):
        try:
            self.A = aslinearoperator(A)
        except TypeError:
            self.A = aslinearoperator(np.asarray(A))

    def __call__(self, transa, m, n, x, y, parm):
        if transa == 'n':
            y[:m] = self.A.matvec(x[:n])
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



def svdp(A, k, which='L', irl_mode=False, kmax=None,
         compute_u=True, compute_v=True, v0=None, full_output=False, tol=0,
         delta=None, eta=None, anorm=0, cgs=False, elr=True, blocksize=1,
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
    which : string (optional)
        which singluar triplets to compute:
        - 'L': compute triplets corresponding to the k largest singular values
        - 'S': compute triplets corresponding to the k smallest singular values
        which='S' requires irl=True.
        Computes largest singular values by default.
    irl_mode : boolean (optional)
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
    cgs : boolean (optional)
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
    if v0 is None:
# bjs debug
        u[:, 0] = np.random.random(m)
        if np.iscomplexobj(np.zeros(0, dtype=typ)):  # complex type
            u[:, 0] += 1j * np.random.random(m)
#        u[:,0] = np.array([0.25662954773130986+0.57062648579991326*1j, 0.20366349962301655+0.17718181865061933*1j, 0.13613447425116998+0.34298286063336469*1j, 0.62030760743804181+0.60769243986903876*1j, 0.31760778194161077+0.70415718758479073*1j, 0.29932099778923460+0.66607309269575221*1j])

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

    u, sigma, v, info, bnd = lansvd( jobu, jobv, m, n, k, kmax, aprod, u, v, tol, doption, ioption, parm)

    neig = len(sigma)

#    return neig, sigma, u[:, :neig], v[:, :neig].conj().T, info

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

def test(flavor, filename, profiling=False):

    if filename:
        s = open(filename, "rb").read()
        root = ElementTree.fromstring(s)

        # The XML is divided into two sections, one for input and the other
        # for (expected) output.
        input_ = util_xml.element_to_dict(root.find("input"))

        if filename == 'data_press_cp0.xml':
            sing0 = np.array([73010.55, 19342.88, 11212.103, 10580.863, 9682.249, 7823.373, 5280.5444, 4496.4272, 3146.9573, 2632.8362, 2087.5728, 1365.9307, 1240.6813, 1133.4723, 1042.3507, 821.9423, 774.10333, 760.31476, 673.0105, 648.33887, 623.1397, 612.5165, 593.5787, 575.90576, 563.35913, 555.03613, 530.4214, 521.7685, 515.6749, 506.13913, 493.8507, 488.3877, 481.42786, 478.7373, 465.66724])
        elif filename == 'data_press_cp5.xml':
            sing0 = np.array([51083.14, 13463.821, 7659.858, 7398.313, 5639.7817, 3891.2114, 1964.6312, 1787.428, 1681.4921, 1416.3357,1140.9312, 987.8655, 884.4408, 760.91614, 740.42426, 692.06335, 689.62585, 668.5781, 623.6184, 584.23596, 579.9215, 570.84375, 565.2316, 560.9096, 554.43567, 548.17957, 532.3629, 528.10956, 516.4973, 507.73767, 504.09927, 496.571, 477.55408, 476.50812, 474.7325, 465.68942, 464.37546, 455.55417, 453.48428, 446.487, 434.9923, 429.20972, 425.54483, 417.49387, 415.32175])
        elif filename == 'data_laser_2010102901.xml':
            sing0 = np.array([275421.64025479293, 16424.18555182939, 2476.5730101379118, 1712.2413005928718, 824.6613384777334, 779.6706625817577, 496.38987966618686, 398.1442243749407, 326.94657869762193, 297.7971078989293, 265.1014109166263, 222.62169658905572, 217.03976473483394, 199.1515093280466, 190.50515431198957, 145.54978846531134, 133.00454060410834, 125.3174060365698, 121.64224345177496, 104.73477901897321])

        signals = input_["signals"]
        A       = signals.copy().astype(np.complex128)
        k       = 30  
        kmax    = k*5
        M       = 256   # empirical
        
        hankel_size = len(A) // 8                   # or hard fix to 512
        if hankel_size > 1024: hankel_size = 1024

        A = A[0:hankel_size].copy()
        N = len(A)
        L = N - M - 1
        X = scipy.linalg.hankel(A[:L + 1], A[L:])

    else:

        m = 200
        n = 100
        k = 40
        kmax = 100
        M = 12
        A = np.random.random((m, n)) + 1j * np.random.random((m, n))
        
        A = np.array(A)
        X= A.astype(np.complex128)


#        # Create a random matrix 
#        if flavor == 's':
#            A = np.random.random((m, n)).astype(np.float32)
#            flavorlabel = 'SINGLE'
#        elif flavor =='d':
#            A = np.random.random((m, n)).astype(np.float64)
#            flavorlabel = 'DOUBLE'
#        elif flavor =='c':
#            A = np.random.random((m, n)) + 1j*np.random.random((m, n))
#            A = A.astype(np.complex64)
#            flavorlabel = 'COMPLEX8'
#        elif flavor =='z':
#            A = np.random.random((m, n)) + 1j*np.random.random((m, n))
#            A = A.astype(np.complex128)
#            flavorlabel = 'COMPLEX16'

#        N = len(A)
#        L = N - M - 1
#        X = scipy.linalg.hankel(A[:L + 1], A[L:])


        # recall abs(np.dot(U0.conjugate().T, U0)) = eye
        #  and   abs(np.dot(V0, V0.conjugate().T)) = eye for complex
        # recall abs(np.dot(U0.T, U0)) = eye
        #  and   abs(np.dot(V0, V0.T)) = eye  for double



    if not profiling:
        reps = 1
    else:
        reps = 50

    print('X.shape = '+str(X.shape))

    for i in range(reps):
        u0, sigma0, v0 = np.linalg.svd(X, full_matrices=True)
        neig0 = len(sigma0)

    reps1 = [0,0,0]

    for i in range(reps):
        result1 = svdp(X, k, kmax=kmax, cgs=True)
        u1, sigma1, v1 = result1[0:3]
        neig1 = len(sigma1)
        if neig1 == k:
            imax = min(neig1, neig0)
            if (max(sigma1[:imax]-sigma0[:imax]) > 0.1):
                reps1[1] += 1
            else:
                reps1[0] += 1
        else:
            reps1[2] += 1

    # print the results
    np.set_printoptions(suppress=True, precision=4, linewidth=120)

    utmp1 = abs(np.dot(u1.conjugate().T,u1))
    utmp1[utmp1<1e-6] = 0.0
    vtmp1 = abs(np.dot(v1, v1.conjugate().T))
    vtmp1[vtmp1<1e-3] = 0.0

    if neig1 == k:
        kmin = min(neig1, neig0)
        print(' Results hlsvdpro_v21, max(abs(sigma1-sigma0)) = ', neig1, max(abs(sigma1[:kmin] - sigma0[:kmin])))
        print(' Results hlsvdpro_v21, good_fit, wrong_fit, bad_fit = ', reps1[0], reps1[1], reps1[2])
        print('sigma1 = ', sigma1)
        print('sigma0 = ', sigma0)
        print('  np.clip(abs(np.dot(u1.conjugate().T, u1))) = ')
        print(utmp1)
        print('  np.clip(abs(np.dot(v1, v1.conjugate().T))) = ')
        print(vtmp1)
    else:
        print(' Results hlsvdpro_v21, max(abs(sigma1-sigma0)) =  null')
        print(' Results hlsvdpro_v21, good_fit, wrong_fit, bad_fit = ',reps1[0], reps1[1], reps1[2])
        print('sigma1 = ', sigma1)
        print('  np.clip(abs(np.dot(u1.conjugate().T, u1))) = ')
        print(utmp1)
        print('  np.clip(abs(np.dot(v1, v1.conjugate().T))) = ')
        print(vtmp1)

    print('')


def test_profiling():
    for flavor in ['z']:  #['s','d','c','z']:
        test(flavor, profiling=True)

if __name__ == "__main__":

    filename = ''
    #filename = "testdata/data_25215695.xml"
    #filename = "testdata/data_press_cp0.xml"
    #filename = "testdata/data_press_cp5.xml"
    #filename = "testdata/data_laser_2010102901.xml"

    if os.path.exists("profile.data"):
        os.remove("profile.data")

    flavor = 'z'
    
    test(flavor, filename)

    # cProfile.run('test_profiling(filename)', 'profile.data')
    # import pstats as ps
    # p = ps.Stats('profile.data')
    # p.strip_dirs().sort_stats('cumulative').print_stats()
