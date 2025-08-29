#
#  (C) Brian J Soher, 2020
#

import numpy as np

from dbsvd import dbdqr
from dbdsdc_cython import dbdsdc


def zritzvec(which, jobu, jobv, m, n, k, dim, D, E, S, U, ldu, V, ldv, kmax):
    """
    DRITZVEC: Compute Ritz-vectors corresponding to the K largest 
              or smallest (depending on the value of WHICH) 
              Ritz-values for A from the Lanczos bidiagonalization 
              A*V_{dim} = U_{dim+1}*B_{dim}.
    
    Parameters:
    
    WHICH: CHARACTER*1. Decides which singular triplets to compute. 
           If WHICH.EQ.'L' then compute triplets corresponding to the K
           largest singular values. 
           If WHICH.EQ.'S' then compute triplets corresponding to the K
           smallest singular values. 
    JOBU: CHARACTER*1. If JOBU.EQ.'Y' then compute the left singular vectors.
          Otherwise the array U is not touched.
    JOBV: CHARACTER*1. If JOBV.EQ.'Y' then compute the right singular 
          vectors. Otherwise the array V is not touched.
    M:    INTEGER. Number of rows of A.
    N:    INTEGER. Number of columns of A.
    K:    INTEGER. Number of desired singular triplets. K <= MIN(DIM,M,N)
    DIM:  INTEGER. Dimension of the Krylov subspace.
    D(DIM): DOUBLE PRECISION array. Contains the diagonal of B.
    E(DIM): DOUBLE PRECISION array. Contains the first sub-diagonal of B.
    S(K): DOUBLE PRECISION array. On return S contains approximation
              to the K largest or smallest (depending on the 
              value of WHICH) singular values of A.
    U(LDU,DIM+1): DOUBLE COMPLEX array. On return the first K columns of U
              will contain approximations to the left singular vectors 
              corresponding to the K largest or smallest (depending on the 
              value of WHICH)  singular values of A.
              On entry the first column of U contains the starting vector
              for the Lanczos bidiagonalization. A random starting vector
              is used if U is zero.
    LDU: INTEGER. Leading dimension of the array U. LDV >= M.
    V(LDV,DIM): DOUBLE PRECISION array. On return the first K columns of V
              will contain approximations to the right singular vectors 
              corresponding to the K largest or smallest (depending on the 
              value of WHICH) singular values of A.
    LDV: INTEGER. Leading dimension of the array V. LDV >= N.
    WORK(LWORK): DOUBLE PRECISION array. Workspace of dimension LWORK.
    IN_LWORK: INTEGER. Dimension of WORK. 
           LWORK should be at least 3*DIM**2 + 
           MAX(3*DIM**2+4*DIM+4, NB*MAX(M,N)), where NB>0 is a block 
           size, which determines how large a fraction of the work in
           setting up the singular vectors is done using fast BLAS-3 
           operations. NB should probably be at least 32 to achieve good
           performance.
    IWORK(8*DIM): INTEGER array. Integer workspace of dimension >= 8*DIM. 
    
    """

    # %----------------------------------------------------------------------%
    # | The bidiagonal SVD is computed in a two-stage procedure:             
    # |                                                   
    # | 1. Compute a QR-factorization M^T*B = [R; 0] of the (k+1)-by-k lower 
    # |    bidiagonal matrix B.
    # | 2. Compute the SVD of the k-by-k upper bidiagonal matrix 
    # |    R = P*S*Q^T. The SVD of B is then (M*P)*S*Q^T.
    # %----------------------------------------------------------------------%

    # %-----------------------------------------%
    # | Set pointers into workspace array
    # %-----------------------------------------%

    # iwrk = 1
    # lwrk = in_lwrk
    # imt  = 1
    # iqt  = imt + (dim + 1) ** 2
    # ip   = iqt + dim ** 2
    # iwrk = ip + dim ** 2
    # lwrk = in_lwrk - iwrk + 1

    # Fortran        Python
    # -------        ------
    # work(imt)   =  workmt
    # work(iqt)   =  workqt
    # work(ip)    =  workp
    # work(iwrk)  =  worki

    # so jold is the input from zlansvd() for the dim param. In zlansvd(),
    # jold = j and j can never be >= lanmax or the 'Krylov subspace dimension'
    # is exceeded and the loop quits. So, max(dim) will always be <= lanmax.
    # I am choosing to make my array dimensions based on dim = lanmax except
    # for iwork, which was explicitly set to 8*kmax in zlansvd().

    lanmax = min(n+1, m+1, kmax)    # as per zlansvd()
    nb = 32                         # number of blocks, as per above, even if not used

    in_lwrk = 3*lanmax**2 + max(3*lanmax**2+4*lanmax+4, nb*max(m,n))     # max needed, as per above

    imt  = 0
    iqt  = imt + (lanmax+1)**2
    ip   = iqt +  lanmax**2
    iwrk = ip  +  lanmax**2
    lwrk = in_lwrk - iwrk+1

    len_workimt   = (lanmax+1)**2
    len_workiqt   = lanmax**2
    len_workip    = lanmax**2
    len_workiwrk  = lwrk
    len_iwork     = 8*kmax
    lzwrk         = m*n

    workimt  = np.zeros((len_workimt,  ), dtype=np.float64)
    workiqt  = np.zeros((len_workiqt,  ), dtype=np.float64)
    workip   = np.zeros((len_workip,   ), dtype=np.float64)
    workiwrk = np.zeros((len_workiwrk, ), dtype=np.float64)

    iwork = np.zeros((len_iwork, ), dtype=np.int)
    zwork = np.zeros((lzwrk,     ), dtype=np.complex128)


    # %-----------------------------------------%
    # | Compute QR-factorization
    # |   B = M * [R; 0]
    # %-----------------------------------------%

    c1 = np.float(0.0)
    c2 = np.float(0.0)

    dbdqr((dim==min(m, n)), jobu, dim, D, E, c1, c2, workimt, dim+1)

    # %-----------------------------------------%
    # | Compute SVD of R:
    # |      R = P * S * Q^T, 
    # | using the Divide-and-conquer SVD
    # %-----------------------------------------%
    
    info = 0
    idd = np.ndarray([1,], dtype='i', order='F')
    ddd = np.ndarray([1,], dtype='d', order='F')

    dbdsdc('u', 'I', dim, D, E, workip, dim, workiqt, dim, ddd, idd, workiwrk, iwork, info)

    # %-----------------------------------------%
    # | Compute left singular vectors for B
    # |    X = P^T * M^T 
    # %-----------------------------------------%

    # dgemm_ovwr('t', dim, dim+1, dim, 1.0, workip, dim, 0.0, workimt, dim+1, workiwrk, lwrk)

    # dgemm_ovwr defined as, compute B <- alpha*op(A)*B + beta*B
    #     but beta=0.0 and op(A) is A' so it's simpler
    #     NB. we have to insert result back into workimt submatrix, then flatten
    workip  = workip.reshape(lanmax, -1).T
    workimt = workimt.reshape(lanmax+1, -1).T
    workimt[:dim,:dim+1] = np.dot(workip[:dim, :dim].T, workimt[:dim, :dim+1])
    workimt = workimt.T.flatten()

    if jobu in ['y','Y']:
        
        # %-----------------------------------------%
        # | Form left Ritz-vectors                  |
        # |   U = U * X^T                           |
        # %-----------------------------------------%
        if which in ['s','S']:
            mstart = dim-k      # dim-k+1
        else:
            mstart = 0          # 1

        # zdgemm_ovwr_left('t', m, k, dim+1, U, ldu, workimt[mstart:], dim+1, zwork, lzwrk)

        # workmt = workimt[mstart:].copy().reshape(lanmax+1, -1)
        # U[:m,:k] = np.dot(U[:m,:dim+1], workmt[:dim+1,:k])  # same as zdgemm_ovrw_left BUT not same as Fortran

        workmt = workimt.copy().reshape(lanmax+1, -1)
        U[:m,:k] = np.dot(U[:m,:dim+1], workmt[:dim+1,mstart:mstart+k])  # same as zdgemm_ovrw_left BUT not same as Fortran


    if jobv in ['y','Y']:

        # %-----------------------------------------%
        # | Form right Ritz-vectors
        # |   V = V * Q
        # %-----------------------------------------%
        
        if which in ['s','S']:
            mstart = dim-k
        else:
            mstart = 0

        # zdgemm_ovwr_left('t', n, k, dim, V, ldv, workiqt[mstart:], dim, zwork, lzwrk)

        # workqt = workiqt[mstart:].copy().reshape(lanmax, -1)
        # V[:m,:k] = np.dot(V[:m,:dim], workqt[:dim,:k])  # same as zdgemm_ovrw_left BUT not same as Fortran

        workqt = workiqt.copy().reshape(lanmax, -1)
        V[:m,:k] = np.dot(V[:m,:dim], workqt[:dim,mstart:mstart+k])  # same as zdgemm_ovrw_left BUT not same as Fortran

