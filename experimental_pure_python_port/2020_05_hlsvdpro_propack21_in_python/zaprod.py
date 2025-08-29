#
#     (C) Brian J Soher, April 3,2020
#

import numpy as np
from scipy.sparse.linalg import aslinearoperator


class _AProd(object):
    """Wrapper class for linear operator

    The call signature of the __call__ method matches the callback of
    the PROPACK routines.
    """
    def __init__(self, A):
        try:
            self.A = aslinearoperator(A)
        except TypeError:
            self.A = aslinearoperator(np.asarray(A))

    def __call__(self, transa, m, n, x, y, sparm, iparm):
        if transa == 'n':
            y[:] = self.A.matvec(x)
        else:
            y[:] = self.A.rmatvec(x)

    @property
    def shape(self):
        return self.A.shape

    @property
    def dtype(self):
        try:
            return self.A.dtype
        except AttributeError:
            return self.A.matvec(np.zeros(self.A.shape[1])).dtype

# c
# c     Note: We assume the matrix is always stored in major 
# c     column order.
# c
# 
#       subroutine zaprod(transa,m,n,x,y,zparm,iparm)
# 
#       implicit none
# 
#       character*1 transa
#       integer m,n,iparm(*)
#       complex*16 x(*),y(*),zparm(*)
# 
#       complex*16, pointer :: pA(:,:)  
#       common /csvdp/ pA     
# 
#       external zgemv
#      
#       call zgemv(transa,m,n,dcmplx(1d0,0d0),pA,m,
#      c     x,1,dcmplx(0d0,0d0),y,1)
# 
#       end


