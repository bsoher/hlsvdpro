April 21, 2020 - Brian J Soher

So, it's been a slog trying to convert PROPACK v2.1 to a 'pure' Python 
version of HLSVDPRO.  Note, as of a couple of years ago based on the Laudadio
paper JMR 157 2002, we migrated to using her workflow, rather than the one
originally set up by the Pijnappfel paper lo- these many years ago. The new
workflow skips any and all FFT calls.  Yay! And so we just need PROPACK's 
zlansvd.f call to return U, Sigma, V (mainly Sigma) SVD values to us. 

SO, the motivation is to get zlansvd() callable from Python, both as a py2f
call to the original Fortran library, and as a pure Python port. We use the
really nice Python call svdp() from Jake Vanderplas' Pypropack Github 
repository as our front end to both.

Step 1 was to compile pypropack for Windows.  Did that (oh ... therein lies a 
story). Bottom line, direct call to py2f from command line creates a single
*.pyd file. But, building the library using the python setup.py ... call 
at the command line creates a stup *.pyd file and some mangled *.dlls that
the *.pyd needs to run. After a bit of research, I think that this has been
a somewhat regular workflow in the numpy/scipy build to try to avoid .DLL/.so
collisions as various packages might contain multiple copies of the same
helper library that could affect one but not the other if there are slight
differences within. End result was a test function that calls the 
pypropack.svdp() function and returns a result.  I set that up against the 
standard hlsvdpro.hlsvd() PyPI version for comparison.

Step 2 was to complete the PROPACK v1.x port. I found a method for using the
lapack_cython library bindings to access the one (?) method I needed, which
was dbdsqr(). Maybe I did dgemm() too, and I know I did dbdsdc() for the
v2.1 pure Python port.  It was cool finding out how to use ctypes to do that.
Those modules don't have much/any error checking.  Use with care. Added 
those to Philip's fine original work, and voila!  Something sort of finished
even though I wanted v2.1 now.  I set that up agains the other two fortran
library call versions and tested against three data sets for 100 calls using
cProfile.

Here's copy of the April 17, 2020 speed runs just in case the original notes
for v1.x get lost.

  FILE (nsv_sought)    LOCAL(native py v1)   PROPACK               HLSVD
                        svDif (nsv) (s/100)   svDif (nsv) (s/100)   (nsv) (s/100)
  press_cp0 (20)     -  ~10-7 (20) (33.732)  ~10-10  (18) (0.815)    (18) (2.165)
  press_cp5 (30)     -  ~10-7 (30) (33.075)  ~10-10  (15) (1.578)    (15) (2.223)
  laser_2010102901 (35) !10-7 (35) (33.272)  ~10-14  (21) (1.525)    (21) (2.292)
  
Step 3 was another epic tale. Begun in earnest right at the beginning of the
pandemic of 2020.  Woo hoo!  I did 6 takes. Tried to splice in a lot of the
code from v1.x that was similar. In the end (around version 4/5/6) I broke 
down and wrote a 'debug' version of the fortran lib, locked down an example
data matrix (and 'random' start values) for both Python and fortran and did a
side by side tweak fest.  Ahhhhhhhhrrrrggg!  But in the end, I've got a pure
port of v2.1 zlansvd.f that works about 98% the same. The SVD values all seem
the same, but the U and V arrays sometimes go south.  I know that the machine
precision values are different a bit between Python and Fortran ... ? But, 
otherwise, I'm stuck without another 2 weeks to devote.  

So I set the native Python version, native_pypropack.svdp(), up against the
np.linalg.svd(A,full_matrices=True) (sort of a cheat, yeah, but ...) and also
the fortran_pypropack.svdpf_aprod() call (built in Aprod funct) but otherwise
just a renamed module from step 2 to differentiate from native_pypropack. 

Here's my speed run for 100 reps, nsv_sought(20-30), hankel_size (512), M (256) 
etc. values as much the same as the April 17 run as I could make them.

  FILE (nsv_sought)    NATIVE (native py v2) PROPACK               NP.LINALG.SVD(full)
                        svDif (nsv) (s/100)   svDif (nsv) (s/100)   (nsv) (s/100)
  press_cp0 (20)     -  ~10-9 (20)  (5.906)  ~10-9  (20) (1.291)    (20) (2.027)
  press_cp0 (30)     -  ~10-9 (30) (10.935)  ~10-9  (30) (1.982)    (30) (1.918)
  press_cp5 (20)     -  ~10-9 (20) (10.669)  ~10-9  (20) (1.845)    (20) (1.898)
  press_cp5 (30)     -  ~10-9 (30) (13.468)  ~10-9  (30) (2.224)    (30) (1.909)
  laser_2010102901 (20) ~10-12 (20) (4.132)  ~10-12 (20) (1.297)    (20) (1.915)
  laser_2010102901 (35) ~10-12 (35) (9.782)  ~10-12 (35) (1.896)    (35) (1.894)
  
  1. I played a bit with the hankel_matrix size (to 384) and kmax (from 200 to 150)
     and the press_cp0 (30) times did not budge. But the SingVals were the same.
  
===========================================================
April 17, 2020  -  Brian J Soher  added these for provenance

This is still notes for converting PROPACK v1.x (?) since I was able to 
use lapack_cython for some of the methods not previously available.

Version3 of this project has been tested. It uses scipy's lapack_cython to 
access the dbdsqr() module. There is very little error checking in the 
dbdsqr_cython.dbdsqr() module. Mainly, should check that the incoming arrays
are numpy ndarrays. 

Found that USE_FORTRAN_RANDOMS True option slows it way down since it has to 
access a file to get values.  I think that using the scipy.random module 
works but occasionally gives glitches. It also is not possible to run the same
code in Python and then in Fortran (via print *,'' statements) to check against
each other if we always have random code. The solution is to use the get-from-
file method, because it always uses the same random numbers.  NB. in v2.1
I use zlanrnv() in both Python and Fortran and they always give the same
random sequence. So I can do a side by side there more easily.

In test.py results for hlsvdpro_local (native python) and hlsvdpro_propack 
(fortran of v2.1 code) are compared to gold standard hlsvdpro.hlsvd() (fortran
v1.x code). In all cases, Hankel matrix size was len(observed) // 8. So that 
was data size = 512 for press_cp0 and _cp5. Not sure for laser. M was 256 
explicitly for hlsvdpro_propack() call. M defaults to len(observed) // 2 for
hlsvdpro_local, so that should have been 256 also.  All in all, the comparison
between methods is more important than the actual fairness, I think for now.

Here are the results for a few data sets:
  FILE (nsv_sought)    LOCAL                 PROPACK               HLSVD
                        svDif (nsv) (s/100)   svDif (nsv) (s/100)   (nsv) (s/100)
  press_cp0 (20)     -  ~10-7 (20) (33.732)  ~10-10  (18) (0.815)    (18) (2.165)
  press_cp5 (30)     -  ~10-7 (30) (33.075)  ~10-10  (15) (1.578)    (15) (2.223)
  laser_2010102901 (35) !10-7 (35) (33.272)  ~10-14  (21) (1.525)    (21) (2.292)
  
  Notes. 
   1) 25215695.xml had too few pts (128) in FID for current setup of hankel size = 256
   2) s/100 for SVD call only if possible, lanczopw(), svdp_aprod(), hlsvd(), respectively
   3) LOCAL and HLSVD had hankel size of 512, PROPACK was 256. When I dropped LOCAL and 
       HLSVD down to 256, LOCAL was 5% faster, HLSVD was 30% faster.  (FYI, PROPACK was at
       hankel size 128 and was also ~35% faster) not sure about accuracy for any of them.

========================================================
These are Philip's notes from way-back - he was converting PROPACK v1x (?) to pure python

This is the latest draft of a port of HLSVDPRO to Python. In keeping with my intention to "port it like an onion", I started by removing the first two layers. http://scion.duhs.duke.edu/vespa/analysis/wiki/HlsvdPorting

All Fortran files have been ported with two exceptions (see below). 

== Caveats ==

1) There's a small block of code at the end of zlansvdw.f that I didn't port. It's only activated if the input param jobv == 'Y' and it's hardcoded to 'N' in this code. It's just a single call to zgemmina() so if I debug that properly, porting this last call shouldn't be difficult.

2) There's no numpy entry point for the LAPACK function dbdsqr() and it's too complex to be worth re-writing in Python. I'm currently using the version inside the HLSVDPRO library but obviously any binary dependency defeats the main purpose of porting to Python. I'm still trying to figure out how to deal with this.

http://www.netlib.org/lapack/double/dbdsqr.f

== Results ==

Results match the Fortran code for all of the press_cpX files. laser_2010102901.xml goes badly. There's still some bugs in the Python to be worked out, apparently.


== Code Organization ==

My Python code is at the top level. The directory 'original' contains source and test files copied from vespa/hlsvd. I've modified the Fortran source to provide a new entry points, a bug fix (see below), and a bazillion print statements to check equivalence with my Python code.

== Performance ==
Performance is bad; I'm not worried about it yet. 

== Quality ==
My code is completely unoptimized, wildly messy by any standard and overrun with FIXMEs.



== Next Steps ==

My next steps are to  --
- Fix the bugs that cause laser_2010102901.xml to fail
- Improve performance which is quite bad
- Figure out how  to eliminate the dependency on dbdsqr()



===========================================================
maybe something promising ... ?

March 20, 2020 - Brian J Soher
Found at: https://gist.github.com/insertinterestingnamehere/b50390fc720e1e554af0


call_dgemm.py
----------------------------
# This shows how to use SciPy's Cython-exposure of the BLAS and LAPACK libraries from within ctypes

Here's the doc entry on BLAS DGEMM

 NAME
      DGEMM - perform one of the matrix-matrix operations   C :=
      alpha*op( A )*op( B ) + beta*C,

 SYNOPSIS
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA,
                       B, LDB, BETA, C, LDC )

          CHARACTER*1  TRANSA, TRANSB
          INTEGER      M, N, K, LDA, LDB, LDC
          DOUBLE       PRECISION ALPHA, BETA
          DOUBLE       PRECISION A( LDA, * ), B( LDB, * ), C( LDC,* )

 PURPOSE
      DGEMM  performs one of the matrix-matrix operations

      where  op( X ) is one of

         op( X ) = X   or   op( X ) = X',

      alpha and beta are scalars, and A, B and C are matrices,
      with op( A ) an m by k matrix,  op( B )  a  k by n matrix
      and  C an m by n matrix.

--------------------------------
import numpy as np
import ctypes as ct
from scipy.linalg import cython_blas

#ct.cdll.LoadLibrary("libopenblas.dll")
#openblas = ct.CDLL("libopenblas.dll")
#dgemm = openblas.dgemm

ct.pythonapi.PyCapsule_GetPointer.restype = ct.c_void_p
ct.pythonapi.PyCapsule_GetPointer.argtypes = [ct.py_object, ct.c_char_p]
ct.pythonapi.PyCapsule_GetName.restype = ct.c_char_p
ct.pythonapi.PyCapsule_GetName.argtypes = [ct.py_object]
ptr_type = ct.CFUNCTYPE(ct.c_void_p, ct.c_char_p, ct.c_char_p,
                        ct.POINTER(ct.c_int), ct.POINTER(ct.c_int),
                        ct.POINTER(ct.c_int), ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_double), ct.POINTER(ct.c_int),
                        ct.POINTER(ct.c_double), ct.POINTER(ct.c_int),
                        ct.POINTER(ct.c_double), ct.POINTER(ct.c_double),
                        ct.POINTER(ct.c_int))
dgemm2 = ptr_type(
             ct.pythonapi.PyCapsule_GetPointer(
                 cython_blas.__pyx_capi__['dgemm'],
                 ct.pythonapi.PyCapsule_GetName(
                     cython_blas.__pyx_capi__['dgemm'])))

a = np.array([[1,2],[3,4]], 'd', order='F')
b = np.array([[5,6],[7,8]], 'd', order='F')
c = np.empty((2,2), order='F')

transa = ct.c_char('N')
transb = ct.c_char('N')
alpha = ct.c_double(1.)
beta = ct.c_double(0.)
lda = ct.c_int(2)
ldb = ct.c_int(2)
ldc = ct.c_int(2)
m = ct.c_int(2)
n = ct.c_int(2)
k = ct.c_int(2)
args = (ct.byref(transa), ct.byref(transb), ct.byref(m), ct.byref(n),
        ct.byref(k), ct.byref(alpha),
        a.ctypes.data_as(ct.POINTER(ct.c_double)), ct.byref(lda),
        b.ctypes.data_as(ct.POINTER(ct.c_double)), ct.byref(ldb),
        ct.byref(beta), c.ctypes.data_as(ct.POINTER(ct.c_double)),
        ct.byref(ldc))
#dgemm(*args)
#print c
c[:] = 0
dgemm2(*args)
print c
print a.dot(b)

#--------------------------------------------
Here is a C example for using DGEMM
#--------------------------------------------

/* C source code is found in dgemm_example.c */

#define min(x,y) (((x) < (y)) ? (x) : (y))

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

int main()
{
    double *A, *B, *C;
    int m, n, k, i, j;
    double alpha, beta;

    printf ("\n This example computes real matrix C=alpha*A*B+beta*C using \n"
            " Intel(R) MKL function dgemm, where A, B, and  C are matrices and \n"
            " alpha and beta are double precision scalars\n\n");

    m = 2000, k = 200, n = 1000;
    printf (" Initializing data for matrix multiplication C=A*B for matrix \n"
            " A(%ix%i) and matrix B(%ix%i)\n\n", m, k, k, n);
    alpha = 1.0; beta = 0.0;

    printf (" Allocating memory for matrices aligned on 64-byte boundary for better \n"
            " performance \n\n");
    A = (double *)mkl_malloc( m*k*sizeof( double ), 64 );
    B = (double *)mkl_malloc( k*n*sizeof( double ), 64 );
    C = (double *)mkl_malloc( m*n*sizeof( double ), 64 );
    if (A == NULL || B == NULL || C == NULL) {
      printf( "\n ERROR: Can't allocate memory for matrices. Aborting... \n\n");
      mkl_free(A);
      mkl_free(B);
      mkl_free(C);
      return 1;
    }

    printf (" Intializing matrix data \n\n");
    for (i = 0; i < (m*k); i++) {
        A[i] = (double)(i+1);
    }

    for (i = 0; i < (k*n); i++) {
        B[i] = (double)(-i-1);
    }

    for (i = 0; i < (m*n); i++) {
        C[i] = 0.0;
    }

    printf (" Computing matrix product using Intel(R) MKL dgemm function via CBLAS interface \n\n");
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 
                m, n, k, alpha, A, k, B, n, beta, C, n);
    printf ("\n Computations completed.\n\n");

    printf (" Top left corner of matrix A: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(k,6); j++) {
        printf ("%12.0f", A[j+i*k]);
      }
      printf ("\n");
    }

    printf ("\n Top left corner of matrix B: \n");
    for (i=0; i<min(k,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.0f", B[j+i*n]);
      }
      printf ("\n");
    }
    
    printf ("\n Top left corner of matrix C: \n");
    for (i=0; i<min(m,6); i++) {
      for (j=0; j<min(n,6); j++) {
        printf ("%12.5G", C[j+i*n]);
      }
      printf ("\n");
    }

    printf ("\n Deallocating memory \n\n");
    mkl_free(A);
    mkl_free(B);
    mkl_free(C);

    printf (" Example completed. \n\n");
    return 0;
}



