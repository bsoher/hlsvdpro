HLSVDPRO 
======

Overview - Black box fitting of Time Domain Signals
------
The HLSVDPRO package provides code to fit a model function (sum of lorentzians) 
to time-domain data via a 'black box' state space approach (see references below). 
One frequent use for this is by the clinical MRS community for residual water 
removal from MRS signals in the time domain.    

Internally, we use the PROPACK library which contains a set of functions for 
computing the singular value decomposition (SVD) of large and sparse matrices. 
The SVD routines are based on the Lanczos bidiagonalization algorithm with 
partial reorthogonalization (BPRO). 

The calculated singular values and column vectors are subsequently processed 
into lists of parameters that describe the sum of lorentzians that approximate 
the data based on the paper by Laudadio (see below). Parameters are numpy 
arrays of: frequencies, damping factors, amplitudes, and phases. 

**Example:**

```
import hlsvdpro
import numpy as np
import matplotlib.pyplot as plt

data = hlsvdpro.get_testdata()
npts = len(data)
indat = hlsvdpro.TESTDATA   		# this is a built-in dict with test data 
dwell = float(indat['step_size'])
nsv_sought = indat['n_singular_values']

result = hlsvdpro.hlsvd(data, nsv_sought, dwell)

nsv_found, singvals, freq, damp, ampl, phas = result

print("np.allclose(freq, indat['freq0']) = ", np.allclose(freq, np.array(indat['freq0'])) )

fid = hlsvdpro.create_hlsvd_fids(result, npts, dwell, sum_results=True, convert=False)

chop = ((((np.arange(len(fid)) + 1) % 2) * 2) - 1)
dat = data * chop
fit = fid * chop
dat[0] *= 0.5
fit[0] *= 0.5

plt.plot(np.fft.fft(dat).real, color='r') 
plt.plot(np.fft.fft(fit).real, color='b') 
plt.plot(np.fft.fft(dat-fit).real, color='g')
plt.show()

```

HLSVDPRO Methods
------

- `hlsvdpro.hlsvdpro(data, nsv_sought, m=None, sparse=True)` - the main method 
for running the hlsvdpro algorithm. It does not require the dwell time of the 
time domain data, but it also does not convert the results to standard units. It
does allow the user to specify the dimensions of the Hankel matrix, and whether
a sparse SVD is performed or not.

- `hlsvdpro.hlsvd(data, nsv_sought, dwell_time)` - provides backwards 
compatibility to the API for HLSVDPRO version 1.x. It calls the hlsvdpro() method
with default values corresponding to the algorithm used in version 1.x. See 
docstring for more information on the default values used.

HLSVDPRO Utility Methods
------
- `hlsvdpro.create_hlsvd_fids(result, npts, dwell, sum_results=False, convert=True)` - 
can be used to create FIDs from the results tuple from either the `hlsvd()` 
or the `hlsvdpro()` methods. It can return either individual FIDs or a sum of 
all FIDs as a result.  

- `hlsvdpro.convert_hlsvd_result(result, dwell)` - uses the dwell time to 
convert the `hlsvdpro()` result tuple to more standard units. Frequencies 
convert to [kHz], and damping factors to [ms]. Phases convert to [degrees]. 
Singular values, amplitudes and row and column matrices are maintained at 
their same values and output tuple locations. Note - the `hlsvd()` method 
automatically calls this internally, so you don't have to convert values
if you use that method.

- `hlsvdpro.get_testdata()` - returns a numpy array of 1024 complex data 
points that represents a real world short TE single voxel PRESS data set.
This function converts the base64 encoded string saved in the TESTDATA dict
into a numpy array for you. Additional information about the data and the 
known values for fitting it via the hlsvd() method can be retrieved from 
the TESTDATA dict.  See 'Example' for more usage information.


Additional Functionality
------

As of version 2.0.0, HLSVDPRO also provides access to internal PROPACK SVD 
functions via the 'propack' module. These have been compiled from PROPACK 
version 2.1, and include:

- `hlsvdpro.propack.lansvd()` - matrix-vector product via user defined Python 
callback 
vector/matrix multiplication
- `hlsvdpro.propack.lansvd_irl()` - matrix-vector product via user defined Python 
callback 
- `hlsvdpro.propack.lansvd_aprod()` - matrix-vector product via fixed internal 
Fortran call 
- `hlsvdpro.propack.lansvd_irl_aprod()` - matrix-vector product via fixed internal 
Fortran call

Technical Overview and References
------

The PROPACK SVD Fortran code was written by Rasmus Munk Larsen. Teresa Laudadio 
and Diana Sima adapted it into HLSVDPRO algorithm. Brian J Soher adapted the 
algorithm for PROPACK v2.1 use and the python wrappers used to access the code.

For complete copyright and license information, see the LICENSE file.

The Python wrapper is compatible with Python 2.7 and 3.x. The binaries
that the Python wrapper calls are 64-bit, so this code is not compatible with
32-bit Python. 


The state space approach is described in: _S.Y. Kung, K.S. Arun and D.V. Bhaskar
Rao, J. Opt. Soc. Am. 73, 1799 (1983)._

HLSVDPRO version 1.0.x made use of code from PROPACK version 1.x and was 
implemented based on the paper by: _W.W.F. Pijnappel, A. van den Boogaart, R. de 
Beer, D. van Ormondt, J. Magn. Reson. 97, 122 (1992)._

HLSVDPRO version 2.0.0 was adaptated to use PROPACK library version 2.1 to 
implement the HLSVDPRO algorithm as described in: _T. Laudadio, N. Mastronardi
L. Vanhamme, P. Van Hecke and S. Van Huffel, "Improved Lanczos algorithms for 
blackbox MRS data quantitation", Journal of Magnetic Resonance, Volume 157, 
pages 292-297, 2002._ 

The Python wrapper used to access the PROPACK library was adapted from the 
pypropack repository create by Jake van der Plas on Github: 
https://github.com/jakevdp/pypropack with extensions by Brian J Soher to 
include the matrix-vector product calculation within the Fortran code.

Acknowlegements and Kudos
-------

We (the HLSVDPRO team) have made some minor modifications to this version of
HLSVDPRO. Most significantly, we have updated the PROPACK library to version 
2.1, added Makefiles and Python wrappers for both the HLSVDPRO and PROPACK 
methods. We are very grateful to all of the above scientists for sharing their
work and agreeing to release it under a BSD license.
