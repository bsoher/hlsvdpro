# Python modules
from __future__ import division
import ctypes
import math
import os
import pprint
import re
pp = pprint.pprint

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.util.misc as util_misc

U_REGEX = "U[(]\s+(\d{1,4})\s+(\d{1,4})\s*[)][:]\s*[(]\s*(.*),\s*(.*)\s*[)]"
U_REGEX = re.compile(U_REGEX)

KMAX = 50
MAX_SINGULAR_VALUES = KMAX
NDPMAX = 8192

# _LIBRARY_NAMES maps platforms to the name of the library on said platform.
_LIBRARY_NAMES = { "osx"       : "libhlsvd-%d.dylib",
                   "windows"   : "hlsvd-%d.dll",
                   "linux"     : "libhlsvd-%d.so",
                 }

def lanczopw(signals, lrow, mcol, kuser):
# Call -- 
    # lanczopw(signals, len(signals), lrow, mcol, kuser, kmax, lsinval, u, v, 
    #          work, lwork, zwork, lzwork, zwork[ilambda], zwork[itrlambda]
    #          zwork[ifvect], ndiv)

    # Output: lsinval, uuu, ndiv (kfit == nsv_found)

# Signature --
    #   subroutine lanczopw(signal,ndp,m,n,k,kmax,sigma,U,V,
    #  c     work,lwrk,zwork,lzwrk,lambda,trlambda,fvect,info)
    #   implicit none

    #   double precision zero, one
    #   double complex zeroc
    #   parameter(zero = 0.0d0, one = 1.0d0, zeroc=(0.0d0,0.0d0))
    #   integer i,m,n,k,kmax,ioption(10),iwork(2*kmax+1),info,lwrk,lzwrk
    #   integer ind1,ndp
    #   integer*8 planF,planB
    #   double precision work(lwrk)
    #   double precision doption(10)
    #   double precision tol,pi,bnd(kmax)
    #   double precision dnrm2,sigma(kmax)
    #   complex*16 lambda(ndp),trlambda(ndp),
    #  c          U(m,kmax+1),V(n,kmax),zwork(lzwrk),
    #  c           fvect(ndp),signal(*)

    libhlsvd = _load_the_library()


    # At this point libhlsvd should be a valid reference to the hlsvd library.
    # I test that assertion here as best as I can by generating a throwaway
    # reference to the function I'm going to call. 
    f = libhlsvd.lanczopw_python_

    # OK, all is well. I create and populate the variables I need as parameters
    # to the function.
    n_data_points = len(signals)

#signal,ndp,m,n,k,kmax,sigma,U,V,
    #  c     work,lwrk,zwork,lzwrk,lambda,trlambda,fvect,info


    signals = signals.tolist()

    # There's a couple of odd loops here that set up fvect. I think the 'f'
    # might stand for 'fased' (phased)? These loops basically shove all the
    # points to the right by ndp/2.
    fvect = np.zeros( (n_data_points, ), np.complex128)

    for i in range(lrow):
        j = i + mcol - 1
        #print "signals[{}] = {:.17E}".format(j, signals[j])
        #fvect.append(signals[j])
        fvect[i] = signals[j]
        #print "fvect[{}] = {:.17E}".format(i, fvect[i])

    # for i, z in enumerate(fvect):
    #     print "fvect[{}] = {:.17E}".format(i + 1, z)

    # fortran: 
    # first i =         2047
    # ind1 =         4096 i =         2047
    # ind1 =         4095 i =         2046
    # ind1 =         4094 i =         2045
    # ind1 =         2052 i =            3
    # ind1 =         2051 i =            2
    # ind1 =         2050 i =            1
    # exited loop, ppoints =         2047 i =            0


    j = n_data_points - 1
    ppoints = 0
    # print "first i = {}".format(mcol - 1)
    for i in range(mcol - 2, -1, -1):
        # print "copying signals[{}] = {} to fvect[{}]".format(i + 1, signals[i], j + 1)
        fvect[j] = signals[i]
        j -= 1
        ppoints += 1

    # print "ppoints = {}, last i = {}".format(ppoints, i)

    # for i, z in enumerate(fvect):
    #     print "fvect[{}] = {:.17E}".format(i + 1, z)

    lam = np.fft.fft(fvect)

    for i, z in enumerate(lam):
        print "lambda[{}] = {}".format(i, z)




    import matplotlib.pyplot as pyplot


    pyplot.plot(lam, color='r')
    pyplot.show()


    kjdfjkd


    InputDoubleArrayType = ctypes.c_double * n_data_points
    OutputDoubleArrayType = ctypes.c_double * MAX_SINGULAR_VALUES

    UuuArrayType = ctypes.c_double * (lrow * (KMAX + 1))

    # In-line comments refer to the names of the corresponding variables
    # in hlsvdpro.f.

    # Input params
    real_signals = InputDoubleArrayType()
    imaginary_signals = InputDoubleArrayType()


    
    for i, signal in enumerate(signals):
        real_signals[i] = signal.real                   # signal_r
        imaginary_signals[i] = signal.imag              # signal_i

    nsv_sought = kuser

    lwrk = lrow + mcol + (13*KMAX) + (8*KMAX**2) + (32 * lrow) + (n_data_points*KMAX)
    lwrk = ctypes.c_long(lwrk)

    lzwrk = lrow + mcol + (32*lrow) + (7*KMAX) + 2 + (2*KMAX**2) + (5 * n_data_points)
    lzwrk = ctypes.c_long(lzwrk)

    n_data_points = ctypes.c_long(n_data_points)        # ndp

    kmax = ctypes.c_long(KMAX)                          # kmax

    lsinval = OutputDoubleArrayType()

    uuu_r = UuuArrayType()
    uuu_i = UuuArrayType()

    # for i, uuu_item in enumerate(uuu):
    #     uuu_r[i] = uuu_item.real
    #     uuu_i[i] = uuu_item.imag

    lrow = ctypes.c_long(lrow)                          # Lrow/m
    mcol = ctypes.c_long(mcol)                          # mcoL/n
    
    # # Output params
    # frequencies = OutputDoubleArrayType()               # freq
    # damping_factors = OutputDoubleArrayType()           # damp
    # amplitudes = OutputDoubleArrayType()                # ampl
    # phases = OutputDoubleArrayType()                    # fase
    # singular_values = OutputDoubleArrayType()           # Lsinval

     #  subroutine lanczopw_python(signal_r,signal_i,ndp,m,n,k,kmax,sigma,
     # c     uuu_r, uuu_i,
     # c     lwrk,lzwrk,info)
    
    done = False
    while not done:
        nsv_sought = ctypes.c_long(nsv_sought)              # kuser/k/kfit
        info = ctypes.c_long()                              # info (0, < 0, > 0)
                                                            # 0 is good!

        # Note - nsv_sought == kuser going in. lanczopw() changes this to
        # the number of singular values it found (kfit).
        libhlsvd.lanczopw_python_(real_signals,
                                  imaginary_signals,
                                  ctypes.pointer(n_data_points),
                                  ctypes.pointer(lrow), 
                                  ctypes.pointer(mcol), 
                                  ctypes.pointer(nsv_sought),
                                  ctypes.pointer(kmax),
                                  ctypes.pointer(lsinval),
                                  ctypes.pointer(uuu_r),
                                  ctypes.pointer(uuu_i),
                                  ctypes.pointer(lwrk), 
                                  ctypes.pointer(lzwrk), 
                                  ctypes.pointer(info),
                                 )

        nsv_sought = nsv_sought.value
        info = info.value

        if (info == -1):
            if nsv_sought > 0:
                # Need to call lanczopw() again.
                print "Python: calling lanczopw again with kuser = {}".format(nsv_sought)
                pass
            else:
                # FIXME - raise an error; PROPACK didn't converge.
                done = True
        else:
            done = True

    nsv_found = nsv_sought

    #print "Python land: nsv_sought = {}".format(nsv_sought)

    mcol = mcol.value
    lrow = lrow.value

    # Turn these into ordinary Python lists
    uuu_r = [value for value in uuu_r]
    uuu_i = [value for value in uuu_i]

    # lines = open("good_u.txt").read().split('\n')
    # lines = [line.strip() for line in lines if line.strip()]

    # d = { }

    # for line in lines:
    #     m = U_REGEX.match(line)
    #     i, j, real, imag = m.groups()
    #     i = int(i)
    #     j = int(j)
    #     real = float(real)
    #     imag = float(imag)
    #     d[ (i, j) ] = (real, imag)
    #     #print "U({},{}): ({:.17E},{:.17E})".format(i, j, real, imag)    

    # import vespa.common.util.math_ as util_math
    # for j in range(lrow):
    #     for i in range(KMAX + 1):
    #         k = (i * lrow) + j
    #         #if uuu_r[k] or uuu_i[k]:

    #         real, imag = d[ (j + 1, i + 1) ]
    #         try:
    #             assert(util_math.eq(real, uuu_r[k]))
    #             assert(util_math.eq(imag, uuu_i[k]))
    #         except AssertionError:
    #             print "AssertionError, k = {}, i = {}, j = {}".format(k, i ,j )
    #             print real, imag, uuu_r[k], uuu_i[k]
    #         except IndexError:
    #             print "IndexError, k = {}, i = {}, j = {}".format(k, i ,j )
    #         # print "U({},{}): ({:.17E},{:.17E})".format(j+1, i+1, uuu_r[k], uuu_i[k])    

    uuu = [complex(uuu_r[i], uuu_i[i]) for i in range(lrow * (KMAX + 1))]

    reorg = np.zeros( (lrow, KMAX + 1), np.complex128 )

    # for i, z in enumerate(uuu):
    #     print "{}: {:.17E}".format(i, z)

    # convert uuu into a numpy array and reorganize from Fortran's column-major
    # order.
    for j in range(lrow):
        for i in range(KMAX + 1):
            k = (i * lrow) + j

            reorg[j][i] = uuu[k]

    uuu = reorg

    # for i in range(lrow):
    #     for j in range(KMAX + 1):
    #         print "U({}, {}) = {:.17E}".format(i, j, uuu[i][j])



    # for i, z in enumerate(uuu):
    #     # U(           1           1 ): ( 9.10583222219728122E-002,-6.27880386519005795E-002)
    #     print "U({}): ({:.17E},{:.17E})".format(i, z.real, z.imag)
                               
    # # I tease the returned variables into Python types before passing them
    # # back to the caller. (Slicing a ctypes array returns a list.) The Fortran
    # # code has already sorted them the way we like (largest singular value 
    # # first).
    # singular_values = singular_values[:nsv_found]
    # frequencies     = frequencies[:nsv_found]
    # damping_factors = damping_factors[:nsv_found]
    # amplitudes      = amplitudes[:nsv_found]
    # phases          = phases[:nsv_found]

    # damping_factors = [1 / df for df in damping_factors]
    # damping_factors = [df * dwell_time for df in damping_factors]

    # frequencies = [frequency / dwell_time for frequency in frequencies]

    # phases = [phase * constants.RADIANS_TO_DEGREES for phase in phases]

    # Convert this from a ctypes thingy to regular Python, and trim.
    lsinval = [x for x in lsinval][:nsv_found]

    return uuu, lsinval, nsv_found


def _load_the_library(libhlsvd=None):
    if not libhlsvd:
        bits = util_misc.get_bit_mode()

        # Construct a name for the library.
        local_path = os.path.dirname(os.path.realpath(__file__))

        local_path = os.path.join(local_path, 'original/dist/bin')

        platform = util_misc.get_platform()

        local_path = os.path.join(local_path, "%s-%d" % (platform, bits))

        libhlsvd = os.path.join(local_path, _LIBRARY_NAMES[platform])

        libhlsvd = libhlsvd % bits

    if util_misc.get_platform() == "windows":
        # Under Windows 7, we need to set the CWD before attempting to load the
        # HLSVD DLL. See http://scion.duhs.duke.edu/vespa/analysis/ticket/51
        restore_cwd = os.getcwd()
        libpath = os.path.dirname(libhlsvd)
        os.chdir(libpath)
    else:
        restore_cwd = None

    try:
        libhlsvd = ctypes.CDLL(libhlsvd)
    except OSError, message:
        raise OSError, "Unable to load '%s'. The OS reports: %s" % \
                                (libhlsvd, message)
    finally:
        if restore_cwd:
            os.chdir(restore_cwd)

    return libhlsvd