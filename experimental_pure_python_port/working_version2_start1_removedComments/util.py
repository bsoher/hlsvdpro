# Python modules
from __future__ import division
import os
import ctypes
import math

# 3rd party modules
import numpy as np

# Our modules
import vespa.common.util.misc as util_misc

# _LIBRARY_NAMES maps platforms to the name of the library on said platform.
_LIBRARY_NAMES = { "osx"       : "libhlsvd-%d.dylib",
                   "windows"   : "hlsvdpro-%d.dll", #"hlsvd-%d.dll",
                   "linux"     : "libhlsvd-%d.so",
                 }

# Variable name translations 
# Fortran    Python             Description
# -------    ------             -------------------------------------------
# kuser      nsv_sought         # of singular values requested by caller
# kfit       nsv_found          # of singular values found
# lsinval    singular_values    Array containing the singular values (floats)
# ndp        n_data_points      
# kmax       MAX_SINGULAR_VALUES  Hardcoded to 50.

def dlapy2(a, b):
    
#    return np.hypot(a,b)
    return math.sqrt((a ** 2) + (b ** 2))


def flat_fortran_to_2d_python(a_list, dim0, dim1):
    """Given a flat (one dimensional) Python list of numbers in Fortran order
    (flattened column-major), reorders them in C order (row major) and shapes
    them according to the two dimensions given.

    Useful for when e.g. a Fortran array like U must be returned as a flat
    list from Fortran and then reshaped in to a numpy array.
    """
    temp = [None] * len(a_list)

    for i, value in enumerate(a_list):
       j = (i // dim0) + ((i % dim0) * (dim1))
       temp[j] = value

    temp = np.array(temp)

    temp.shape = (dim0, dim1)

    return temp



def load_the_library(libhlsvd=None):
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

