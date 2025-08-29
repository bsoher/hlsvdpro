# Python modules
from __future__ import division
import ctypes
import pdb

# 3rd party modules
import scipy
import scipy.linalg.lapack
import scipy.linalg

print "scipy version is {}".format(scipy.__version__)
print "clapack.__file__ is {}".format(scipy.linalg.lapack.clapack.__file__)
print "flapack.__file__ is {}".format(scipy.linalg.lapack.flapack.__file__)

name = 'gelss'

raw_names = [name, 'z' + name, 'd' + name, ]

names =  raw_names[:]
names += ['_' + name       for name in raw_names]
names += [      name + '_' for name in raw_names]
names += ['_' + name + '_' for name in raw_names]


lib = ctypes.CDLL(scipy.linalg.lapack.clapack.__file__)

for name in names:
    found = "PRESENT  " if hasattr(lib, name) else "NOT FOUND"

    print "{} in clapack: {}".format(found, name)

lib = ctypes.CDLL(scipy.linalg.lapack.flapack.__file__)

for name in names:
    found = "PRESENT  " if hasattr(lib, name) else "NOT FOUND"

    print "{} in flapack: {}".format(found, name)

# try:
#     lib.dbdsqr
#     print "dbdsqr is PRESENT in clapack"
# except AttributeError:
#     print "dbdsqr NOT FOUND in clapack"

# lib = ctypes.CDLL(scipy.linalg.lapack.flapack.__file__)

# try:
#     lib.dbdsqr
#     print "dbdsqr is PRESENT in flapack"
# except AttributeError:
#     print "dbdsqr NOT FOUND in flapack"
