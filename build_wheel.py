# Python modules
from __future__ import division
from __future__ import print_function

import sys
import os
import os.path
import shutil
import distutils.util

def build_wheel(which='ugh',test=False):

    base_path = os.path.dirname(os.path.abspath(__file__))
    dstdir = os.path.join(base_path, 'hlsvdpro')

    pyver2 = str(sys.version_info.major) + str(sys.version_info.minor)
    if pyver2 != which:
        raise ValueError('Error, pytag does not match Python version run from cmd line = '+str(which))
    
    if which == '27':
        pytag = 'py27'
    elif which == '36':
        pytag = 'py36'
    elif which == '37':
        pytag = 'py37'
    elif which == '38':
        pytag = 'py38'
    elif which == '39':
        pytag = 'py39'
    elif which == '310':
        pytag = 'py310'
    elif which == '311':
        pytag = 'py311'
    elif which == '312':
        pytag = 'py312'
    elif which == '313':
        pytag = 'py313'
    else:
        raise ValueError('Error, python tag not recognized')
        
    plat_name = distutils.util.get_platform().lower()
    if 'win-amd64' in plat_name:
        plat_name = 'win_amd64'
        libdir = '_win_'+pytag
        libnames = ['_cpropack.pyd','_dpropack.pyd','_spropack.pyd','_zpropack.pyd']
    elif 'linux' in plat_name:
        plat_name = 'manylinux2014_x86_64'
        libdir = '_linux_'+pytag
        libnames = ['_cpropack.so','_dpropack.so','_spropack.so','_zpropack.so','libgfortran.so.3']
    elif 'macosx' in plat_name:
        plat_name = 'macosx_10_9_x86_64'
        libdir = '_osx_'+pytag
        libnames = ['_cpropack.so','_dpropack.so','_spropack.so','_zpropack.so']
    else:
        raise ValueError('Error, platform name not recognized')

    # 1. remove build/ and egg-info/ dirs from hlsvdpro/
    
    if not test:
        shutil.rmtree(os.path.join(base_path,'build'), ignore_errors=True)
        shutil.rmtree(os.path.join(base_path,'hlsvdpro.egg-info'), ignore_errors=True)
    else:
        print("\n1. remove build/ and egg_info dirs from hlsvdpro/ \n")
        print("shutil.rmtree("+os.path.join(base_path,'build')+", ignore_errors=True)")
        print("shutil.rmtree("+os.path.join(base_path,'hlsvdpro.egg-info')+", ignore_errors=True)")
    
    # 2. clean lib files from hlsvdpro/hlsvdpro

    for libname in libnames:
        if not test:
            fname = os.path.join(base_path,'hlsvdpro',libname)
            if os.path.exists(fname):
                os.remove(fname)
        else:
            print("\n2. clean lib files from ~hlsvdpro/hlsvdpro \n")
            print("os.remove("+os.path.join(base_path,'hlsvdpro',libname)+")")
    
    # 3. copy python version appropriate libs to hlsvdpro/hlsvdpro

    for libname in libnames:
        srcdir = os.path.join(base_path, 'hlsvdpro', libdir, libname)
        if not test:
            shutil.copy2(srcdir, dstdir)
        else:
            print("\n3. copy python version appropriate libs to ~hlsvdpro/hlsvdpro \n")
            print("shutil.copy2("+srcdir+","+dstdir+")")

    # 4. run setup.py

    cmd = "python setup.py bdist_wheel --python-tag "+pytag+" --plat-name "+plat_name
    if not test:
        os.system(cmd)
    else:  
        print("\n4. run setup.py \n")
        print("os.system("+cmd+")")
        print("\ndone - build_wheel")
    
#------------------------------------------------------------

if __name__ == "__main__":

    if len(sys.argv) < 2:
        print( "Please supply the Python version to be used (ex. '39')")
    elif len(sys.argv) < 3:
        pyver = str(sys.argv[1])
        test = False
    else:
        pyver = str(sys.argv[1])
        test = True
    
    build_wheel(which=pyver, test=test)
