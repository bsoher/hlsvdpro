#!/bin/bash 
# This batch file builds the hlsvdpro wheels for distribution
#   It switches between Python versions using 'conda activate xxx'
#   and then calls the 'build_wheel.py' script
source ~/miniconda2/etc/profile.d/conda.sh
clear
echo Building wheel for - Python 36
conda activate python36
python -V
cd ./pypropack
make clean
make
make dist
cd ..
python build_wheel.py 36 
echo ---- done -----------------------------------------------------------------
echo Building wheel for - Python 37
conda activate python37
python -V
cd ./pypropack
make clean
make
make dist
cd ..
python build_wheel.py 37
echo ---- done -----------------------------------------------------------------
echo Building wheel for - Python 38
conda activate python38
python -V
cd ./pypropack
make clean
make
make dist
cd ..
python build_wheel.py 38
echo ---- done -----------------------------------------------------------------
echo Building wheel for - Python 27
conda activate python27
python -V
cd ./pypropack
make clean
make
make dist
cd ..
python build_wheel.py 27 
echo ---- done -----------------------------------------------------------------
echo Done building wheels
echo See ya ... !
