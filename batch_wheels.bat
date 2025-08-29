@ECHO OFF
:: This batch file builds the hlsvdpro wheels for distribution
::   It switches between Python versions using 'conda activate xxx'
::   and then calls the 'build_wheel.py' script
ECHO Building wheel for - Python 36
call conda activate python36
python -V
cd ./pypropack
call mingw32-make clean
call mingw32-make
call mingw32-make dist
cd ..
call python build_wheel.py 36 
ECHO ---- done -----------------------------------------------------------------
::PAUSE
ECHO Building wheel for - Python 37
call conda activate python37
python -V
cd ./pypropack
call mingw32-make clean
call mingw32-make
call mingw32-make dist
cd ..
call python build_wheel.py 37
ECHO ---- done -----------------------------------------------------------------
::PAUSE
ECHO Building wheel for - Python 38
call conda activate python38
python -V
cd ./pypropack
call mingw32-make clean
call mingw32-make
call mingw32-make dist
cd ..
python build_wheel.py 38
ECHO ---- done -----------------------------------------------------------------
::PAUSE
ECHO Building wheel for - Python 27
call conda activate python27
python -V
cd ./pypropack
call mingw32-make clean
call mingw32-make
call mingw32-make dist
cd ..
python build_wheel.py 27 
ECHO ---- done -----------------------------------------------------------------
ECHO Done building wheels
ECHO See ya ... !
