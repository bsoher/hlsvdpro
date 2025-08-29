
This version of the pypropack library was originally named: hlsvdpro_make_v04_worksAll

I've just updated the make files to dump output into the 'hlsvdpro' top level directory to better work with the 'setup.py bdist_wheel' command

--------------------------------

June 14th, 2020 - Next big issue, Python 3 versions will compile, but won't import. Saw:

ImportError: DLL load failed: A dynamic link library (DLL) initialization routine failed. 

Tried a bunch of things, including checking dependencies with DependenciesGui.exe program found on github (very nice) and nothing obvious jumped out. Found one stackoverflow that might be relevant at: https://stackoverflow.com/questions/46505778/f2py-with-openmp-gives-import-error-in-python and their solution was to use a different mingw-w64 compiler.  I had been using TDF-gcc version from http://tdm-gcc.tdragon.net/ but now I'm going to try a new install from https://sourceforge.net/projects/mingw-w64/  fingers crossed.

--------------------------------

At this point, May 1st 2020, I have tried:

1) Just use f2py from cmd line OR from setup.py (numpy.disutils.misc_util.xxx) and in both cases I could not get them to statically link libgfortran or libgcc. So, when I tested out the *.pyd libraries the mingw gfortran.so (gfortran.so.5 ??) needed to be in my PATH or the import to python failed.  OH, and also I could not get the setup.py method of f2py to just create *.pyd files, only *.pyd stubs and *.dll pairs. 

2) However, both the above methods compiled for Win10 and gave good result from the example_propack.py test.

3) BUT, on Linux the pypropack_setup_v05_linux eventually compiled, but only the SINGLE and DOUBLE shared objects worked. And they needed accces to libgfortran.so.5 to run. And it had to be from some 'more typical' directory. I couldn't put them in /usr/local/lib64 and put that at the start of the PATH var and have that shared lib found??? Annoying.  

4) Still, like for Win10, the library needs to be statically linked to libgfortan.a I think.  Thus this next round of Makefile low level, control freak, OCD, mucking about in this hlsvdpro_make_v01_addAproInt_winWorks project.

5) This project started from the pypropack_setup_v04_addAprodInt as did the semi-defunct pypropack_setup_v05_linux project.

6) I am leaving the notes below for provenance.

---------------------

How I made Makefile work for me.

As in all things, start with the work that Philip did.  

I wanted to add in the f2py interface to our original ctypes work. This would allow me to forego translating the complex numpy arrays to twin float and then reassembling them in Fortran to COMPLEX. I did this by grabbing the project dependent _dpropack-f2pywrappers.f and _dpropackmodule.c files, and the project independent fortranobject.c and fortranobject.h files, that f2py generated from each *.pyf file we ran in the pypropack_setup_v04_addAprodInt.  I copied these to both an archive directory at top level and into each, respective, subdirectory for single, double, complex8 and complex16 libraries. These files go into a sub-sub-directroy called f2py_Util/ to simplify creating their Makefile rules. I also dumped all the additional files (mostly lapack and blas files missing from the original PROPACK 2.1 zip) into the Lapack_Util/ sub-sub-directory. 

I adapted the Makefile.win to have compile rules for the f2py_Util *.f and *.c files with different flags than those for the propack *.F and Lapack_Util/*.f file compile rules.  I also figured out how from within the Makefile to find the path to the python installation against which the library would be linked. Finally, I set up the linker flags to statically link to the libgfortran.a library under the D:\TDM-GCC-64 directory.  Just FYI, when I renamed that directory, temporarily, all the other windows *.pyd that I'd created various ways all failed to work ... |<:  Sigh. Change the name back and they all import and work again.

Other FYI things - Had to remove the OpenMP preprocessor sections from dblasext.F (and other versions of this).  Made sure that all files were listed in the Makefile as needed.

Here are some important snippets from the Makefile for the DOUBLE project, NB. NOT a full and working Makefile:

# These steps will locate the first python installation in your PATH variable.
# I use conda and activate the python that I want before doing my compile, and
# conda ensures that the right path is first in my PATH.
PYDIR := $(shell where python | head -n1)
PYDIR := $(subst python.exe,,$(PYDIR))
PYDIR := $(subst \,/,$(PYDIR))

F2PY_IFLAGS = -I$(PYDIR)lib/site-packages/numpy/core/include -I$(PYDIR)include
F2PY_CFLAGS = -g -DDEBUG -DMS_WIN64 -O0 -Wall -Wstrict-prototypes $(F2PY_IFLAGS)
F2PY_FFLAGS = -Wall -g -ffixed-form -fno-second-underscore -O3 -funroll-loops $(F2PY_IFLAGS)

# Pattern rules for turning f2py_Util src files into object files
$(BIN_DIR)/%.o: f2py_Util/%.c
    gcc $(F2PY_CFLAGS) -c $< -o $@

$(BIN_DIR)/%.o: f2py_Util/%.f
    gcc $(F2PY_FFLAGS) -c $< -o $@

 This builds a target that statically links in GFortran's libraries so
# we don't have to require mingw32 gfortran be installed for this to work.
$(TARGET_LIB_NAME): $(BIN_DIR) $(OBJS) lapack 
    gfortran -static-libgcc     \
            -shared                                 \
            -Wall $(ARCH)                           \
            -o $(TARGET_LIB_NAME)                   \
            $(OBJS)                                 \
            -L$(PYDIR)libs                          \
            -lpython27                              \
            $(BIN_DIR)/lapack.a                     \
            -Wl,-Bstatic                            \
            -lgfortran
            
Now in Linux ----

https://gcc.gnu.org/faq.html#multiple has comment on findind ld with multiple gcc installations

            

======================================================================================
Everything below here is from the pypropack_setup_v04_addAprodInt original notes. Pre May 1st 2020 stuff.
- main issue with all this work is that I could NOT get any of these build methods to statically link libgfortran and libgcc to the _xpropack.pyd modules. 
- So that limits my distribution a lot.
- Still, it was fun learning stuff.


This was a source that helped me figure out how to run f2py for PROPACK
----------------------------------------------------------------------------------
Compile fortran module with f2py and Python 3.6 on Windows 10
https://stackoverflow.com/questions/48826283/compile-fortran-module-with-f2py-and-python-3-6-on-windows-10

List of all files to be compiled:
-----------------------------------------------
SINGLE     slasq2.f slasq3.f slasq4.f slasq5.f slasq6.f slasr.f slasrt.f slassq.f slasv2.f smgs.pentium.F smgs.risc.F sreorth.F sritzvec.F ssafescal.F stat.h xerbla.f ieeeck.f ilaenv.f lsame.f printstat.F sbdsdc.f sbdsqr.f sblasext.F sbsvd.F second.F sgemm_ovwr.F sgetu0.F slacpy.f slaed6.f slamch.f slamrg.f slanbpro.F slanst.f slansvd.F slansvd_irl.F slapy2.f slarnv.f slartg.f slaruv.f slas2.f slascl.f slasd0.f slasd1.f slasd2.f slasd3.f slasd4.f slasd5.f slasd6.f slasd7.f slasd8.f slasda.f slasdq.f slasdt.f slaset.f slasq1.f snrm2.f srot.f sscal.f sswap.f saxpy.f scopy.f sgemm.f sgemv.f sdot.f saprod.f slansvd_aprod.F
DOUBLE     dlassq.f dlasv2.f zmgs.risc.F dmgs.pentium.F dreorth.F dritzvec.F dsafescal.F ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dgemm_ovwr.F dgetu0.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanbpro.F dlanst.f dlansvd.F dlansvd_irl.F dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dcabs1.f dnrm2.f dscal.f daxpy.f drot.f dswap.f ddot.f dznrm2.f dgemm.f dcopy.f dgemv.f daprod.f dlansvd_aprod.F
COMPLEX8   critzvec.F csafescal.F ieeeck.f ilaenv.f lsame.f printstat.F sbdsdc.f sbdsqr.f sblasext.F sbsvd.F second.F sgemm_ovwr.F slacpy.f slaed6.f slamch.f slamrg.f slanst.f slapy2.f slarnv.f slartg.f slaruv.f slas2.f slascl.f slasd0.f slasd1.f slasd2.f slasd3.f slasd4.f slasd5.f slasd6.f slasd7.f slasd8.f slasda.f slasdq.f slasdt.f slaset.f slasq1.f slasq2.f slasq3.f slasq4.f slasq5.f slasq6.f slasr.f slasrt.f slassq.f slasv2.f stat.h xerbla.f cblasext.F cgemm_ovwr.F cgetu0.F clanbpro.F clansvd.F clansvd_irl.F clarnv.f clascl.f cmgs.pentium.F cmgs.risc.F creorth.F saxpy.f scnrm2.f scopy.f sdot.f sgemm.f snrm2.f srot.f sscal.f sswap.f caxpy.f ccopy.f cdotc.f cdotu.f cscal.f csscal.f cgemv.f caprod.f clansvd_aprod.F
COMPLEX16  dgemm_ovwr.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanst.f dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dlassq.f dlasv2.f dmgs.pentium.F dnrm2.f drot.f dscal.f dswap.f dznrm2.f ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f zaxpy.f zblasext.F zcopy.f zdotc.f zdotu.f zdscal.f zgemm_ovwr.F zgemv.f zgetu0.F zlanbpro.F zlansvd.F zlansvd_irl.F zlarnv.f zlascl.f zmgs.pentium.F zreorth.F zritzvec.F zsafescal.F zscal.f daxpy.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dcabs1.f dcopy.f ddot.f dgemm.f zaprod.f zlansvd_aprod.F


What I am trying to do here:  I want to have the aprod() call be 'internal' to the zlansvd() call to see if that speeds up the overall run time for finding k singular values of A.  To do this, I have to have a available inside aprod() because it is not passed in as a parameter to the call. I think that this was to make it more flexible for users to declare as they see fit? Whatever. My problem was that I could not just put the array A into a common block because it was not declared to be of fixed dimensions. PROPACK's example does this by having a max array size A set at run time that the data can be read into.  

Found a solution that seemed to work using a pointer to point at A and putting the pointer into the common block.  Fingers crossed for zlansvd_aprod.f implementation ...


Create Signature File - this is a two step process:
---------------------------------------------------
1. Create a brand new *.pyf file for one of the subdirectories using the command below. Modify the result as directed and save. This only adds the *_aprod() functions to the signature files.

f2py -m _zpropack -h zpropack.pyf zlansvd_aprod.F zlansvd_irl_aprod.F zlansvd.F zlansvd_irl.F --overwrite-signature

NB. need to switch names for the 'z', 'c', 's', 'd' prefixes

-NB NB NB !!! two wrinkles - both are caused by the auto-magic of the signature file creation step. Both can be solved manually and the *.pyf file saved for (many re-) use.
  1) I had to manually remove the lines that refer to pointer 'pa' and common block 'common /csvdp/ pa' to make f2py work and create the fortranobject.c file properly. Some issue with complex pointers?  Still a mystery.
  2) I had to manually change lines with intent(inout) to be intent(in,out) for the results to return to Python correctly. I can not just put them into the Fortran code as intent(in,out) either, the signature file creations balks at that.

2. Now we add in the lansvd() and lansvd_irl() functions from the original PROPACK repos.  I don't know how they created these signature files. They have all the needed bells and whistles, but none of the decorators (ie. intent(out), dimension() etc.) are in the fortran files.  So when I just include the lansvd.F and lansvd_irl.F file to the call in Step 1, I get the wrong signatures.  Sigh.  So, we copy and paste the subroutine chunks from zpropack.pyf.Orig into the zpropack.pyf file we created in Step 1.  NB. don't forget to copy the Aprod section, too.



Create Notes for SINGLE library
--------------------------------------------------
Comand run at cmd line:
f2py -c spropack.pyf --compiler=mingw32 --fcompiler=gnu95 --build-dir build3 -include"<setjmpex.h>" -m _spropack slasq2.f slasq3.f slasq4.f slasq5.f slasq6.f slasr.f slasrt.f slassq.f slasv2.f smgs.pentium.F sreorth.F sritzvec.F ssafescal.F stat.h xerbla.f ieeeck.f ilaenv.f lsame.f printstat.F sbdsdc.f sbdsqr.f sblasext.F sbsvd.F second.F sgemm_ovwr.F sgetu0.F slacpy.f slaed6.f slamch.f slamrg.f slanbpro.F slanst.f slansvd.F slansvd_irl.F slapy2.f slarnv.f slartg.f slaruv.f slas2.f slascl.f slasd0.f slasd1.f slasd2.f slasd3.f slasd4.f slasd5.f slasd6.f slasd7.f slasd8.f slasda.f slasdq.f slasdt.f slaset.f slasq1.f snrm2.f srot.f sscal.f sswap.f saxpy.f scopy.f sgemm.f sgemv.f sdot.f saprod.f slansvd_irl_aprod.F slansvd_aprod.F


Create Notes for DOUBLE library
--------------------------------------------------
Comand run at cmd line:
f2py -c dpropack.pyf --compiler=mingw32 --fcompiler=gnu95 -include"<setjmpex.h>" -m _dpropack dlassq.f dlasv2.f dmgs.pentium.F dreorth.F dritzvec.F dsafescal.F ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dgemm_ovwr.F dgetu0.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanbpro.F dlanst.f dlansvd.F dlansvd_irl.F dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dcabs1.f dnrm2.f dscal.f daxpy.f drot.f dswap.f ddot.f dznrm2.f dgemm.f dcopy.f dgemv.f daprod.F dlansvd_irl_aprod.F dlansvd_aprod.F

(optional - for local build directory)
f2py -c dpropack.pyf --compiler=mingw32 --fcompiler=gnu95 --build-dir build3 -include"<setjmpex.h>" -m _dpropack dlassq.f dlasv2.f dmgs.pentium.F dreorth.F dritzvec.F dsafescal.F ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dgemm_ovwr.F dgetu0.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanbpro.F dlanst.f dlansvd.F dlansvd_irl.F dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dcabs1.f dnrm2.f dscal.f daxpy.f drot.f dswap.f ddot.f dznrm2.f dgemm.f dcopy.f dgemv.f daprod.f dlansvd_irl_aprod.F dlansvd_aprod.F


Create Notes for COMPLEX8 library
--------------------------------------------------
Comand run at cmd line:
f2py -c cpropack.pyf --compiler=mingw32 --fcompiler=gnu95 --build-dir build3 -include"<setjmpex.h>" -m _cpropack critzvec.F csafescal.F ieeeck.f ilaenv.f lsame.f printstat.F sbdsdc.f sbdsqr.f sblasext.F sbsvd.F second.F sgemm_ovwr.F slacpy.f slaed6.f slamch.f slamrg.f slanst.f slapy2.f slarnv.f slartg.f slaruv.f slas2.f slascl.f slasd0.f slasd1.f slasd2.f slasd3.f slasd4.f slasd5.f slasd6.f slasd7.f slasd8.f slasda.f slasdq.f slasdt.f slaset.f slasq1.f slasq2.f slasq3.f slasq4.f slasq5.f slasq6.f slasr.f slasrt.f slassq.f slasv2.f stat.h xerbla.f cblasext.F cgemm_ovwr.F cgetu0.F clanbpro.F clansvd.F clansvd_irl.F clarnv.f clascl.f cmgs.pentium.F creorth.F saxpy.f scnrm2.f scopy.f sdot.f sgemm.f snrm2.f srot.f sscal.f sswap.f caxpy.f ccopy.f cdotc.f cdotu.f cscal.f csscal.f cgemv.f caprod.f clansvd_irl_aprod.F clansvd_aprod.F


Create Notes for COMPLEX16 library
--------------------------------------------------
Comand run at cmd line:
f2py -c zpropack.pyf --compiler=mingw32 --fcompiler=gnu95 -include"<setjmpex.h>" -m _zpropack dgemm_ovwr.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanst.f dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dlassq.f dlasv2.f dmgs.pentium.F dnrm2.f drot.f dscal.f dswap.f dznrm2.f ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f zaxpy.f zblasext.F zcopy.f zdotc.f zdotu.f zdscal.f zgemm_ovwr.F zgemv.f zgetu0.F zlanbpro.F zlansvd.F zlansvd_irl.F zlarnv.f zlascl.f zmgs.pentium.F zreorth.F zritzvec.F zsafescal.F zscal.f daxpy.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dcabs1.f dcopy.f ddot.f dgemm.f zaprod.f zlansvd_irl_aprod.F zlansvd_aprod.F

(optional - for local build directory)
f2py -c zpropack.pyf --compiler=mingw32 --fcompiler=gnu95 --build-dir build3 -include"<setjmpex.h>" -m _zpropack dgemm_ovwr.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanst.f dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dlassq.f dlasv2.f dmgs.pentium.F dnrm2.f drot.f dscal.f dswap.f dznrm2.f ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f zaxpy.f zblasext.F zcopy.f zdotc.f zdotu.f zdscal.f zgemm_ovwr.F zgemv.f zgetu0.F zlanbpro.F zlansvd.F zlansvd_irl.F zlarnv.f zlascl.f zmgs.pentium.F zreorth.F zritzvec.F zsafescal.F zscal.f daxpy.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dcabs1.f dcopy.f ddot.f dgemm.f zaprod.f zlansvd_irl_aprod.F zlansvd_aprod.F



Global Error Notes
--------------------------------------------
NB. dmgs.risc.F  and  dmgs.pentium.F are NOT both needed, pick one
NB. zmgs.risc.F  and  zmgs.pentium.F are NOT both needed, pick one




Error due to lack of <setjmpex.h> include file 
----------------------------------------------------------------------

1. Original fix for COMPLEX16 -> Manual edit of _zpropackmodule.c and manual re-compile/link from midway thru f2py steps.

    I saw:
    C:/TDM-GCC-64/bin/../lib/gcc/x86_64-w64-mingw32/9.2.0/../../../../x86_64-w64-mingw32/bin/ld.exe: d:\users\bsoher\appdata\local\temp\tmp4jx2uy\Release\users\bsoher\appdata\local\temp\tmp4jx2uy\src.win-amd64-2.7\_dpropackmodule.o: in function `f2py_rout__dpropack_dlansvd':
    d:/users/bsoher/appdata/local/temp/tmp4jx2uy/src.win-amd64-2.7/_dpropackmodule.c:950: undefined reference to `__intrinsic_setjmpex'
    Found suggested fix, but could only implement manually ...
    d:/users/bsoher/appdata/local/temp/tmphvg_t5/src.win-amd64-2.7/_zpropackmodule.c
    - changed //#include <setjmp.h>
          #include <setjmpex.h>
    - reran this line: (careful, in the cmd window it was concatenated with a similar next line.  Look for '.ogcc ' to separate)
        gcc -g -DDEBUG -DMS_WIN64 -O0 -Wall -Wstrict-prototypes -DNPY_MINGW_USE_CUSTOM_MSVCR -D__MSVCRT_VERSION__=0x1500 -Id:\users\bsoher\appdata\local\temp\tmphvg_t5\src.win-amd64-2.7 -ID:\Users\bsoher\miniconda2\envs\python27_vespa\lib\site-packages\numpy\core\include -ID:\Users\bsoher\miniconda2\envs\python27_vespa\include -ID:\Users\bsoher\miniconda2\envs\python27_vespa\PC -c d:\users\bsoher\appdata\local\temp\tmphvg_t5\src.win-amd64-2.7\_zpropackmodule.c -o d:\users\bsoher\appdata\local\temp\tmphvg_t5\Release\users\bsoher\appdata\local\temp\tmphvg_t5\src.win-amd64-2.7\_zpropackmodule.o
    and it worked.

2. NB. Saw the same error on DOUBLE compile, so did same fix and it worked

3. Found way to AUTOMATE step 1 using -include option in f2py re:
https://docs.scipy.org/doc/numpy/f2py/usage.html#command-f2py







