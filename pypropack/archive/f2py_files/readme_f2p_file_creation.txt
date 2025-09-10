This is a synopsis of how these files were created for use in compiling the hlsvdpro package.

Background
----------
Between hlsvdpro version 1.0.x and 1.1.x we made considerable changes to the underlying code. But, with the itention of achieving comprable results. 
 - We still use the PROPACK library, written by Rasmus Munk, to calculate a subset of SVD values, but we moved from PROPACK v1.1 to v2.1 the most recent version. 
 - We are now using f2py style access to the PROPACK fortran code because it more easily handles parameters that are complex number arrays. Using f2py added some complexity, but in a 'one and done' sort of way. The additional files are fixed for the current release and don't have to be recreated at compile time. 
 - Between version 1.0.x and 1.1.x we have reworked our hlsvd() algorithm to more closely match the workflow in the Laudadio paper, than the one in the original Pijnapple paper. The original algo used an FFT to massage the data going into the SVD calc. The Laudadio paper just arranges it into a Hankel matrix and off to the races we go.
 - We have also moved more of the algo steps into Python/Scipy (1.1.x) rather than doing it all in fortran (1.0.x). This simplifies the compile, icreases transparency and adds only a small time hit to processing.
 - In addition to the hlsvd() method we also provide direct access to the PROPACK SVD calculation routines via an svdp() method.
 - Finally, we have written a native Python version of the PROPACK complex16 SVD routine for situations were using the Fortan library might not be feasible. Final time tests remain to be taken, but for now it seems to be 4x-5x slower.
 
Create Signature File - this is a two step process:
---------------------------------------------------
1. Create a brand new *.pyf file for one of the subdirectories using the command below. Modify the result as directed and save. This only adds the *_aprod() functions to the signature files.

f2py -m _zpropack -h zpropack.pyf zlansvd_aprod.F zlansvd_irl_aprod.F zlansvd.F zlansvd_irl.F --overwrite-signature

NB. need to switch names for the 'z', 'c', 's', 'd' prefixes to create in different subdirectories of PROPACK

-NB NB NB !!! two wrinkles - both are caused by the auto-magic of the signature file creation step. Both can be solved manually and the *.pyf file saved for (many re-) use.
  1) I had to manually remove the lines that refer to pointer 'pa' and common block 'common /csvdp/ pa' to make f2py work and create the fortranobject.c file properly. Some issue with complex pointers?  Still a mystery.
  2) (python 3.8 and maybe 3.9 and 3.10 and older - definitely for numpy 1.xx) I had to manually change lines with intent(inout) to be intent(in,out) for the results to return to Python correctly. I can not just put them into the Fortran code as intent(in,out) either, the signature file creations balks at that.
  2) (python 3.11 and newer - definitely numpy 2.x - newer f2py call) lines with intent(in,out) were OK, but not automatically added to bottom two functions clansvd() and clansvd_irl(). Same issue for intent(in) and intent(out). Had to bring up OLD version of signature file and WinMerge the two files to determine where to modify.

2. Now we add in the lansvd() and lansvd_irl() functions from the original PROPACK repos.  I don't know how they created these signature files. They have all the needed bells and whistles, but none of the decorators (ie. intent(out), dimension() etc.) are in the fortran files.  So when I just include the lansvd.F and lansvd_irl.F file to the call in Step 1, I get the wrong signatures.  Sigh.  So, we copy and paste the subroutine chunks from zpropack.pyf.Orig into the zpropack.pyf file we created in Step 1.  NB. don't forget to copy the Aprod section, too.



Compile Notes for SINGLE library - NB. bjs 2025 - build3 local dir option did not work for SINGLE ONLY, mytemp option did work .. weird.
--------------------------------------------------
Comand run at cmd line:
f2py -c spropack.pyf --compiler=mingw32 --fcompiler=gnu95 --build-dir mytemp -include"<setjmpex.h>" -m _spropack slasq2.f slasq3.f slasq4.f slasq5.f slasq6.f slasr.f slasrt.f slassq.f slasv2.f smgs.pentium.F sreorth.F sritzvec.F ssafescal.F stat.h xerbla.f ieeeck.f ilaenv.f lsame.f printstat.F sbdsdc.f sbdsqr.f sblasext.F sbsvd.F second.F sgemm_ovwr.F sgetu0.F slacpy.f slaed6.f slamch.f slamrg.f slanbpro.F slanst.f slansvd.F slansvd_irl.F slapy2.f slarnv.f slartg.f slaruv.f slas2.f slascl.f slasd0.f slasd1.f slasd2.f slasd3.f slasd4.f slasd5.f slasd6.f slasd7.f slasd8.f slasda.f slasdq.f slasdt.f slaset.f slasq1.f snrm2.f srot.f sscal.f sswap.f saxpy.f scopy.f sgemm.f sgemv.f sdot.f saprod.f slansvd_irl_aprod.F slansvd_aprod.F

(optional - for local build directory)
f2py -c spropack.pyf --compiler=mingw32 --fcompiler=gnu95 --build-dir build3 --build-dir mytemp -include"<setjmpex.h>" -m _spropack slasq2.f slasq3.f slasq4.f slasq5.f slasq6.f slasr.f slasrt.f slassq.f slasv2.f smgs.pentium.F sreorth.F sritzvec.F ssafescal.F stat.h xerbla.f ieeeck.f ilaenv.f lsame.f printstat.F sbdsdc.f sbdsqr.f sblasext.F sbsvd.F second.F sgemm_ovwr.F sgetu0.F slacpy.f slaed6.f slamch.f slamrg.f slanbpro.F slanst.f slansvd.F slansvd_irl.F slapy2.f slarnv.f slartg.f slaruv.f slas2.f slascl.f slasd0.f slasd1.f slasd2.f slasd3.f slasd4.f slasd5.f slasd6.f slasd7.f slasd8.f slasda.f slasdq.f slasdt.f slaset.f slasq1.f snrm2.f srot.f sscal.f sswap.f saxpy.f scopy.f sgemm.f sgemv.f sdot.f saprod.f slansvd_irl_aprod.F slansvd_aprod.F

Create Notes for DOUBLE library
--------------------------------------------------
Comand run at cmd line:
f2py -c dpropack.pyf --compiler=mingw32 --fcompiler=gnu95 -include"<setjmpex.h>" -m _dpropack dlassq.f dlasv2.f dmgs.pentium.F dreorth.F dritzvec.F dsafescal.F ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dgemm_ovwr.F dgetu0.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanbpro.F dlanst.f dlansvd.F dlansvd_irl.F dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dcabs1.f dnrm2.f dscal.f daxpy.f drot.f dswap.f ddot.f dznrm2.f dgemm.f dcopy.f dgemv.f daprod.F dlansvd_irl_aprod.F dlansvd_aprod.F

(optional - for local build directory)
f2py -c dpropack.pyf --compiler=mingw32 --fcompiler=gnu95 --build-dir build3 -include"<setjmpex.h>" -m _dpropack dlassq.f dlasv2.f dmgs.pentium.F dreorth.F dritzvec.F dsafescal.F ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dgemm_ovwr.F dgetu0.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanbpro.F dlanst.f dlansvd.F dlansvd_irl.F dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dcabs1.f dnrm2.f dscal.f daxpy.f drot.f dswap.f ddot.f dznrm2.f dgemm.f dcopy.f dgemv.f daprod.f dlansvd_irl_aprod.F dlansvd_aprod.F


Create Notes for COMPLEX8 library
--------------------------------------------------
Comand run at cmd line:
f2py -c cpropack.pyf --compiler=mingw32 --fcompiler=gnu95 -include"<setjmpex.h>" -m _cpropack critzvec.F csafescal.F ieeeck.f ilaenv.f lsame.f printstat.F sbdsdc.f sbdsqr.f sblasext.F sbsvd.F second.F sgemm_ovwr.F slacpy.f slaed6.f slamch.f slamrg.f slanst.f slapy2.f slarnv.f slartg.f slaruv.f slas2.f slascl.f slasd0.f slasd1.f slasd2.f slasd3.f slasd4.f slasd5.f slasd6.f slasd7.f slasd8.f slasda.f slasdq.f slasdt.f slaset.f slasq1.f slasq2.f slasq3.f slasq4.f slasq5.f slasq6.f slasr.f slasrt.f slassq.f slasv2.f stat.h xerbla.f cblasext.F cgemm_ovwr.F cgetu0.F clanbpro.F clansvd.F clansvd_irl.F clarnv.f clascl.f cmgs.pentium.F creorth.F saxpy.f scnrm2.f scopy.f sdot.f sgemm.f snrm2.f srot.f sscal.f sswap.f caxpy.f ccopy.f cdotc.f cdotu.f cscal.f csscal.f cgemv.f caprod.f clansvd_irl_aprod.F clansvd_aprod.F

(optional - for local build directory)
f2py -c cpropack.pyf --compiler=mingw32 --fcompiler=gnu95 --build-dir build3 -include"<setjmpex.h>" -m _cpropack critzvec.F csafescal.F ieeeck.f ilaenv.f lsame.f printstat.F sbdsdc.f sbdsqr.f sblasext.F sbsvd.F second.F sgemm_ovwr.F slacpy.f slaed6.f slamch.f slamrg.f slanst.f slapy2.f slarnv.f slartg.f slaruv.f slas2.f slascl.f slasd0.f slasd1.f slasd2.f slasd3.f slasd4.f slasd5.f slasd6.f slasd7.f slasd8.f slasda.f slasdq.f slasdt.f slaset.f slasq1.f slasq2.f slasq3.f slasq4.f slasq5.f slasq6.f slasr.f slasrt.f slassq.f slasv2.f stat.h xerbla.f cblasext.F cgemm_ovwr.F cgetu0.F clanbpro.F clansvd.F clansvd_irl.F clarnv.f clascl.f cmgs.pentium.F creorth.F saxpy.f scnrm2.f scopy.f sdot.f sgemm.f snrm2.f srot.f sscal.f sswap.f caxpy.f ccopy.f cdotc.f cdotu.f cscal.f csscal.f cgemv.f caprod.f clansvd_irl_aprod.F clansvd_aprod.F

Create Notes for COMPLEX16 library
--------------------------------------------------
Comand run at cmd line:
f2py -c zpropack.pyf --compiler=mingw32 --fcompiler=gnu95 -include"<setjmpex.h>" -m _zpropack dgemm_ovwr.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanst.f dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dlassq.f dlasv2.f dmgs.pentium.F dnrm2.f drot.f dscal.f dswap.f dznrm2.f ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f zaxpy.f zblasext.F zcopy.f zdotc.f zdotu.f zdscal.f zgemm_ovwr.F zgemv.f zgetu0.F zlanbpro.F zlansvd.F zlansvd_irl.F zlarnv.f zlascl.f zmgs.pentium.F zreorth.F zritzvec.F zsafescal.F zscal.f daxpy.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dcabs1.f dcopy.f ddot.f dgemm.f zaprod.f zlansvd_irl_aprod.F zlansvd_aprod.F

(optional - for local build directory)
f2py -c zpropack.pyf --compiler=mingw32 --fcompiler=gnu95 --build-dir build3 -include"<setjmpex.h>" -m _zpropack dgemm_ovwr.F dlacpy.f dlaed6.f dlamch.f dlamrg.f dlanst.f dlapy2.f dlarnv.f dlartg.f dlaruv.f dlas2.f dlascl.f dlasd0.f dlasd1.f dlasd2.f dlasd3.f dlasd4.f dlasd5.f dlasd6.f dlasd7.f dlasd8.f dlasda.f dlasdq.f dlasdt.f dlaset.f dlasq1.f dlasq2.f dlasq3.f dlasq4.f dlasq5.f dlasq6.f dlasr.f dlasrt.f dlassq.f dlasv2.f dmgs.pentium.F dnrm2.f drot.f dscal.f dswap.f dznrm2.f ieeeck.f ilaenv.f lsame.f printstat.F second.F stat.h xerbla.f zaxpy.f zblasext.F zcopy.f zdotc.f zdotu.f zdscal.f zgemm_ovwr.F zgemv.f zgetu0.F zlanbpro.F zlansvd.F zlansvd_irl.F zlarnv.f zlascl.f zmgs.pentium.F zreorth.F zritzvec.F zsafescal.F zscal.f daxpy.f dbdsdc.f dbdsqr.f dblasext.F dbsvd.F dcabs1.f dcopy.f ddot.f dgemm.f zaprod.f zlansvd_irl_aprod.F zlansvd_aprod.F

Global Error Notes
--------------------------------------------
NB. (z,c,d,s)mgs.risc.F  and  (z,c,d,s)mgs.pentium.F are NOT both needed, pick one



Python 2 vs Python 3 Notes
------------------------------------

When I ran the above f2py commands, it created some temp files in a system dir (or the local build3 dir) that are wrappers for C and Fortran that allow my Fortran code to talk to Python. These had names like : _cpropack-f2pywrappers.f  _cpropackmodule.c  _zpropack-f2pywrappers.f  _zpropackmodule.c and so on for 's' and 'd'.  These differ by input/output variable datatypes. There were also two other files, one C header and one C code, that appear to be more generic 'helper' code that were generated. All copies of these were identical for the four data types. They were called fortranobject.c and fortranobject.h

When I compare running f2py in Python 2 and then in Python 3, all of these files (except for fortranobject.h - but lets just say all files for simplicity) they all differed.  I've run the f2py call manually for each data type and kept copies of each of these wrapper files. I store them in pypropack/archive/f2py_files in either the py2_f2py or py3_f2py subdirectories.  These are also located in each of the complex16, complex8, single and double subdirectories in the f2py_Util sub-subdirectory. 

