
These are results for: working_version2_start1_removedComments
Using data set:  original/dist/press_cp0.xml
- fun fact, code is a bit (10%) faster for press_cp5.xml
Ran both hlsvdpro_local and hlsvd 10 times for same data
Plots look very similar

Code Changes (best as I can remember)
--------------------------------------
- removed a bunch of comments to read more clearly
- moved all the dlapy2 (util.py) calls to inline - half their time spent in overhead (0.15 sec or so)
- maybe others, can't recall


pydev debugger: starting (pid: 26324)
hankel_size len() // 8 = 512
nsv_sought             = 20
Mon Mar 30 10:24:20 2020    profile.data

         734817 function calls (734807 primitive calls) in 4.977 seconds

   Ordered by: cumulative time

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000    4.977    4.977 <string>:1(<module>)
        1    0.007    0.007    4.976    4.976 test.py:37(test)
       10    0.027    0.003    4.079    0.408 hlsvdpro_local.py:38(hlsvdpro)
       10    0.007    0.001    3.679    0.368 lanczopw.py:29(lanczopw)
       20    0.020    0.001    3.668    0.183 zlansvdw.py:61(zlansvdw)
       76    0.404    0.005    3.516    0.046 zlanbprow.py:75(zlanbprow)
     1804    1.195    0.001    2.058    0.001 zreorth.py:192(zreorth2)
       10    0.858    0.086    0.864    0.086 hlsvd.py:21(hlsvd)
     1912    0.064    0.000    0.593    0.000 aprodw.py:13(aprodw)
     1912    0.056    0.000    0.529    0.000 aprodw.py:53(ztmultz)
    49842    0.155    0.000    0.465    0.000 {method 'sum' of 'numpy.ndarray' objects}
     3844    0.039    0.000    0.404    0.000 fftpack.py:55(_raw_fft)
    49842    0.025    0.000    0.309    0.000 _methods.py:34(_sum)
    47592    0.298    0.000    0.298    0.000 {method 'conjugate' of 'numpy.ndarray' objects}
    50866    0.290    0.000    0.290    0.000 {method 'reduce' of 'numpy.ufunc' objects}
       10    0.265    0.026    0.271    0.027 hlsvdpro_local.py:136(_vanmon)
     1912    0.037    0.000    0.242    0.000 fftpack.py:212(ifft)
     1932    0.016    0.000    0.234    0.000 fftpack.py:115(fft)
    92914    0.154    0.000    0.154    0.000 {method 'conjugate' of 'numpy.generic' objects}
      936    0.129    0.000    0.145    0.000 zlanbprow.py:694(dupdate_mu)
      936    0.124    0.000    0.140    0.000 zlanbprow.py:750(dupdate_nu)
     1932    0.138    0.000    0.138    0.000 {numpy.fft.fftpack_lite.cfftf}
     1912    0.131    0.000    0.131    0.000 {numpy.fft.fftpack_lite.cfftb}
       10    0.045    0.004    0.099    0.010 hlsvdpro_local.py:155(_zcalc)
       40    0.034    0.001    0.067    0.002 zget0w.py:50(zgetu0w)
     1014    0.008    0.000    0.054    0.000 misc.py:19(norm)
     3844    0.034    0.000    0.051    0.000 helper.py:259(put_twiddle_factors)
     3844    0.018    0.000    0.039    0.000 helper.py:285(pop_twiddle_factors)
     1014    0.035    0.000    0.035    0.000 {numpy.vdot}
       76    0.019    0.000    0.033    0.000 dbdsqr.py:180(dbdsqr_cython)
     9724    0.029    0.000    0.031    0.000 {numpy.array}
     1014    0.017    0.000    0.027    0.000 function_base.py:434(asarray_chkfinite)
     7687    0.017    0.000    0.025    0.000 collections.py:149(pop)
        1    0.000    0.000    0.025    0.025 xml_.py:355(element_to_dict)
        6    0.000    0.000    0.025    0.004 xml_.py:260(_element_to_any_type)
        1    0.000    0.000    0.025    0.025 xml_.py:496(element_to_numpy_array)
        1    0.000    0.000    0.025    0.025 xml_.py:542(decode_numeric_list)
     1892    0.024    0.000    0.024    0.000 zsafescal.py:11(zsafescal)
        1    0.000    0.000    0.024    0.024 fileio.py:181(decode_xdr)
        1    0.012    0.012    0.023    0.023 xdrlib.py:240(unpack_farray)
     1014    0.004    0.000    0.019    0.000 blas.py:343(get_blas_funcs)
      879    0.014    0.000    0.017    0.000 zlanbprow.py:621(dcompute_int)
     1804    0.012    0.000    0.017    0.000 zlanbprow.py:605(dset_mu)
       10    0.016    0.002    0.017    0.002 dbdsqr.py:35(dbdsqr)
     1014    0.010    0.000    0.015    0.000 blas.py:297(_get_funcs)
     2445    0.014    0.000    0.015    0.000 {method 'copy' of 'numpy.ndarray' objects}
    91142    0.012    0.000    0.012    0.000 {math.sqrt}
     7256    0.005    0.000    0.011    0.000 numeric.py:469(asarray)
      456    0.001    0.000    0.011    0.000 _internal.py:289(__init__)
       10    0.000    0.000    0.010    0.001 defmatrix.py:217(__mul__)
       10    0.010    0.001    0.010    0.001 {numpy.dot}
     8855    0.010    0.000    0.010    0.000 {range}
    52306    0.010    0.000    0.010    0.000 {max}
     3709    0.004    0.000    0.010    0.000 zlanbprow.py:55(fortran_range)
     8192    0.008    0.000    0.009    0.000 xdrlib.py:194(unpack_float)
      456    0.003    0.000    0.009    0.000 _internal.py:272(_get_void_ptr)
       86    0.008    0.000    0.009    0.000 zlansvdw.py:429(dbdqr)
       76    0.009    0.000    0.009    0.000 zlansvdw.py:508(refinebounds)
     1024    0.002    0.000    0.009    0.000 {method 'all' of 'numpy.ndarray' objects}
     1050    0.008    0.000    0.008    0.000 {numpy.zeros}
     3843    0.007    0.000    0.008    0.000 collections.py:81(__delitem__)
    43994    0.007    0.000    0.007    0.000 {math.copysign}
       60    0.007    0.000    0.007    0.000 {zip}
       10    0.000    0.000    0.007    0.001 hlsvd.py:146(_load_the_library)
     1024    0.001    0.000    0.007    0.000 _methods.py:45(_all)
     3844    0.006    0.000    0.006    0.000 collections.py:71(__setitem__)
       30    0.000    0.000    0.006    0.000 defmatrix.py:117(__new__)
    22073    0.006    0.000    0.006    0.000 {method 'append' of 'list' objects}
    71987    0.006    0.000    0.006    0.000 {abs}
       20    0.003    0.000    0.003    0.000 {nt.chdir}
     3844    0.003    0.000    0.003    0.000 helper.py:313(_prune_cache)
       10    0.003    0.000    0.003    0.000 linalg.py:1182(eig)
     1952    0.003    0.000    0.003    0.000 {method 'astype' of 'numpy.ndarray' objects}
       10    0.000    0.000    0.002    0.000 __init__.py:350(__init__)
     3843    0.002    0.000    0.002    0.000 {method 'pop' of 'list' objects}
       80    0.002    0.000    0.002    0.000 {method 'uniform' of 'mtrand.RandomState' objects}
      456    0.002    0.000    0.002    0.000 _internal.py:260(__array_interface__)
     1014    0.002    0.000    0.002    0.000 blas.py:236(find_best_blas_type)
       10    0.002    0.000    0.002    0.000 {_ctypes.LoadLibrary}
      912    0.002    0.000    0.002    0.000 __init__.py:502(cast)
     2038    0.002    0.000    0.002    0.000 {getattr}
     3854    0.001    0.000    0.001    0.000 {numpy.core._multiarray_umath.normalize_axis_index}
     3844    0.001    0.000    0.001    0.000 fftpack.py:104(_unitary)
      456    0.001    0.000    0.001    0.000 _internal.py:308(data_as)
        1    0.001    0.001    0.001    0.001 fileio.py:124(collapse_complexes)
    12593    0.001    0.000    0.001    0.000 {len}
    20/10    0.001    0.000    0.001    0.000 numeric.py:1401(roll)
     8192    0.001    0.000    0.001    0.000 {_struct.unpack}
     3995    0.001    0.000    0.001    0.000 {method 'pop' of 'dict' objects}
      686    0.001    0.000    0.001    0.000 {_ctypes.pointer}
     6928    0.001    0.000    0.001    0.000 {method 'item' of 'numpy.ndarray' objects}
     1963    0.000    0.000    0.000    0.000 {method 'lower' of 'str' objects}
     2048    0.000    0.000    0.000    0.000 {method 'get' of 'dict' objects}
        1    0.000    0.000    0.000    0.000 {zlib.decompress}
       10    0.000    0.000    0.000    0.000 ntpath.py:483(abspath)
        1    0.000    0.000    0.000    0.000 <string>:122(XML)
        1    0.000    0.000    0.000    0.000 {built-in method feed}
     1156    0.000    0.000    0.000    0.000 {isinstance}
       76    0.000    0.000    0.000    0.000 __init__.py:75(CFUNCTYPE)
      456    0.000    0.000    0.000    0.000 {method 'from_buffer' of '_ctypes.PyCArrayType' objects}
     1488    0.000    0.000    0.000    0.000 {min}
       10    0.000    0.000    0.000    0.000 ntpath.py:415(normpath)
        1    0.000    0.000    0.000    0.000 base64.py:60(b64decode)
        1    0.000    0.000    0.000    0.000 {binascii.a2b_base64}
       20    0.000    0.000    0.000    0.000 ntpath.py:213(dirname)
     1520    0.000    0.000    0.000    0.000 {_ctypes.POINTER}
        1    0.000    0.000    0.000    0.000 {numpy.fromiter}
        1    0.000    0.000    0.000    0.000 {numpy.fft.fftpack_lite.cffti}
       12    0.000    0.000    0.000    0.000 __init__.py:376(__getattr__)
      117    0.000    0.000    0.000    0.000 {method 'reshape' of 'numpy.ndarray' objects}
       10    0.000    0.000    0.000    0.000 linalg.py:144(_commonType)
       20    0.000    0.000    0.000    0.000 ntpath.py:174(split)
      456    0.000    0.000    0.000    0.000 _internal.py:257(__init__)
       10    0.000    0.000    0.000    0.000 linalg.py:215(_assertFinite)
       20    0.000    0.000    0.000    0.000 ntpath.py:63(join)
       70    0.000    0.000    0.000    0.000 ntpath.py:96(splitdrive)
       30    0.000    0.000    0.000    0.000 {_warnings.warn}
       12    0.000    0.000    0.000    0.000 __init__.py:383(__getitem__)
       20    0.000    0.000    0.000    0.000 {method 'view' of 'numpy.ndarray' objects}
       10    0.000    0.000    0.000    0.000 {sorted}
        1    0.000    0.000    0.000    0.000 {open}
       10    0.000    0.000    0.000    0.000 numeric.py:1557(normalize_axis_tuple)
       10    0.000    0.000    0.000    0.000 {nt._getfullpathname}
        1    0.000    0.000    0.000    0.000 {method 'read' of 'file' objects}
       60    0.000    0.000    0.000    0.000 defmatrix.py:169(__array_finalize__)
      684    0.000    0.000    0.000    0.000 {_ctypes.byref}
       10    0.000    0.000    0.000    0.000 defmatrix.py:38(asmatrix)
       10    0.000    0.000    0.000    0.000 {method 'flatten' of 'numpy.ndarray' objects}
       10    0.000    0.000    0.000    0.000 linalg.py:116(_makearray)
       10    0.000    0.000    0.000    0.000 {method 'tolist' of 'numpy.ndarray' objects}
       10    0.000    0.000    0.000    0.000 {numpy.concatenate}
       20    0.000    0.000    0.000    0.000 numeric.py:541(asanyarray)
       10    0.000    0.000    0.000    0.000 {nt.getcwd}
      172    0.000    0.000    0.000    0.000 {method 'upper' of 'str' objects}
        1    0.000    0.000    0.000    0.000 test.py:191(convert_result)
       30    0.000    0.000    0.000    0.000 linalg.py:121(isComplexType)
       10    0.000    0.000    0.000    0.000 linalg.py:111(get_linalg_error_extobj)
       42    0.000    0.000    0.000    0.000 {method 'startswith' of 'str' objects}
       80    0.000    0.000    0.000    0.000 {method 'replace' of 'str' objects}
       40    0.000    0.000    0.000    0.000 zreorth.py:38(zreorth)
        1    0.000    0.000    0.000    0.000 {_elementtree.XMLParser}
       10    0.000    0.000    0.000    0.000 {numpy.empty_like}
       40    0.000    0.000    0.000    0.000 {issubclass}
       10    0.000    0.000    0.000    0.000 linalg.py:209(_assertNdSquareness)
       12    0.000    0.000    0.000    0.000 {method 'split' of 'str' objects}
       10    0.000    0.000    0.000    0.000 linalg.py:134(_realType)
        1    0.000    0.000    0.000    0.000 xdrlib.py:138(__init__)
       10    0.000    0.000    0.000    0.000 {method 'ravel' of 'numpy.ndarray' objects}
       10    0.000    0.000    0.000    0.000 linalg.py:203(_assertRankAtLeast2)
       10    0.000    0.000    0.000    0.000 __init__.py:360(_FuncPtr)
        2    0.000    0.000    0.000    0.000 constants.py:467(any_type_to_internal)
       10    0.000    0.000    0.000    0.000 linalg.py:137(_complexType)
       10    0.000    0.000    0.000    0.000 {method 'transpose' of 'numpy.ndarray' objects}
       10    0.000    0.000    0.000    0.000 numeric.py:1471(<dictcomp>)
       10    0.000    0.000    0.000    0.000 {method 'join' of 'str' objects}
       20    0.000    0.000    0.000    0.000 {built-in method get}
        1    0.000    0.000    0.000    0.000 {built-in method find}
       12    0.000    0.000    0.000    0.000 {setattr}
       10    0.000    0.000    0.000    0.000 {method 'lstrip' of 'str' objects}
       10    0.000    0.000    0.000    0.000 {method 'items' of 'dict' objects}
        1    0.000    0.000    0.000    0.000 constants.py:485(any_type_to_numpy)
        1    0.000    0.000    0.000    0.000 xdrlib.py:154(done)
       10    0.000    0.000    0.000    0.000 {operator.index}
       10    0.000    0.000    0.000    0.000 {method '__array_prepare__' of 'numpy.ndarray' objects}
        2    0.000    0.000    0.000    0.000 constants.py:462(is_complex)
        1    0.000    0.000    0.000    0.000 xdrlib.py:141(reset)
        1    0.000    0.000    0.000    0.000 __init__.py:101(CFunctionType)
        1    0.000    0.000    0.000    0.000 {built-in method close}
        1    0.000    0.000    0.000    0.000 {built-in method getchildren}
        1    0.000    0.000    0.000    0.000 {iter}
        1    0.000    0.000    0.000    0.000 {method 'strip' of 'str' objects}
        1    0.000    0.000    0.000    0.000 {hasattr}
        1    0.000    0.000    0.000    0.000 {method 'reverse' of 'list' objects}
        1    0.000    0.000    0.000    0.000 {method 'disable' of '_lsprof.Profiler' objects}


