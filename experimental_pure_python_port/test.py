#!/usr/bin/env python

# Python modules
from __future__ import division
import sys
import xml.etree.cElementTree as ElementTree
import cProfile
import os
import pprint
from __builtin__ import False
from pickle import FALSE
pp=pprint.pprint

# 3rd party modules
import numpy as np

# Vespa modules
import vespa.common.util.xml_ as util_xml
import vespa.common.constants as vespa_constants
import hlsvdpro_local
import hlsvdpro_bjs
import hlsvdpro_scipy
import hlsvdpro_propack

import hlsvdpro


#SHOW_PLOT = False
SHOW_PLOT = True
DTOR = np.pi / 180.0

"""
This performs a simple test of HLSVD to see that it's working and
producing sane results.
"""



def test(filename):
    s = open(filename, "rb").read()
    root = ElementTree.fromstring(s)

    # The XML is divided into two sections, one for input and the other
    # for (expected) output.
    input_ = util_xml.element_to_dict(root.find("input"))

    signals     = input_["signals"]
    observed    = signals.copy()
    step_size   = input_["step_size"]
    nsv_sought  = input_["n_singular_values"]
    hankel_size = len(observed) // 8            # or hard fix to 512
    dwell_time  = float(step_size)

    observed   = signals[0:hankel_size].copy()

    print 'hankel_size len() // 8 = ' + str(hankel_size)
    print 'nsv_sought             = ' + str(nsv_sought)

    reps = 1         # 4  for time  trial

    result1 = []
    result2 = []
    result3 = []
    result4 = []

    for i in range(reps):
        result1 = hlsvdpro_propack.hlsvdpro_propack(observed, nsv_sought, 256)
# #        result1 = hlsvdpro_scipy.hlsvdpro_scipy(signals, nsv_sought, hankel_size)
    result1 = convert_result(result1, dwell_time)
       
    # for i in range(reps):
    #     result2 = hlsvdpro_bjs.hlsvdpro_bjs(observed, nsv_sought, 512)
    # result2 = convert_result(result2, dwell_time)
         
    # for i in range(reps):
    #     result3 = hlsvdpro_local.hlsvdpro(observed, nsv_sought, step_size)
    # result3 = convert_result(result3, dwell_time)

    for i in range(reps):
        # this is fortran call in site-packages module - NO param convert needed
        result4 = hlsvdpro.hlsvd(observed, nsv_sought, step_size)
    # result4 converted within the Fortran call

    if result1: pp_result_compare(result1, result4, 1, 4)
    if result2: pp_result_compare(result2, result4, 2, 4)
    if result1 and result2: pp_result_compare(result1, result2, 1, 2)
    if result3: pp_result_compare(result3, result4, 3, 4)


    if SHOW_PLOT:

        estimateds = []
        for indx,res in enumerate([result1, result2, result3, result4]):
    
            if res:
                nsv_found, singular_values, frequencies, damping_factors, amplitudes, phases = res
        
                estimated = _create_hlsvd_sum_fid(frequencies, damping_factors, amplitudes, 
                                                  phases, len(observed), step_size)
                estimateds.append(estimated)
        
                # We're only concerned with the real portion of the data
                obs = np.real(np.array(observed))
                est = np.real(estimated)
            
                rmse = np.sqrt(np.sum((obs - est) ** 2) / len(obs))
                nrmse = rmse / max(obs) - min(obs)
            
                print str(indx)+" - RMSE = %.2f, normalized RMSE = %.2f%%" % (rmse, nrmse / 100)

        import matplotlib.pyplot as plt

        nest = len(estimateds)

        fig, axes = plt.subplots(nrows=nest+1, ncols=2)

        phase = np.exp(1j*DTOR*(-110))     # was -110 for cp0 -200 for cp5
        obst  = np.array(observed)
        obsf  = np.fft.fft(_chop(obst*phase))

        for i,est in enumerate(estimateds):
            
            estt = np.array(est)
            estf = np.fft.fft(_chop(estt*phase))

            # The observed values are always noisier, so plotting them first allows
            # the cleaner estimation to be displayed on top. 
            axes[i,0]
            axes[i,0].plot(obst, color='b')
            axes[i,0].plot(estt, color='r')
            axes[i,0].plot(obst-estt, color='g')
    
            axes[i,1].plot(obsf, color='b')
            axes[i,1].plot(estf, color='r')
            axes[i,1].plot(obsf-estf, color='g')

        esta = estimateds[-2]
        estb = estimateds[-1] 
        estt = esta - estb
        estaf = np.fft.fft(_chop(esta*phase))
        estbf = np.fft.fft(_chop(estb*phase))
        estf = estaf - estbf

        axes[nest,0].plot(esta, color='b')
        axes[nest,0].plot(estb, color='r')
        axes[nest,0].plot(estt, color='g')
        axes[nest,1].plot(estaf, color='b')
        axes[nest,1].plot(estbf, color='r')
        axes[nest,1].plot(estf,  color='g')
            
        plt.show()
        bob = 1
        bob += 1




# swiped from functors/funct_hlsvd_create_fids.py, modified for this code
def _create_hlsvd_sum_fid(frequencies, decays, areas, phases, 
                           ndata_points, dwell_time):

    """ 
    Construct time domain signal from the estimated parameters.
    
    """
    result = np.zeros((len(frequencies), ndata_points), dtype=complex)
    t = np.arange(ndata_points) * dwell_time 
    
    # K is an invariant for the loop below
    K = 1j * 2 * np.pi

    for i, decay in enumerate(decays):
        if decays[i]:
            # We frequently force the exp() function here into very small 
            # values that raise the following exception:
            #    FloatingPointError: underflow encountered in exp
            # Handling those exceptions with try/except is not an option 
            # because raising the exception interrupts the calculation.
            # Instead, we temporarily disable numpy's error reporting,
            # allow it to silently generate NaNs, and then change those NaNs
            # to 0.0.
            old_settings = np.seterr(all='ignore')
            
            line_result = areas[i] * np.exp((t/decays[i]) + K * (frequencies[i]*t+phases[i]/360.0))
            
            zeros = np.zeros_like(line_result)

            result[i,:] = np.where(np.isnan(line_result), zeros, line_result)

            np.seterr(**old_settings)
        else:
            result[i,:] = result[i,:] * 0

    return np.sum(result, axis=0)


def convert_result(result, dwell_time):

    nsv_found, singular_values, frequencies, damping_factors, amplitudes, phases, _, _ = result
    
    damping_factors = [1 / df for df in damping_factors]
    damping_factors = [df * dwell_time for df in damping_factors]
 
    frequencies = [frequency / dwell_time for frequency in frequencies]
 
    phases = [phase * vespa_constants.RADIANS_TO_DEGREES for phase in phases]   
    
    return  nsv_found, singular_values, frequencies, damping_factors, amplitudes, phases


def _chop(data):
    # numpy broadcasting deals with (N,dim0) sized data arrays
    chop_array = ((((np.arange(len(data)) + 1) % 2) * 2) - 1)
    return data * chop_array


def pp_result_compare(res1, res2, id1, id2):
    
    labl = ['nsv_found', 'singular_values', 'frequencies', 'damping_terms', 'amplitudes', 'phases']

    # NB. The following comparisons come out in different order, thus the need to 
    #     do the np.argsort() gymnastics.
    
    _, _, frequencies1, _, _, _ = res1
    _, _, frequencies2, _, _, _ = res2
    
    indices1 = np.argsort(frequencies1)
    indices2 = np.argsort(frequencies2)
    
    for labl, item1, item2 in zip(labl, res1, res2):
        
        if labl == 'nsv_found':
            print('')
            print(labl, str(item1), str(item2))
        elif labl == 'singular_values':
            print('')
            pp(labl+str(id1))
            pp(np.array(item1))
            pp(labl+str(id2))
            pp(np.array(item2))
            pp(labl+ ' delta')
            pp(np.array([val1-val2 for val1,val2 in zip(np.array(item1),np.array(item2))]))    
        else:
            print('')
            # pp(labl+str(id1))
            # pp(np.array(item1)[indices1])
            # pp(labl+str(id2))
            # pp(np.array(item2)[indices2])
            # pp(labl+ ' delta')
            # pp(np.array([val1-val2 for val1,val2 in zip(np.array(item1)[indices1],np.array(item2)[indices2])]))


 


if __name__ == "__main__":
    if len(sys.argv) < 2:
        #filename = "original/dist/25215695.xml"
        #filename = "original/dist/laser_2010102901.xml"
        filename = "original/dist/press_cp0.xml"
        #filename = "original/dist/press_cp5.xml"
        # print "Please supply the name of the XML data file."
    else:
        filename = sys.argv[1]

    if os.path.exists("profile.data"):
        os.remove("profile.data")

    test(filename)
     
    # cProfile.run('test(filename)', 'profile.data')
    # import pstats as ps
    # p = ps.Stats('profile.data')
    # p.strip_dirs().sort_stats('cumulative').print_stats()


