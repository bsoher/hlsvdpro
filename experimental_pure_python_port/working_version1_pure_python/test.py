#!/usr/bin/env python

# Python modules
from __future__ import division
import sys
import xml.etree.cElementTree as ElementTree
import cProfile
import os
import pprint
pp=pprint.pprint

# 3rd party modules
import numpy as np

# Vespa modules
import vespa.common.util.xml_ as util_xml
import vespa.common.constants as vespa_constants
import hlsvdpro_local
import hlsvdpro_bjs
import hlsvdpro_scipy

import hlsvdpro


SHOW_PLOT = True


"""This performs a simple test of HLSVD to see that it's working and
producing sane results.
"""



def test(filename):
    s = open(filename, "rb").read()
    root = ElementTree.fromstring(s)

    # The XML is divided into two sections, one for input and the other
    # for (expected) output.
    input_ = util_xml.element_to_dict(root.find("input"))

    observed    = input_["signals"]
    step_size   = input_["step_size"]
    nsv_sought  = input_["n_singular_values"]
    hankel_size = len(observed) // 8

    print 'hankel_size // 8 = ' + str(hankel_size)
    print 'nsv_sought       = ' + str(nsv_sought)

    reps = 1        # 4  for time  trial

#     for i in range(reps):
#         results1 = hlsvdpro_scipy.hlsvdpro_scipy(observed, nsv_sought, hankel_size)
#  
#     for i in range(reps):
#         results2 = hlsvdpro_bjs.hlsvdpro_bjs(observed, input_["n_singular_values"], 512)
         
    for i in range(reps):
        results = hlsvdpro_local.hlsvdpro(observed, nsv_sought, step_size)
 
    for i in range(reps):
        results3 = hlsvdpro.hlsvd(observed, nsv_sought, step_size)


#    return

    nsv_found, singular_values, frequencies, damping_factors, amplitudes, phases = results



    dwell_time = step_size

    damping_factors = [1 / df for df in damping_factors]
    damping_factors = [df * dwell_time for df in damping_factors]
 
    frequencies = [frequency / dwell_time for frequency in frequencies]
 
    phases = [phase * vespa_constants.RADIANS_TO_DEGREES for phase in phases]

    # print nsv_found
    # pp(singular_values)
    # pp(frequencies)
    # pp(damping_factors)
    # pp(amplitudes)
    # pp(phases)

    nsv_found3, singular_values3, frequencies3, damping_factors3, amplitudes3, phases3 = results3
    
    # NB. The following comparisons come out in different order, thus the need to 
    #     do the np.argsort() gymnastics.
    
    print('Differences between local and fortran hlsvd module results')
    print('')
    print('nsv_found',str(nsv_found), str(nsv_found3))
    print('singular_values ')
    print(str(singular_values))
    print('singular_values3 ')
    print(str(singular_values3))
    print('singular_values delta')
    print(str([item1-item2 for item1,item2 in zip(singular_values,singular_values3)]))
    
    indices  = np.argsort(frequencies)
    indices3 = np.argsort(frequencies3)
    
    print('')
    print('frequencies ')
    print(np.array(frequencies)[indices])
    print('frequencies3 ')
    print(np.array(frequencies3)[indices3])
    print('frequencies delta ')
    print(str([item1-item2 for item1,item2 in zip(np.array(frequencies)[indices],np.array(frequencies3)[indices3])]))
    
    print('')
    print('damping_factors ')
    print(np.array(damping_factors)[indices])
    print('damping_factors3 ')
    print(np.array(damping_factors3)[indices3])
    print('damping_factors delta ')
    print(str([item1-item2 for item1,item2 in zip(np.array(damping_factors)[indices],np.array(damping_factors3)[indices3])]))

    print('')
    print('amplitudes ')
    print(np.array(amplitudes)[indices])
    print('amplitudes3 ')
    print(np.array(amplitudes3)[indices3])
    print('amplitudes delta ')
    print(str([item1-item2 for item1,item2 in zip(np.array(amplitudes)[indices],np.array(amplitudes3)[indices3])]))
    
    print('')
    print('phases ')
    print(np.array(phases)[indices])
    print('phases3 ')
    print(np.array(phases3)[indices3])
    print('phases delta ')
    print(str([item1-item2 for item1,item2 in zip(np.array(phases)[indices],np.array(phases3)[indices3])]))


    estimated = _create_hlsvd_sum_fid(frequencies, damping_factors, amplitudes, 
                                      phases, len(observed), step_size)

    # We're only concerned with the real portion of the data
    observed = np.real(np.array(observed))
    estimated = np.real(estimated)

    rmse = np.sqrt(np.sum((observed - estimated) ** 2) / len(observed))

    nrmse = rmse / max(observed) - min(observed)

    print "RMSE = %.2f, normalized RMSE = %.2f%%" % (rmse, nrmse / 100)

    if SHOW_PLOT:
        import matplotlib.pyplot as pyplot


        # The observed values are always noisier, so plotting them first allows
        # the cleaner estimation to be displayed on top. 
        pyplot.plot(observed, color='r')
        pyplot.plot(estimated, color='b')
        pyplot.plot(observed-estimated, color='g')
        pyplot.show()


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
            
            line_result = areas[i] * np.exp((t/decays[i]) + \
                          K * (frequencies[i]*t+phases[i]/360.0))
            
            zeros = np.zeros_like(line_result)

            result[i,:] = np.where(np.isnan(line_result), zeros, line_result)

            np.seterr(**old_settings)
        else:
            result[i,:] = result[i,:] * 0

    return np.sum(result, axis=0)
 


if __name__ == "__main__":
    if len(sys.argv) < 2:
        #filename = "original/dist/25215695.xml"
        filename = "original/dist/press_cp0.xml"
        # print "Please supply the name of the XML data file."
    else:
        filename = sys.argv[1]

    if os.path.exists("profile.data"):
        os.remove("profile.data")

    test(filename)
    
#     cProfile.run('test(filename)', 'profile.data')
#     import pstats as ps
#     p = ps.Stats('profile.data')
#     p.strip_dirs().sort_stats('cumulative').print_stats()


