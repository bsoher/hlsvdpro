# Python modules
from __future__ import division
import sys
import xml.etree.cElementTree as ElementTree

# 3rd party modules
import numpy as np

# Vespa modules
import vespa.common.util.xml_ as util_xml
import vespa.hlsvd.dist.hlsvd as hlsvd

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

    observed = input_["signals"]
    step_size = input_["step_size"]

    # Call hlsvd()
    results = hlsvd.hlsvd(observed, input_["n_singular_values"], step_size)


    nsv_found, singular_values, frequencies, \
        damping_factors, amplitudes, phases = results

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
        print "Please supply the name of the XML data file."
    else:
        test(sys.argv[1])


