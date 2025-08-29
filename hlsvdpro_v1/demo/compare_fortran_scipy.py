# Python modules
from __future__ import division
import json
import sys

# 3rd party modules
import numpy as np
import matplotlib.pyplot as pyplot

# This project's modules
import hlsvdpro
import hlsvdpro_scipy

"""This performs a simple test of HLSVD to see that it's working and
producing sane results. It should show a plot of the observed values
in red and the values estimated by HLSVDPRO in blue.
"""

def tests(data_filename=None):
    with open(data_filename) as f:
        input_ = json.load(f)

    # Signals (input data) are stored as correlated arrays of real/imaginary
    # pairs. Here I knit them into a single complex list.
    input_['signals'] = [complex(r, i) for r, i in zip(input_['signals_r'],
                                                       input_['signals_i'])]
    del input_['signals_r']
    del input_['signals_i']

    signals     = np.asarray(input_["signals"], complex)
    observed    = signals.copy()
    step_size   = input_["step_size"]
    nsv_sought  = input_["n_singular_values"]
    hankel_size = 256 #int(len(signals) // 8)
    
    print( 'hankel_size // 8 = ' + str(hankel_size) )
    print( 'nsv_sought       = ' + str(nsv_sought) )

    reps = 1  #4

    for i in range(reps):
        # Call fortran hlsvd()
        results1 = hlsvdpro.hlsvd(observed, nsv_sought, step_size)
        
    for i in range(reps):
        results2 = hlsvdpro_scipy.hlsvdpro_scipy(observed, nsv_sought, hankel_size)

    nsv_found1, singular_values1, frequencies1, \
        damping_factors1, amplitudes1, phases1 = results1

    nsv_found2, singular_values2, frequencies2, \
        damping_factors2, amplitudes2, phases2 = results2

    estimated1 = _create_hlsvd_sum_fid(frequencies1, damping_factors1, amplitudes1,
                                      phases1, len(observed), step_size)

    estimated2 = _create_hlsvd_sum_fid(frequencies2, damping_factors2, amplitudes2,
                                      phases2, len(observed), step_size)

    # We're only concerned with the real portion of the data
    observed = np.real(np.array(observed))
    estimated1 = np.real(estimated1)
    estimated2 = np.real(estimated2)

    rmse1 = np.sqrt(np.sum((observed - estimated1) ** 2) / len(observed))
    nrmse1 = rmse1 / max(observed) - min(observed)

    rmse2 = np.sqrt(np.sum((observed - estimated2) ** 2) / len(observed))
    nrmse2 = rmse2 / max(observed) - min(observed)

    print("RMSE1 = %.2f, normalized RMSE1 = %.2f%%" % (rmse1, nrmse1 / 100))
    print("RMSE2 = %.2f, normalized RMSE2 = %.2f%%" % (rmse2, nrmse2 / 100))

    # The observed values are always noisier, so plotting them first allows
    # the cleaner estimation to be displayed on top.
    pyplot.plot(observed, color='r')
    pyplot.plot(estimated1, color='b')
    pyplot.plot(estimated2, color='g')
    pyplot.show()




# This function written by Brian Soher for Vespa:
# https://scion.duhs.duke.edu/vespa/
# Modified by me (Philip Semanchuk); any errors are mine.
def _create_hlsvd_sum_fid(frequencies, decays, areas, phases,
                          ndata_points, dwell_time):
    """Construct time domain signal from the estimated parameters."""
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

            line_result = areas[i] * np.exp((t/decays[i]) +
                          K * (frequencies[i]*t+phases[i]/360.0))

            zeros = np.zeros_like(line_result)

            result[i,:] = np.where(np.isnan(line_result), zeros, line_result)

            np.seterr(**old_settings)
        else:
            result[i,:] = result[i,:] * 0

    return np.sum(result, axis=0)


def run_tests():
    data_filename = 'press_cp0.json'
    tests(data_filename)



if __name__ == "__main__":
    run_tests()

#    import cProfile
#    cProfile.run('run_tests()')
 
    # to time this call user
    #>python -m cProfile -s cumulative compare_fortran_scipy.py