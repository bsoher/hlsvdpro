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
import hlsvdpro_local
import hlsvdpro_propack
import hlsvdpro


SHOW_PLOT = True
#SHOW_PLOT = True
DTOR = np.pi / 180.0

"""This performs a simple test of HLSVD to see that it's working and
producing sane results.
"""



def do_test(signals, k, sing0, U0, V0, reps=1):

    observed    = signals.copy()
    step_size   = 0.001
    nsv_sought  = k
    hankel_size = len(observed) // 1            # or hard fix to 512
    dwell_time  = float(step_size)

    observed   = signals[0:hankel_size].copy()

    print 'hankel_size len() // 1 = ' + str(hankel_size)
    print 'nsv_sought             = ' + str(nsv_sought)


    reps1 = [0,0,0]     # good_fit, good_wrong (right # singval, wrong values, bad_fit (wrong # singval)
    reps3 = [0,0,0]
    reps4 = [0,0,0]

    for i in range(reps):
        result1 = hlsvdpro_propack.hlsvdpro_propack(observed, nsv_sought, 4)
        if result1[0] == nsv_sought:
            if (max(result1[1]-sing0) > 0.01):
                reps1[1] += 1
            else:
                reps1[0] += 1
        else:
            reps1[2] += 1

    # for i in range(reps):
    #     result3 = hlsvdpro_local.hlsvdpro(observed, nsv_sought, step_size)
    #     if result3[0] == nsv_sought:
    #         if (max(result3[1]-sing0) > 0.01):
    #             reps3[1] += 1
    #         else:
    #             reps3[0] += 1
    #     else:
    #         reps3[2] += 1
    #
    #
    # for i in range(reps):
    #     # this is fortran call in site-packages module - NO param convert needed
    #     result4 = hlsvdpro.hlsvd(observed, nsv_sought, step_size)
    #     if result4[0] == nsv_sought:
    #         if (max(result4[1]-sing0) > 0.01):
    #             reps4[1] += 1
    #         else:
    #             reps4[0] += 1
    #     else:
    #         reps4[2] += 1

    # nsing1, sing1 = result1[0:2]
    # nsing4, sing4 = result4[0:2]

    u1 = result1[6]
    v1 = result1[7]
    utmp = abs(np.dot(u1.conjugate().T,u1))
    utmp[utmp<1e-6] = 0.0
    vtmp = abs(np.dot(v1, v1.conjugate().T))
    vtmp[vtmp<1e-6] = 0.0

    print(' Results hlsvdpro_propack, max(abs(sing1-sing0)) = ', max(abs(result1[1]-sing0)))
    print('     np.clip(abs(np.dot(u1.conjugate().T, u1))) = ', utmp)
    print('     np.clip(abs(np.dot(v1, v1.conjugate().T))) = ', vtmp)
    print(' Results hlsvdpro_propack, good_fit, wrong_fit, bad_fit = ',reps1[0], reps1[1], reps1[2])




    # print(' Results hlsvdpro_local  , max(sing3-sing0) = ', max(result3[1]-sing0))
    # print(' Results hlsvdpro_local  , good_fit, wrong_fit, bad_fit = ',reps3[0], reps3[1], reps3[2])
    # print(' Results hlsvdpro,hlsvd  , max(sing4-sing0) = ', max(result4[1]-sing0))
    # print(' Results hlsvdpro.hlsvd  , good_fit, wrong_fit, bad_fit = ',reps4[0], reps4[1], reps4[2])


    bob = 1
    bob += 1




def test():

    m = 20
    n = 10

    A = [0.54881350 + 0.31179588j, 0.71518937 + 0.69634349j, 0.60276338 + 0.37775184j, 0.54488318 + 0.17960368j, 0.42365480 + 0.02467873j, 0.64589411 + 0.06724963j, 0.43758721 + 0.67939277j, 0.89177300 + 0.45369684j, 0.96366276 + 0.53657921j, 0.38344152 + 0.89667129j]

    # A = np.random.random((m, n)) + 1j * np.random.random((m, n))
    A = np.array(A)
    A = A.astype(np.complex128)
    flavorlabel = 'COMPLEX16'

#    sing0 = np.array([10.21100239,  2.55836133,  2.33302783,  2.14057459,  1.850737  ])
#   hlsvdpro_local         = [3.844931797977525, 1.0606943021737258, 0.8539664440754421, 0.6082379328434804, 0.47619445709816216]
#   hlsvd original fortran   [3.844931797977525, 1.0606943021737258, 0.8539664440754421, 0.6082379328434804, 0.47619445709816216]
#   hlsvdpro_propack         [3.8449318          1.0606943           0.85396644          0.60823793          0.47619446]

    sing0 = np.array([3.844931797977525, 1.0606943021737258, 0.8539664440754421, 0.6082379328434804, 0.47619445709816216], dtype=float)

    k = 5

    # u, sigma, v = np.linalg.svd(A, full_matrices=False)

    do_test(A, k, sing0, U0, V0, reps=1)


    # # print the results
    # np.set_printoptions(suppress=True, precision=8)
    #
    # print(A)
    # print('m = ',m,'  n = ',n, 'k = ',k)
    #
    # print('Flavor Label = ' + flavorlabel)
    # print('SingVals sing0      = ', sing0)
    # print('SingVals do_test    = ', sigma1)
    # print('SingVal Diffs (sing0 vs svdp  ) = ', sing0 - sigma1)
    print('')




if __name__ == "__main__":

    test()
