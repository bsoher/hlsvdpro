#!/usr/bin/env python

# Python modules
from __future__ import division
import sys
import pstats

p = pstats.Stats("profile.data")

p.sort_stats('cumulative').print_stats(30)


