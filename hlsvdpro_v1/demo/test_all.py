# Python modules
from __future__ import division
import os

# Local modules
import demo

filenames = [filename for filename in os.listdir('.') if filename.endswith('.json')]


for filename in filenames:
    print "\nTesting %s..." % filename
    demo.main(filename)
