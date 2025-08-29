# Python modules
from __future__ import division

# Vespa modules
import vespa.common.util.dir_list as util_dir_list 

# Local modules
import test

test.SHOW_PLOT = False

filenames = util_dir_list.DirList('.', util_dir_list.FilterEndsWith('xml')).files


for filename in filenames:
    print "\nTesting %s..." % filename
    test.test(filename)
