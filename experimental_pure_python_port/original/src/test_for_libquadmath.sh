#!/bin/bash

# This script tests to see if the dylib passed on the command line (usually
# libhlsvd.dylib) depends on libquadmath.dylib. If it does, the script
# prints a warning/reminder to stdout.

# Created by the Vespa team
# Copyright & license info for this script (not HLSVD itself) is in the 
# Vespa LICENSE file.


GREP_OUT=`otool -L $1|grep libquadmath`

if [ -z "$GREP_OUT" ] ; then
    :  # No-op. All is well.
else
    echo 
    echo "      *************   SOMETHING SMELLS FISHY   *************"
    echo "$1 requires libquadmath.dylib at runtime. This is"
    echo "probably not what you want. Did you remember to hide libquadmath.dylib"
    echo "from the linker?"
    echo "See here for more info:"
    echo "http://scion.duhs.duke.edu/vespa/project/wiki/HlsvdBuilding"
    echo
fi
