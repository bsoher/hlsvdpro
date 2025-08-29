# Makefile for HLSVD
# Created by Brian J. Soher
# Copyright & license info for this Makefile (not HLSVD itself) is in the
# Vespa LICENSE file.

# BITS can be 32 or 64 and controls the bit mode of the binary.
BITS = 64
ARCH = -m$(BITS)
OPTIMIZE = -O
DISABLED_WARNINGS = -Wno-unused
WARNINGS = -Wall $(DISABLED_WARNINGS)
FFLAGS = $(OPTIMIZE) $(WARNINGS) $(ARCH)
# When building LAPACK we force ffpe-trap to be blank (disabled) because
# LAPACK's make file says, "LAPACK is designed to check for and handle these
# cases internally and enabling these traps will likely cause LAPACK to crash."
# ref: https://icl.cs.utk.edu/svn/lapack-dev/lapack/trunk/CMAKE/CheckLAPACKCompilerFlags.cmake
LAPACK_FFLAGS =    -ffpe-trap=    $(FFLAGS)

# If GNU msys is installed, make's default shell will (might?) become sh.exe
# which is Unix-y. That doesn't play well with the DOS-specific commands for
# recreating BIN_DIR.
SHELL=cmd

# These steps will locate the first python installation in your PATH variable.
# I use conda and activate the python that I want before doing my compile, and
# conda ensures that the right path is first in my PATH.
PYDIR := $(shell python -c "from __future__ import print_function;import sys; print(sys.executable)")
PYDIR := $(subst python.exe,,$(PYDIR))
PYDIR := $(subst \,/,$(PYDIR))
PYMIN := $(shell python -c "from __future__ import print_function;import sys; print(sys.version_info.minor)")
PYMAJ := $(shell python -c "from __future__ import print_function;import sys; print(sys.version_info.major)")
PYVER := $(PYMAJ).$(PYMIN)
PYVER2 := $(PYMAJ)$(PYMIN)


F2PY_IFLAGS = -I$(PYDIR)lib/site-packages/numpy/core/include -I$(PYDIR)include -I$(PYDIR)lib/site-packages/numpy/_core/include
F2PY_CFLAGS = -g -DDEBUG -DMS_WIN64 -O0 -Wall -Wstrict-prototypes $(F2PY_IFLAGS)
F2PY_FFLAGS = -Wall -g -ffixed-form -fno-second-underscore -O3 -funroll-loops $(F2PY_IFLAGS)


BIN_DIR = ./bin
DIST_TARGET_DIR = ../../hlsvdpro/
TARGET_LIB_NAME = $(BIN_DIR)/_spropack.pyd


all: $(TARGET_LIB_NAME)


# Create the BIN_DIR if it doesn't exist
$(BIN_DIR): force
	if NOT EXIST $(BIN_DIR) mkdir $(subst /,\,$(BIN_DIR))

# Pattern rules for turning f2py_Util src files into object files
$(BIN_DIR)/%.o: f2py_Util/%.c
	gcc $(F2PY_CFLAGS) -c $< -o $@

$(BIN_DIR)/%.o: f2py_Util/%.f
	gcc $(F2PY_FFLAGS) -c $< -o $@

# Pattern rule for turning Lapack_Util/*.f into object files
$(BIN_DIR)/%.o: Lapack_Util/%.f
	gfortran $(LAPACK_FFLAGS) -c $< -o $@
    
# Pattern rule for turning HLSVD *.f into object files
$(BIN_DIR)/%.o: %.F
	gfortran $(FFLAGS) -c $< -o $@


# Variables to define lapack object files in the BIN_DIR
LAPACK_SRCS := $(wildcard Lapack_Util/*.f)
LAPACK_OBJS := $(addsuffix .o,$(basename $(notdir $(LAPACK_SRCS))))
LAPACK_OBJS := $(addprefix $(BIN_DIR)/,$(LAPACK_OBJS))

F2PY_SRCS := $(wildcard f2py_Util/*.f)
F2PY_OBJS := $(addsuffix .o,$(basename $(notdir $(F2PY_SRCS))))
F2PY_OBJS := $(addprefix $(BIN_DIR)/,$(F2PY_OBJS))

F2PYC_SRCS := $(wildcard f2py_Util/*.c)
F2PYC_OBJS := $(addsuffix .o,$(basename $(notdir $(F2PYC_SRCS))))
F2PYC_OBJS := $(addprefix $(BIN_DIR)/,$(F2PYC_OBJS))


# Rule for building lapack.a
$(BIN_DIR)/lapack.a: $(BIN_DIR) $(LAPACK_OBJS)
	ar -r $(BIN_DIR)/lapack.a $(LAPACK_OBJS)

lapack: $(BIN_DIR)/lapack.a

# SRCS lists the source (Fortran) files for the HLSVD code.
# OBJS are the compiled versions of these files that land in the BIN_DIR
SRCS := _spropack-f2pywrappers.f _spropackmodule.c fortranobject.c \
        printstat.F saprod.F sblasext.F sbsvd.F second.F sgemm_ovwr.F \
        sgetu0.F slanbpro.F slansvd.F slansvd_aprod.F slansvd_irl.F \
        slansvd_irl_aprod.F smgs.pentium.F sreorth.F sritzvec.F \
        ssafescal.F 

OBJS := $(addprefix $(BIN_DIR)/, $(addsuffix .o,$(basename $(SRCS))))


# This builds a target that statically links in GFortran's libraries so
# we don't have to require mingw32 gfortran be installed for this to work.
$(TARGET_LIB_NAME): $(BIN_DIR) $(OBJS) lapack 
	gfortran -static-libgcc 	                    \
			-shared  								\
			-Wall $(ARCH)                          	\
			-o $(TARGET_LIB_NAME)                  	\
			$(OBJS)				  					\
			-L$(PYDIR)libs							\
			-lpython$(PYVER2)	    				\
			$(BIN_DIR)/lapack.a						\
			-Wl,-Bstatic							\
			-lgfortran



## In Vespa we only use the library, not an executable. But the executable
## is useful for stepping through with gdb so I've used it occasionally.
## See hlsvdmain in hlsvdpro.f for details.
#exe: $(BIN_DIR) $(OBJS) lapack
#	gfortran -static-libgcc -static-libgfortran   	\
#	         -Wall $(ARCH)                      	\
#	         -o hlsvdpro.exe						\
#	    	 -Wl,-Bdynamic							\
#	    	 -L.									\
#	    	 -lfftw3-3 								\
#	    	 -Wl,-static 							\
#	   	 	 $(OBJS)								\
#		 	 $(BIN_DIR)/lapack.a


# depcheck lists the runtime libraries that the target relies on.
# It relies on grep which isn't on Windows by default.
depcheck: $(TARGET_LIB_NAME)
	objdump -p $(TARGET_LIB_NAME)|grep DLL


# Make passes these commands to cmd.exe, so they must be DOS-compatible
# (backslashes and all).
dist: 
	copy $(subst /,\,$(TARGET_LIB_NAME)) $(subst /,\,$(DIST_TARGET_DIR))
	copy $(subst /,\,$(TARGET_LIB_NAME)) $(subst /,\,$(DIST_TARGET_DIR))_win_py$(PYVER2)

clean:
	del $(subst /,\,$(BIN_DIR))\**


force:

