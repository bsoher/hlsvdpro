# Makefile for HLSVD
# Created by the Vespa team
# Copyright & license info for this Makefile (not HLSVD itself) is in the
# Vespa LICENSE file.

# BITS can be 32 or 64 and controls the bit mode of the binary.
BITS = 32
ARCH = -m$(BITS)
OPTIMIZE = -O
DISABLED_WARNINGS = -Wno-unused
WARNINGS = -Wall $(DISABLED_WARNINGS)
FFLAGS = $(OPTIMIZE) $(WARNINGS) $(ARCH) 
# When building LAPACK we force ffpe-trap to be blank (disabled) because 
# LAPACK's make file says, "LAPACK is designed to check for and handle these 
# cases internally and enabling these traps will likely cause LAPACK to crash."
# ref: https://icl.cs.utk.edu/svn/lapack-dev/lapack/trunk/CMAKE/CheckLAPACKCompilerFlags.cmake
LAPACK_FFLAGS =   -ffpe-trap=  $(FFLAGS)


# BIN_DIR is OS- and architecture-dependent
BIN_DIR := ./bin/osx-$(BITS)

TARGET_LIB_NAME = $(BIN_DIR)/libhlsvd-$(BITS).dylib

TARGET_DIR = ../dist/$(BIN_DIR)


all: $(TARGET_LIB_NAME)


# Create the BIN_DIR if it doesn't exist
$(BIN_DIR): force
	@if (test ! -d $(BIN_DIR)) then \
        mkdir -p $(BIN_DIR); \
	fi


# Pattern rule for turning lapack/*.f into object files
$(BIN_DIR)/%.o: lapack/%.f
	gfortran $(LAPACK_FFLAGS) -c $< -o $@


# Pattern rule for turning HLSVD *.f into object files
$(BIN_DIR)/%.o: %.f
	gfortran $(FFLAGS) -c $< -o $@


# Variables to define lapack object files in the BIN_DIR
LAPACK_SRCS := $(wildcard lapack/*.f)
LAPACK_OBJS := $(addsuffix .o,$(basename $(notdir $(LAPACK_SRCS))))
LAPACK_OBJS := $(addprefix $(BIN_DIR)/,$(LAPACK_OBJS))


# Rule for building lapack.a
$(BIN_DIR)/lapack.a: $(BIN_DIR) $(LAPACK_OBJS) 
	ar -r $(BIN_DIR)/lapack.a $(LAPACK_OBJS)

lapack: $(BIN_DIR)/lapack.a


# SRCS lists the source (Fortran) files for the HLSVD code.
# OBJS are the compiled versions of these files that land in the BIN_DIR
SRCS := aprodw.f blasext.f hlsvdpro.f lanczopw.f printstat.f \
		zsafescal.f zgemm_ovwr.f zget0w.f zlanbprow.f zlansvdw.f zreorth.f

OBJS := $(addprefix $(BIN_DIR)/, $(addsuffix .o,$(basename $(SRCS))))

$(TARGET_LIB_NAME): $(BIN_DIR) $(OBJS) lapack
	gfortran -static-libgcc -static-libgfortran -dynamiclib     \
             -Wall $(ARCH)                                      \
             -o $(TARGET_LIB_NAME)                              \
             $(BIN_DIR)/lapack.a                                \
             -lgfortran                                         \
             -lfftw3-$(BITS)									\
             $(OBJS)
	./test_for_libquadmath.sh $(TARGET_LIB_NAME)
	

# In Vespa we only use the library, not an executable. But the executable
# is useful for stepping through with gdb so I've used it occasionally. 
# See hlsvdmain in hlsvdpro.f for details.
exe: $(BIN_DIR) $(OBJS) lapack
	gfortran -static-libgcc -static-libgfortran   	\
	         -Wall $(ARCH)                      	\
	         -o hlsvdpro							\
		 	 $(BIN_DIR)/lapack.a					\
		 	 -lgfortran								\
		 	 -lfftw3-$(BITS)						\
	   	 	 $(OBJS)
	chmod +x hlsvdpro


# depcheck lists the runtime libraries that the target relies on. 
depcheck: $(TARGET_LIB_NAME)
	otool -L $(TARGET_LIB_NAME)


dist: $(TARGET_LIB_NAME)
	cp $(TARGET_LIB_NAME) $(TARGET_DIR)


clean:
	rm -rf $(BIN_DIR)

force:
