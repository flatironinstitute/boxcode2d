# Makefile for boxcode2d
# This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 


# compiler, and linking from C, fortran
CC=gcc
CXX=g++
FC=gfortran


# set compiler flags for c and fortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy -w
FFLAGSTAB = -fPIC -std=legacy -w
CFLAGS= -std=c99 
CFLAGS+= $(FFLAGS) 
CXXFLAGS= -std=c++11 -DSCTL_PROFILE=-1
CXXFLAGS+=$(FFLAGS)

# set linking libraries
CLIBS = -lgfortran -lm -ldl 
LIBS = -lm

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

# Python Exetucable
PYTHON=python


# flags for MATLAB MEX compilation..
MFLAGS=-compatibleArrayDims -DMWF77_UNDERSCORE1 
MWFLAGS=-c99complex 
MOMPFLAGS = -D_OPENMP

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran

FMM_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	FMM_INSTALL_DIR = ${HOME}/lib
endif

DYLIBS = $(LIBS)

LIBNAME=libboxcode2d
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LLINKLIB = -lboxcode2d


# For your OS, override the above by placing make variables in make.inc
-include make.inc


# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)

CFLAGS += $(OMPFLAGS)
FFLAGS += $(OMPFLAGS)
MFLAGS += $(MOMPFLAGS)
LIBS += $(OMPLIBS)
DYLIBS += $(OMPLIBS)
MEXLIBS += $(OMPLIBS)

endif


# folders
SRC = ./src
TAB = ./tab
GEN = ./gen
TEST = ./test

# objects to compile
#
# Common objects

OBJS = $(SRC)/lrtree2d.o $(SRC)/compositegrids2d.o \
	$(SRC)/legetens.o $(SRC)/legeexps.o $(SRC)/poissbox2d.o \
	$(SRC)/poissbox2drouts.o $(SRC)/l2dpwrouts.o $(SRC)/boxcode2drouts.o \
	$(SRC)/pbox2dreftabs.o $(SRC)/voltab2d.o $(SRC)/loadsyms2d.o $(SRC)/lwtsexp_ext.o \
	$(SRC)/l2dterms_basic.o $(SRC)/prini.o $(SRC)/dlaran.o  \
	$(SRC)/pyplot.o \
	$(SRC)/hkrand.o \
	$(SRC)/lrt2d_compositegrid2d_pyplot.o 

TOBJS = $(GEN)/dcuhre.o $(TEST)/laprouts2d.o 

TABLES = $(patsubst %.f,%.o,$(wildcard $(TAB)/pbox2dtab*.f)) 

.PHONY: usage install lib test all clean

default: usage

usage:
	@echo "Makefile for boxcode2d. Specify what to make:"
	@echo "  make install - compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR) - compile and install the main library at custom location given by PREFIX"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test - compile and run validation tests (will take a couple of mins)"
	@echo "  make objclean - remove most object files, preserving tables and libs"
	@echo "  make deepobjclean - remove all object files, preserving libs"
	@echo "  make clean - also remove libs and executables"
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=OFF' for single-threaded"


# implicit rules for objects (note -o ensures writes to correct dir)

tab/%.o: tab/%.f 
	$(FC) -c $(FFLAGSTAB) $< -o $@

%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@



# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif


install: $(STATICLIB) $(DYNAMICLIB)
	echo $(FMM_INSTALL_DIR)
	mkdir -p $(FMM_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(FMM_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(FMM_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp lib/$(LIMPLIB) $(FMM_INSTALL_DIR)/
	@echo "Make sure to include " $(FMM_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(FMM_INSTALL_DIR) " -lfmm2d"


$(STATICLIB): $(OBJS) $(TABLES)
	ar rcs $(STATICLIB) $(OBJS) $(TABLES)
	mv $(STATICLIB) lib-static/
$(DYNAMICLIB): $(OBJS) $(TABLES)
	$(FC) -shared -fPIC $(OBJS) $(TABLES) -o $(DYNAMICLIB) $(DYLIBS) 
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/


# testing routines
#


test/l2dpwrouts:
	$(FC) $(FFLAGS) test/test_l2dpwrouts.f $(TOBJS) $(OBJS) $(TABLES) -o test/int2-test_l2dpwrouts $(LIBS)

test/lrt2d:
	$(FC) $(FFLAGS) test/test_lrt2d.f $(TOBJS) $(OBJS) $(TABLES) -o test/int2-test_lrt2d $(LIBS)

test/pbox2drouts:
	$(FC) $(FFLAGS) test/test_pbox2drouts.f $(TOBJS) $(OBJS) $(TABLES) -o test/int2-test_pbox2drouts $(LIBS)


test: $(STATICLIB) $(TOBJS) test/l2dpwrouts test/lrt2d test/pbox2drouts
	rm -f print_testres.txt
	(cd test/; ./run_unittests.sh)
	cat print_testres.txt
	rm print_testres.txt

objclean:
	rm -f $(OBJS) $(TOBJS)

deepobjclean:
	rm -f $(OBJS) $(TOBJS) $(TABLES)

clean:
	rm -f $(OBJS) $(TOBJS)
	rm -f $(TABLES)
	rm -f lib/$(DYNAMICLIB)
	rm -f lib-static/$(STATICLIB)
	rm -f test/int2*
	rm -f gen/int2*


