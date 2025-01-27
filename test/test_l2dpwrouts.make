# Makefile for test_lrt2d

TESTNAME = test_l2dpwrouts
TESTFILE = $(TESTNAME).f

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy
FFLAGS = -g -std=legacy

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

LBLAS = -lblas -llapack

LIBS = -lm
DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

# For your OS, override the above by placing make variables in make.inc
-include make.inc


# objects to compile
#
# Common objects
SRC = ../src

L2DPWOBJS = $(SRC)/legetens.o $(SRC)/laprouts2d.o $(SRC)/legeexps.o \
	$(SRC)/l2dterms.o $(SRC)/poissbox2drouts.o $(SRC)/boxcode2drouts.o \
	$(SRC)/l2dpwrouts.o  $(SRC)/lrtree2d.o $(SRC)/lwtsexp_ext.o $(SRC)/pyplot.o

TESTOBJS = $(SRC)/hkrand.o $(SRC)/dlaran.o $(SRC)/prini_new.o

OBJS = $(L2DPWOBJS) $(TESTOBJS)

.PHONY: test 

default: test


#
# implicit rules for objects (note -o ensures writes to correct dir)
#
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

test: $(OBJS) $(TESTOBJS)
	$(FC) $(FFLAGS) $(TESTFILE) -o int2-$(TESTNAME) $(OBJS)
	./int2-$(TESTNAME)
#
# housekeeping routines
#
clean: objclean
	rm -f int2-$(TESTNAME)

objclean: 
	rm -f $(OBJS)
	rm -f *.o
