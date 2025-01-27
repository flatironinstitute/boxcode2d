# Makefile for test_lrt2d

TESTNAME = test_lrt2d
TESTFILE = $(TESTNAME).f

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy

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

TREEOBJS = $(SRC)/lrtree2d.o
TESTOBJS = $(SRC)/hkrand.o $(SRC)/dlaran.o $(SRC)/prini_new.o $(SRC)/pyplot.o

OBJS = $(TREEOBJS) $(TESTOBJS)

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
