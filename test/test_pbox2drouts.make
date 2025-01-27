# Makefile for test_pbox2drouts

TESTNAME = test_pbox2drouts
TESTFILE = $(TESTNAME).f

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy
FFLAGSTAB = -fPIC -std=legacy

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
TAB = ../tab
GEN = ../gen

FMMOBJS = $(SRC)/lrtree2d.o $(SRC)/compositegrids2d.o \
	$(SRC)/legetens.o $(SRC)/legeexps.o $(SRC)/poissbox2d.o \
	$(SRC)/poissbox2drouts.o $(SRC)/l2dpwrouts.o $(SRC)/boxcode2drouts.o \
	$(SRC)/pbox2dreftabs.o $(SRC)/voltab2d.o $(SRC)/loadsyms2d.o $(SRC)/lwtsexp_ext.o \
	$(SRC)/l2dterms.o $(SRC)/pyplot.o 

# tables take a long time to compile, keep separate from standard
# clean, etc. 
TABLES = $(patsubst %.f,%.o,$(wildcard $(TAB)/pbox2dtab*.f)) 

TESTOBJS = $(SRC)/hkrand.o $(SRC)/dlaran.o $(SRC)/prini_new.o \
	$(GEN)/dcuhre.o

NONTABLEOBJS = $(FMMOBJS) $(TESTOBJS)

OBJS = $(NONTABLEOBJS) $(TABLES)

.PHONY: test 

default: test

#
# implicit rules for objects (note -o ensures writes to correct dir)
#

$(TAB)/%.o: $(TAB)/%.f 
	$(FC) -c $(FFLAGSTAB) $< -o $@

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
clean:
	rm -f $(NONTABLEOBJS)
	rm -f int2-$(TESTNAME)

deepclean: 
	rm -f $(OBJS)
	rm -f int2-$(TESTNAME)
