# Makefile for test_lrt2d

NAME = gentabs_pbox2d
DRIVERFILE = $(NAME).f

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
#OMPFLAGS = 
OMPLIBS =-lgomp 

FFLAGS += $(OMPFLAGS)

LBLAS = -lblas -llapack

LIBS = -lm
DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

LIBS = $(OMPLIBS)

# For your OS, override the above by placing make variables in make.inc
-include make.inc


# objects to compile
#
# Common objects
SRC = ../src
GEN = .

L2DPWOBJS = $(SRC)/legetens.o $(SRC)/legeexps.o \
	$(GEN)/pbox2d_gentab.o $(GEN)/dcuhre.o

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
	$(FC) $(FFLAGS) $(DRIVERFILE) -o int2-$(NAME) $(OBJS) $(LIBS)
	./int2-$(NAME)
#
# housekeeping routines
#
clean: objclean
	rm -f int2-$(NAME)

objclean: 
	rm -f $(OBJS)
	rm -f *.o
