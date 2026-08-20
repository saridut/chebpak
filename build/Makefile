#Makefile for gfortran compiler

SHELL = /bin/bash

DEBUG = false#{true, false}  (NOTE: no whitespace at end)

FC = gfortran

#*************************************************

ifeq ($(DEBUG), true)
FCFLAGS = -c -O0 -Wall -fcheck=all -fbacktrace -std=f2018 \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
else
FCFLAGS = -c -O3 -std=f2008ts \
		  -funroll-loops \
		  -flto=auto     \
		  -fno-unsafe-math-optimizations \
		  -frounding-math -fsignaling-nans
endif

SRC_DIR = ../src 
BUILD_DIR = .
LIB_DIR = ../lib
LIB_NAME = $(LIB_DIR)/libchebpak.a

vpath %.f90 $(SRC_DIR)

#Redefine the location of the interface files if necessary
MKL_INTRFC = $(HOME)/soft/mkl

INCS = \
	-I$(SRC_DIR) \
	-I$(MKL_INTRFC)/include

LIBS = \
	$(MKL_INTRFC)/lib/libmkl_blas95_lp64.a \
	$(MKL_INTRFC)/lib/libmkl_lapack95_lp64.a \
	-L$(MKLROOT)/lib -Wl,-rpath,$(MKLROOT)/lib \
	-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread \
	-lz -ldl -lm

OBJ_NAMES = \
		  constants_m    \
		  strings_m      \
		  utils_math_m   \
		  chebyshev_m

OBJECTS = $(addsuffix .o,$(OBJ_NAMES))

.PHONY : lib
lib : $(LIB_NAME)

.PHONY : doc
doc :
	ford ../docs/ford.md

$(LIB_NAME) : $(OBJECTS)
	if [ ! -d $(LIB_DIR) ]; then mkdir -p $(LIB_DIR); fi
	ar -rcs $@ $^

$(OBJECTS) : %.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@

.PHONY : clean clobber
clean:
	-rm -rf *.o *.mod
clobber:
	-rm -rf $(LIB_NAME) *.o *.mod
