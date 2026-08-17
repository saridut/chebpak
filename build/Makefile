#Makefile for Taylor with GCC-10 compilers

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

SRC_DIR = $(HOME)/workspace/npl/pplib 
BUILD_DIR = .
INSTALL_DIR = .
EXECUTABLE = calc_disp

vpath %.f90 $(SRC_DIR)
vpath $(EXECUTABLE) $(INSTALL_DIR)

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
		  M_CLI2         \
		  strings_m      \
		  utils_math_m   \
		  qsort_m        \
		  vector_m       \
		  table_m        \
		  trajectory_reader_m   \
		  simbox_m       \
		  specs_m        \
		  atmcfg_m       \
		  config_io_m    \
		  chebyshev_m    \
		  atc_mesh_m     \
		  atc_meshgen_m  \
		  atc_vtk_m      \
		  atc_atm2mesh_m \
		  calc_disp

OBJECTS = $(addsuffix .o,$(OBJ_NAMES))

$(EXECUTABLE) : $(OBJECTS)
	$(FC) $(INCS) $^ $(LIBS) -o $@

$(OBJECTS) : %.o : %.f90
	$(FC) $(INCS) $(FCFLAGS) $^ -o $@

.PHONY : clean clobber
clean:
	-rm -rf *.o *.mod
clobber:
	-rm -rf $(INSTALL_DIR)/$(EXECUTABLE) *.o *.mod
