#!/bin/bash
# there are several option flag to controll the compiling process, 
# you can activate them or not according the running environment 

#_PMLib_ =1
ifdef _PMLib_ 
PMLIB_HOME = /Users/bamboo/Libs/PMlib
PMLIB_INC  = -I${PMLIB_HOME}/include 
LDFLAGS_PM = -L${PMLIB_HOME}/lib -lPMmpi 
OPTFLAGS_PM= -fopenmp -D_PMLib_  
endif

_HDF5_=1
ifdef _HDF5_
HDF5_HOME  = /Users/bamboo/Libs/hdf5
HDF5_INC   = -I${HDF5_HOME}/include
LDFLAGS_H5 = -L${HDF5_HOME}/lib -lhdf5 
OPTFLAGS_H5= -D _HDF5_   
endif

#_HEAT_ = 1
ifdef _HEAT_
TARGET = heat_alpha.exe
HEAT_INC = -Isrc/heat
HEAT_FSRC = $(wildcard src/heat/*.f90) 
HEAT_CSRC = $(wildcard src/heat/*.cpp) 
HEAT_OPT = -D _TEST_HEAT_   
endif

_PFM_ = 1
ifdef _PFM_
TARGET = pfm_alpha.exe
PFM_INC = -Isrc/pfm
PFM_FSRC = $(wildcard src/pfm/*.f90) 
PFM_CSRC = $(wildcard src/pfm/*.cpp) 
PFM_OPT = -D _TEST_PFM_   
endif

INCLUDES= -Isrc/include ${PMLIB_INC} ${HDF5_INC} ${HEAT_INC} ${PFM_INC}

.FUFFIXES: .o .cpp .c .f90

FSRC = $(wildcard src/core/*.f90)  \
	   $(wildcard src/utils/*.f90) \
       $(HEAT_FSRC) \
       $(PFM_FSRC)

CSRC = $(wildcard src/core/*.cpp) \
	   $(wildcard src/utils/*.cpp) \
	   $(wildcard src/*.cpp) \
	   $(HEAT_CSRC) \
	   $(PFM_CSRC) 

CSRC2 = $(wildcard src/utils/*.c) 

	   
SRCS = $(FSRC) $(CSRC)

FOBJS = $(FSRC:.f90=.o)
COBJS = $(CSRC:.cpp=.o) $(CSRC2:.c=.o)
OBJS  = $(FOBJS) $(COBJS)

CXX=mpic++
FC=mpif90


OPTFLAGS= -O3 ${OPTFLAGS_PM} ${OPTFLAGS_H5} $(HEAT_OPT) $(PFM_OPT)
CXXFLAGS= -std=c++11 $(OPTFLAGS) $(INCLUDES)
FFLAGS  = $(OPTFLAGS) $(INCLUDES) -fdefault-real-8 -fdefault-double-8
LDFLAGS = ${LDFLAGS_PM} $(LDFLAGS_H5) 


RM = rm

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $^  
%.o: %.c
	$(CXX) $(CXXFLAGS) -o $@ -c $^  
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $^ 

.PHONY: clean
clean:
	$(RM) $(OBJS)

depend:
	makedepend -Y *.cpp *.h
