#!/bin/bash

TARGET = pfm_alpha.exe

ifdef _PMLib_ 
PMLIB_HOME = /Users/bamboo/Libs/PMlib
PMLIB_INC  = -I${PMLIB_HOME}/include 
LDFLAGS_PM = -L${PMLIB_HOME}/lib -lPMmpi 
OPTFLAGS_PM= -fopenmp -D_PMLib_  
endif

ifdef _HDF5_
HDF5_HOME  = /Users/bamboo/Libs/hdf5
HDF5_INC   = -I${HDF5_HOME}/include
LDFLAGS_H5 = -L${HDF5_HOME}/lib -lhdf5 
OPTFLAGS_H5= -D _HDF5_   
endif

INCLUDES= -Isrc/include -Isrc/heat ${PMLIB_INC} ${HDF5_INC}

.FUFFIXES: .o .cpp .c .f90

FSRC = $(wildcard src/core/*.f90)  \
	   $(wildcard src/utils/*.f90) \
	   $(wildcard src/heat/*.f90) 

CSRC = $(wildcard src/core/*.cpp) \
	   $(wildcard src/heat/*.cpp) \
	   $(wildcard src/utils/*.cpp) \
	   $(wildcard src/*.cpp) 
CSRC2 = $(wildcard src/utils/*.c) 

	   
SRCS = $(FSRC) $(CSRC)

FOBJS = $(FSRC:.f90=.o)
COBJS = $(CSRC:.cpp=.o) $(CSRC2:.c=.o)
OBJS  = $(FOBJS) $(COBJS)

CXX=mpic++
FC=mpif90


OPTFLAGS= -O3 ${OPTFLAGS_PM} ${OPTFLAGS_H5} 
CXXFLAGS= $(OPTFLAGS) $(INCLUDES)
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
