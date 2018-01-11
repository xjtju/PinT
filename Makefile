#!/bin/bash

TARGET = pfm_alpha.exe

PMLIB_HOME = /Users/bamboo/Libs/PMlib
INCLUDES= -Isrc/include -Isrc/heat -I${PMLIB_HOME}/include

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
OPTFLAGS= -O2   
#OPTFLAGS= -fopenmp -O2 -D_PMLib_  
CXXFLAGS= $(OPTFLAGS) $(INCLUDES)
FFLAGS  = $(OPTFLAGS) $(INCLUDES) -fdefault-real-8 -fdefault-double-8
LDFLAGS = -L${PMLIB_HOME}/lib -lPMmpi

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
