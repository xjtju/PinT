#!/bin/bash

TARGET = pfm_alpha.exe 

CSRC = common.cpp blas.cpp PBiCGStab.cpp HeatCG.cpp main.cpp 
SRCS = $(CSRC)

.SUFFIXES: .o .cpp

COBJS = $(CSRC:.cpp=.o)
OBJS  = $(COBJS)

CXX=mpic++
FC=mpif90
OPTFLAGS=-O2 
CXXFLAGS= $(OPTFLAGS)

LDFLAGS=

RM = rm

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(@) $(OBJS) $(LDFLAGS)

.cpp.o:
	$(CXX) $(CXXFLAGS) -o $(@) -c $<

clean:
	$(RM) $(OBJS)

depend:
	makedepend -Y *.cpp *.h
