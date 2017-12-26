#!/bin/bash

TARGET = pfm_alpha.exe 

CSRC = common.cpp blas.cpp Grid.cpp HeatGrid.cpp PBiCGStab.cpp HeatCG.cpp main.cpp 
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
# DO NOT DELETE

Grid.o: Grid.h common.h
HeatCG.o: HeatCG.h PBiCGStab.h common.h blas.h Grid.h
HeatGrid.o: HeatGrid.h Grid.h common.h
PBiCGStab.o: PBiCGStab.h common.h blas.h Grid.h
blas.o: blas.h
common.o: common.h
main.o: common.h PinT.h HeatGrid.h Grid.h HeatCG.h PBiCGStab.h blas.h
Grid.o: common.h
HeatCG.o: PBiCGStab.h common.h blas.h Grid.h
HeatGrid.o: Grid.h common.h
PBiCGStab.o: common.h blas.h Grid.h
