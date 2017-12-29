#!/bin/bash

TARGET = pfm_alpha.exe
INCLUDES= -Isrc/core -Isrc/utils -Isrc/heat

CSRC = $(wildcard src/core/*.cpp) \
	   $(wildcard src/heat/*.cpp) \
	   $(wildcard src/*.cpp) 
CSRC2 = $(wildcard src/utils/*.c) 

FSRC = 
	   
SRCS = $(CSRC)

COBJS = $(CSRC:.cpp=.o) $(CSRC2:.c=.o)
OBJS  = $(COBJS)

CXX=mpic++
FC=mpif90
OPTFLAGS= -O2 
CXXFLAGS= $(OPTFLAGS) $(INCLUDES)

LDFLAGS=

RM = rm

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS) 

.PHONY: clean
clean:
	$(RM) $(OBJS)

depend:
	makedepend -Y *.cpp *.h
