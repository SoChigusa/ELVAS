CXX := g++
CXXFLAGS := -g -Wall -std=c++14 -O2
INCLUDE = -I/Users/sochigusa/works/ELVAS/src/include -I/Users/sochigusa/works/RGEFlow/
RGEFLOW_PATH = /Users/sochigusa/works/RGEFlow

all: elvas.o sm.out

elvas.o: include/elvas.h elvas.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDE) -o elvas.o -c elvas.cpp

sm.out: elvas.o $(addprefix /Users/SoChigusa/works/RGEFlow/, *.o) sm.cpp
		$(CXX) $(CXXFLAGS) $(INCLUDE) -o sm.out sm.cpp elvas.o $(addprefix /Users/SoChigusa/works/RGEFlow/, *.o)