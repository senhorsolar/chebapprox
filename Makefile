CC := g++
CXXFLAGS := -std=c++20

INCLUDES := -Iinclude # local path
SRCFILES := $(shell ls src/*.cpp)
OFILES = $(SRCFILES:.cpp=.o)

all: example1D

%.o: %.cpp
	$(CC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

example1D: example1D.cpp $(OFILES)
	$(CC) $(CXXFLAGS) $(INCLUDES) -o $@ $^

clean:
	rm -f example1D

.PHONY: all clean
