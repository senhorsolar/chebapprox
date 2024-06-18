CC := g++
CXXFLAGS := -std=c++20

INCLUDES := -Iinclude # local path
SRCFILES := $(shell ls src/*.cpp)
OFILES = $(SRCFILES:.cpp=.o)

all: example1D exampleND

%.o: %.cpp
	$(CC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

exampleND: exampleND.cpp $(OFILES)
	$(CC) $(CXXFLAGS) $(INCLUDES) -o $@ $^

example1D: example1D.cpp $(OFILES)
	$(CC) $(CXXFLAGS) $(INCLUDES) -o $@ $^

clean:
	rm -f example1D exampleND
	rm -f $(OFILES)

.PHONY: all clean
