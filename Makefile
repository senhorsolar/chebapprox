CC := g++-13
CXXFLAGS := -std=c++23
INCLUDES := -Iinclude # local path
SRCFILES := $(shell ls src/*.cpp)
OFILES = $(patsubst src/%.cpp,obj/%.o, $(SRCFILES))

EIGEN_PATH ?= /usr/include/eigen3
INCLUDES += -I$(EIGEN_PATH)

all: example1D exampleND

obj:
	mkdir -p $@

obj/%.o: src/%.cpp obj
	$(CC) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

exampleND: exampleND.cpp $(OFILES)
	$(CC) $(CXXFLAGS) $(INCLUDES) -o $@ $^

example1D: example1D.cpp $(OFILES)
	$(CC) $(CXXFLAGS) $(INCLUDES) -o $@ $^

clean:
	rm -f example1D exampleND
	rm -rf obj

.PHONY: all clean
