CC = g++
CXXFLAGS = -std=c++11
CXXDEBUG_FLAGS = $(CXXFLAGS) -O0 -g
CXXRELEASE_FLAGS = $(CXXFLAGS) -O2 -DNDEBUG
OMPFLAGS=-fopenmp

#BOOST_ROOT=
ifdef BOOST_ROOT
BOOST_FLAGS = -I$(BOOST_ROOT)/include -DHAVE_BOOST
else
BOOST_FLAGS =
endif

binaries = test_141 test_6114 test_input 

.PHONY: all
all: $(binaries)
		
test_141: %: %.cpp
	$(CC) $(CXXDEBUG_FLAGS) -o $@ $<

test_6114: %: %.cpp
	$(CC) $(CXXRELEASE_FLAGS) $(OMPFLAGS) -o $@ $<

test_input: %: %.cpp
	$(CC) $(CXXRELEASE_FLAGS) $(OMPFLAGS) $(BOOST_FLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -f $(binaries)
