CC = g++
CXXFLAGS = -std=c++11
CXXRELEASE_FLAGS = $(CXXFLAGS) -O2 -DNDEBUG

#BOOST_ROOT=
ifdef BOOST_ROOT
BOOST_FLAGS = -I$(BOOST_ROOT)/include -DHAVE_BOOST
else
BOOST_FLAGS =
endif

binaries = test_alg1

.PHONY: all
all: $(binaries)
		
test_alg1: %: %.cpp
	$(CC) $(CXXRELEASE_FLAGS) $(BOOST_FLAGS) -o $@ $<

.PHONY: clean
clean:
	rm -f $(binaries)
