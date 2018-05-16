DEBUG = 0

ifeq ($(shell pkg-config --exists fftw3 && echo 1), 1)
HAS_FFTW = 1
else
HAS_FFTW = 0
endif

DEMO_TARGETS = bin/benchmark_ff bin/benchmark_ff_omp
ifeq ($(HAS_FFTW), 1)
	 DEMO_TARGETS += bin/benchmark_fftw
endif

CXX = g++
LD = $(CXX)
CXXFLAGS = -std=c++14 -Wall -Wextra -Wno-sign-compare -Isrc/
LDFLAGS = -lm
ifeq ($(HAS_FFTW), 1)
CFLAGS += $(shell pkg-config --clfags fftw3)
LDFLAGS += $(shell pkg-config --libs fftw3)
endif

ifeq ($(DEBUG), 1)
CXXFLAGS += -g -DEBUG
ifeq ($(HAS_FFTW), 1)
CXXFLAGS += -DCHECK_WITH_FFTW
endif
else
CXXFLAGS += -O3 -DNDEBUG
endif

LIB_SRCS = $(wildcard src/*.cpp)
LIBS_HEADERS = $(wildcard src/*.hpp)
DEMO_SRCS = $(wildcard demo/*.cpp)

.PHONY: all clean

all: demo

demo: $(DEMO_TARGETS)

bin/benchmark_ff: $(LIBS_HEADERS) $(LIB_SRCS) $(DEMO_SRCS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -pthread $(LIB_SRCS) $(DEMO_SRCS) -o $@ $(LDFLAGS)

bin/benchmark_ff_omp: $(LIBS_HEADERS) $(LIB_SRCS) $(DEMO_SRCS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -DOMP $(LIB_SRCS) $(DEMO_SRCS) -o $@ $(LDFLAGS) -fopenmp

bin/benchmark_fftw: $(LIBS_HEADERS) $(LIB_SRCS) $(DEMO_SRCS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -DUSE_FFTW $(DEMO_SRCS) -o $@ $(LDFLAGS)

clean:
	$(RM) bin/*
