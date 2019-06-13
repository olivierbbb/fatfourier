DEBUG = 0

PACKAGE = fat_fourier
LIB_TARGET = build/lib$(PACKAGE).so
DEMO_TARGETS = build/benchmark

CXX = g++
LD = $(CXX)
CXXFLAGS = -std=c++14 -Wall -Wextra -Iinclude
LDFLAGS =

ifeq ($(DEBUG), 1)
CXXFLAGS += -g -DEBUG
else
CXXFLAGS += -O2 -DNDEBUG
endif

LIB_CXXFLAGS = -fPIC
LIB_LDFLAGS = -lm

DEMO_CXXFLAGS += $(shell pkg-config --cflags fftw3)
DEMO_LDFLAGS += $(shell pkg-config --libs fftw3)

LIB_SRCS = $(wildcard src/*.cpp)
LIB_OBJS = $(patsubst src/%.cpp, build/lib/%.o, $(LIB_SRCS))
LIB_DEPS = $(wildcard build/deps/lib/*.d)
DEMO_SRCS = $(wildcard demo/*.cpp)
DEMO_DEPS = $(wildcard build/deps/demo/*.d)

.PHONY: all lib demo clean

all: lib demo

lib: $(LIB_TARGET)

demo: $(DEMO_TARGETS)

build/lib$(PACKAGE).so: $(LIB_OBJS)
	@mkdir -p $(@D)
	$(LD) -shared $^ -o $@ $(LDFLAGS) $(LIB_LDFLAGS)

build/%: build/demo/%.o $(LIB_OBJS)
	@mkdir -p $(@D)
	$(LD) -o $@ $^ $(LDFLAGS) $(DEMO_LDFLAGS)

build/lib/%.o: src/%.cpp
	@mkdir -p $(@D) build/deps/lib
	$(CXX) $(CXXFLAGS) $(LIB_CXXFLAGS) -MMD -MF build/deps/lib/$*.d -c -o $@ $<

build/demo/%.o: demo/%.cpp $(LIB_OBJS)
	@mkdir -p $(@D) build/deps/demo
	$(CXX) $(CXXFLAGS) $(DEMO_CXXFLAGS) -MMD -MF build/deps/demo/$*.d -c -o $@ $<

ifneq ($(MAKECMDGOALS), clean)
-include $(LIB_DEPS)
-include $(DEMO_DEPS)
endif

clean:
	$(RM) -r build/*
