CXX = g++
LD = $(CXX)
CXXFLAGS = -std=c++14 -Wall -Wextra -Wno-sign-compare
CXXFLAGS_DBG = -g -DEBUG -DCHECK_WITH_FFTW
CXXFLAGS_RLS = -O3 -DNDEBUG
CXXFLAGS += $(CXXFLAGS_DBG)
LDFLAGS = -lm -lfftw3

SRC_DIR = src
SRC_FILES = $(wildcard $(SRC_DIR)/*.cpp)
HEADER_FILES = $(wildcard $(SRC_DIR)/*.hpp)
BIN_DIR = bin

.PHONY: all clean

all: $(BIN_DIR)/fatfourier \
$(BIN_DIR)/fatfourier_omp $(BIN_DIR)/fatfourier_threaded_2 $(BIN_DIR)/fatfourier_threaded_4 \
$(BIN_DIR)/fatfourier_threaded_auto $(BIN_DIR)/fatfourier_fftw

$(BIN_DIR)/fatfourier: $(HEADER_FILES) $(SRC_FILES)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(SRC_FILES) -o $@ $(LDFLAGS)

$(BIN_DIR)/fatfourier_omp: $(HEADER_FILES) $(SRC_FILES)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -DOMP $(SRC_FILES) -o $@ $(LDFLAGS) -fopenmp

$(BIN_DIR)/fatfourier_threaded_2: $(HEADER_FILES) $(SRC_FILES)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -DTHREADED -DTHREADS_COUNT=2 -pthread $(SRC_FILES) -o $@ $(LDFLAGS)

$(BIN_DIR)/fatfourier_threaded_4: $(HEADER_FILES) $(SRC_FILES)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -DTHREADED -DTHREADS_COUNT=4 -pthread $(SRC_FILES) -o $@ $(LDFLAGS)

$(BIN_DIR)/fatfourier_threaded_8: $(HEADER_FILES) $(SRC_FILES)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -DTHREADED -DTHREADS_COUNT=4 -pthread $(SRC_FILES) -o $@ $(LDFLAGS)

$(BIN_DIR)/fatfourier_threaded_auto: $(HEADER_FILES) $(SRC_FILES)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -DTHREADED -pthread $(SRC_FILES) -o $@ $(LDFLAGS) 

$(BIN_DIR)/fatfourier_fftw: $(HEADER_FILES) $(SRC_FILES)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -DFFTW_SEQ $(SRC_FILES) -o $@ $(LDFLAGS)

clean:
	$(RM) $(BIN_DIR)/*
