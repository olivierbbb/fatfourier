#pragma once

#include <cstdint>
#include <fstream>
#include <string>

class WaveReader {
public:
    WaveReader(std::string filename);
    uint64_t samples_count() { return samples_count_; }
    void read_samples(double* buffer, uint64_t buffer_size);

private:
    std::ifstream ifs_;
    uint64_t data_end_;
    uint64_t samples_count_;
};
