#include "wave_reader.hpp"
#include <cassert>
#include <iostream>

using namespace std;

#define WAVE_FORMAT_PCM 1

struct WaveHeader {
    char riff_id[4]; // "RIFF"
    uint32_t riff_size;
    char wave_id[4]; // "WAVE"

    char fmt_id[4]; // "fmt "
    uint32_t fmt_size;
    uint16_t fmt_tag; // WAVE_FORMAT_PCM
    uint16_t channels_count; // 1
    uint32_t sample_rate;
    uint32_t byte_rate;
    uint16_t block_align; // 1 * 16 / 8
    uint16_t bit_depth; // 16

    char data_id[4]; // "data"
    uint32_t data_size;
};

WaveReader::WaveReader(string filename) :
    ifs_(filename, ios::in | ios::binary) {
    if (!ifs_.is_open()) {
        cerr << "cannot open file \"" << filename << "\"" << endl;
        exit(1);
    }

    WaveHeader header;
    ifs_.read((char*) &header.riff_id, sizeof(header.riff_id));
    ifs_.read((char*) &header.riff_size, sizeof(header.riff_size));
    ifs_.read((char*) &header.wave_id, sizeof(header.wave_id));

    ifs_.read((char*) &header.fmt_id, sizeof(header.fmt_id));
    ifs_.read((char*) &header.fmt_size, sizeof(header.fmt_size));

    unsigned int count = 0;
    while (header.fmt_id[0] != 'f' && count < 5) { // skip until fmt
        ifs_.seekg((streamoff) header.fmt_size, ios::cur);
        ifs_.read((char*) &header.fmt_id, sizeof(header.fmt_id));
        ifs_.read((char*) &header.fmt_size, sizeof(header.fmt_size));
        count++;
    }

    ifs_.read((char*) &header.fmt_tag, sizeof(header.fmt_tag));
    assert(header.fmt_tag == WAVE_FORMAT_PCM);
    ifs_.read((char*) &header.channels_count, sizeof(header.channels_count));
    assert(header.channels_count == 1);
    ifs_.read((char*) &header.sample_rate, sizeof(header.sample_rate));
    ifs_.read((char*) &header.byte_rate, sizeof(header.byte_rate));
    ifs_.read((char*) &header.block_align, sizeof(header.block_align));
    assert(header.block_align == 1 * 16 / 8);
    ifs_.read((char*) &header.bit_depth, sizeof(header.bit_depth));
    assert(header.bit_depth == 16);

    ifs_.read((char*) &header.data_id, sizeof(header.data_id));
    ifs_.read((char*) &header.data_size, sizeof(header.data_size));

    count = 0;
    while (header.data_id[0] != 'd' && count < 5) { // skip until data chunk
        ifs_.seekg((streamoff) header.data_size, ios::cur);
        ifs_.read((char*) &header.data_id, sizeof(header.data_id));
        ifs_.read((char*) &header.data_size, sizeof(header.data_size));
        count++;
    }

    data_end_ = (uint64_t) ifs_.tellg() + header.data_size;
    samples_count_ = (uint64_t) header.data_size * 8 / header.bit_depth;
}

void WaveReader::read_samples(double* buffer, uint64_t buffer_size) {
    assert(ifs_.is_open());

    int16_t* int_buffer = new int16_t[buffer_size];
    streamsize read_size = 2 * buffer_size;
    assert(ifs_.tellg() + (streamsize) read_size <= data_end_);

    ifs_.read((char*) int_buffer, read_size);
    assert(ifs_.good());

    for (uint64_t i = 0; i < buffer_size; i++) {
        buffer[i] = (double) int_buffer[i] / INT16_MAX;
    }

    delete[] int_buffer;
}
