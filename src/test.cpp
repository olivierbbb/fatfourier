#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <iostream>

#ifdef FFTW_SEQ
#include <fftw3.h>
#else
#include "fft.hpp"
#ifdef CHECK_WITH_FFTW
#include <fftw3.h>
#endif
#endif

#include "wave_reader.hpp"

using namespace std;

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "usage: " << argv[0] << " <audio_in.wav> <samples_count>" << endl;
        return 1;
    }

    uint64_t buffer_size = strtoll(argv[2], NULL, 10);
    if (buffer_size != pow(2, log2(buffer_size))) {
        cerr << "<samples_count> must be a power of 2" << endl;
        return 1;
    }

    WaveReader wave(argv[1]);
    if (buffer_size >= wave.samples_count()) {
        cerr << "<samples_count> is bigger than file" << endl;
        return 1;
    }

    double* in_buffer = new double[buffer_size];
    double* out_buffer = new double[buffer_size];

    wave.read_samples(in_buffer, buffer_size);

// Compute values with fftw in order to check our own results later
#ifdef CHECK_WITH_FFTW
    // TODO use fftw malloc
    double* check_buffer = new double[buffer_size];
    copy(in_buffer, in_buffer + buffer_size, check_buffer);

    fftw_plan check_plan = fftw_plan_r2r_1d(
        buffer_size,
        check_buffer,
        check_buffer,
        FFTW_R2HC,
        FFTW_ESTIMATE);

    fftw_execute(check_plan);
#endif

    chrono::time_point<chrono::system_clock> start, end;
    start = chrono::system_clock::now();

// compute with FFTW or ourselves
#ifdef FFTW_SEQ
    fftw_plan plan = fftw_plan_r2r_1d(
        buffer_size,
        in_buffer,
        out_buffer,
        FFTW_R2HC,
        FFTW_ESTIMATE);
    fftw_execute(plan);
#else
    FFT fft(in_buffer, out_buffer, buffer_size);
    fft.compute();
#endif

    end = chrono::system_clock::now();

    unsigned int ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "fft of " << buffer_size << " samples computed in " << ms << "ms" << endl;

#ifdef CHECK_WITH_FFTW
    uint64_t freqs_count = buffer_size / 2 + 1;
    for (uint64_t i = 0; i < freqs_count; i++) {
        int precision = 100000;
        double value = precision * trunc(out_buffer[i]) / precision;
        double fftw_value = precision * trunc(check_buffer[i]) / precision;
        if (value != fftw_value) {
            cerr << "Wrong result" << value << fftw_value << endl;
            return 1;
        }
    }
#endif

#ifdef FFTW_SEQ
    fftw_destroy_plan(plan);
#endif
    delete[] in_buffer;
    delete[] out_buffer;

#ifdef CHECK_WITH_FFTW
    fftw_destroy_plan(check_plan);
    delete[] check_buffer;
#endif
}
