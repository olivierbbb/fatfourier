#include <chrono>
#include <cmath>
#include <cstdint>
#include <getopt.h>
#include <iostream>

#ifdef USE_FFTW
#include <fftw3.h>
#else
#include "fat_fourier.hpp"
#ifdef CHECK_WITH_FFTW
#include <fftw3.h>
#endif
#endif

#define IS_POWER_OF_2(x) (((x) & ((x) -1)) == 0)

using namespace std;
#if !defined(OMP) && !defined(USE_FFTW)
    #define THREADS_COUNT_OPT
#endif

void usage(string exec_name) {
    #ifdef THREADS_COUNT_OPT
    cerr << "Usage: " << exec_name << " -l <length>\n";
    #else
    cerr << "Usage: " << exec_name << " -l <length> [-t <nb_threads>]\n";
    #endif
    cerr << "Required:\n";
    cerr << "  -l <length>\t\tpower of two of the computed sequence length\n";
    #ifdef THREADS_COUNT_OPT
    cerr << "Options:\n";
    cerr << "  -t <nb_threads>\t\tnumber of threads to spawn [default: auto]\n";
    #endif
}

int main(int argc, char* argv[]) {
    int seq_length = -1;
    #ifdef THREADS_COUNT_OPT
        int threads_count = -1;
    #endif

    char opt;
    #ifdef THREADS_COUNT_OPT
    char const* opt_list = "l:t:";
    #else
    char const* opt_list = "l:";
    #endif

    while ((opt = getopt(argc, argv, opt_list)) != -1) {
        switch (opt) {
        case 'l':
            seq_length = atoi(optarg);
            break;
        #ifdef THREADS_COUNT_OPT
        case 't':
            threads_count = atoi(optarg);
            break;
        #endif
        default:
            usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    if (seq_length <= 0) {
        cerr << "Sequence length missing\n";
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

#ifdef THREADS_COUNT_OPT
    if (threads_count <= 0) {
        cerr << "Number of threads missing\n";
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }
#endif

    uint64_t buffer_size = pow(2, seq_length);
    double* in_buffer = new double[buffer_size];
    double* out_buffer = new double[buffer_size];

    // fill buffer with random vlaues
    for (uint64_t i = 0; i < buffer_size; i++) {
        in_buffer[i] = rand() / RAND_MAX;
    }

// Compute values with fftw in order to check results later
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

// compute with fatfourier or FFTW
#ifdef USE_FFTW
    fftw_plan plan = fftw_plan_r2r_1d(
        buffer_size,
        in_buffer,
        out_buffer,
        FFTW_R2HC,
        FFTW_ESTIMATE);
    fftw_execute(plan);
#else
    FatFourier ff(in_buffer, out_buffer, buffer_size);
#ifdef OMP
    ff.compute();
#else
    ff.compute(threads_count);
#endif
#endif

    end = chrono::system_clock::now();

    unsigned int ms = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "FFT of " << buffer_size << " samples computed in " << ms << "ms\n";

#ifdef CHECK_WITH_FFTW
    uint64_t freqs_count = buffer_size / 2 + 1;
    for (uint64_t i = 0; i < freqs_count; i++) {
        int precision = 100000;
        double value = precision * trunc(out_buffer[i]) / precision;
        double fftw_value = precision * trunc(check_buffer[i]) / precision;
        if (value != fftw_value) {
            cerr << "Incorect result: \"" << value << "\", should be: \"" << fftw_value << "\"\n";
            break;
        }
    }
#endif

#ifdef USE_FFTW
    fftw_destroy_plan(plan);
#endif
    delete[] in_buffer;
    delete[] out_buffer;

#ifdef CHECK_WITH_FFTW
    fftw_destroy_plan(check_plan);
    delete[] check_buffer;
#endif
}
