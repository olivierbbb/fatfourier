#include "fatfourier.hpp"
#include <chrono>
#include <cmath>
#include <complex>
#include <fftw3.h>
#include <getopt.h>
#include <iostream>
#include <utility>

#define IS_POWER_OF_2(x) (((x) & ((x) -1)) == 0)

using std::complex;
using std::string;
using std::tuple;
using std::vector;
using namespace fatfourier;

const int FFT_SIZE_DEFAULT = 4096;
const int NB_ITERS_DEFAULT = 1000;

double truncate_precision(double value, int precision = 1e9)
{
    return precision * trunc(value) / precision;
}

void usage(string exec_name)
{
    std::cerr << "Usage: " << exec_name << " [options]\n";
    std::cerr << "Options:\n";
    std::cerr << "  -l <fft-size>\tLength of FFT to compute [default: " << FFT_SIZE_DEFAULT << "]\n";
    std::cerr << "  -n <nb-iters>\tIterations to run [default: " << NB_ITERS_DEFAULT << "]\n";
    std::cerr << "  -c\tCompute complex FFt instead of real\n";
}

bool check_values(complex<double>* complex_ff, complex<double>* complex_fftw, size_t size)
{
    double* real_ff = reinterpret_cast<double*>(complex_ff);
    double* real_fftw = reinterpret_cast<double*>(complex_fftw);

    for (size_t i = 0; i < size * 2; i++) {
        if (truncate_precision(real_ff[i]) != truncate_precision(real_fftw[i])) {
            return false;
        }
    }
    return true;
}

tuple<int, int> benchmark_real(int fft_size, int nb_iters)
{
    int in_count = fft_size;
    int out_count = fft_size / 2 + 1;

    vector<double> in_real_fftw(in_count);
    vector<double> in_real_ff(in_count);
    vector<complex<double>> out_complex_fftw(out_count);
    vector<complex<double>> out_complex_ff(out_count);

    fftw_plan plan = fftw_plan_dft_r2c_1d(
        in_count,
        in_real_fftw.data(),
        reinterpret_cast<fftw_complex*>(out_complex_fftw.data()),
        FFTW_ESTIMATE);

    auto ff_real_to_complex = FatFourier::real_to_complex(in_count);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    unsigned int time_fftw = 0;
    unsigned int time_ff = 0;
    for (int i = 0; i < nb_iters; i++) {
        // fill input buffer with random values
        for (int i = 0; i < in_count; i++) {
            double real = rand() / RAND_MAX;
            in_real_fftw[i] = real;
        }

        start = std::chrono::system_clock::now();
        fftw_execute(plan);
        end = std::chrono::system_clock::now();

        time_fftw += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        start = std::chrono::system_clock::now();
        ff_real_to_complex(in_real_ff.data(), out_complex_ff.data());
        end = std::chrono::system_clock::now();

        time_ff += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        if (!check_values(out_complex_ff.data(), out_complex_fftw.data(), out_count)) {
            std::cerr << "Incorrect result found\n";
            fftw_destroy_plan(plan);
            exit(EXIT_FAILURE);
        }
    }

    fftw_destroy_plan(plan);
    return { time_fftw, time_ff };
}

tuple<int, int> benchmark_complex(int fft_size, int nb_iters)
{
    int in_count = fft_size;
    int out_count = in_count + 1;

    vector<complex<double>> in_complex_fftw(in_count);
    vector<complex<double>> in_complex_ff(in_count);
    vector<complex<double>> out_complex_fftw(out_count);
    vector<complex<double>> out_complex_ff(out_count);

    fftw_plan plan = fftw_plan_dft_1d(
        in_count,
        reinterpret_cast<fftw_complex*>(in_complex_fftw.data()),
        reinterpret_cast<fftw_complex*>(out_complex_fftw.data()),
        FFTW_FORWARD,
        FFTW_ESTIMATE);

    auto ff_complex_to_complex = FatFourier::complex_to_complex(in_count);

    std::chrono::time_point<std::chrono::system_clock> start, end;
    unsigned int time_fftw = 0;
    unsigned int time_ff = 0;
    for (int i = 0; i < nb_iters; i++) {
        // fill input buffer with random values
        for (int i = 0; i < in_count; i++) {
            double real = rand() / RAND_MAX;
            double imag = rand() / RAND_MAX;
            in_complex_fftw[i].real(real);
            in_complex_ff[i].imag(imag);
        }

        start = std::chrono::system_clock::now();
        fftw_execute(plan);
        end = std::chrono::system_clock::now();

        time_fftw += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        start = std::chrono::system_clock::now();
        ff_complex_to_complex(in_complex_ff.data(), out_complex_ff.data());
        end = std::chrono::system_clock::now();

        time_ff += std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();

        if (!check_values(out_complex_ff.data(), out_complex_fftw.data(), out_count)) {
            std::cerr << "Incorrect result found\n";
            fftw_destroy_plan(plan);
            exit(EXIT_FAILURE);
        }
    }

    fftw_destroy_plan(plan);
    return { time_fftw, time_ff };
}

int main(int argc, char* argv[])
{
    int fft_size = FFT_SIZE_DEFAULT;
    int nb_iters = NB_ITERS_DEFAULT;
    bool complex_fft = false;

    char opt;
    while ((opt = getopt(argc, argv, "cl:n:")) != -1) {
        switch (opt) {
        case 'l':
            fft_size = atoi(optarg);
            if (!IS_POWER_OF_2(fft_size)) {
                std::cerr << "FFT must be power of 2\n";
                usage(argv[0]);
                exit(EXIT_FAILURE);
            }
            break;

        case 'n':
            nb_iters = atoi(optarg);
            if (nb_iters <= 0) {
                std::cerr << "Number of iterations must be greater than 0\n";
                usage(argv[0]);
                exit(EXIT_FAILURE);
            }
            break;

        case 'c':
            complex_fft = true;
            break;

        default:
            usage(argv[0]);
            exit(EXIT_FAILURE);
        }
    }

    unsigned int time_fftw;
    unsigned int time_ff;

    if (complex_fft) {
        std::tie(time_fftw, time_ff) = benchmark_complex(fft_size, nb_iters);
    } else {
        std::tie(time_fftw, time_ff) = benchmark_real(fft_size, nb_iters);
    }

    std::cout << "Computed " << nb_iters << " FFT(s) of " << fft_size << " samples\n";
    std::cout << "FFTW:\t\t" << time_fftw << " μs\n";
    std::cout << "FatFourier:\t" << time_ff << " μs\n";
}
