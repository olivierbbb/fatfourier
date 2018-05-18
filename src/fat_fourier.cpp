#include <cassert>
#include <iostream>

#include <functional>
#include <stack>
#include <vector>
#ifdef OMP
#include <omp.h>
#else
#include <thread>
#endif

#include "fat_fourier.hpp"

#define IS_POWER_OF_2(x) (((x) & ((x) -1)) == 0)

using namespace std;

FatFourier::FatFourier(double* in_buffer, double* out_buffer, uint64_t buffer_size)
    : in_buffer_(in_buffer)
    , out_buffer_(out_buffer)
    , buffer_size_(buffer_size)
    // cast real array to half-size complex array (they are binary-compatible)
    , complexes_(reinterpret_cast<complex<double>*>(in_buffer))
    , complexes_count_(buffer_size / 2) {
    if (!IS_POWER_OF_2(buffer_size)) {
        cerr << "Buffer size must be power of 2" << endl;
        exit(EXIT_FAILURE);
    }

    twiddles_ = new complex<double>[buffer_size / 2];

#ifdef OMP
#pragma omp parallel for
#endif
    for (uint64_t k = 0; k < buffer_size / 2; k++) {
        twiddles_[k] = polar(1.0, -2.0 * M_PI * k / buffer_size);
    }
}

FatFourier::~FatFourier() {
    delete[] twiddles_;
}

/**
@see http://www.engineeringproductivitytools.com/stuff/T0001/PT10.HTM#Head574
About k == 0 et k == count / 2, @see http://processors.wiki.ti.com/index.php/Efficient_FFT_Computation_of_Real_Input
*/
#ifdef OMP
void FatFourier::compute() {
    unsigned int procs_count = omp_get_num_procs();
    unsigned int threads_count = pow(2, log2(procs_count));
#else
void FatFourier::compute(unsigned int threads_count) {
    // auto threads count
    if (threads_count == 0) {
        unsigned int procs_count = thread::hardware_concurrency();
        // floor to power of 2
        threads_count = pow(2, log2(procs_count));
    } else if (!IS_POWER_OF_2(threads_count)) {
        cerr << "Number of threads must be power of 2" << endl;
        exit(EXIT_FAILURE);
    }
#endif

    bit_reverse_permute_buffer_();

    if (threads_count == 1) {
        fft_(complexes_, complexes_count_);
    } else {
        parallel_fft_(threads_count);
    }

    // Recompute real values from complex values
    // (we are not interested in the second half of real values)
    complex<double> even, odd, twiddle, result;
    // k == 0
    even = complexes_[0] + conj(complexes_[0]);
    odd = -1i * (complexes_[0] - conj(complexes_[0]));
    twiddle = polar(1.0, -2.0 * M_PI);
    result = 0.5 * (even + twiddle * odd);
    out_buffer_[0] = result.real();

    // k == 1 .. count / 2 - 1
    // TODO parallel for?
    for (uint64_t k = 1; k < buffer_size_ / 2; k++) {
        even = complexes_[k] + conj(complexes_[complexes_count_ - k]);
        odd = -1i * (complexes_[k] - conj(complexes_[complexes_count_ - k]));
        twiddle = twiddles_[k];

        result = 0.5 * (even + twiddle * odd);
        out_buffer_[k] = result.real();
    }

    // k == count / 2
    result = 0.5 * (complexes_[0] + conj(complexes_[0]) + 1i * (complexes_[0] - conj(complexes_[0])));
    out_buffer_[buffer_size_ / 2] = result.real();
}

/**
Compute fft of @p values serie of length @p count, storing results in-place in @p values
@p depth: recursivity depth (used to retrieve twiddle factor)
*/
void FatFourier::fft_(complex<double>* values, uint64_t count, unsigned int depth) {
    assert(IS_POWER_OF_2(count));

    if (count < 2) {
        return;
    }

    const uint64_t half_count = count / 2;
    complex<double>* evens = values;
    complex<double>* odds = values + half_count;

    fft_(evens, half_count, depth + 1);
    fft_(odds, half_count, depth + 1);

    // locating twiddle factor with depth
    unsigned int twiddle_index_multiplier = 1 << depth;
    complex<double> even, odd, twiddle;
    for (uint64_t k = 0; k < half_count; k++) {
        even = evens[k];
        odd = odds[k];
        twiddle = twiddles_[k * twiddle_index_multiplier];

        evens[k] = even + twiddle * odd;
        odds[k] = even - twiddle * odd;
    }
}

/**
@returns @p bit_length first bits of @p value mirror-reversed
@see https://graphics.stanford.edu/~seander/bithacks.html#BitReverseObvious
*/
static uint64_t reverse_bits(uint64_t value, unsigned int bit_length) {
    uint64_t reversed = value & 1;
    unsigned int shift = bit_length - 1;
    for (value >>= 1; value; value >>= 1) {
        reversed <<= 1;
        reversed |= value & 1;
        shift--;
    }
    return reversed << shift;
}

/** Performs in place bit reversal permutation on @p buffer */
void FatFourier::bit_reverse_permute_buffer_() {
    const unsigned int bit_length = log2(complexes_count_);

#ifdef OMP
#pragma omp parallel for
#endif
    for (uint64_t k = 0; k < complexes_count_; k++) {
        complex<double> swap;
        uint64_t reversed_k;
        reversed_k = reverse_bits(k, bit_length);
        if (reversed_k > k) { // dont swap twice
            swap = complexes_[k];
            complexes_[k] = complexes_[reversed_k];
            complexes_[reversed_k] = swap;
        }
    }
}

void FatFourier::parallel_fft_(unsigned int threads_count) {
#ifdef OMP
    vector<function<void()>> parallel_ffts;
#else
    vector<thread> parallel_ffts;
#endif
    // input values for async (pre-parallel) ffts
    stack<tuple<complex<double>*, uint64_t, unsigned int>> async_ffts;
    // butterfly closures for async (pre-parallel) ffts
    vector<function<void()>> async_butterflies;
    async_ffts.emplace(complexes_, complexes_count_, 1);

    complex<double>* values;
    uint64_t count;
    unsigned int depth;
    while (!async_ffts.empty()) {
        tie(values, count, depth) = async_ffts.top();
        async_ffts.pop();

        const uint64_t half_count = count / 2;
        complex<double>* evens = values;
        complex<double>* odds = values + half_count;

        if (pow(2, depth) < threads_count) {
            async_ffts.emplace(evens, half_count, depth + 1);
            async_ffts.emplace(odds, half_count, depth + 1);
        } else {
#ifdef OMP
            parallel_ffts.emplace_back([this, evens, half_count, depth] {
                this->fft_(evens, half_count, depth + 1);
            });
            parallel_ffts.emplace_back([this, odds, half_count, depth] {
                this->fft_(odds, half_count, depth + 1);
            });
#else
            parallel_ffts.emplace_back(&FatFourier::fft_, this, evens, half_count, depth + 1);
            parallel_ffts.emplace_back(&FatFourier::fft_, this, odds, half_count, depth + 1);
#endif
        }

        // stack butterflies computations to be done into closures, with all needed context,
        // to handle them when threads are done
        async_butterflies.emplace_back([this, evens, odds, half_count, depth] {
            unsigned int twiddle_index_multiplier = 1 << depth;
            complex<double> even, odd, twiddle;

            for (uint64_t k = 0; k < half_count; k++) {
                even = evens[k];
                odd = odds[k];
                twiddle = twiddles_[k * twiddle_index_multiplier];

                evens[k] = even + twiddle * odd;
                odds[k] = even - twiddle * odd;
            }
        });
    }

#ifdef OMP
#pragma omp parallel for
    for (int i = 0; i < parallel_ffts.size(); i++) {
        parallel_ffts[i]();
    }
#else
    for (int i = 0; i < parallel_ffts.size(); i++) {
        parallel_ffts[i].join();
    }
#endif

    // perform butterflies in inner to outer order
    while (!async_butterflies.empty()) {
        // TODO parallelize each pair?
        async_butterflies.back()();
        async_butterflies.pop_back();
    }
}
