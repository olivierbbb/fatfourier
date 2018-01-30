#include <cassert>
#if defined (OMP)
    #include <omp.h>
#endif
#include "fft.hpp"
#include <iostream>
#define IS_POWER_OF_2(x) (((x) & ((x) - 1)) == 0)

using namespace std;

FFT::FFT(double* buffer, uint64_t buffer_size) :
    reals_(buffer),
    reals_count_(buffer_size),
    // cast real array to half-size complex array (they are binary-compatible)
    complexes_(reinterpret_cast<complex<double>*>(buffer)),
    complexes_count_(buffer_size / 2)
{
    assert(IS_POWER_OF_2(buffer_size));

    twiddles_ = new complex<double>[buffer_size / 2];
    for (uint64_t k = 0; k < buffer_size / 2; k++) {
        twiddles_[k] = polar(1.0, -2.0 * M_PI * k / buffer_size);
    }

    #ifdef OMP
        omp_set_num_threads(2);
    #elif defined THREAD
        #ifdef THREADS_COUNT
            threads_count_ = THREADS_COUNT;
        #else
            threads_count_ = thread::hardware_concurrency();
            threads_count_ = pow(2, log2(threads_count_));
        #endif
        assert(IS_POWER_OF_2(threads_count_));
    #endif
}

FFT::~FFT() {
    delete[] twiddles_;
}

/**
 * Cree une nouvelle serie de complexes a partir du buffer d'input
 * faite des valeurs aux index pairs en partie reelme et des valeurs aux index
 * impairs en partie imaginaire, et calcule la FFT de cette serie.
 * @see http://www.engineeringproductivitytools.com/stuff/T0001/PT10.HTM#Head574
 * About k == 0 et k == count / 2, @see http://processors.wiki.ti.com/index.php/Efficient_FFT_Computation_of_Real_Input
 */
void FFT::compute() {
    bit_reverse_permute_buffer_();

    #if defined(OMP)
        omp_fft_(complexes_, complexes_count_);
    #elif defined(THREAD)
        threaded_fft_(complexes_, complexes_count_);
    #else
        fft_(complexes_, complexes_count_);
    #endif
    
    // Recompute real values from complex values
    // (we are not interested in the second half of real values)
    complex<double> even, odd, twiddle, result;
    // k == 0
    even = complexes_[0] + conj(complexes_[0]);
    odd = -1i * (complexes_[0] - conj(complexes_[0]));
    twiddle = polar(1.0, -2.0 * M_PI);
    result = 0.5 * (even + twiddle * odd);
    reals_[0] = result.real();

    // k == 1 .. count / 2 - 1
    for (uint64_t k = 1; k < reals_count_ / 2; k++) {
        even = complexes_[k] + conj(complexes_[complexes_count_ - k]);
        odd = -1i * (complexes_[k] - conj(complexes_[complexes_count_ - k]));
        twiddle = twiddles_[k];

        result = 0.5 * (even + twiddle * odd);
        reals_[k] = result.real();
    }

    // k == count / 2
    result = 0.5 * (complexes_[0] + conj(complexes_[0]) + 1i * (complexes_[0] - conj(complexes_[0])));
    reals_[reals_count_ / 2] = result.real();
}

/**
 * Compute fft of @p values serie of length @p count, storing results in-place in @p values
 * @p depth: recursivity depth (used to retrieve twiddle factor)
 */
void FFT::fft_(complex<double>* values, uint64_t count, unsigned int depth) {
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
 * @returns @p bit_length first bits of @p value mirror-reversed
 * @see https://graphics.stanford.edu/~seander/bithacks.html#BitReverseObvious
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

/**
 * Performs in place bit reversal permutation on @p buffer
 */
void FFT::bit_reverse_permute_buffer_() {
    complex<double> swap;
    uint64_t reversed_k;
    const unsigned int bit_length = log2(complexes_count_);
    for (uint64_t k = 0; k < complexes_count_; k++) {
        reversed_k = reverse_bits(k, bit_length);
        if (reversed_k > k) { // dont swap twice
            swap = complexes_[k];
            complexes_[k] = complexes_[reversed_k];
            complexes_[reversed_k] = swap;
        }
    }
}

#ifdef OMP

void FFT::omp_fft_(complex<double>* values, uint64_t count, unsigned int depth) {
    assert(IS_POWER_OF_2(count));

    if (count < 2) {
        return;
    }

    const uint64_t half_count = count / 2;
    complex<double>* evens = values;
    complex<double>* odds = values + half_count;

    // un thread pour chaque sous-fft
    #pragma omp parallel num_threads(2)
    {
        if (omp_get_thread_num() % 2 == 0) {
            fft_(evens, half_count, depth + 1);
        } else {
            fft_(odds, half_count, depth + 1);
        }
    }

    // un seul thread pour les butterflies
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

#elif defined THREAD

// TODO use iterative algorithm instead of recursivity
/*void FFT::threaded_fft_() {
    const uint64_t half_count = complexes_count_ / 2;
    complex<double>* evens = complexes_;
    complex<double>* odds = complexes_ + half_count;

    for (int i = 0; i < ; i++) {
        // stacks butterflies computations to be done into closures, with all needed context
        async_stack_.emplace_back([this, evens, odds, half_count, depth] {
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

        evens_threads_.emplace_back(&FFT::fft_, this, evens, half_count, depth + 1);
        odds_threads_.emplace_back(&FFT::fft_, this, odds, half_count, depth + 1);
    }

    for (int i = 0; i < evens_threads_.size(); i++) {
        evens_threads_[i].join();
        odds_threads_[i].join();
    }
    evens_threads_.clear();
    odds_threads_.clear();

    while (async_stack_.size() > 0) {
        async_stack_.back()();
        async_stack_.pop_back();
    }
}*/


void FFT::threaded_fft_(complex<double>* values, uint64_t count, unsigned int depth) {
    assert(IS_POWER_OF_2(count));

    const uint64_t half_count = count / 2;
    complex<double>* evens = values;
    complex<double>* odds = values + half_count;

    for (unsigned int i = 0; i < threads_count_; i++) {
    if (pow(2, depth) == threads_count_) {
        async_stack_.emplace_back([this, evens, odds, half_count, depth] {
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

        evens_threads_.emplace_back(&FFT::fft_, this, evens, half_count, depth + 1);
        odds_threads_.emplace_back(&FFT::fft_, this, odds, half_count, depth + 1);
    } else {
        // stacks butterflies computations to be done into closures, with all needed context
        // NB: must be stacked BEFORE recursive calls...
        async_stack_.emplace_back([this, evens, odds, half_count, depth] {
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

        threaded_fft_(evens, half_count, depth + 1);
        threaded_fft_(odds, half_count, depth + 1);
    }

    if (depth == 1) {
        assert(evens_threads_.size() == odds_threads_.size());
        for (int i = 0; i < evens_threads_.size(); i++) {
            evens_threads_[i].join();
            odds_threads_[i].join();
        }
        evens_threads_.clear();
        odds_threads_.clear();

        while (async_stack_.size() > 0) {
            async_stack_.back()();
            async_stack_.pop_back();
        }
    }
}

#endif
