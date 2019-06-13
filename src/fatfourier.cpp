#include "fatfourier.hpp"
#include <cassert>
#include <cstdint>

#define IS_POWER_OF_2(x) (((x) & ((x) -1)) == 0)

namespace fatfourier {

using std::complex;
using std::vector;
using namespace std::complex_literals;

static uint64_t reverse_bits(uint64_t value, unsigned int bit_length);
template <class T>
void bit_reverse_permute(T* buffer, size_t size);

FatFourier::FatFourier(size_t size)
    : in_count(size),
      out_count(size + 1),
      twiddles(size)
{
    if (!IS_POWER_OF_2(size)) {
        throw std::invalid_argument("FFT size must be power of 2");
    }

    for (size_t k = 0; k < in_count; k++) {
        twiddles[k] = std::polar(1.0, -1.0 * M_PI * (double) k / (double) size);
    }
}

RealToComplex FatFourier::real_to_complex(size_t in_count)
{
    FatFourier ff = FatFourier(in_count / 2);
    return [ff](double* in_real, std::complex<double>* out_complex) {
        ff.compute(in_real, out_complex);
    };
}

ComplexToComplex FatFourier::complex_to_complex(size_t in_count)
{
    FatFourier ff = FatFourier(in_count);
    return [ff](std::complex<double>* in_complex, std::complex<double>* out_complex) {
        ff.compute(in_complex, out_complex);
    };
}

/**
@see http://www.engineeringproductivitytools.com/stuff/T0001/PT10.HTM#Head574
About k == 0 and k == count / 2, @see http://processors.wiki.ti.com/index.php/Efficient_FFT_Computation_of_Real_Input
*/
void FatFourier::compute(double* in_real, complex<double>* out_complex) const
{
    // cast real array to half-size complex array (they are binary-compatible)
    complex<double>* in_complex = reinterpret_cast<complex<double>*>(in_real);

    bit_reverse_permute(in_complex, in_count);
    fft_recursive(in_complex, in_count);

    complex<double> in_conj, even, odd, twiddle, result;
    // k == 0
    in_conj = std::conj(in_complex[0]);
    even = in_complex[0] + in_conj;
    odd = -1i * (in_complex[0] - in_conj);
    twiddle = std::polar(1.0, -2.0 * M_PI);
    result = 0.5 * (even + twiddle * odd);
    out_complex[0] = result;

    // k == 1 .. N - 1
    for (size_t k = 1; k < in_count; k++) {
        in_conj = std::conj(in_complex[in_count - k]);
        even = in_complex[k] + in_conj;
        odd = -1i * (in_complex[k] - in_conj);
        twiddle = twiddles[k];

        result = 0.5 * (even + twiddle * odd);
        out_complex[k] = result;
    }

    // k == N
    in_conj = std::conj(in_complex[0]);
    even = in_complex[0] + in_conj;
    odd = 1i * (in_complex[0] - in_conj);
    result = 0.5 * (even + odd);
    out_complex[out_count] = result;
}

void FatFourier::compute(complex<double>* in_complex, complex<double>* out_complex) const
{
    bit_reverse_permute(in_complex, in_count);
    fft_recursive(in_complex, in_count);
}
/**
Compute fft of @p values serie of length @p count, storing results in-place in @p values
@p depth: recursivity depth (used to retrieve twiddle factor)
*/
void FatFourier::fft_recursive(complex<double>* values, size_t count, unsigned int depth) const
{
    assert(IS_POWER_OF_2(count));

    if (count < 2) {
        return;
    }

    const size_t half_count = count / 2;
    complex<double>* evens = values;
    complex<double>* odds = values + half_count;

    fft_recursive(evens, half_count, depth + 1);
    fft_recursive(odds, half_count, depth + 1);

    // locating twiddle factor with depth
    const unsigned int twiddle_index_factor = 1 << depth;
    complex<double> even, odd, twiddle;
    for (size_t k = 0; k < half_count; k++) {
        even = evens[k];
        odd = odds[k];
        twiddle = twiddles[k * twiddle_index_factor];

        evens[k] = even + twiddle * odd;
        odds[k] = even - twiddle * odd;
    }
}

/** Performs in place bit reversal permutation on input buffer */
template <class T>
void bit_reverse_permute(T* buffer, size_t size)
{
    const unsigned int bit_length = std::log2(size);

    T swap;
    uint64_t reversed_k;
    for (size_t k = 0; k < size; k++) {
        reversed_k = reverse_bits(k, bit_length);
        if (reversed_k > k) { // dont swap twice
            swap = buffer[k];
            buffer[k] = buffer[reversed_k];
            buffer[reversed_k] = swap;
        }
    }
}

/**
@returns @p bit_length first bits of @p value mirror-reversed
@see https://graphics.stanford.edu/~seander/bithacks.html#BitReverseObvious
*/
uint64_t reverse_bits(uint64_t value, unsigned int bit_length)
{
    uint64_t reversed = value & 1;
    unsigned int shift = bit_length - 1;
    for (value >>= 1; value; value >>= 1) {
        reversed <<= 1;
        reversed |= value & 1;
        shift--;
    }
    return reversed << shift;
}
}
