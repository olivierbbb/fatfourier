#pragma once
#include <complex>
#include <functional>
#include <vector>

namespace fatfourier {

using RealToComplex = std::function<void(double* in_real, std::complex<double>* out_complex)>;
using ComplexToComplex = std::function<void(std::complex<double>* in_complex, std::complex<double>* out_complex)>;

class FatFourier
{
public:
    static RealToComplex real_to_complex(size_t in_count);
    static ComplexToComplex complex_to_complex(size_t in_count);

private:
    const size_t in_count;
    const size_t out_count;
    /**
    Twiddle factors lookup tabble
    @see http://www.engineeringproductivitytools.com/stuff/T0001/PT04.HTM#Head363
    */
    std::vector<std::complex<double>> twiddles;

    FatFourier(size_t in_count);
    /** Real to complex computation */
    void compute(double* in_real, std::complex<double>* out_complex) const;
    /** Complex to complex computation */
    void compute(std::complex<double>* in_complex, std::complex<double>* out_complex) const;
    void fft_recursive(std::complex<double>* values, size_t count, unsigned int depth = 1) const;
};
}
