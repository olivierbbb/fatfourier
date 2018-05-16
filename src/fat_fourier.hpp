#pragma once

#include <complex>
#include <cstdint>

class FatFourier {
public:
    FatFourier(double* in_buffer, double* out_buffer, uint64_t buffer_size);
    ~FatFourier();
#ifdef OMP
    void compute();
#else
    /** @p threads_count: Number of threads to spawn (not including main thread).
    If equal to 0, will be set automatically */
    void compute(unsigned int threads_count = 0);
#endif

private:
    /** Input and output buffers as arrays of real values */
    double* in_buffer_;
    double* out_buffer_;
    uint64_t buffer_size_;

    /** Pointer to input buffer as array of complex values */
    std::complex<double>* complexes_;
    uint64_t complexes_count_;

    /**
    Twiddle factors lookup tabble
    @see http://www.engineeringproductivitytools.com/stuff/T0001/PT04.HTM#Head363
    */
    std::complex<double>* twiddles_;

    void bit_reverse_permute_buffer_();
    void fft_(std::complex<double>* values, uint64_t count, unsigned int depth = 1);
    void parallel_fft_(unsigned int threads_count);
};
