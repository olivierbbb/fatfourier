#pragma once

#include <cstdint>
#include <complex>

class FFT {
public:
    FFT(double* in_buffer, double* out_buffer, uint64_t buffer_size);
    ~FFT();
    void compute();

private:
    /**Input and output buffers as arrays of real values */
    double* in_buffer_;
    double* out_buffer_;
    uint64_t buffer_size_;

    /** Pointer to input buffer as array of complex values */
    std::complex<double>* complexes_;
    uint64_t complexes_count_;
    
    /**
     * Twiddle factors lookup tabble
     * @see http://www.engineeringproductivitytools.com/stuff/T0001/PT04.HTM#Head363
     * */
    std::complex<double>* twiddles_;

    #if defined(OMP) || defined(THREADED)
        /** Number of threads to spawn (not including main thread) or number of threads in parallel sections */
        unsigned int threads_count_;
    #endif

    void fft_(std::complex<double>* values, uint64_t count, unsigned int depth = 1);
    void half_real_fft_();
    void bit_reverse_permute_buffer_();

    #if defined(OMP) || defined(THREADED)
        void parallel_fft_();
    #endif
};
