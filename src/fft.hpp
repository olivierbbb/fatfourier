#pragma once

#include <cstdint>
#include <complex>
#ifdef THREAD
    #include <vector>
    #include <functional>
    #include <thread>
#endif

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

    #ifdef THREAD
        /** Number of threads to spawn (not including main thread) */
        unsigned int threads_count_;
        /** Post-ffts bufferflies to compute once all threads are done */
        std::vector<std::function<void()>> async_stack_;
        /** Threads handling ffts of even-indexed (complex) values */
        std::vector<std::thread> evens_threads_;
        /** Threads handling ffts of odd-indexed (complex) values */
        std::vector<std::thread> odds_threads_;
    #endif

    void fft_(std::complex<double>* values, uint64_t count, unsigned int depth = 1);
    void half_real_fft_();
    void bit_reverse_permute_buffer_();

    #ifdef OMP
        void omp_fft_(std::complex<double>* values, uint64_t count, unsigned int depth = 1);
    #elif defined THREAD
        void threaded_fft_(std::complex<double>* values, uint64_t count, unsigned int depth = 1);
    #endif
};
