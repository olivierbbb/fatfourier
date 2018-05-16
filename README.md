## Fast Fourier Transform

The Fourier transform maps values in time to their frequency counterparts. Applied to a discreet signal of *N* real values, Discrete Fourier Transform returns a sequence of *N* complex values representing the energy magnitude (real part) and phase (imaginary part) at frequency *n/N*. The longer the signal, the finer the spectral resolution.

While a naive implementation may have a complexity of *O(n²)*, the Cooley-Tukey Fast Fourier Transform is *O(n log n)*. It uses a divide-and-conquer approach by diving the input sequence in odd-index and even-index values and calling itself recursively on both half-length sequences. The partial results obtained are then combined with a [butterfly][1] computation to produce the full results.

[1]: https://en.wikipedia.org/wiki/Butterfly_diagram

## Sequential optimizations

Before parallelizing the DFT, there are some low-hanging fruits that can be worked on for the sequential routine (see [The FFT Demystified][2]):

[2]: http://www.engineeringproductivitytools.com/stuff/T0001/index.html

### Half real FFT

If we assume that our input signal has no imaginary part (no phase information), then we can take advantage of the symmetry in the results: the second half represents the magnitude for negative frequencies and it is symmetrical to the first half. By building a new sequence of complex values half the length of the original input sequence, with even-index values as real parts and odd-index values as imaginary parts, we may obtain the same information, but computing only an half-size FFT with a slightly altered butterfly computation.

### Twiddle factor table

The butterfly phase, performed in each FFT and the computation of the two nested FFTs, involves some twiddle factors for each value in the sequence. Theses factors are expensive to compute but they are actually shared between all recursive FFTs: the twiddle factor at index *n/4* of a FFT of depth *d* is identical to the one at index *n/2* of its enclosing FFT of depth *d-1*. Precomputing and storing beforehand the twiddle factors therefor allows for another performance gain.

### Bit-reversal sort

The sorting of the sequence between even and odd index values has to be performed at each recursion stage and leads to expensive copy operations. A better solution is to use [bit-reversal][3] permutation, which will rearrange once and for all the input sequence in a way that, at each recursion stage, values that were initial even-indexed will be located in the first half and thoses odd-indexed in the second half. Once the bit-reversal permutation is done, the input sequence remains as-is in memory, each recursive FFT only needing to know the position and length of the values it will handle.

[3]: https://graphics.stanford.edu/~seander/bithacks.html#BitReverseObvious

## Parallelization

The divide-and-conquer nature of the Cooley-Tukey algorithm makes it an easy candidate for parallel execution. The volume of the input and output data makes us prefer a shared memory solution such as C++ threads or OpenMP over a distributed-memory system like MPI.

As a first step, an easy but limited solution consists in parallelizing only the root FFT by sending the each of the two nested FFTs on a different process and then combining the results on the main process. Generalizing the solution to a higher number of processes (constrained to a power of two), two parallelized FFTs have to been started at each recursive FFT of depth equal to the square of the number of processes. The processes are not spawned anymore from the same function call, but rather from different function calls at the same recursion depth.

The tricky part lays in the butterfly phases that are left to perform once all the nested parallelized FFTs have been computed. Theses computations require some context that was available in the local scopes of each recursive call,but theses scopes are gone because the processes are spawned in a non-blocking way and their termination is waited for in the root FFT call.

The solution is to move from a recursive approach to an iterative one for all "pre-parallel" FFTs performed on the main process. A FIFO queue stores the input parameters of each of theses FFTs. An iterative loop pops them and pushes back on the queue new input parameters for each pair of nested FFT, unless the parallelization depth has been reached, in which case two threads are spawned to compute them. A LIFO stack stores the closures that will perform the post-nested-FFTs computations, each closure embedding all the context needed. When all the sub-FFT threads are done, the closures on the LIFO stack are popped and executed:

```python
PARALLEL_FFT(ALL_DATA):
    D_QUEUE ← {ALL_DATA}
    B_STACK ← {}
    while D_QUEUE ≠ ∅
        DATA ← POP_FRONT(D_QUEUE)
        EVEN_HALF, ODD_HALF ← SPLIT(DATA)
        if DEPTH² < NB_THREADS
            D_QUEUE ← D_QUEUE ∪ {EVEN_HALF, ODD_HALF}
            BUTTERFLY ← INIT_BUTTERFLY(EVEN_HALF, ODD_HALF)
            B_STACK ← B_STACK ∪ {BUTTERFLY}
        else
            THREAD(FFT(EVEN_HALF))
            THREAD(FFT(ODD_HALF))

    WAIT_FOR_THREADS()
    while B_STACK ≠ ∅
        BUTTERFLY ← POP_BACK(B_STACK)
        COMPUTE_BUTTERFLY(BUTTERFLY)

    return ALL_DATA
```

The bit-reversal permutation greatly helps with the parallelization of the Cooley-Tukey algorithm. Once sorted, the memory layout of all the input data will not move. Results of each FFT can be stored in-place instead of the input data, without risking data race.

## Results

TODO plot results for sequential, C++ threads, OpenMP and FFTW.
SPOILER: it's faster (almost twice with 4 threads), but not as fast as a non-parallel FFT with [FFTW](http://www.fftw.org/)!
