## Sluggishly computing your FFTs since 2017

This SSE-less, SSE-less, AVX-less, AVX2-less FFT will happily compute the frequency representation of whatever time series you input it, be it complex or real, as long as it has a size of 2ⁿ. It uses a radix-2 Cooley–Tukey algorithm and it is thread-safe. Because sluggish doesn't mean lazy, this FFT tries its best at implementing some of the optimizations described in [The FFT Demystified][2].

[2]: https://web.archive.org/web/20180307115735/http://www.engineeringproductivitytools.com/stuff/T0001/index.html

### Half real FFT

If our input signal is real-only (as in the case of an digital audio signal) then we can take advantage of the symmetry in the results: the second half represents the magnitude for negative frequencies and it is symmetrical to the first half. By building a new sequence of complex values half the length of the original input sequence, with even-index values as real parts and odd-index values as imaginary parts, we can obtain the same information, but computing only an half-size FFT with a slightly altered butterfly computation.

### Twiddle factor table

The butterfly phase, performed in each FFT at each recursion stage, involves some twiddle factors for each value in the sequence. These factors are expensive to compute but they are actually shared between all recursive FFTs: the twiddle factor at index *n/4* of a FFT of depth *d* is identical to the one at index *n/2* of its enclosing FFT of depth *d-1*. Precomputing and storing beforehand the twiddle factors hence allows for another easy performance gain.

### Bit-reversal sort

The sorting of the sequence between even and odd index values has to be performed at each recursion stage and leads to expensive copy operations. A better solution is to use [bit-reversal][3] permutation, which will rearrange once and for all the input sequence in a way that, at each recursion stage, values that were initial even-indexed will be located in the first half and thoses odd-indexed in the second half. Once the bit-reversal permutation is done, the input sequence remains as-is in memory, each recursive FFT only needing to know the position and length of the values it will handle.

[3]: https://graphics.stanford.edu/~seander/bithacks.html#BitReverseObvious
