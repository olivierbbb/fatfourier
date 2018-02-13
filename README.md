## How it works

The Fourier transform maps values in time to their frequency counterparts. Applied to a discreet signal of *N* real values, a DFT (Discrete Fourier Transform) returns a serie of *N* complex values representings the energy magnitude (real part) and phase (imaginary part) at frequency *n/N*.

While a naive implemention may have a complexity of *O(n^2)*, the Cooley-Tukey Fast Fourier Transform is *O(n log n)*. It splits the series between odd-indexed and even-indexed values, before calling itself recursively on both. The partial results are then combined with a [butterfly][1] computation to form the full results series.

[1]: https://en.wikipedia.org/wiki/Butterfly_diagram

## Sequential optimisations

Before parallelising the DFT, there are some low-hanging fruits that can be worked on for the sequential routine:

* __half real FFT__.

* __bit reversal sort__

## Parallelisation

## Results



