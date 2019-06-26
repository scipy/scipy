pypocketfft
===========

This package provides Fast Fourier and Hartley transforms with a simple
Python interface.

The central algorithms are derived from Paul Swarztrauber's FFTPACK code
(http://www.netlib.org/fftpack).

Features
--------
- supports fully complex and half-complex (i.e. complex-to-real and
  real-to-complex) FFTs
- supports multidimensional arrays and selection of the axes to be transformed.
- supports single and double precision
- makes use of CPU vector instructions when performing 2D and higher-dimensional
  transforms
- does not have persistent transform plans, which makes the interface simpler
- supports prime-length transforms without degrading to O(N**2) performance
- Has optional OpenMP support for multidimensional transforms
