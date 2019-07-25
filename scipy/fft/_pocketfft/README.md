pypocketfft
===========

This package provides Fast Fourier, trigonometric and Hartley transforms with a
simple Python interface.

The central algorithms are derived from Paul Swarztrauber's FFTPACK code
(http://www.netlib.org/fftpack).

Features
--------
- supports fully complex and half-complex (i.e. complex-to-real and
  real-to-complex) FFTs, discrete sine/cosine transforms and Hartley transforms
- achieves very high accuracy for all transforms
- supports multidimensional arrays and selection of the axes to be transformed
- supports single, double, and long double precision
- makes use of CPU vector instructions when performing 2D and higher-dimensional
  transforms
- supports prime-length transforms without degrading to O(N**2) performance
- has optional OpenMP support for multidimensional transforms
