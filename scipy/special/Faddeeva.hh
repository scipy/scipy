/* Copyright (c) 2012 Massachusetts Institute of Technology
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE. 
 */

/* Available at: http://ab-initio.mit.edu/Faddeeva

   Header file for Faddeeva.cc; see that file for more information. */

#ifndef FADDEEVA_HH
#define FADDEEVA_HH 1

#include <complex>

namespace Faddeeva {

// compute w(z) = exp(-z^2) erfc(-iz) [ Faddeeva / scaled complex error func ]
extern std::complex<double> w(std::complex<double> z,double relerr=0);
extern double w_im(double x); // special-case code for Im[w(x)] of real x

// Various functions that we can compute with the help of w(z)

// compute erfcx(z) = exp(z^2) erfc(z)
extern std::complex<double> erfcx(std::complex<double> z, double relerr=0);
extern double erfcx(double x); // special case for real x

// compute erf(z), the error function of complex arguments
extern std::complex<double> erf(std::complex<double> z, double relerr=0);
extern double erf(double x); // special case for real x

// compute erfi(z) = -i erf(iz), the imaginary error function
extern std::complex<double> erfi(std::complex<double> z, double relerr=0);
extern double erfi(double x); // special case for real x

// compute erfc(z) = 1 - erf(z), the complementary error function
extern std::complex<double> erfc(std::complex<double> z, double relerr=0);
extern double erfc(double x); // special case for real x

// compute Dawson(z) = sqrt(pi)/2  *  exp(-z^2) * erfi(z)
extern std::complex<double> Dawson(std::complex<double> z, double relerr=0);
extern double Dawson(double x); // special case for real x

} // namespace Faddeeva

#endif // FADDEEVA_HH
