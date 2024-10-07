/*
 *
 * This file accompanied with the translation unit cdflib.c is a C rewrite of
 * the Fortran code with the following as its original description:
 *
 * Cumulative distribution functions, inverses and parameters for
 * Beta, Binomial, Chi-square, noncentral Chi-square, F, noncentral F, Gamma,
 * negative Binomial, Normal, Poisson, Student's t distributions.
 * It uses various TOMS algorithms and Abramowitz & Stegun, also Bus Dekker
 * zero-finding algorithm.
 *
 * The original Fortran code can be found at Netlib
 * https://www.netlib.org/random/
 *
 *
 * References
 * ----------
 *
 *  J. C. P. Bus, T. J. Dekker, Two Efficient Algorithms with Guaranteed
 *  Convergence for Finding a Zero of a Function, ACM Trans. Math. Software 1:4
 *  (1975) 330-345, DOI:10.1145/355656.355659
 *
 *  M. Abramowitz and I. A. Stegun (Eds.) (1964) Handbook of Mathematical
 *  Functions with Formulas, Graphs, and Mathematical Tables. National Bureau
 *  of Standards Applied Mathematics Series, U.S. Government Printing Office,
 *  Washington, D.C..
 *
 */

/*
 *
 * Copyright (C) 2024 SciPy developers
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * a. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * b. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * c. Names of the SciPy Developers may not be used to endorse or promote
 *    products derived from this software without specific prior written
 *    permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 * OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef CDFLIB_H
#define CDFLIB_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include <math.h>

struct TupleDD
{
    double d1;
    double d2;
};

struct TupleDI
{
    double d1;
    int i1;
};

struct TupleDDI
{
    double d1;
    double d2;
    int i1;
};

struct TupleDID
{
    double d1;
    int i1;
    double d2;
};

struct TupleDDID
{
    double d1;
    double d2;
    int i1;
    double d3;
};


typedef struct DinvrState DinvrState;
typedef struct DzrorState DzrorState;


struct TupleDDID cdfbet_which1( double, double, double, double);
struct TupleDDID cdfbet_which2(double, double, double, double);
struct TupleDID cdfbet_which3(double, double, double, double, double);
struct TupleDID cdfbet_which4(double, double, double, double, double);
struct TupleDDID cdfbin_which1(double, double, double, double);
struct TupleDID cdfbin_which2(double, double, double, double, double);
struct TupleDID cdfbin_which3(double, double, double, double, double);
struct TupleDDID cdfbin_which4(double, double, double, double);
struct TupleDDID cdfchi_which1(double, double);
struct TupleDID cdfchi_which2(double, double, double);
struct TupleDID cdfchi_which3(double, double, double);
struct TupleDDID cdfchn_which1(double, double, double);
struct TupleDID cdfchn_which2(double, double, double);
struct TupleDID cdfchn_which3(double, double, double);
struct TupleDID cdfchn_which4(double, double, double);
struct TupleDDID cdff_which1(double, double, double);
struct TupleDID cdff_which2(double, double, double, double);
struct TupleDID cdff_which3(double, double, double, double);
struct TupleDID cdff_which4(double, double, double, double);
struct TupleDDID cdffnc_which1(double, double, double, double);
struct TupleDID cdffnc_which2(double, double, double, double, double);
struct TupleDID cdffnc_which3(double, double, double, double, double);
struct TupleDID cdffnc_which4(double, double, double, double, double);
struct TupleDID cdffnc_which5(double, double, double, double, double);
struct TupleDDID cdfgam_which1(double, double, double);
struct TupleDID cdfgam_which2(double, double, double, double);
struct TupleDID cdfgam_which3(double, double, double, double);
struct TupleDID cdfgam_which4(double, double, double, double);
struct TupleDDID cdfnbn_which1(double, double, double, double);
struct TupleDID cdfnbn_which2(double, double, double, double, double);
struct TupleDID cdfnbn_which3(double, double, double, double, double);
struct TupleDDID cdfnbn_which4(double, double, double, double);
struct TupleDDID cdfnor_which1(double, double, double);
struct TupleDID cdfnor_which2(double, double, double, double);
struct TupleDID cdfnor_which3(double, double, double, double);
struct TupleDID cdfnor_which4(double, double, double, double);
struct TupleDDID cdfpoi_which1(double, double);
struct TupleDID cdfpoi_which2(double, double, double);
struct TupleDID cdfpoi_which3(double, double, double);
struct TupleDDID cdft_which1(double, double);
struct TupleDID cdft_which2(double, double, double);
struct TupleDID cdft_which3(double, double, double);
struct TupleDDID cdftnc_which1(double, double, double);
struct TupleDID cdftnc_which2(double, double, double, double);
struct TupleDID cdftnc_which3(double, double, double, double);
struct TupleDID cdftnc_which4(double, double, double, double);

#ifdef __cplusplus
}      /* extern "C" */
#endif /* __cplusplus */

#endif /* ifndef */
