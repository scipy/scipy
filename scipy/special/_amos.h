/*
 *
 * This file accompanied with the implementation file _amos.c is a
 * C translation of the Fortran code written by D.E. Amos with the
 * following original description:
 *
 *
 * A Portable Package for Bessel Functions of a Complex Argument
 * and Nonnegative Order
 * 
 * This algorithm is a package of subroutines for computing Bessel
 * functions and Airy functions.  The routines are updated
 * versions of those routines found in TOMS algorithm 644.
 *
 * Disclaimer:
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                 ISSUED BY SANDIA LABORATORIES,
 *                   A PRIME CONTRACTOR TO THE
 *               UNITED STATES DEPARTMENT OF ENERGY
 * * * * * * * * * * * * * *  NOTICE   * * * * * * * * * * * * * * *
 * THIS REPORT WAS PREPARED AS AN ACCOUNT OF WORK SPONSORED BY THE
 * UNITED STATES GOVERNMENT.  NEITHER THE UNITED STATES NOR THE
 * UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
 * EMPLOYEES, NOR ANY OF THEIR CONTRACTORS, SUBCONTRACTORS, OR THEIR
 * EMPLOYEES, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 * LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS
 * OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT OR PROCESS
 * DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT INFRINGE
 * PRIVATELY OWNED RIGHTS.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * THIS CODE HAS BEEN APPROVED FOR UNLIMITED RELEASE.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *
 *
 * The original Fortran code can be found at https://www.netlib.org/amos/
 *
 * References:
 * 
 * [1]: Abramowitz, M. and Stegun, I. A., Handbook of Mathematical
 *      Functions, NBS Applied Math Series 55, U.S. Dept. of Commerce,
 *      Washington, D.C., 1955
 * 
 * [2]: Amos, D. E., Algorithm 644, A Portable Package For Bessel
 *      Functions of a Complex Argument and Nonnegative Order, ACM
 *      Transactions on Mathematical Software, Vol. 12, No. 3,
 *      September 1986, Pages 265-273, DOI:10.1145/7921.214331
 * 
 * [3]: Amos, D. E., Remark on Algorithm 644, ACM Transactions on
 *      Mathematical Software, Vol. 16, No. 4, December 1990, Page
 *      404, DOI:10.1145/98267.98299
 * 
 * [4]: Amos, D. E., A remark on Algorithm 644: "A portable package
 *      for Bessel functions of a complex argument and nonnegative order",
 *      ACM Transactions on Mathematical Software, Vol. 21, No. 4,
 *      December 1995, Pages 388-393, DOI:10.1145/212066.212078
 * 
 * [5]: Cody, W. J., Algorithm 665, MACHAR: A Subroutine to
 *      Dynamically Determine Machine Parameters, ACM Transactions on
 *      Mathematical Software, Vol. 14, No. 4, December 1988, Pages
 *      303-311, DOI:10.1145/50063.51907
 *
 */

/*
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
 */
 
#ifndef _AMOS_H
#define _AMOS_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include <math.h>
#include <complex.h>

double complex amos_airy(double complex, int, int, int *, int *);
int amos_besh(double complex, double, int, int, int, double complex *, int *);
int amos_besi(double complex, double, int, int, double complex *, int *);
int amos_besj(double complex, double, int, int, double complex *, int *);
int amos_besk(double complex, double, int, int, double complex *, int *);
int amos_besy(double complex, double, int, int, double complex *, int *);
double complex amos_biry(double complex,int, int, int *);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* ifndef */
