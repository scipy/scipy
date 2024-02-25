/*
 *
 * This file accompanied with the implementation file specfun.c is
 * a partial C translation of the Fortran code by Zhang and Jin following
 * original description:
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *       COMPUTATION OF SPECIAL FUNCTIONS
 *
 *          Shanjie Zhang and Jianming Jin
 *
 *       Copyrighted but permission granted to use code in programs.
 *       Buy their book:
 *
 *          Shanjie Zhang, Jianming Jin,
 *          Computation of Special Functions,
 *          Wiley, 1996,
 *          ISBN: 0-471-11963-6,
 *          LC: QA351.C45.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *       Scipy changes:
 *       - Compiled into a single source file and changed REAL To DBLE throughout.
 *       - Changed according to ERRATA.
 *       - Changed GAMMA to GAMMA2 and PSI to PSI_SPEC to avoid potential conflicts.
 *       - Made functions return sf_error codes in ISFER variables instead
 *         of printing warnings. The codes are
 *         - SF_ERROR_OK        = 0: no error
 *         - SF_ERROR_SINGULAR  = 1: singularity encountered
 *         - SF_ERROR_UNDERFLOW = 2: floating point underflow
 *         - SF_ERROR_OVERFLOW  = 3: floating point overflow
 *         - SF_ERROR_SLOW      = 4: too many iterations required
 *         - SF_ERROR_LOSS      = 5: loss of precision
 *         - SF_ERROR_NO_RESULT = 6: no result obtained
 *         - SF_ERROR_DOMAIN    = 7: out of domain
 *         - SF_ERROR_ARG       = 8: invalid input parameter
 *         - SF_ERROR_OTHER     = 9: unclassified error
 *       - Improved initial guesses for roots in JYZO.
 *
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

#ifndef SPECFUN_H
#define SPECFUN_H

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

#include <stdlib.h>
#include <math.h>
#include <complex.h>

void specfun_airyzo(int, int, double *, double *, double *, double *);
void specfun_aswfa(double, int, int, double, int, double, double *, double *);
void specfun_bernob(int, double *);
double complex specfun_cerror(double complex);
void specfun_cerzo(int, double complex *);
double complex specfun_cchg(double, double, double complex);
void specfun_cfc(double complex, double complex *, double complex *);
void specfun_cfs(double complex, double complex *, double complex *);
double complex specfun_cgama(double complex, int);
double specfun_chgm(double, double, double);
double specfun_chgu(double, double, double, int *, int *);
void specfun_clpmn(double complex, int, int, int, double complex *, double complex *);
void specfun_clpn(int, double complex, double complex *, double complex *);
void specfun_clqmn(double complex, int, int, double complex *, double complex *);
void specfun_clqn(int, double complex, double complex *, double complex *);
void specfun_cpbdn(int, double complex, double complex *, double complex *);
double specfun_cva2(int, int, double);
void specfun_cyzo(int, int, int, double complex*, double complex *);
double specfun_eix(double);
double specfun_e1xb(double);
double complex specfun_eixz(double complex);
double complex specfun_e1z(double complex);
void specfun_eulerb(int, double *);
void specfun_fcoef(int, int, double, double, double *);
void specfun_fcszo(int, int, double complex *);
void specfun_ffk(int, double, double *, double *, double *, double *, double *, double *, double *, double *);
double complex specfun_hygfz(double, double, double, double complex, int*);
void specfun_itairy(double, double *, double *, double *, double *);
void specfun_itika(double, double *, double *);
void specfun_itjya(double, double *, double *);
double specfun_itsh0(double);
double specfun_itsl0(double);
double specfun_itth0(double);
void specfun_ittika(double, double *, double *);
void specfun_ittjya(double, double *, double *);
void specfun_jdzo(int, double *, int *, int *, int *);
void specfun_jyzo(int, int, double *, double *, double *, double *);
void specfun_klvna(double, double *, double *, double *, double *, double *, double *, double *, double *);
void specfun_klvnzo(int, int, double *);
void specfun_lamn(int, double, int *, double *, double *);
void specfun_lamv(double, double, double *, double *, double *);
void specfun_lpmn(int, int, double, double *, double *);
double specfun_lpmv(double, int, double);
void specfun_lpn(int, double, double *, double *);
void specfun_lqmn(double, int, int, double *, double *);
void specfun_lqnb(int, double, double *, double *);
void specfun_mtu0(int, int, double, double, double *, double *);
void specfun_mtu12(int, int, int, double, double, double *, double *, double *, double *);
void specfun_pbdv(double, double, double *, double *, double *, double *);
void specfun_pbvv(double, double, double *, double *, double *, double *);
void specfun_pbwa(double, double, double *, double *, double *, double *);
void specfun_rctj(int, double, int *, double *, double *);
void specfun_rcty(int, double, int *, double *, double *);
void specfun_rswfp(int, int, double, double, double, int, double *, double *, double *, double *);
void specfun_rswfo(int, int, double, double, double, int, double *, double *, double *, double *);
void specfun_sdmn(int, int, double, double, int, double *);
void specfun_segv(int, int, double, int, double *, double *);



#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* ifndef */
