# Copyright (c) 2012, Jaydeep P. Bardhan
# Copyright (c) 2012, Matthew G. Knepley
# Copyright (c) 2014, Janani Padmanabhan
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without modification,
# are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
# OF SUCH DAMAGE.

import cython

from . cimport sf_error

from libc.math cimport sqrt, fabs, pow
from libc.stdlib cimport malloc, free
from numpy.math cimport NAN, PI

cdef extern from "lapack_defs.h":
    void c_dstevr(char *jobz, char *range, int *n, double *d, double *e,
                  double *vl, double *vu, int *il, int *iu, double *abstol,
                  int *m, double *w, double *z, int *ldz, int *isuppz,
                  double *work, int *lwork, int *iwork, int *liwork, int *info) nogil


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef inline double* lame_coefficients(double h2, double k2, int n, int p,
                                      void **bufferp, double signm,
                                      double signn) nogil:

    # Ensure that the caller can safely call free(*bufferp) even if an
    # invalid argument is found in the following validation code.
    bufferp[0] = NULL

    if n < 0:
        sf_error.error("ellip_harm", sf_error.ARG, "invalid value for n")
        return NULL

    if p < 1 or p > 2*n + 1:
        sf_error.error("ellip_harm", sf_error.ARG, "invalid value for p")
        return NULL

    if fabs(signm) != 1 or fabs(signn) != 1:
        sf_error.error("ellip_harm", sf_error.ARG, "invalid signm or signn")
        return NULL

    cdef double s2, alpha, beta, gamma, lamba_romain, pp, psi, t1, tol, vl, vu
    cdef int r, tp, j, size, i, info, lwork, liwork, c, iu
    cdef Py_UNICODE t

    r = n/2
    alpha = h2
    beta = k2 - h2
    gamma = alpha - beta

    if p - 1 < r + 1:
        t, tp, size = 'K', p, r + 1
    elif p - 1 < (n - r) + (r + 1):
        t, tp, size = 'L', p - (r + 1), n - r
    elif p - 1 < (n - r) + (n - r) + (r + 1):
        t, tp, size = 'M', p - (n - r) - (r + 1), n - r
    elif p - 1 < 2*n + 1:
        t, tp, size = 'N', p - (n - r) - (n - r) - (r + 1), r

    lwork = 60*size
    liwork = 30*size
    tol = 0.0
    vl = 0
    vu = 0

    cdef void *buffer = malloc((sizeof(double)*(7*size + lwork))
                               + (sizeof(int)*(2*size + liwork)))
    bufferp[0] = buffer
    if not buffer:
        sf_error.error("ellip_harm", sf_error.NO_RESULT, "failed to allocate memory")
        return NULL

    cdef double *g = <double *>buffer
    cdef double *d = g + size
    cdef double *f = d + size
    cdef double *ss =  f + size
    cdef double *w =  ss + size
    cdef double *dd = w + size
    cdef double *eigv = dd + size
    cdef double *work = eigv + size

    cdef int *iwork = <int *>(work + lwork)
    cdef int *isuppz = iwork + liwork

    if t == 'K':
        for j in range(0, r + 1):
           g[j] = (-(2*j + 2)*(2*j + 1)*beta)
           if n%2:
               f[j] = (-alpha*(2*(r- (j + 1)) + 2)*(2*((j + 1) + r) + 1))
               d[j] = ((2*r + 1)*(2*r + 2) - 4*j*j)*alpha + (2*j + 1)*(2*j + 1)*beta
           else:
               f[j] = (-alpha*(2*(r - (j + 1)) + 2)*(2*(r + (j + 1)) - 1))
               d[j] = 2*r*(2*r + 1)*alpha - 4*j*j*gamma
		
    elif t == 'L':
        for j in range(0, n - r):
           g[j] = (-(2*j + 2)*(2*j + 3)*beta)
           if n%2:
               f[j] = (-alpha*(2*(r- (j + 1)) + 2)*(2*((j + 1) + r) + 1))
               d[j] = (2*r + 1)*(2*r + 2)*alpha - (2*j + 1)*(2*j + 1)*gamma
           else:
               f[j] = (-alpha*(2*(r - (j + 1)))*(2*(r+(j + 1)) + 1))
               d[j] = (2*r*(2*r + 1) - (2*j + 1)*(2*j + 1))*alpha + (2*j + 2)*(2*j + 2)*beta
		
    elif t == 'M':
        for j in range(0, n - r):
           g[j] = (-(2*j + 2)*(2*j + 1)*beta)
           if n%2:
               f[j] = (-alpha*(2*(r - (j + 1)) + 2)*(2*((j + 1) + r) + 1))
               d[j] = ((2*r + 1)*(2*r + 2) - (2*j + 1)*(2*j + 1))*alpha + 4*j*j*beta
           else:
               f[j] = (-alpha*(2*(r - (j + 1)))*(2*(r+(j + 1)) + 1))
               d[j] = 2*r*(2*r + 1)*alpha - (2*j + 1)*(2*j + 1)*gamma	

    elif t == 'N':
        for j in range(0, r):
           g[j] = (-(2*j + 2)*(2*j + 3)*beta)
           if n%2:
               f[j] = (-alpha*(2*(r- (j + 1)))*(2*((j + 1) + r) + 3))
               d[j] = (2*r + 1)*(2*r + 2)*alpha - (2*j + 2)*(2*j + 2)*gamma	
           else:
               f[j] = (-alpha*(2*(r - (j + 1)))*(2*(r+(j + 1)) + 1))
               d[j] = 2*r*(2*r + 1)*alpha - (2*j + 2)*(2*j +2)*alpha + (2*j + 1)*(2*j + 1)*beta

    for i in range(0, size):
        if i == 0:
           ss[i] = 1
        else:
           ss[i] = sqrt(g[i - 1]/f[i - 1])*ss[i - 1]

    for i in range(0, size-1):
        dd[i] = g[i]*ss[i]/ss[i+1]

    c_dstevr("V", "I", &size, d, dd, &vl, &vu, &tp, &tp, &tol, &c, w, eigv,
             &size, isuppz, work, &lwork, iwork, &liwork, &info)

    if info != 0:
        sf_error.error("ellip_harm", sf_error.NO_RESULT, "failed to allocate memory")
        return NULL

    for i in range(0, size):
        eigv[i] /= ss[i]

    for i in range(0, size):
        eigv[i] = eigv[i]/(eigv[size - 1]/pow(-h2, size - 1))
    return eigv


@cython.wraparound(False)
@cython.boundscheck(False)
@cython.cdivision(True)
cdef inline double ellip_harm_eval(double h2, double k2, int n, int p,
                                   double s, double *eigv, double signm,
                                   double signn) nogil:
    cdef int size, tp, r, j
    cdef double s2, pp, lambda_romain, psi
    s2 = s*s
    r = n/2
    if p - 1 < r + 1:
        size, psi = r + 1, pow(s, n - 2*r)
    elif p - 1 < (n - r) + (r + 1):
        size, psi = n - r, pow(s, 1 - n + 2*r)*signm*sqrt(fabs(s2 - h2))
    elif p - 1 < (n - r) + (n - r) + (r + 1):
        size, psi = n - r, pow(s, 1 - n + 2*r)*signn*sqrt(fabs(s2 - k2))
    elif p - 1 < 2*n + 1:
        size, psi = r, pow(s,  n - 2*r)*signm*signn*sqrt(fabs((s2 - h2)*(s2 - k2)))
    lambda_romain = 1.0 - <double>s2/<double>h2
    pp = eigv[size - 1]

    for j in range(size - 2, -1, -1):
        pp = pp*lambda_romain + eigv[j]
    pp = pp*psi
    return pp


cdef inline double ellip_harmonic(double h2, double k2, int n, int p, double s,
                                  double signm, double signn) nogil:
    cdef double result
    cdef double *eigv
    cdef void *bufferp
    eigv = lame_coefficients(h2, k2, n, p, &bufferp, signm, signn)
    if not eigv:
        free(bufferp)
        return NAN
    result = ellip_harm_eval(h2, k2, n, p, s, eigv, signm, signn)
    free(bufferp)
    return result
