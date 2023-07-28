/* Copyright (C) 2003-2005 Peter J. Verveer
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 *    products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "ni_support.h"
#include "ni_fourier.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include <numpy/npy_math.h>
#include "npy_2_complexcompat.h"

#if !defined(M_PI)
#define M_PI 3.14159265358979323846
#endif

#define _NI_GAUSSIAN 0
#define _NI_UNIFORM 1
#define _NI_ELLIPSOID 2

static double polevl(double x, const double coef[], int N)
{
    double ans;
    const double *p = coef;
    int i = N;

    ans = *p++;
    do
        ans = ans * x + *p++;
    while(--i);

    return ans ;
}

double p1evl(double x, const double coef[], int N)
{
    double ans;
    const double *p = coef;
    int i = N - 1;

    ans = x + *p++;
    do
        ans = ans * x + *p++;
    while(--i);

    return ans;
}

#define THPIO4 2.35619449019234492885
#define SQ2OPI .79788456080286535588
#define Z1 1.46819706421238932572E1
#define Z2 4.92184563216946036703E1

static double _bessel_j1(double x)
{
    double w, z, p, q, xn;
    const double RP[4] = {
        -8.99971225705559398224E8,
        4.52228297998194034323E11,
        -7.27494245221818276015E13,
        3.68295732863852883286E15,
    };
    const double RQ[8] = {
        6.20836478118054335476E2,
        2.56987256757748830383E5,
        8.35146791431949253037E7,
        2.21511595479792499675E10,
        4.74914122079991414898E12,
        7.84369607876235854894E14,
        8.95222336184627338078E16,
        5.32278620332680085395E18,
    };
    const double PP[7] = {
        7.62125616208173112003E-4,
        7.31397056940917570436E-2,
        1.12719608129684925192E0,
        5.11207951146807644818E0,
        8.42404590141772420927E0,
        5.21451598682361504063E0,
        1.00000000000000000254E0,
    };
    const double PQ[7] = {
        5.71323128072548699714E-4,
        6.88455908754495404082E-2,
        1.10514232634061696926E0,
        5.07386386128601488557E0,
        8.39985554327604159757E0,
        5.20982848682361821619E0,
        9.99999999999999997461E-1,
    };
    const double QP[8] = {
        5.10862594750176621635E-2,
        4.98213872951233449420E0,
        7.58238284132545283818E1,
        3.66779609360150777800E2,
        7.10856304998926107277E2,
        5.97489612400613639965E2,
        2.11688757100572135698E2,
        2.52070205858023719784E1,
    };
    const double QQ[7] = {
        7.42373277035675149943E1,
        1.05644886038262816351E3,
        4.98641058337653607651E3,
        9.56231892404756170795E3,
        7.99704160447350683650E3,
        2.82619278517639096600E3,
        3.36093607810698293419E2,
    };

    w = x;
    if (x < 0)
        w = -x;

    if (w <= 5.0) {
        z = x * x;
        w = polevl(z, RP, 3) / p1evl(z, RQ, 8);
        w = w * x * (z - Z1) * (z - Z2);
        return w ;
    }

    w = 5.0 / x;
    z = w * w;
    p = polevl(z, PP, 6) / polevl(z, PQ, 6);
    q = polevl(z, QP, 7) / p1evl(z, QQ, 7);
    xn = x - THPIO4;
    p = p * cos(xn) - w * q * sin(xn);
    return p * SQ2OPI / sqrt(x);
}

#define CASE_FOURIER_OUT_RR(_TYPE, _type, _po, _tmp) \
case _TYPE:                                          \
    *(_type *)_po = _tmp;                            \
    break

#define CASE_FOURIER_OUT_RC(_TYPE, _type, _T, _po, _tmp) \
case _TYPE:                                          \
    NPY_CSETREAL##_T((_type *)_po, tmp);                      \
    NPY_CSETIMAG##_T((_type *)_po, 0.0);                      \
    break

#define CASE_FOURIER_OUT_CC(_TYPE, _type, _T, _po, _tmp_r, _tmp_i) \
case _TYPE:                                                    \
    NPY_CSETREAL##_T((_type *)_po, _tmp_r);                             \
    NPY_CSETIMAG##_T((_type *)_po, _tmp_i);                             \
    break

#define CASE_FOURIER_FILTER_RC(_TYPE, _type, _t, _pi, _tmp, _tmp_r, _tmp_i) \
case _TYPE:                                                             \
    _tmp_r = npy_creal##_t(*((_type *)_pi)) * _tmp;                               \
    _tmp_i = npy_cimag##_t(*((_type *)_pi)) * _tmp;                               \
    break

#define CASE_FOURIER_FILTER_RR(_TYPE, _type, _pi, _tmp) \
case _TYPE:                                             \
    _tmp *= *(_type *)_pi;                              \
    break

int NI_FourierFilter(PyArrayObject *input, PyArrayObject* parameter_array,
                     npy_intp n, int axis, PyArrayObject* output,
                     int filter_type)
{
    NI_Iterator ii, io;
    char *pi, *po;
    double *parameters = NULL, **params = NULL;
    npy_intp kk, hh, size;
    npy_double *iparameters = (void *)PyArray_DATA(parameter_array);
    NPY_BEGIN_THREADS_DEF;

    /* precalculate the parameters: */
    parameters = malloc(PyArray_NDIM(input) * sizeof(double));
    if (!parameters) {
        PyErr_NoMemory();
        goto exit;
    }
    for (kk = 0; kk < PyArray_NDIM(input); kk++) {
        /* along the direction of the real transform we must use the given
             length of that dimensons, unless a complex transform is assumed
             (n < 0): */
        int shape = kk == axis ?
                (n < 0 ? PyArray_DIM(input, kk) : n) : PyArray_DIM(input, kk);
        switch (filter_type) {
            case _NI_GAUSSIAN:
                parameters[kk] = *iparameters++ * M_PI / (double)shape;
                parameters[kk] = -2.0 * parameters[kk] * parameters[kk];
                break;
            case _NI_ELLIPSOID:
            case _NI_UNIFORM:
                parameters[kk] = *iparameters++;
                break;
        }
    }
    /* allocate memory for tables: */
    params = malloc(PyArray_NDIM(input) * sizeof(double*));
    if (!params) {
        PyErr_NoMemory();
        goto exit;
    }
    for (kk = 0; kk < PyArray_NDIM(input); kk++) {
        params[kk] = NULL;
    }
    for (kk = 0; kk < PyArray_NDIM(input); kk++) {
        if (PyArray_DIM(input, kk) > 1 || filter_type == _NI_ELLIPSOID) {
            params[kk] = malloc(PyArray_DIM(input, kk) * sizeof(double));
            if (!params[kk]) {
                PyErr_NoMemory();
                goto exit;
            }
        }
    }

    NPY_BEGIN_THREADS;

    switch (filter_type) {
        case _NI_GAUSSIAN:
            /* calculate the tables of exponentials: */
            for (hh = 0; hh < PyArray_NDIM(input); hh++) {
                if (params[hh]) {
                    if (hh == axis && n >= 0) {
                        for (kk = 0; kk < PyArray_DIM(input, hh); kk++) {
                            double tmp = parameters[hh] * kk * kk;
                            params[hh][kk] = fabs(tmp) > 50.0 ? 0.0 : exp(tmp);
                        }
                    }
                    else {
                        const npy_intp dim = PyArray_DIM(input, hh);
                        int jj = 0;
                        for (kk = 0; kk < (dim + 1) / 2; kk++) {
                            double tmp = parameters[hh] * kk * kk;
                            params[hh][jj++] =
                                            fabs(tmp) > 50.0 ? 0.0 : exp(tmp);
                        }
                        for (kk = -(dim / 2); kk < 0; kk++) {
                            double tmp = parameters[hh] * kk * kk;
                            params[hh][jj++] =
                                            fabs(tmp) > 50.0 ? 0.0 : exp(tmp);
                        }
                    }
                }
            }
            break;
        case _NI_UNIFORM:
            /* calculate the tables of parameters: */
            for (hh = 0; hh < PyArray_NDIM(input); hh++) {
                if (params[hh]) {
                    params[hh][0] = 1.0;
                    if (hh == axis && n >= 0) {
                        double tmp = M_PI * parameters[hh] / n;
                        for (kk = 1; kk < PyArray_DIM(input, hh); kk++) {
                            params[hh][kk] = tmp > 0.0 ?
                                             sin(tmp * kk) / (tmp * kk) : 0.0;
                        }
                    }
                    else {
                        const npy_intp dim = PyArray_DIM(input, hh);
                        double tmp = M_PI * parameters[hh] / dim;
                        int jj = 1;
                        for (kk = 1; kk < (dim + 1) / 2; kk++) {
                            params[hh][jj++] = tmp > 0.0 ?
                                            sin(tmp * kk) / (tmp * kk) : 0.0;
                        }
                        for (kk = -(dim / 2); kk < 0; kk++) {
                            params[hh][jj++] = tmp > 0.0 ?
                                            sin(tmp * kk) / (tmp * kk) : 0.0;
                        }
                    }
                }
            }
            break;
        case _NI_ELLIPSOID:
            /* calculate the tables of parameters: */
            for (hh = 0; hh < PyArray_NDIM(input); hh++) {
                if (params[hh]) {
                    params[hh][0] = 1.0;
                    if (hh == axis && n >= 0) {
                        double tmp = M_PI * parameters[hh] / n;
                        for (kk = 0; kk < PyArray_DIM(input, hh); kk++) {
                            params[hh][kk] = (double)kk * tmp;
                        }
                    }
                    else {
                        const npy_intp dim = PyArray_DIM(input, hh);
                        double tmp = M_PI * parameters[hh] / dim;
                        int jj = 0;
                        for(kk = 0; kk < (dim + 1) / 2; kk++) {
                            params[hh][jj++] = (double)kk * tmp;
                        }
                        for(kk = -(dim / 2); kk < 0; kk++) {
                            params[hh][jj++] = (double)kk * tmp;
                        }
                    }
                }
                else if (PyArray_DIM(input, hh) > 0) {
                    params[hh][0] = 1.0;
                }
            }
            if (PyArray_NDIM(input) > 1)
                for(hh = 0; hh < PyArray_NDIM(input); hh++)
                    for(kk = 0; kk < PyArray_DIM(input, hh); kk++) {
                        params[hh][kk] *= params[hh][kk];
                    }
            break;
        default:
            break;
    }
    /* initialize input element iterator: */
    if (!NI_InitPointIterator(input, &ii))
        goto exit;
    /* initialize output element iterator: */
    if (!NI_InitPointIterator(output, &io))
        goto exit;
    pi = (void *)PyArray_DATA(input);
    po = (void *)PyArray_DATA(output);
    size = PyArray_SIZE(input);
    /* iterator over the elements: */
    for(hh = 0; hh < size; hh++) {
        double tmp = 1.0;
        switch (filter_type) {
        case _NI_GAUSSIAN:
        case _NI_UNIFORM:
            for (kk = 0; kk < PyArray_NDIM(input); kk++) {
                if (params[kk])
                    tmp *= params[kk][ii.coordinates[kk]];
            }
            break;
        case _NI_ELLIPSOID:
            switch (PyArray_NDIM(input)) {
            case 1:
                tmp = params[0][ii.coordinates[0]];
                tmp = tmp != 0 ? sin(tmp) / (tmp) : 1.0;
                break;
            case 2:
                tmp = 0.0;
                for(kk = 0; kk < 2; kk++)
                    tmp += params[kk][ii.coordinates[kk]];
                tmp = sqrt(tmp);
                tmp = tmp > 0.0 ? 2.0 * _bessel_j1(tmp) / tmp : 1.0;
                break;
            case 3:
                {
                    double r = 0.0;
                    for(kk = 0; kk < 3; kk++)
                        r += params[kk][ii.coordinates[kk]];
                    r = sqrt(r);
                    if (r > 0.0) {
                        tmp = 3.0 * (sin(r) - r * cos(r));
                        tmp /= r * r * r;
                    } else {
                        tmp = 1.0;
                    }
                }
                break;
            }
            break;
        default:
            break;
        }
        if (PyArray_TYPE(input) == NPY_CFLOAT ||
                PyArray_TYPE(input) == NPY_CDOUBLE) {
            double tmp_r = 0.0, tmp_i = 0.0;
            switch (PyArray_TYPE(input)) {
                CASE_FOURIER_FILTER_RC(NPY_CFLOAT, npy_cfloat, f,
                                       pi, tmp, tmp_r, tmp_i);
                CASE_FOURIER_FILTER_RC(NPY_CDOUBLE, npy_cdouble,,
                                       pi, tmp, tmp_r, tmp_i);
            default:
                NPY_END_THREADS;
                PyErr_SetString(PyExc_RuntimeError, "data type not supported");
                goto exit;
            }
            switch (PyArray_TYPE(output)) {
                CASE_FOURIER_OUT_CC(NPY_CFLOAT, npy_cfloat, F, po, tmp_r, tmp_i);
                CASE_FOURIER_OUT_CC(NPY_CDOUBLE, npy_cdouble,, po, tmp_r, tmp_i);
            default:
                NPY_END_THREADS;
                PyErr_SetString(PyExc_RuntimeError, "data type not supported");
                goto exit;
            }
        }
        else {
            switch (PyArray_TYPE(input)) {
                CASE_FOURIER_FILTER_RR(NPY_BOOL, npy_bool, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_UBYTE, npy_ubyte, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_USHORT, npy_ushort, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_UINT, npy_uint, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_ULONG, npy_ulong, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_ULONGLONG, npy_ulonglong, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_BYTE, npy_byte, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_SHORT, npy_short, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_INT, npy_int, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_LONG, npy_long, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_LONGLONG, npy_longlong, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_FLOAT, npy_float, pi, tmp);
                CASE_FOURIER_FILTER_RR(NPY_DOUBLE, npy_double, pi, tmp);
            default:
                NPY_END_THREADS;
                PyErr_SetString(PyExc_RuntimeError, "data type not supported");
                goto exit;
            }
            switch (PyArray_TYPE(output)) {
                CASE_FOURIER_OUT_RR(NPY_FLOAT, npy_float, po, tmp);
                CASE_FOURIER_OUT_RR(NPY_DOUBLE, npy_double, po, tmp);
                CASE_FOURIER_OUT_RC(NPY_CFLOAT, npy_cfloat, F, po, tmp);
                CASE_FOURIER_OUT_RC(NPY_CDOUBLE, npy_cdouble,, po, tmp);
            default:
                NPY_END_THREADS;
                PyErr_SetString(PyExc_RuntimeError, "data type not supported");
                goto exit;
            }
        }
        NI_ITERATOR_NEXT2(ii, io, pi, po);
    }

 exit:
    NPY_END_THREADS;
    free(parameters);
    if (params) {
        for (kk = 0; kk < PyArray_NDIM(input); kk++) {
            free(params[kk]);
        }
        free(params);
    }
    return PyErr_Occurred() ? 0 : 1;
}

#define CASE_FOURIER_SHIFT_R(_TYPE, _type, _pi, _tmp, _r, _i, _cost, _sint) \
case _TYPE:                                                                 \
    _tmp = *(_type *)_pi;                                                   \
    _r = _tmp * _cost;                                                      \
    _i = _tmp * _sint;                                                      \
    break

#define CASE_FOURIER_SHIFT_C(_TYPE, _type, _t, _pi, _r, _i, _cost, _sint) \
case _TYPE:                                                           \
    _r = npy_creal##_t(*((_type *)_pi)) * _cost - npy_cimag##_t(*((_type *)_pi)) * _sint; \
    _i = npy_creal##_t(*((_type *)_pi)) * _sint + npy_cimag##_t(*((_type *)_pi)) * _cost; \
    break

int NI_FourierShift(PyArrayObject *input, PyArrayObject* shift_array,
            npy_intp n, int axis, PyArrayObject* output)
{
    NI_Iterator ii, io;
    char *pi, *po;
    double *shifts = NULL, **params = NULL;
    npy_intp kk, hh, size;
    npy_double *ishifts = (void *)PyArray_DATA(shift_array);
    NPY_BEGIN_THREADS_DEF;

    /* precalculate the shifts: */
    shifts = malloc(PyArray_NDIM(input) * sizeof(double));
    if (!shifts) {
        PyErr_NoMemory();
        goto exit;
    }
    for (kk = 0; kk < PyArray_NDIM(input); kk++) {
        /* along the direction of the real transform we must use the given
             length of that dimensons, unless a complex transform is assumed
             (n < 0): */
        int shape = kk == axis ?
                (n < 0 ? PyArray_DIM(input, kk) : n) : PyArray_DIM(input, kk);
        shifts[kk] = -2.0 * M_PI * *ishifts++ / (double)shape;
    }
    /* allocate memory for tables: */
    params = malloc(PyArray_NDIM(input) * sizeof(double*));
    if (!params) {
        PyErr_NoMemory();
        goto exit;
    }
    for (kk = 0; kk < PyArray_NDIM(input); kk++) {
        params[kk] = NULL;
    }
    for (kk = 0; kk < PyArray_NDIM(input); kk++) {
        if (PyArray_DIM(input, kk) > 1) {
            params[kk] = malloc(PyArray_DIM(input, kk) * sizeof(double));
            if (!params[kk]) {
                PyErr_NoMemory();
                goto exit;
            }
        }
    }

    NPY_BEGIN_THREADS;

    for (hh = 0; hh < PyArray_NDIM(input); hh++) {
        if (params[hh]) {
            if (hh == axis && n >= 0) {
                for (kk = 0; kk < PyArray_DIM(input, hh); kk++) {
                    params[hh][kk] = shifts[hh] * kk;
                }
            }
            else {
                int jj = 0;
                for (kk = 0; kk < (PyArray_DIM(input, hh) + 1) / 2; kk++) {
                    params[hh][jj++] = shifts[hh] * kk;
                }
                for (kk = -(PyArray_DIM(input, hh) / 2); kk < 0; kk++) {
                    params[hh][jj++] = shifts[hh] * kk;
                }
            }
        }
    }
    /* initialize input element iterator: */
    if (!NI_InitPointIterator(input, &ii))
        goto exit;
    /* initialize output element iterator: */
    if (!NI_InitPointIterator(output, &io))
        goto exit;
    pi = (void *)PyArray_DATA(input);
    po = (void *)PyArray_DATA(output);
    size = PyArray_SIZE(input);
    /* iterator over the elements: */
    for(hh = 0; hh < size; hh++) {
        double tmp = 0.0, sint, cost, r = 0.0, i = 0.0;
        for (kk = 0; kk < PyArray_NDIM(input); kk++) {
            if (params[kk])
                tmp += params[kk][ii.coordinates[kk]];
        }
        sint = sin(tmp);
        cost = cos(tmp);
        switch (PyArray_TYPE(input)) {
            CASE_FOURIER_SHIFT_R(NPY_BOOL, npy_bool,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_UBYTE, npy_ubyte,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_USHORT, npy_ushort,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_UINT, npy_uint,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_ULONG, npy_ulong,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_ULONGLONG, npy_ulonglong,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_BYTE, npy_byte,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_SHORT, npy_short,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_INT, npy_int,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_LONG, npy_long,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_LONGLONG, npy_longlong,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_FLOAT, npy_float,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_R(NPY_DOUBLE, npy_double,
                                 pi, tmp, r, i, cost, sint);
            CASE_FOURIER_SHIFT_C(NPY_CFLOAT, npy_cfloat, f,
                                 pi, r, i, cost, sint);
            CASE_FOURIER_SHIFT_C(NPY_CDOUBLE, npy_cdouble,,
                                 pi, r, i, cost, sint);
        default:
            NPY_END_THREADS;
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        switch (PyArray_TYPE(output)) {
            CASE_FOURIER_OUT_CC(NPY_CFLOAT, npy_cfloat, F, po, r, i);
            CASE_FOURIER_OUT_CC(NPY_CDOUBLE, npy_cdouble,, po, r, i);
        default:
            NPY_END_THREADS;
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        NI_ITERATOR_NEXT2(ii, io, pi, po);
    }

 exit:
    NPY_END_THREADS;
    free(shifts);
    if (params) {
        for (kk = 0; kk < PyArray_NDIM(input); kk++) {
            free(params[kk]);
        }
        free(params);
    }
    return PyErr_Occurred() ? 0 : 1;
}
