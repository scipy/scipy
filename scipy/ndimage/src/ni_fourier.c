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
#include <stdlib.h>
#include <math.h>
#include <assert.h>

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

#define CASE_FOURIER_OUT_RR(_po, _tmp, _type) \
case t ## _type:                              \
    *(_type*)_po = _tmp;                        \
    break

#define CASE_FOURIER_OUT_RC(_po, _tmp, _type) \
case t ## _type:                              \
    (*(_type*)_po).real = tmp;                  \
    (*(_type*)_po).imag = 0.0;                  \
    break

#define CASE_FOURIER_OUT_CC(_po, _tmp_r, _tmp_i, _type) \
case t ## _type:                                        \
    (*(_type*)_po).real = _tmp_r;                         \
    (*(_type*)_po).imag = _tmp_i;                         \
    break

#define CASE_FOURIER_FILTER_RC(_pi, _tmp, _tmp_r, _tmp_i, _type) \
case t ## _type:                                                 \
    _tmp_r = (*(_type*)_pi).real * _tmp;                           \
    _tmp_i = (*(_type*)_pi).imag * _tmp;                           \
    break;

#define CASE_FOURIER_FILTER_RR(_pi, _tmp, _type) \
case t ## _type:                                 \
    _tmp *= *(_type*)_pi;                          \
    break;

int NI_FourierFilter(PyArrayObject *input, PyArrayObject* parameter_array,
                     npy_intp n, int axis, PyArrayObject* output,
                     int filter_type)
{
    NI_Iterator ii, io;
    char *pi, *po;
    double *parameters = NULL, **params = NULL;
    npy_intp kk, hh, size;
    Float64 *iparameters = (void *)PyArray_DATA(parameter_array);
    int ll;

    /* precalculate the parameters: */
    parameters = (double*)malloc(input->nd * sizeof(double));
    if (!parameters) {
        PyErr_NoMemory();
        goto exit;
    }
    for(kk = 0; kk < input->nd; kk++) {
        /* along the direction of the real transform we must use the given
             length of that dimensons, unless a complex transform is assumed
             (n < 0): */
        int shape = kk == axis ?
                        (n < 0 ? input->dimensions[kk] : n) : input->dimensions[kk];
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
    params = (double**) malloc(input->nd * sizeof(double*));
    if (!params) {
        PyErr_NoMemory();
        goto exit;
    }
    for(kk = 0; kk < input->nd; kk++)
        params[kk] = NULL;
    for(kk = 0; kk < input->nd; kk++) {
        if (input->dimensions[kk] > 1) {
            params[kk] = (double*)malloc(input->dimensions[kk] * sizeof(double));
            if (!params[kk]) {
                PyErr_NoMemory();
                goto exit;
            }
        }
    }
    switch (filter_type) {
        case _NI_GAUSSIAN:
            /* calculate the tables of exponentials: */
            for (hh = 0; hh < input->nd; hh++) {
                if (params[hh]) {
                    if (hh == axis && n >= 0) {
                        for(kk = 0; kk < input->dimensions[hh]; kk++) {
                            double tmp = parameters[hh] * kk * kk;
                            params[hh][kk] = fabs(tmp) > 50.0 ? 0.0 : exp(tmp);
                        }
                    } else {
                        int jj = 0;
                        for(kk = 0; kk < (input->dimensions[hh] + 1) / 2; kk++) {
                            double tmp = parameters[hh] * kk * kk;
                            params[hh][jj++] = fabs(tmp) > 50.0 ? 0.0 : exp(tmp);
                        }
                        for(kk = -(input->dimensions[hh] / 2); kk < 0; kk++) {
                            double tmp = parameters[hh] * kk * kk;
                            params[hh][jj++] = fabs(tmp) > 50.0 ? 0.0 : exp(tmp);
                        }
                    }
                }
            }
            break;
        case _NI_UNIFORM:
            /* calculate the tables of parameters: */
            for (hh = 0; hh < input->nd; hh++) {
                if (params[hh]) {
                    params[hh][0] = 1.0;
                    if (hh == axis && n >= 0) {
                        double tmp = M_PI * parameters[hh] / n;
                        for(kk = 1; kk < input->dimensions[hh]; kk++)
                            params[hh][kk] = tmp > 0.0 ?
                                                                            sin(tmp * kk) / (tmp * kk) : 0.0;
                    } else {
                        double tmp = M_PI * parameters[hh] / input->dimensions[hh];
                        int jj = 1;
                        for(kk = 1; kk < (input->dimensions[hh] + 1) / 2; kk++)
                            params[hh][jj++] = tmp > 0.0 ?
                                                                            sin(tmp * kk) / (tmp * kk) : 0.0;
                        for(kk = -(input->dimensions[hh] / 2); kk < 0; kk++)
                            params[hh][jj++] = tmp > 0.0 ?
                                                                            sin(tmp * kk) / (tmp * kk) : 0.0;
                    }
                }
            }
            break;
        case _NI_ELLIPSOID:
            /* calculate the tables of parameters: */
            for (hh = 0; hh < input->nd; hh++) {
                if (params[hh]) {
                    params[hh][0] = 1.0;
                    if (hh == axis && n >= 0) {
                        double tmp = M_PI * parameters[hh] / n;
                        for(kk = 0; kk < input->dimensions[hh]; kk++)
                            params[hh][kk] = (double)kk * tmp;
                    } else {
                        double tmp = M_PI * parameters[hh] / input->dimensions[hh];
                        int jj = 0;
                        for(kk = 0; kk < (input->dimensions[hh] + 1) / 2; kk++)
                            params[hh][jj++] = (double)kk * tmp;
                        for(kk = -(input->dimensions[hh] / 2); kk < 0; kk++)
                            params[hh][jj++] = (double)kk * tmp;
                    }
                } else if (input->dimensions[hh] > 0) {
                    params[hh][0] = 1.0;
                }
            }
            if (input->nd > 1)
                for(hh = 0; hh < input->nd; hh++)
                    for(kk = 0; kk < input->dimensions[hh]; kk++)
                        params[hh][kk] = params[hh][kk] * params[hh][kk];
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
    size = 1;
    for(ll = 0; ll < input->nd; ll++)
        size *= input->dimensions[ll];
    /* iterator over the elements: */
    for(hh = 0; hh < size; hh++) {
        double tmp = 1.0;
        switch (filter_type) {
        case _NI_GAUSSIAN:
        case _NI_UNIFORM:
            for(kk = 0; kk < input->nd; kk++)
                if (params[kk])
                    tmp *= params[kk][ii.coordinates[kk]];
            break;
        case _NI_ELLIPSOID:
            switch (input->nd) {
            case 1:
                tmp = params[0][ii.coordinates[0]];
                tmp = tmp > 0.0 ? sin(tmp) / (tmp) : 1.0;
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
        if (input->descr->type_num == tComplex64 ||
                input->descr->type_num == tComplex128) {
            double tmp_r = 0.0, tmp_i = 0.0;
            switch (input->descr->type_num) {
                CASE_FOURIER_FILTER_RC(pi, tmp, tmp_r, tmp_i, Complex64);
                CASE_FOURIER_FILTER_RC(pi, tmp, tmp_r, tmp_i, Complex128);
            default:
                PyErr_SetString(PyExc_RuntimeError, "data type not supported");
                goto exit;
            }
            switch (output->descr->type_num) {
                CASE_FOURIER_OUT_CC(po, tmp_r, tmp_i, Complex64);
                CASE_FOURIER_OUT_CC(po, tmp_r, tmp_i, Complex128);
            default:
                PyErr_SetString(PyExc_RuntimeError, "data type not supported");
                goto exit;
            }
        } else {
            switch (input->descr->type_num) {
                CASE_FOURIER_FILTER_RR(pi, tmp, Bool)
                CASE_FOURIER_FILTER_RR(pi, tmp, UInt8)
                CASE_FOURIER_FILTER_RR(pi, tmp, UInt16)
                CASE_FOURIER_FILTER_RR(pi, tmp, UInt32)
#if HAS_UINT64
                CASE_FOURIER_FILTER_RR(pi, tmp, UInt64)
#endif
                CASE_FOURIER_FILTER_RR(pi, tmp, Int8)
                CASE_FOURIER_FILTER_RR(pi, tmp, Int16)
                CASE_FOURIER_FILTER_RR(pi, tmp, Int32)
                CASE_FOURIER_FILTER_RR(pi, tmp, Int64)
                CASE_FOURIER_FILTER_RR(pi, tmp, Float32)
                CASE_FOURIER_FILTER_RR(pi, tmp, Float64)
            default:
                PyErr_SetString(PyExc_RuntimeError, "data type not supported");
                goto exit;
            }
            switch (output->descr->type_num) {
                CASE_FOURIER_OUT_RR(po, tmp, Float32);
                CASE_FOURIER_OUT_RR(po, tmp, Float64);
                CASE_FOURIER_OUT_RC(po, tmp, Complex64);
                CASE_FOURIER_OUT_RC(po, tmp, Complex128);
            default:
                PyErr_SetString(PyExc_RuntimeError, "data type not supported");
                goto exit;
            }
        }
        NI_ITERATOR_NEXT2(ii, io, pi, po);
    }

 exit:
    if (parameters) free(parameters);
    if (params) {
        for(kk = 0; kk < input->nd; kk++)
            if (params[kk]) free(params[kk]);
        free(params);
    }
    return PyErr_Occurred() ? 0 : 1;
}

#define CASE_FOURIER_SHIFT_R(_pi, _tmp, _r, _i, _cost, _sint, _type) \
case t ## _type:                                                     \
    _tmp = *(_type*)_pi;                                               \
    _r = _tmp * _cost;                                                 \
    _i = _tmp * _sint;                                                 \
    break;

#define CASE_FOURIER_SHIFT_C(_pi, _r, _i, _cost, _sint, _type)     \
case t ## _type:                                                   \
    _r = (*(_type*)_pi).real * _cost - (*(_type*)_pi).imag * _sint;  \
    _i = (*(_type*)_pi).real * _sint + (*(_type*)_pi).imag * _cost;  \
    break;

int NI_FourierShift(PyArrayObject *input, PyArrayObject* shift_array,
            npy_intp n, int axis, PyArrayObject* output)
{
    NI_Iterator ii, io;
    char *pi, *po;
    double *shifts = NULL, **params = NULL;
    npy_intp kk, hh, size;
    Float64 *ishifts = (void *)PyArray_DATA(shift_array);
    int ll;

    /* precalculate the shifts: */
    shifts = (double*)malloc(input->nd * sizeof(double));
    if (!shifts) {
        PyErr_NoMemory();
        goto exit;
    }
    for(kk = 0; kk < input->nd; kk++) {
        /* along the direction of the real transform we must use the given
             length of that dimensons, unless a complex transform is assumed
             (n < 0): */
        int shape = kk == axis ?
                        (n < 0 ? input->dimensions[kk] : n) : input->dimensions[kk];
        shifts[kk] = -2.0 * M_PI * *ishifts++ / (double)shape;
    }
    /* allocate memory for tables: */
    params = (double**) malloc(input->nd * sizeof(double*));
    if (!params) {
        PyErr_NoMemory();
        goto exit;
    }
    for(kk = 0; kk < input->nd; kk++)
        params[kk] = NULL;
    for(kk = 0; kk < input->nd; kk++) {
        if (input->dimensions[kk] > 1) {
            params[kk] = (double*)malloc(input->dimensions[kk] * sizeof(double));
            if (!params[kk]) {
                PyErr_NoMemory();
                goto exit;
            }
        }
    }
    for (hh = 0; hh < input->nd; hh++) {
        if (params[hh]) {
            if (hh == axis && n >= 0) {
                for(kk = 0; kk < input->dimensions[hh]; kk++)
                    params[hh][kk] = shifts[hh] * kk;
            } else {
                int jj = 0;
                for(kk = 0; kk < (input->dimensions[hh] + 1) / 2; kk++) {
                    params[hh][jj++] = shifts[hh] * kk;
                }
                for(kk = -(input->dimensions[hh] / 2); kk < 0; kk++) {
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
    size = 1;
    for(ll = 0; ll < input->nd; ll++)
        size *= input->dimensions[ll];
    /* iterator over the elements: */
    for(hh = 0; hh < size; hh++) {
        double tmp = 0.0, sint, cost, r = 0.0, i = 0.0;
        for(kk = 0; kk < input->nd; kk++)
            if (params[kk])
                tmp += params[kk][ii.coordinates[kk]];
        sint = sin(tmp);
        cost = cos(tmp);
        switch (input->descr->type_num) {
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, Bool)
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, UInt8)
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, UInt16)
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, UInt32)
#if HAS_UINT64
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, UInt64)
#endif
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, Int8)
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, Int16)
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, Int32)
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, Int64)
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, Float32)
            CASE_FOURIER_SHIFT_R(pi, tmp, r, i, cost, sint, Float64)
            CASE_FOURIER_SHIFT_C(pi, r, i, cost, sint, Complex64)
            CASE_FOURIER_SHIFT_C(pi, r, i, cost, sint, Complex128)
        default:
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        switch (output->descr->type_num) {
            CASE_FOURIER_OUT_CC(po, r, i, Complex64);
            CASE_FOURIER_OUT_CC(po, r, i, Complex128);
        default:
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        NI_ITERATOR_NEXT2(ii, io, pi, po);
    }

 exit:
    if (shifts) free(shifts);
    if (params) {
        for(kk = 0; kk < input->nd; kk++)
            if (params[kk]) free(params[kk]);
        free(params);
    }
    return PyErr_Occurred() ? 0 : 1;
}
