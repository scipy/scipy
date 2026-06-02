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
#include "ni_filters.h"
#include <math.h>

#define BUFFER_SIZE 256000

int NI_Correlate1D(PyArrayObject *input, PyArrayObject *weights,
                   int axis, PyArrayObject *output, NI_ExtendMode mode,
                   double cval, npy_intp origin)
{
    int symmetric = 0, more;
    npy_intp ii, jj, ll, lines, length, size1, size2, filter_size;
    double *ibuffer = NULL, *obuffer = NULL;
    npy_double *fw;
    NI_LineBuffer iline_buffer, oline_buffer;
    NPY_BEGIN_THREADS_DEF;

    /* test for symmetry or anti-symmetry: */
    filter_size = PyArray_SIZE(weights);
    size1 = filter_size / 2;
    size2 = filter_size - size1 - 1;
    fw = (void *)PyArray_DATA(weights);
    if (filter_size & 0x1) {
        symmetric = 1;
        for(ii = 1; ii <= filter_size / 2; ii++) {
            if (fabs(fw[ii + size1] - fw[size1 - ii]) > DBL_EPSILON) {
                symmetric = 0;
                break;
            }
        }
        if (symmetric == 0) {
            symmetric = -1;
            for(ii = 1; ii <= filter_size / 2; ii++) {
                if (fabs(fw[size1 + ii] + fw[size1 - ii]) > DBL_EPSILON) {
                    symmetric = 0;
                    break;
                }
            }
        }
    }
    /* allocate and initialize the line buffers: */
    lines = -1;
    if (!NI_AllocateLineBuffer(input, axis, size1 + origin, size2 - origin,
                                                         &lines, BUFFER_SIZE, &ibuffer))
        goto exit;
    if (!NI_AllocateLineBuffer(output, axis, 0, 0, &lines, BUFFER_SIZE,
                                                         &obuffer))
        goto exit;
    if (!NI_InitLineBuffer(input, axis, size1 + origin, size2 - origin,
                                                            lines, ibuffer, mode, cval, &iline_buffer))
        goto exit;
    if (!NI_InitLineBuffer(output, axis, 0, 0, lines, obuffer, mode, 0.0,
                                                 &oline_buffer))
        goto exit;

    NPY_BEGIN_THREADS;
    length = PyArray_NDIM(input) > 0 ? PyArray_DIM(input, axis) : 1;
    fw += size1;
    /* iterate over all the array lines: */
    do {
        /* copy lines from array to buffer: */
        if (!NI_ArrayToLineBuffer(&iline_buffer, &lines, &more)) {
            goto exit;
        }
        /* iterate over the lines in the buffers: */
        for(ii = 0; ii < lines; ii++) {
            /* get lines: */
            double *iline = NI_GET_LINE(iline_buffer, ii) + size1;
            double *oline = NI_GET_LINE(oline_buffer, ii);
            /* the correlation calculation: */
            if (symmetric > 0) {
                for(ll = 0; ll < length; ll++) {
                    oline[ll] = iline[0] * fw[0];
                    for(jj = -size1 ; jj < 0; jj++)
                        oline[ll] += (iline[jj] + iline[-jj]) * fw[jj];
                    ++iline;
                }
            } else if (symmetric < 0) {
                for(ll = 0; ll < length; ll++) {
                    oline[ll] = iline[0] * fw[0];
                    for(jj = -size1 ; jj < 0; jj++)
                        oline[ll] += (iline[jj] - iline[-jj]) * fw[jj];
                    ++iline;
                }
            } else {
                for(ll = 0; ll < length; ll++) {
                    oline[ll] = iline[size2] * fw[size2];
                    for(jj = -size1; jj < size2; jj++)
                        oline[ll] += iline[jj] * fw[jj];
                    ++iline;
                }
            }
        }
        /* copy lines from buffer to array: */
        if (!NI_LineBufferToArray(&oline_buffer)) {
            goto exit;
        }
    } while(more);
exit:
    NPY_END_THREADS;
    free(ibuffer);
    free(obuffer);
    return PyErr_Occurred() ? 0 : 1;
}

#define CASE_CORRELATE_POINT(_TYPE, _type, _pi, _weights, _offsets,        \
                             _filter_size, _cvalue, _res, _mv)             \
case _TYPE:                                                                \
{                                                                          \
    npy_intp _ii, _offset;                                                 \
    for (_ii = 0; _ii < _filter_size; ++_ii) {                             \
        _offset = _offsets[_ii];                                           \
        if (_offset == _mv) {                                              \
            _res += _weights[_ii] * _cvalue;                               \
        }                                                                  \
        else {                                                             \
            _res += _weights[_ii] * (double)(*((_type *)(_pi + _offset))); \
        }                                                                  \
    }                                                                      \
}                                                                          \
break

#define CASE_FILTER_OUT(_TYPE, _type, _po, _tmp) \
case _TYPE:                                      \
    *(_type *)_po = (_type)_tmp;                 \
    break

/* Avoid undefined behaviour of float -> unsigned conversions. */
#define CASE_FILTER_OUT_SAFE(_TYPE, _type, _po, _tmp)                  \
case _TYPE:                                                            \
    *(_type *)_po = (_tmp) > -1. ? (_type)(_tmp) : -(_type)(-_tmp);    \
    break

int NI_Correlate(PyArrayObject* input, PyArrayObject* weights,
                 PyArrayObject* output, NI_ExtendMode mode,
                 double cvalue, npy_intp *origins)
{
    npy_bool *pf = NULL;
    npy_intp fsize, jj, kk, filter_size = 0, border_flag_value;
    npy_intp *offsets = NULL, *oo, size;
    NI_FilterIterator fi;
    NI_Iterator ii, io;
    char *pi, *po;
    npy_double *pw;
    npy_double *ww = NULL;
    int err = 0;
    NPY_BEGIN_THREADS_DEF;

    /* get the footprint: */
    fsize = PyArray_SIZE(weights);
    pw = (npy_double*)PyArray_DATA(weights);
    pf = malloc(fsize * sizeof(npy_bool));
    if (!pf) {
        PyErr_NoMemory();
        goto exit;
    }
    for(jj = 0; jj < fsize; jj++) {
        if (fabs(pw[jj]) > DBL_EPSILON) {
            pf[jj] = 1;
            ++filter_size;
        } else {
            pf[jj] = 0;
        }
    }
    /* copy the weights to contiguous memory: */
    ww = malloc(filter_size * sizeof(npy_double));
    if (!ww) {
        PyErr_NoMemory();
        goto exit;
    }
    jj = 0;
    for(kk = 0; kk < fsize; kk++) {
        if (pf[kk]) {
            ww[jj++] = pw[kk];
        }
    }
    /* initialize filter offsets: */
    if (!NI_InitFilterOffsets(input, pf, PyArray_DIMS(weights), origins,
                              mode, &offsets, &border_flag_value, NULL)) {
        goto exit;
    }
    /* initialize filter iterator: */
    if (!NI_InitFilterIterator(PyArray_NDIM(input), PyArray_DIMS(weights),
                               filter_size, PyArray_DIMS(input), origins,
                               &fi)) {
        goto exit;
    }
    /* initialize input element iterator: */
    if (!NI_InitPointIterator(input, &ii))
        goto exit;
    /* initialize output element iterator: */
    if (!NI_InitPointIterator(output, &io))
        goto exit;

    NPY_BEGIN_THREADS;
    /* get data pointers an array size: */
    pi = (void *)PyArray_DATA(input);
    po = (void *)PyArray_DATA(output);
    size = PyArray_SIZE(input);
    /* iterator over the elements: */
    oo = offsets;
    for(jj = 0; jj < size; jj++) {
        double tmp = 0.0;
        switch (PyArray_TYPE(input)) {
            CASE_CORRELATE_POINT(NPY_BOOL, npy_bool,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_UBYTE, npy_ubyte,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_USHORT, npy_ushort,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_UINT, npy_uint,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_ULONG, npy_ulong,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_ULONGLONG, npy_ulonglong,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_BYTE, npy_byte,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_SHORT, npy_short,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_INT, npy_int,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_LONG, npy_long,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_LONGLONG, npy_longlong,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_FLOAT, npy_float,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            CASE_CORRELATE_POINT(NPY_DOUBLE, npy_double,
                                 pi, ww, oo, filter_size, cvalue, tmp,
                                 border_flag_value);
            default:
                err = 1;
                goto exit;
        }
        switch (PyArray_TYPE(output)) {
            CASE_FILTER_OUT_SAFE(NPY_BOOL, npy_bool, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_UBYTE, npy_ubyte, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_USHORT, npy_ushort, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_UINT, npy_uint, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_ULONG, npy_ulong, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_ULONGLONG, npy_ulonglong, po, tmp);
            CASE_FILTER_OUT(NPY_BYTE, npy_byte, po, tmp);
            CASE_FILTER_OUT(NPY_SHORT, npy_short, po, tmp);
            CASE_FILTER_OUT(NPY_INT, npy_int, po, tmp);
            CASE_FILTER_OUT(NPY_LONG, npy_long, po, tmp);
            CASE_FILTER_OUT(NPY_LONGLONG, npy_longlong, po, tmp);
            CASE_FILTER_OUT(NPY_FLOAT, npy_float, po, tmp);
            CASE_FILTER_OUT(NPY_DOUBLE, npy_double, po, tmp);
            default:
                err = 1;
                goto exit;
        }
        NI_FILTER_NEXT2(fi, ii, io, oo, pi, po);
    }
exit:
    NPY_END_THREADS;
    if (err == 1) {
        PyErr_SetString(PyExc_RuntimeError, "array type not supported");
    }
    free(offsets);
    free(ww);
    free(pf);
    return PyErr_Occurred() ? 0 : 1;
}

int
NI_UniformFilter1D(PyArrayObject *input, npy_intp filter_size,
                   int axis, PyArrayObject *output, NI_ExtendMode mode,
                   double cval, npy_intp origin)
{
    npy_intp lines, kk, ll, length, size1, size2;
    int more;
    double *ibuffer = NULL, *obuffer = NULL;
    NI_LineBuffer iline_buffer, oline_buffer;
    NPY_BEGIN_THREADS_DEF;

    size1 = filter_size / 2;
    size2 = filter_size - size1 - 1;
    /* allocate and initialize the line buffers: */
    lines = -1;
    if (!NI_AllocateLineBuffer(input, axis, size1 + origin, size2 - origin,
                                                         &lines, BUFFER_SIZE, &ibuffer))
        goto exit;
    if (!NI_AllocateLineBuffer(output, axis, 0, 0, &lines, BUFFER_SIZE,
                                                         &obuffer))
        goto exit;
    if (!NI_InitLineBuffer(input, axis, size1 + origin, size2 - origin,
                                                            lines, ibuffer, mode, cval, &iline_buffer))
        goto exit;
    if (!NI_InitLineBuffer(output, axis, 0, 0, lines, obuffer, mode, 0.0,
                                                 &oline_buffer))
        goto exit;
    NPY_BEGIN_THREADS;
    length = PyArray_NDIM(input) > 0 ? PyArray_DIM(input, axis) : 1;

    /* iterate over all the array lines: */
    do {
        /* copy lines from array to buffer: */
        if (!NI_ArrayToLineBuffer(&iline_buffer, &lines, &more)) {
            goto exit;
        }
        /* iterate over the lines in the buffers: */
        for(kk = 0; kk < lines; kk++) {
            /* get lines: */
            double *iline = NI_GET_LINE(iline_buffer, kk);
            double *oline = NI_GET_LINE(oline_buffer, kk);
            /* do the uniform filter: */
            double tmp = 0.0;
            double *l1 = iline;
            double *l2 = iline + filter_size;
            for (ll = 0; ll < filter_size; ++ll) {
                tmp += iline[ll];
            }
            oline[0] = tmp / filter_size;
            for (ll = 1; ll < length; ++ll) {
                tmp += *l2++ - *l1++;
                oline[ll] = tmp / filter_size;
            }
        }
        /* copy lines from buffer to array: */
        if (!NI_LineBufferToArray(&oline_buffer)) {
            goto exit;
        }
    } while(more);

 exit:
    NPY_END_THREADS;
    free(ibuffer);
    free(obuffer);
    return PyErr_Occurred() ? 0 : 1;
}

#define INCREASE_RING_PTR(ptr) \
    (ptr)++;                   \
    if ((ptr) >= end) {        \
        (ptr) = ring;          \
    }

#define DECREASE_RING_PTR(ptr) \
    if ((ptr) == ring) {       \
        (ptr) = end;           \
    }                          \
    (ptr)--;

int
NI_MinOrMaxFilter1D(PyArrayObject *input, npy_intp filter_size,
                    int axis, PyArrayObject *output, NI_ExtendMode mode,
                    double cval, npy_intp origin, int minimum)
{
    npy_intp lines, kk, ll, length, size1, size2;
    int more;
    double *ibuffer = NULL, *obuffer = NULL;
    NI_LineBuffer iline_buffer, oline_buffer;

    struct pairs {
        double value;
        npy_intp death;
    } *ring = NULL, *minpair, *end, *last;

    NPY_BEGIN_THREADS_DEF;

    size1 = filter_size / 2;
    size2 = filter_size - size1 - 1;
    /* allocate and initialize the line buffers: */
    lines = -1;
    if (!NI_AllocateLineBuffer(input, axis, size1 + origin, size2 - origin,
                                            &lines, BUFFER_SIZE, &ibuffer))
        goto exit;
    if (!NI_AllocateLineBuffer(output, axis, 0, 0, &lines, BUFFER_SIZE,
                                                                &obuffer))
        goto exit;
    if (!NI_InitLineBuffer(input, axis, size1 + origin, size2 - origin,
                                lines, ibuffer, mode, cval, &iline_buffer))
        goto exit;
    if (!NI_InitLineBuffer(output, axis, 0, 0, lines, obuffer, mode, 0.0,
                                                            &oline_buffer))
        goto exit;

    NPY_BEGIN_THREADS;
    length = PyArray_NDIM(input) > 0 ? PyArray_DIM(input, axis) : 1;

    /* ring is a dequeue of pairs implemented as a circular array */
    ring = malloc(filter_size * sizeof(struct pairs));
    if (!ring) {
        goto exit;
    }
    end = ring + filter_size;

    /* iterate over all the array lines: */
    do {
        /* copy lines from array to buffer: */
        if (!NI_ArrayToLineBuffer(&iline_buffer, &lines, &more)) {
            goto exit;
        }
        /* iterate over the lines in the buffers: */
        for(kk = 0; kk < lines; kk++) {
            /* get lines: */
            double *iline = NI_GET_LINE(iline_buffer, kk);
            double *oline = NI_GET_LINE(oline_buffer, kk);

            /* This check could be moved out to the Python wrapper */
            if (filter_size == 1) {
                memcpy(oline, iline, sizeof(double) * length);
            }
            else {
                /*
                 * Original code by Richard Harter, adapted from:
                 * http://www.richardhartersworld.com/cri/2001/slidingmin.html
                 */
                minpair = ring;
                minpair->value = *iline++;
                minpair->death = filter_size;
                last = ring;

                for (ll = 1; ll < filter_size + length - 1; ll++) {
                    double val = *iline++;
                    if (minpair->death == ll) {
                        INCREASE_RING_PTR(minpair)
                    }
                    if ((minimum && val <= minpair->value) ||
                        (!minimum && val >= minpair->value)) {
                        minpair->value = val;
                        minpair->death = ll + filter_size;
                        last = minpair;
                    }
                    else {
                        while ((minimum && last->value >= val) ||
                               (!minimum && last->value <= val)) {
                            DECREASE_RING_PTR(last)
                        }
                        INCREASE_RING_PTR(last)
                        last->value = val;
                        last->death = ll + filter_size;
                    }
                    if (ll >= filter_size - 1) {
                        *oline++ = minpair->value;
                    }
                }
            }
        }
        /* copy lines from buffer to array: */
        if (!NI_LineBufferToArray(&oline_buffer)) {
            goto exit;
        }
    } while(more);

 exit:
    NPY_END_THREADS;
    free(ibuffer);
    free(obuffer);
    free(ring);
    return PyErr_Occurred() ? 0 : 1;
}

#undef DECREASE_RING_PTR
#undef INCREASE_RING_PTR


#define CASE_MIN_OR_MAX_POINT(_TYPE, _type, _pi, _offsets, _filter_size, \
                              _cval, _minimum, _res, _mv, _ss)           \
case _TYPE:                                                              \
{                                                                        \
    npy_intp _ii;                                                        \
    npy_intp _oo = _offsets[0];                                          \
    _type _tmp;                                                          \
    _type _cv = (_type)_cval;                                            \
    _res = _oo == _mv ? _cv : *(_type *)(_pi + _oo);                     \
    if (_ss != NULL) {                                                   \
        _res += _ss[0];                                                  \
    }                                                                    \
    for (_ii = 1; _ii < _filter_size; ++_ii) {                           \
        _oo = _offsets[_ii];                                             \
        _tmp = _oo == _mv ? _cv : *(_type *)(_pi + _oo);                 \
        if (_ss != NULL) {                                               \
            _tmp += (_type) _ss[_ii];                                    \
        }                                                                \
        if (_minimum) {                                                  \
            if (_tmp < _res) {                                           \
                _res = _tmp;                                             \
            }                                                            \
        }                                                                \
        else {                                                           \
            if (_tmp > _res) {                                           \
                _res = _tmp;                                             \
            }                                                            \
        }                                                                \
    }                                                                    \
}                                                                        \
break

int NI_MinOrMaxFilter(PyArrayObject* input, PyArrayObject* footprint,
                      PyArrayObject* structure, PyArrayObject* output,
                      NI_ExtendMode mode, double cvalue, npy_intp *origins,
                      int minimum)
{
    npy_bool *pf = NULL;
    npy_intp fsize, jj, kk, filter_size = 0, border_flag_value;
    npy_intp *offsets = NULL, *oo, size;
    NI_FilterIterator fi;
    NI_Iterator ii, io;
    char *pi, *po;
    int err = 0;
    double *ss = NULL;
    npy_double *ps;
    NPY_BEGIN_THREADS_DEF;

    /* get the footprint: */
    fsize = PyArray_SIZE(footprint);
    pf = (npy_bool*)PyArray_DATA(footprint);
    for(jj = 0; jj < fsize; jj++) {
        if (pf[jj]) {
            ++filter_size;
        }
    }
    /* get the structure: */
    if (structure) {
        ss = malloc(filter_size * sizeof(double));
        if (!ss) {
            PyErr_NoMemory();
            goto exit;
        }
        /* copy the weights to contiguous memory: */
        ps = (npy_double*)PyArray_DATA(structure);
        jj = 0;
        for(kk = 0; kk < fsize; kk++)
            if (pf[kk])
                ss[jj++] = minimum ? -ps[kk] : ps[kk];
    }
    /* initialize filter offsets: */
    if (!NI_InitFilterOffsets(input, pf, PyArray_DIMS(footprint), origins,
                              mode, &offsets, &border_flag_value, NULL)) {
        goto exit;
    }
    /* initialize filter iterator: */
    if (!NI_InitFilterIterator(PyArray_NDIM(input), PyArray_DIMS(footprint),
                               filter_size, PyArray_DIMS(input), origins,
                               &fi)) {
        goto exit;
    }
    /* initialize input element iterator: */
    if (!NI_InitPointIterator(input, &ii))
        goto exit;
    /* initialize output element iterator: */
    if (!NI_InitPointIterator(output, &io))
        goto exit;

    NPY_BEGIN_THREADS;

    /* get data pointers an array size: */
    pi = (void *)PyArray_DATA(input);
    po = (void *)PyArray_DATA(output);
    size = PyArray_SIZE(input);
    /* iterator over the elements: */
    oo = offsets;
    for(jj = 0; jj < size; jj++) {
        double tmp = 0.0;
        switch (PyArray_TYPE(input)) {
            CASE_MIN_OR_MAX_POINT(NPY_BOOL, npy_bool,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_UBYTE, npy_ubyte,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_USHORT, npy_ushort,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_UINT, npy_uint,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_ULONG, npy_ulong,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_ULONGLONG, npy_ulonglong,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_BYTE, npy_byte,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_SHORT, npy_short,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_INT, npy_int,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_LONG, npy_long,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_LONGLONG, npy_longlong,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_FLOAT, npy_float,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            CASE_MIN_OR_MAX_POINT(NPY_DOUBLE, npy_double,
                                  pi, oo, filter_size, cvalue, minimum, tmp,
                                  border_flag_value, ss);
            default:
                err = 1;
                goto exit;
        }
        switch (PyArray_TYPE(output)) {
            CASE_FILTER_OUT_SAFE(NPY_BOOL, npy_bool, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_UBYTE, npy_ubyte, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_USHORT, npy_ushort, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_UINT, npy_uint, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_ULONG, npy_ulong, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_ULONGLONG, npy_ulonglong, po, tmp);
            CASE_FILTER_OUT(NPY_BYTE, npy_byte, po, tmp);
            CASE_FILTER_OUT(NPY_SHORT, npy_short, po, tmp);
            CASE_FILTER_OUT(NPY_INT, npy_int, po, tmp);
            CASE_FILTER_OUT(NPY_LONG, npy_long, po, tmp);
            CASE_FILTER_OUT(NPY_LONGLONG, npy_longlong, po, tmp);
            CASE_FILTER_OUT(NPY_FLOAT, npy_float, po, tmp);
            CASE_FILTER_OUT(NPY_DOUBLE, npy_double, po, tmp);
            default:
                err = 1;
                goto exit;
        }
        NI_FILTER_NEXT2(fi, ii, io, oo, pi, po);
    }
exit:
    NPY_END_THREADS;
    if (err == 1) {
        PyErr_SetString(PyExc_RuntimeError, "array type not supported");
    }
    free(offsets);
    free(ss);
    return PyErr_Occurred() ? 0 : 1;
}

static double NI_Select(double *buffer, npy_intp min,
                        npy_intp max, npy_intp rank)
{
    npy_intp ii, jj;
    double x, t;

    if (min == max)
        return buffer[min];

    x = buffer[min];
    ii = min - 1;
    jj = max + 1;
    for(;;) {
        do
            jj--;
        while(buffer[jj] > x);
        do
            ii++;
        while(buffer[ii] < x);
        if (ii < jj) {
            t = buffer[ii];
            buffer[ii] = buffer[jj];
            buffer[jj] = t;
        } else {
            break;
        }
    }

    ii = jj - min + 1;
    if (rank < ii)
        return NI_Select(buffer, min, jj, rank);
    else
        return NI_Select(buffer, jj + 1, max, rank - ii);
}

#define CASE_RANK_POINT(_TYPE, _type, _pi, _offsets, _filter_size, \
                        _cval, _rank, _buffer, _res, _mv)          \
case _TYPE:                                                        \
{                                                                  \
    npy_intp _ii;                                                  \
    for (_ii = 0; _ii < _filter_size; ++_ii) {                     \
        npy_intp _offset = _offsets[_ii];                          \
        if (_offset == _mv) {                                      \
            _buffer[_ii] = (_type)_cval;                           \
        }                                                          \
        else {                                                     \
            _buffer[_ii] = *(_type *)(_pi + _offset);              \
        }                                                          \
    }                                                              \
    _res = (_type)NI_Select(_buffer, 0, _filter_size - 1, _rank);  \
}                                                                  \
break



int NI_RankFilter(PyArrayObject* input, int rank,
                  PyArrayObject* footprint, PyArrayObject* output,
                  NI_ExtendMode mode, double cvalue, npy_intp *origins)
{
    npy_intp fsize, jj, filter_size = 0, border_flag_value;
    npy_intp *offsets = NULL, *oo, size;
    NI_FilterIterator fi;
    NI_Iterator ii, io;
    char *pi, *po;
    npy_bool *pf = NULL;
    double *buffer = NULL;
    int err = 0;
    NPY_BEGIN_THREADS_DEF;

    /* get the footprint: */
    fsize = PyArray_SIZE(footprint);
    pf = (npy_bool*)PyArray_DATA(footprint);
    for(jj = 0; jj < fsize; jj++) {
        if (pf[jj]) {
            ++filter_size;
        }
    }
    /* buffer for rank calculation: */
    buffer = malloc(filter_size * sizeof(double));
    if (!buffer) {
        PyErr_NoMemory();
        goto exit;
    }
    /* iterator over the elements: */
    oo = offsets;
    /* initialize filter offsets: */
    if (!NI_InitFilterOffsets(input, pf, PyArray_DIMS(footprint), origins,
                              mode, &offsets, &border_flag_value, NULL)) {
        goto exit;
    }
    /* initialize filter iterator: */
    if (!NI_InitFilterIterator(PyArray_NDIM(input), PyArray_DIMS(footprint),
                               filter_size, PyArray_DIMS(input), origins,
                               &fi)) {
        goto exit;
    }
    /* initialize input element iterator: */
    if (!NI_InitPointIterator(input, &ii))
        goto exit;
    /* initialize output element iterator: */
    if (!NI_InitPointIterator(output, &io))
        goto exit;

    NPY_BEGIN_THREADS;
    /* get data pointers an array size: */
    pi = (void *)PyArray_DATA(input);
    po = (void *)PyArray_DATA(output);
    size = PyArray_SIZE(input);
    /* iterator over the elements: */
    oo = offsets;
    for(jj = 0; jj < size; jj++) {
        double tmp = 0.0;
        switch (PyArray_TYPE(input)) {
            CASE_RANK_POINT(NPY_BOOL, npy_bool,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_UBYTE, npy_ubyte,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_USHORT, npy_ushort,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_UINT, npy_uint,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_ULONG, npy_ulong,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_ULONGLONG, npy_ulonglong,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_BYTE, npy_byte,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_SHORT, npy_short,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_INT, npy_int,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_LONG, npy_long,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_LONGLONG, npy_longlong,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_FLOAT, npy_float,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            CASE_RANK_POINT(NPY_DOUBLE, npy_double,
                            pi, oo, filter_size, cvalue, rank, buffer, tmp,
                            border_flag_value);
            default:
                err = 1;
                goto exit;
        }
        switch (PyArray_TYPE(output)) {
            CASE_FILTER_OUT_SAFE(NPY_BOOL, npy_bool, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_UBYTE, npy_ubyte, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_USHORT, npy_ushort, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_UINT, npy_uint, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_ULONG, npy_ulong, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_ULONGLONG, npy_ulonglong, po, tmp);
            CASE_FILTER_OUT(NPY_BYTE, npy_byte, po, tmp);
            CASE_FILTER_OUT(NPY_SHORT, npy_short, po, tmp);
            CASE_FILTER_OUT(NPY_INT, npy_int, po, tmp);
            CASE_FILTER_OUT(NPY_LONG, npy_long, po, tmp);
            CASE_FILTER_OUT(NPY_LONGLONG, npy_longlong, po, tmp);
            CASE_FILTER_OUT(NPY_FLOAT, npy_float, po, tmp);
            CASE_FILTER_OUT(NPY_DOUBLE, npy_double, po, tmp);
            default:
                err = 1;
                goto exit;
        }
        NI_FILTER_NEXT2(fi, ii, io, oo, pi, po);
    }
exit:
    NPY_END_THREADS;
    if (err == 1) {
        PyErr_SetString(PyExc_RuntimeError, "array type not supported");
    }
    free(offsets);
    free(buffer);
    return PyErr_Occurred() ? 0 : 1;
}

int NI_GenericFilter1D(PyArrayObject *input,
            int (*function)(double*, npy_intp, double*, npy_intp, void*),
            void* data, npy_intp filter_size, int axis, PyArrayObject *output,
            NI_ExtendMode mode, double cval, npy_intp origin)
{
    int more;
    npy_intp ii, lines, length, size1, size2;
    double *ibuffer = NULL, *obuffer = NULL;
    NI_LineBuffer iline_buffer, oline_buffer;

    /* allocate and initialize the line buffers: */
    size1 = filter_size / 2;
    size2 = filter_size - size1 - 1;
    lines = -1;
    if (!NI_AllocateLineBuffer(input, axis, size1 + origin, size2 - origin,
                                                         &lines, BUFFER_SIZE, &ibuffer))
        goto exit;
    if (!NI_AllocateLineBuffer(output, axis, 0, 0, &lines, BUFFER_SIZE,
                                                         &obuffer))
        goto exit;
    if (!NI_InitLineBuffer(input, axis, size1 + origin, size2 - origin,
                                                            lines, ibuffer, mode, cval, &iline_buffer))
        goto exit;
    if (!NI_InitLineBuffer(output, axis, 0, 0, lines, obuffer, mode, 0.0,
                                                 &oline_buffer))
        goto exit;
    length = PyArray_NDIM(input) > 0 ? PyArray_DIM(input, axis) : 1;
    /* iterate over all the array lines: */
    do {
        /* copy lines from array to buffer: */
        if (!NI_ArrayToLineBuffer(&iline_buffer, &lines, &more)) {
            goto exit;
        }
        /* iterate over the lines in the buffers: */
        for(ii = 0; ii < lines; ii++) {
            /* get lines: */
            double *iline = NI_GET_LINE(iline_buffer, ii);
            double *oline = NI_GET_LINE(oline_buffer, ii);
            if (!function(iline, length + size1 + size2, oline, length, data)) {
                if (!PyErr_Occurred()) {
                    PyErr_SetString(PyExc_RuntimeError,
                                "unknown error in line processing function");
                }
                goto exit;
            }
        }
        /* copy lines from buffer to array: */
        if (!NI_LineBufferToArray(&oline_buffer)) {
            goto exit;
        }
    } while(more);
exit:
    free(ibuffer);
    free(obuffer);
    return PyErr_Occurred() ? 0 : 1;
}


#define CASE_FILTER_POINT(_TYPE, _type, _pi, _offsets, _filter_size, _cvalue, \
                          _res, _mv, _function, _data, _buffer)               \
case _TYPE:                                                                   \
{                                                                             \
    npy_intp _ii;                                                             \
    for (_ii = 0; _ii < _filter_size; ++_ii) {                                \
        const npy_intp _offset = _offsets[_ii];                               \
        if (_offset == _mv) {                                                 \
            _buffer[_ii] = (double)_cvalue;                                   \
        }                                                                     \
        else {                                                                \
            _buffer[_ii] = (double)(*(_type*)(_pi + _offset));                \
        }                                                                     \
    }                                                                         \
    if (!_function(_buffer, _filter_size, &_res, _data)) {                    \
        if (!PyErr_Occurred()) {                                              \
            PyErr_SetString(PyExc_RuntimeError,                               \
                            "unknown error in filter function");              \
            goto exit;                                                        \
        }                                                                     \
    }                                                                         \
}                                                                             \
break

int NI_GenericFilter(PyArrayObject* input,
            int (*function)(double*, npy_intp, double*, void*), void *data,
            PyArrayObject* footprint, PyArrayObject* output,
            NI_ExtendMode mode, double cvalue, npy_intp *origins)
{
    npy_bool *pf = NULL;
    npy_intp fsize, jj, filter_size = 0, border_flag_value;
    npy_intp *offsets = NULL, *oo, size;
    NI_FilterIterator fi;
    NI_Iterator ii, io;
    char *pi, *po;
    double *buffer = NULL;

    /* get the footprint: */
    fsize = PyArray_SIZE(footprint);
    pf = (npy_bool*)PyArray_DATA(footprint);
    for(jj = 0; jj < fsize; jj++) {
        if (pf[jj])
            ++filter_size;
    }
    /* initialize filter offsets: */
    if (!NI_InitFilterOffsets(input, pf, PyArray_DIMS(footprint), origins,
                              mode, &offsets, &border_flag_value, NULL)) {
        goto exit;
    }
    /* initialize filter iterator: */
    if (!NI_InitFilterIterator(PyArray_NDIM(input), PyArray_DIMS(footprint),
                               filter_size, PyArray_DIMS(input), origins,
                               &fi)) {
        goto exit;
    }
    /* initialize input element iterator: */
    if (!NI_InitPointIterator(input, &ii))
        goto exit;
    /* initialize output element iterator: */
    if (!NI_InitPointIterator(output, &io))
        goto exit;
    /* get data pointers an array size: */
    pi = (void *)PyArray_DATA(input);
    po = (void *)PyArray_DATA(output);
    size = PyArray_SIZE(input);
    /* buffer for filter calculation: */
    buffer = malloc(filter_size * sizeof(double));
    if (!buffer) {
        PyErr_NoMemory();
        goto exit;
    }
    /* iterate over the elements: */
    oo = offsets;
    for(jj = 0; jj < size; jj++) {
        double tmp = 0.0;
        switch (PyArray_TYPE(input)) {
            CASE_FILTER_POINT(NPY_BOOL, npy_bool,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_UBYTE, npy_ubyte,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_USHORT, npy_ushort,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_UINT, npy_uint,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_ULONG, npy_ulong,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_ULONGLONG, npy_ulonglong,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_BYTE, npy_byte,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_SHORT, npy_short,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_INT, npy_int,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_LONG, npy_long,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_LONGLONG, npy_longlong,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_FLOAT, npy_float,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            CASE_FILTER_POINT(NPY_DOUBLE, npy_double,
                              pi, oo, filter_size, cvalue, tmp,
                              border_flag_value, function, data, buffer);
            default:
                PyErr_SetString(PyExc_RuntimeError,
                                "array type not supported");
                goto exit;
        }
        switch (PyArray_TYPE(output)) {
            CASE_FILTER_OUT_SAFE(NPY_BOOL, npy_bool, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_UBYTE, npy_ubyte, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_USHORT, npy_ushort, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_UINT, npy_uint, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_ULONG, npy_ulong, po, tmp);
            CASE_FILTER_OUT_SAFE(NPY_ULONGLONG, npy_ulonglong, po, tmp);
            CASE_FILTER_OUT(NPY_BYTE, npy_byte, po, tmp);
            CASE_FILTER_OUT(NPY_SHORT, npy_short, po, tmp);
            CASE_FILTER_OUT(NPY_INT, npy_int, po, tmp);
            CASE_FILTER_OUT(NPY_LONG, npy_long, po, tmp);
            CASE_FILTER_OUT(NPY_LONGLONG, npy_longlong, po, tmp);
            CASE_FILTER_OUT(NPY_FLOAT, npy_float, po, tmp);
            CASE_FILTER_OUT(NPY_DOUBLE, npy_double, po, tmp);
            default:
                PyErr_SetString(PyExc_RuntimeError,
                                "array type not supported");
                goto exit;
        }
        NI_FILTER_NEXT2(fi, ii, io, oo, pi, po);
    }
exit:
    free(offsets);
    free(buffer);
    return PyErr_Occurred() ? 0 : 1;
}
