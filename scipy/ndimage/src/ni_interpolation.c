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
#include "ni_interpolation.h"
#include "ni_splines.h"
#include <stdlib.h>
#include <math.h>


/* map a coordinate outside the borders, according to the requested
     boundary condition: */
static double
map_coordinate(double in, npy_intp len, int mode)
{
    if (in < 0) {
        switch (mode) {
        case NI_EXTEND_MIRROR:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz2 = 2 * len - 2;
                in = sz2 * (npy_intp)(-in / sz2) + in;
                in = in <= 1 - len ? in + sz2 : -in;
            }
            break;
        case NI_EXTEND_REFLECT:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz2 = 2 * len;
                if (in < -sz2)
                    in = sz2 * (npy_intp)(-in / sz2) + in;
                // -1e-15 check to avoid possibility that: (-in - 1) == -1
                in = in < -len ? in + sz2 : (in > -1e-15 ? 1e-15 : -in) - 1;
            }
            break;
        case NI_EXTEND_WRAP:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz = len - 1;
                // Integer division of -in/sz gives (-in mod sz)
                // Note that 'in' is negative
                in += sz * ((npy_intp)(-in / sz) + 1);
            }
            break;
        case NI_EXTEND_GRID_WRAP:
            if (len <= 1) {
                in = 0;
            } else {
                // in = len - 1 + fmod(in + 1, len);
                in += len * ((npy_intp)((-1 - in) / len) + 1);
            }
            break;
        case NI_EXTEND_NEAREST:
            in = 0;
            break;
        case NI_EXTEND_CONSTANT:
            in = -1;
            break;
        }
    } else if (in > len-1) {
        switch (mode) {
        case NI_EXTEND_MIRROR:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz2 = 2 * len - 2;
                in -= sz2 * (npy_intp)(in / sz2);
                if (in >= len)
                    in = sz2 - in;
            }
            break;
        case NI_EXTEND_REFLECT:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz2 = 2 * len;
                in -= sz2 * (npy_intp)(in / sz2);
                if (in >= len)
                    in = sz2 - in - 1;
            }
            break;
        case NI_EXTEND_WRAP:
            if (len <= 1) {
                in = 0;
            } else {
                npy_intp sz = len - 1;
                in -= sz * (npy_intp)(in / sz);
            }
            break;
        case NI_EXTEND_GRID_WRAP:
            if (len <= 1) {
                in = 0;
            } else {
                in -= len * (npy_intp)(in / len);
            }
            break;
        case NI_EXTEND_NEAREST:
            in = len - 1;
            break;
        case NI_EXTEND_CONSTANT:
            in = -1;
            break;
        }
    }

    return in;
}

#define BUFFER_SIZE 256000
#define TOLERANCE 1e-15


/* one-dimensional spline filter: */
int NI_SplineFilter1D(PyArrayObject *input, int order, int axis,
                      NI_ExtendMode mode, PyArrayObject *output)
{
    int npoles = 0, more;
    npy_intp kk, lines, len;
    double *buffer = NULL, poles[MAX_SPLINE_FILTER_POLES];
    NI_LineBuffer iline_buffer, oline_buffer;
    NPY_BEGIN_THREADS_DEF;

    len = PyArray_NDIM(input) > 0 ? PyArray_DIM(input, axis) : 1;
    if (len < 1)
        goto exit;

    /* these are used in the spline filter calculation below: */
    if (get_filter_poles(order, &npoles, poles)) {
        goto exit;
    }

    /* allocate an initialize the line buffer, only a single one is used,
         because the calculation is in-place: */
    lines = -1;
    if (!NI_AllocateLineBuffer(input, axis, 0, 0, &lines, BUFFER_SIZE,
                               &buffer)) {
        goto exit;
    }
    if (!NI_InitLineBuffer(input, axis, 0, 0, lines, buffer,
                           NI_EXTEND_DEFAULT, 0.0, &iline_buffer)) {
        goto exit;
    }
    if (!NI_InitLineBuffer(output, axis, 0, 0, lines, buffer,
                           NI_EXTEND_DEFAULT, 0.0, &oline_buffer)) {
        goto exit;
    }
    NPY_BEGIN_THREADS;

    /* iterate over all the array lines: */
    do {
        /* copy lines from array to buffer: */
        if (!NI_ArrayToLineBuffer(&iline_buffer, &lines, &more)) {
            goto exit;
        }
        /* iterate over the lines in the buffer: */
        for(kk = 0; kk < lines; kk++) {
            /* get line: */
            double *ln = NI_GET_LINE(iline_buffer, kk);
            /* spline filter: */
            if (len > 1) {
                apply_filter(ln, len, poles, npoles, mode);
            }
        }

        /* copy lines from buffer to array: */
        if (!NI_LineBufferToArray(&oline_buffer)) {
            goto exit;
        }
    } while(more);

 exit:
    NPY_END_THREADS;
    free(buffer);
    return PyErr_Occurred() ? 0 : 1;
}

/* copy row of coordinate array from location at _p to _coor */
#define CASE_MAP_COORDINATES(_TYPE, _type, _p, _coor, _rank, _stride) \
case _TYPE:                                                           \
{                                                                     \
    npy_intp _hh;                                                     \
    for (_hh = 0; _hh < _rank; ++_hh) {                               \
        _coor[_hh] = *(_type *)_p;                                    \
        _p += _stride;                                                \
    }                                                                 \
}                                                                     \
break

#define CASE_INTERP_COEFF(_TYPE, _type, _coeff, _pi, _idx) \
case _TYPE:                                                \
    _coeff = *(_type *)(_pi + _idx);                       \
    break

#define CASE_INTERP_OUT(_TYPE, _type, _po, _t) \
case _TYPE:                                    \
    *(_type *)_po = (_type)_t;                 \
    break

#define CASE_INTERP_OUT_UINT(_TYPE, _type, _po, _t)  \
case NPY_##_TYPE:                                    \
    _t = _t > 0 ? _t + 0.5 : 0;                      \
    _t = _t > NPY_MAX_##_TYPE ? NPY_MAX_##_TYPE : t; \
    _t = _t < 0 ? 0 : t;                             \
    *(_type *)_po = (_type)_t;                       \
    break

#define CASE_INTERP_OUT_INT(_TYPE, _type, _po, _t)   \
case NPY_##_TYPE:                                    \
    _t = _t > 0 ? _t + 0.5 : _t - 0.5;               \
    _t = _t > NPY_MAX_##_TYPE ? NPY_MAX_##_TYPE : t; \
    _t = _t < NPY_MIN_##_TYPE ? NPY_MIN_##_TYPE : t; \
    *(_type *)_po = (_type)_t;                       \
    break

int _get_spline_boundary_mode(int mode)
{
    if ((mode == NI_EXTEND_CONSTANT) || (mode == NI_EXTEND_WRAP))
        // Modes without an analytic prefilter or explicit prepadding use
        // mirror extension.
        return NI_EXTEND_MIRROR;
    return mode;
}

int
NI_GeometricTransform(PyArrayObject *input, int (*map)(npy_intp*, double*,
                int, int, void*), void* map_data, PyArrayObject* matrix_ar,
                PyArrayObject* shift_ar, PyArrayObject *coordinates,
                PyArrayObject *output, int order, int mode, double cval,
                int nprepad)
{
    char *po, *pi, *pc = NULL;
    npy_intp **edge_offsets = NULL, **data_offsets = NULL, filter_size;
    char **edge_grid_const = NULL;
    npy_intp ftmp[NPY_MAXDIMS], *fcoordinates = NULL, *foffsets = NULL;
    npy_intp cstride = 0, kk, hh, ll, jj;
    npy_intp size;
    double **splvals = NULL, icoor[NPY_MAXDIMS], tmp;
    npy_intp idimensions[NPY_MAXDIMS], istrides[NPY_MAXDIMS];
    NI_Iterator io, ic;
    npy_double *matrix = matrix_ar ? (npy_double*)PyArray_DATA(matrix_ar) : NULL;
    npy_double *shift = shift_ar ? (npy_double*)PyArray_DATA(shift_ar) : NULL;
    int irank = 0, orank, spline_mode;
    NPY_BEGIN_THREADS_DEF;

    NPY_BEGIN_THREADS;

    for(kk = 0; kk < PyArray_NDIM(input); kk++) {
        idimensions[kk] = PyArray_DIM(input, kk);
        istrides[kk] = PyArray_STRIDE(input, kk);
    }
    irank = PyArray_NDIM(input);
    orank = PyArray_NDIM(output);

    /* if the mapping is from array coordinates: */
    if (coordinates) {
        /* initialize a line iterator along the first axis: */
        if (!NI_InitPointIterator(coordinates, &ic))
            goto exit;
        cstride = ic.strides[0];
        if (!NI_LineIterator(&ic, 0))
            goto exit;
        pc = (void *)(PyArray_DATA(coordinates));
    }

    /* offsets used at the borders: */
    edge_offsets = malloc(irank * sizeof(npy_intp*));
    data_offsets = malloc(irank * sizeof(npy_intp*));
    if (NPY_UNLIKELY(!edge_offsets || !data_offsets)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }

    if (mode == NI_EXTEND_GRID_CONSTANT) {
        // boolean indicating if the current point in the filter footprint is
        // outside the bounds
        edge_grid_const = malloc(irank * sizeof(char*));
        if (NPY_UNLIKELY(!edge_grid_const)) {
            NPY_END_THREADS;
            PyErr_NoMemory();
            goto exit;
        }
        for(jj = 0; jj < irank; jj++)
            edge_grid_const[jj] = NULL;
        for(jj = 0; jj < irank; jj++) {
            edge_grid_const[jj] = malloc((order + 1) * sizeof(char));
            if (NPY_UNLIKELY(!edge_grid_const[jj])) {
                NPY_END_THREADS;
                PyErr_NoMemory();
                goto exit;
            }
        }
    }

    for(jj = 0; jj < irank; jj++)
        data_offsets[jj] = NULL;
    for(jj = 0; jj < irank; jj++) {
        data_offsets[jj] = malloc((order + 1) * sizeof(npy_intp));
        if (NPY_UNLIKELY(!data_offsets[jj])) {
            NPY_END_THREADS;
            PyErr_NoMemory();
            goto exit;
        }
    }
    /* will hold the spline coefficients: */
    splvals = malloc(irank * sizeof(double*));
    if (NPY_UNLIKELY(!splvals)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }
    for(jj = 0; jj < irank; jj++)
        splvals[jj] = NULL;
    for(jj = 0; jj < irank; jj++) {
        splvals[jj] = malloc((order + 1) * sizeof(double));
        if (NPY_UNLIKELY(!splvals[jj])) {
            NPY_END_THREADS;
            PyErr_NoMemory();
            goto exit;
        }
    }

    filter_size = 1;
    for(jj = 0; jj < irank; jj++)
        filter_size *= order + 1;

    /* initialize output iterator: */
    if (!NI_InitPointIterator(output, &io))
        goto exit;

    /* get data pointers: */
    pi = (void *)PyArray_DATA(input);
    po = (void *)PyArray_DATA(output);

    /* make a table of all possible coordinates within the spline filter: */
    fcoordinates = malloc(irank * filter_size * sizeof(npy_intp));
    /* make a table of all offsets within the spline filter: */
    foffsets = malloc(filter_size * sizeof(npy_intp));
    if (NPY_UNLIKELY(!fcoordinates || !foffsets)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }
    for(jj = 0; jj < irank; jj++)
        ftmp[jj] = 0;
    kk = 0;
    for(hh = 0; hh < filter_size; hh++) {
        for(jj = 0; jj < irank; jj++)
            fcoordinates[jj + hh * irank] = ftmp[jj];
        foffsets[hh] = kk;
        for(jj = irank - 1; jj >= 0; jj--) {
            if (ftmp[jj] < order) {
                ftmp[jj]++;
                kk += istrides[jj];
                break;
            } else {
                ftmp[jj] = 0;
                kk -= istrides[jj] * order;
            }
        }
    }

    spline_mode = _get_spline_boundary_mode(mode);

    size = PyArray_SIZE(output);
    for(kk = 0; kk < size; kk++) {
        double t = 0.0;
        int constant = 0, edge = 0;
        npy_intp offset = 0;
        if (mode == NI_EXTEND_GRID_CONSTANT) {
            // reset edge flags for each location in the filter footprint
            for (hh = 0; hh < irank; hh++) {
                for(ll = 0; ll <= order; ll++) {
                    edge_grid_const[hh][ll] = 0;
                }
            }
        }
        if (map) {
            NPY_END_THREADS;
            /* call mappint functions: */
            if (!map(io.coordinates, icoor, orank, irank, map_data)) {
                if (!PyErr_Occurred())
                    PyErr_SetString(PyExc_RuntimeError,
                                    "unknown error in mapping function");
                goto exit;
            }
            NPY_BEGIN_THREADS;
        } else if (matrix) {
            /* do an affine transformation: */
            npy_double *p = matrix;
            for(hh = 0; hh < irank; hh++) {
                tmp = shift[hh];
                ll = 0;
                for (; ll + 1 < orank; ll += 2) {
                    tmp += io.coordinates[ll] * *p++;
                    tmp += io.coordinates[ll + 1] * *p++;
                }
                if (ll < orank) {
                    tmp += io.coordinates[ll] * *p++;
                }
                icoor[hh] = tmp;
            }
        } else if (coordinates) {
            /* mapping is from an coordinates array: */
            char *p = pc;
            switch (PyArray_TYPE(coordinates)) {
                CASE_MAP_COORDINATES(NPY_BOOL, npy_bool,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_UBYTE, npy_ubyte,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_USHORT, npy_ushort,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_UINT, npy_uint,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_ULONG, npy_ulong,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_ULONGLONG, npy_ulonglong,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_BYTE, npy_byte,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_SHORT, npy_short,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_INT, npy_int,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_LONG, npy_long,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_LONGLONG, npy_longlong,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_FLOAT, npy_float,
                                     p, icoor, irank, cstride);
                CASE_MAP_COORDINATES(NPY_DOUBLE, npy_double,
                                     p, icoor, irank, cstride);
            default:
                NPY_END_THREADS;
                PyErr_SetString(PyExc_RuntimeError,
                                "coordinate array data type not supported");
                goto exit;
            }
        }

        /* iterate over axes: */
        for(hh = 0; hh < irank; hh++) {
            double cc = icoor[hh] + nprepad;
            if ((mode != NI_EXTEND_GRID_CONSTANT) && (mode != NI_EXTEND_NEAREST)) {
                /* if the input coordinate is outside the borders, map it: */
                cc = map_coordinate(cc, idimensions[hh], mode);
            }
            if (cc > -1.0 || mode == NI_EXTEND_GRID_CONSTANT || mode == NI_EXTEND_NEAREST) {
                /* find the filter location along this axis: */
                npy_intp start;
                if (order & 1) {
                    start = (npy_intp)floor(cc) - order / 2;
                } else {
                    start = (npy_intp)floor(cc + 0.5) - order / 2;
                }
                /* get the offset to the start of the filter: */
                offset += istrides[hh] * start;
                npy_intp idx = 0;

                if (mode == NI_EXTEND_GRID_CONSTANT) {
                    // Determine locations in the filter footprint that are
                    // outside the range.
                    for(ll = 0; ll <= order; ll++) {
                        idx = start + ll;
                        edge_grid_const[hh][ll] = (idx < 0 || idx >= idimensions[hh]);
                    }
                } else {

                    if (start < 0 || start + order >= idimensions[hh]) {
                        /* implement border mapping, if outside border: */
                        edge = 1;
                        edge_offsets[hh] = data_offsets[hh];

                        for(ll = 0; ll <= order; ll++) {
                            idx = start + ll;
                            idx = (npy_intp)map_coordinate(idx, idimensions[hh], spline_mode);

                            /* calculate and store the offsets at this edge: */
                            edge_offsets[hh][ll] = istrides[hh] * (idx - start);
                        }
                    } else {
                        /* we are not at the border, use precalculated offsets: */
                        edge_offsets[hh] = NULL;
                    }
                }
                if(order!=0){
                    get_spline_interpolation_weights(cc, order, splvals[hh]);
                }
            } else {
                /* we use the constant border condition: */
                constant = 1;
                break;
            }
        }

        if (!constant) {
            npy_intp *ff = fcoordinates;
            const int type_num = PyArray_TYPE(input);
            t = 0.0;
            for(hh = 0; hh < filter_size; hh++) {
                double coeff = 0.0;
                npy_intp idx = 0;
                char is_cval = 0;
                if (mode == NI_EXTEND_GRID_CONSTANT) {
                    for(ll = 0; ll < irank; ll++) {
                        if (edge_grid_const[ll][ff[ll]]) {
                            is_cval = 1;
                        }
                    }
                }
                if (is_cval) {
                    coeff = cval;
                } else {
                    if (NPY_UNLIKELY(edge)) {
                        for(ll = 0; ll < irank; ll++) {
                            if (edge_offsets[ll])
                                idx += edge_offsets[ll][ff[ll]];
                            else
                                idx += ff[ll] * istrides[ll];
                        }
                    } else {
                        idx = foffsets[hh];
                    }
                    idx += offset;
                    switch (type_num) {
                        CASE_INTERP_COEFF(NPY_BOOL, npy_bool,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_UBYTE, npy_ubyte,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_USHORT, npy_ushort,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_UINT, npy_uint,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_ULONG, npy_ulong,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_ULONGLONG, npy_ulonglong,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_BYTE, npy_byte,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_SHORT, npy_short,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_INT, npy_int,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_LONG, npy_long,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_LONGLONG, npy_longlong,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_FLOAT, npy_float,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_DOUBLE, npy_double,
                                          coeff, pi, idx);
                    default:
                        NPY_END_THREADS;
                        PyErr_SetString(PyExc_RuntimeError,
                                        "data type not supported");
                        goto exit;
                    }
                }
                /* calculate the interpolated value: */
                for(ll = 0; ll < irank; ll++)
                    if (order > 0)
                        coeff *= splvals[ll][ff[ll]];
                t += coeff;
                ff += irank;
            }
        } else {
            t = cval;
        }
        /* store output value: */
        switch (PyArray_TYPE(output)) {
            CASE_INTERP_OUT(NPY_BOOL, npy_bool, po, t);
            CASE_INTERP_OUT_UINT(UBYTE, npy_ubyte, po, t);
            CASE_INTERP_OUT_UINT(USHORT, npy_ushort, po, t);
            CASE_INTERP_OUT_UINT(UINT, npy_uint, po, t);
            CASE_INTERP_OUT_UINT(ULONG, npy_ulong, po, t);
            CASE_INTERP_OUT_UINT(ULONGLONG, npy_ulonglong, po, t);
            CASE_INTERP_OUT_INT(BYTE, npy_byte, po, t);
            CASE_INTERP_OUT_INT(SHORT, npy_short, po, t);
            CASE_INTERP_OUT_INT(INT, npy_int, po, t);
            CASE_INTERP_OUT_INT(LONG, npy_long, po, t);
            CASE_INTERP_OUT_INT(LONGLONG, npy_longlong, po, t);
            CASE_INTERP_OUT(NPY_FLOAT, npy_float, po, t);
            CASE_INTERP_OUT(NPY_DOUBLE, npy_double, po, t);
        default:
            NPY_END_THREADS;
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        if (coordinates) {
            NI_ITERATOR_NEXT2(io, ic, po, pc);
        } else {
            NI_ITERATOR_NEXT(io, po);
        }
    }

 exit:
    NPY_END_THREADS;
    free(edge_offsets);
    if (edge_grid_const) {
        for(jj = 0; jj < irank; jj++)
            free(edge_grid_const[jj]);
        free(edge_grid_const);
    }
    if (data_offsets) {
        for(jj = 0; jj < irank; jj++)
            free(data_offsets[jj]);
        free(data_offsets);
    }
    if (splvals) {
        for(jj = 0; jj < irank; jj++)
            free(splvals[jj]);
        free(splvals);
    }
    free(foffsets);
    free(fcoordinates);
    return PyErr_Occurred() ? 0 : 1;
}

int NI_ZoomShift(PyArrayObject *input, PyArrayObject* zoom_ar,
                 PyArrayObject* shift_ar, PyArrayObject *output,
                 int order, int mode, double cval, int nprepad, int grid_mode)
{
    char *po, *pi;
    npy_intp **zeros = NULL, **offsets = NULL, ***edge_offsets = NULL;
    npy_intp ftmp[NPY_MAXDIMS], *fcoordinates = NULL, *foffsets = NULL;
    npy_intp jj, hh, kk, filter_size, odimensions[NPY_MAXDIMS];
    npy_intp idimensions[NPY_MAXDIMS], istrides[NPY_MAXDIMS];
    npy_intp size;
    char ***edge_grid_const = NULL;
    double ***splvals = NULL;
    NI_Iterator io;
    npy_double *zooms = zoom_ar ? (npy_double*)PyArray_DATA(zoom_ar) : NULL;
    npy_double *shifts = shift_ar ? (npy_double*)PyArray_DATA(shift_ar) : NULL;
    int rank = 0;
    NPY_BEGIN_THREADS_DEF;

    NPY_BEGIN_THREADS;

    for (kk = 0; kk < PyArray_NDIM(input); kk++) {
        idimensions[kk] = PyArray_DIM(input, kk);
        istrides[kk] = PyArray_STRIDE(input, kk);
        odimensions[kk] = PyArray_DIM(output, kk);
    }
    rank = PyArray_NDIM(input);

    /* if the mode is 'constant' we need some temps later: */
    if (mode == NI_EXTEND_CONSTANT) {
        zeros = malloc(rank * sizeof(npy_intp*));
        if (NPY_UNLIKELY(!zeros)) {
            NPY_END_THREADS;
            PyErr_NoMemory();
            goto exit;
        }
        for(jj = 0; jj < rank; jj++)
            zeros[jj] = NULL;
        for(jj = 0; jj < rank; jj++) {
            zeros[jj] = malloc(odimensions[jj] * sizeof(npy_intp));
            if (NPY_UNLIKELY(!zeros[jj])) {
                NPY_END_THREADS;
                PyErr_NoMemory();
                goto exit;
            }
        }
    } else if (mode == NI_EXTEND_GRID_CONSTANT) {
        // boolean indicating if the current point in the filter footprint is
        // outside the bounds
        edge_grid_const = malloc(rank * sizeof(char*));
        if (NPY_UNLIKELY(!edge_grid_const)) {
            NPY_END_THREADS;
            PyErr_NoMemory();
            goto exit;
        }
        for(jj = 0; jj < rank; jj++)
            edge_grid_const[jj] = NULL;
        for(jj = 0; jj < rank; jj++) {
            edge_grid_const[jj] = malloc(odimensions[jj] * sizeof(char*));
            if (NPY_UNLIKELY(!edge_grid_const[jj])) {
                NPY_END_THREADS;
                PyErr_NoMemory();
                goto exit;
            }
            for(hh = 0; hh < odimensions[jj]; hh++) {
                edge_grid_const[jj][hh] = NULL;
            }
        }
    }
    /* store offsets, along each axis: */
    offsets = malloc(rank * sizeof(npy_intp*));
    /* store spline coefficients, along each axis: */
    splvals = malloc(rank * sizeof(double**));
    /* store offsets at all edges: */

    if (NPY_UNLIKELY(!offsets || !splvals)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }
    for(jj = 0; jj < rank; jj++) {
        offsets[jj] = NULL;
        splvals[jj] = NULL;
    }
    for(jj = 0; jj < rank; jj++) {
        offsets[jj] = malloc(odimensions[jj] * sizeof(npy_intp));
        splvals[jj] = malloc(odimensions[jj] * sizeof(double*));
        if (NPY_UNLIKELY(!offsets[jj] || !splvals[jj])) {
            NPY_END_THREADS;
            PyErr_NoMemory();
            goto exit;
        }
        for(hh = 0; hh < odimensions[jj]; hh++) {
            splvals[jj][hh] = NULL;
        }
    }

    if (mode != NI_EXTEND_GRID_CONSTANT){
        edge_offsets = malloc(rank * sizeof(npy_intp**));
        if (NPY_UNLIKELY(!edge_offsets)) {
            NPY_END_THREADS;
            PyErr_NoMemory();
            goto exit;
        }
        for(jj = 0; jj < rank; jj++) {
            edge_offsets[jj] = NULL;
        }
        for(jj = 0; jj < rank; jj++) {
            edge_offsets[jj] = malloc(odimensions[jj] * sizeof(npy_intp*));
            if (NPY_UNLIKELY(!edge_offsets[jj])) {
                NPY_END_THREADS;
                PyErr_NoMemory();
                goto exit;
            }
            for(hh = 0; hh < odimensions[jj]; hh++) {
                edge_offsets[jj][hh] = NULL;
            }
        }
    }

    int spline_mode = _get_spline_boundary_mode(mode);

    for(jj = 0; jj < rank; jj++) {
        double shift = 0.0, zoom = 0.0;
        if (shifts)
            shift = shifts[jj];
        if (zooms)
            zoom = zooms[jj];
        for(kk = 0; kk < odimensions[jj]; kk++) {
            double cc = (double)kk;
            if (shifts)
                cc += shift;
            if (zooms)
            {
                if (grid_mode)
                {
                    cc += 0.5;
                    cc *= zoom;
                    cc -= 0.5;
                } else {
                    cc *= zoom;
                }
            }
            cc += (double)nprepad;
            if ((mode != NI_EXTEND_GRID_CONSTANT) && (mode != NI_EXTEND_NEAREST)) {
                /* if the input coordinate is outside the borders, map it: */
                cc = map_coordinate(cc, idimensions[jj], mode);
            }
            if (cc > -1.0 || mode == NI_EXTEND_GRID_CONSTANT || mode == NI_EXTEND_NEAREST) {
                npy_intp start;
                if (zeros && zeros[jj])
                    zeros[jj][kk] = 0;
                if (order & 1) {
                    start = (npy_intp)floor(cc) - order / 2;
                } else {
                    start = (npy_intp)floor(cc + 0.5) - order / 2;
                }
                offsets[jj][kk] = istrides[jj] * start;
                if (start < 0 || start + order >= idimensions[jj]) {
                    npy_intp idx = 0;

                    if (mode == NI_EXTEND_GRID_CONSTANT) {
                        edge_grid_const[jj][kk] = malloc((order + 1) * sizeof(char));
                        if (NPY_UNLIKELY(!edge_grid_const[jj][kk])) {
                            NPY_END_THREADS;
                            PyErr_NoMemory();
                            goto exit;
                        }
                        for(hh = 0; hh <= order; hh++) {
                            idx = start + hh;
                            edge_grid_const[jj][kk][hh] = (idx < 0 || idx >= idimensions[jj]);
                        }
                    } else {
                        edge_offsets[jj][kk] = malloc((order + 1) * sizeof(npy_intp));
                        if (NPY_UNLIKELY(!edge_offsets[jj][kk])) {
                            NPY_END_THREADS;
                            PyErr_NoMemory();
                            goto exit;
                        }
                        for(hh = 0; hh <= order; hh++) {
                            idx = start + hh;
                            idx = (npy_intp)map_coordinate(idx, idimensions[jj], spline_mode);
                            edge_offsets[jj][kk][hh] = istrides[jj] * (idx - start);
                        }
                    }

                }
                if (order > 0) {
                    splvals[jj][kk] = malloc((order + 1) * sizeof(double));
                    if (NPY_UNLIKELY(!splvals[jj][kk])) {
                        NPY_END_THREADS;
                        PyErr_NoMemory();
                        goto exit;
                    }
                    get_spline_interpolation_weights(cc, order, splvals[jj][kk]);
                }
            } else {
                zeros[jj][kk] = 1;
            }
        }
    }

    filter_size = 1;
    for(jj = 0; jj < rank; jj++)
        filter_size *= order + 1;

    if (!NI_InitPointIterator(output, &io))
        goto exit;

    pi = (void *)PyArray_DATA(input);
    po = (void *)PyArray_DATA(output);

    /* store all coordinates and offsets with filter: */
    fcoordinates = malloc(rank * filter_size * sizeof(npy_intp));
    foffsets = malloc(filter_size * sizeof(npy_intp));
    if (NPY_UNLIKELY(!fcoordinates || !foffsets)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }

    for(jj = 0; jj < rank; jj++)
        ftmp[jj] = 0;
    kk = 0;
    for(hh = 0; hh < filter_size; hh++) {
        for(jj = 0; jj < rank; jj++)
            fcoordinates[jj + hh * rank] = ftmp[jj];
        foffsets[hh] = kk;
        for(jj = rank - 1; jj >= 0; jj--) {
            if (ftmp[jj] < order) {
                ftmp[jj]++;
                kk += istrides[jj];
                break;
            } else {
                ftmp[jj] = 0;
                kk -= istrides[jj] * order;
            }
        }
    }
    size = PyArray_SIZE(output);
    for(kk = 0; kk < size; kk++) {
        double t = 0.0;
        npy_intp edge = 0, oo = 0, zero = 0;

        for(hh = 0; hh < rank; hh++) {
            if (zeros && zeros[hh][io.coordinates[hh]]) {
                /* we use constant border condition */
                zero = 1;
                break;
            }
            oo += offsets[hh][io.coordinates[hh]];
            if (mode != NI_EXTEND_GRID_CONSTANT) {
                if (edge_offsets[hh][io.coordinates[hh]])
                    edge = 1;
            }         }

        if (!zero) {
            npy_intp *ff = fcoordinates;
            const int type_num = PyArray_TYPE(input);
            t = 0.0;
            for(hh = 0; hh < filter_size; hh++) {
                npy_intp idx = 0;
                double coeff = 0.0;
                int is_cval = 0;
                if (mode == NI_EXTEND_GRID_CONSTANT)
                {
                    for(jj = 0; jj < rank; jj++) {
                        if (edge_grid_const[jj][io.coordinates[jj]])
                        {
                            if (edge_grid_const[jj][io.coordinates[jj]][ff[jj]])
                                is_cval = 1;
                        }
                    }
                }
                if (is_cval) {
                    coeff = cval;
                } else {
                    if (NPY_UNLIKELY(edge)) {
                        /* use precalculated edge offsets: */
                        for(jj = 0; jj < rank; jj++) {
                            if (edge_offsets[jj][io.coordinates[jj]])
                                idx += edge_offsets[jj][io.coordinates[jj]][ff[jj]];
                            else
                                idx += ff[jj] * istrides[jj];
                        }
                        idx += oo;
                    } else {
                        /* use normal offsets: */
                        idx += oo + foffsets[hh];
                    }
                    switch (type_num) {
                        CASE_INTERP_COEFF(NPY_BOOL, npy_bool,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_UBYTE, npy_ubyte,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_USHORT, npy_ushort,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_UINT, npy_uint,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_ULONG, npy_ulong,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_ULONGLONG, npy_ulonglong,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_BYTE, npy_byte,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_SHORT, npy_short,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_INT, npy_int,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_LONG, npy_long,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_LONGLONG, npy_longlong,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_FLOAT, npy_float,
                                          coeff, pi, idx);
                        CASE_INTERP_COEFF(NPY_DOUBLE, npy_double,
                                          coeff, pi, idx);
                    default:
                        NPY_END_THREADS;
                        PyErr_SetString(PyExc_RuntimeError,
                                        "data type not supported");
                        goto exit;
                    }
                }
                /* calculate interpolated value: */
                for(jj = 0; jj < rank; jj++)
                    if (order > 0)
                        coeff *= splvals[jj][io.coordinates[jj]][ff[jj]];
                t += coeff;
                ff += rank;
            }
        } else {
            t = cval;
        }
        /* store output: */
        switch (PyArray_TYPE(output)) {
            CASE_INTERP_OUT(NPY_BOOL, npy_bool, po, t);
            CASE_INTERP_OUT_UINT(UBYTE, npy_ubyte, po, t);
            CASE_INTERP_OUT_UINT(USHORT, npy_ushort, po, t);
            CASE_INTERP_OUT_UINT(UINT, npy_uint, po, t);
            CASE_INTERP_OUT_UINT(ULONG, npy_ulong, po, t);
            CASE_INTERP_OUT_UINT(ULONGLONG, npy_ulonglong, po, t);
            CASE_INTERP_OUT_INT(BYTE, npy_byte, po, t);
            CASE_INTERP_OUT_INT(SHORT, npy_short, po, t);
            CASE_INTERP_OUT_INT(INT, npy_int, po, t);
            CASE_INTERP_OUT_INT(LONG, npy_long, po, t);
            CASE_INTERP_OUT_INT(LONGLONG, npy_longlong, po, t);
            CASE_INTERP_OUT(NPY_FLOAT, npy_float, po, t);
            CASE_INTERP_OUT(NPY_DOUBLE, npy_double, po, t);
        default:
            NPY_END_THREADS;
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        NI_ITERATOR_NEXT(io, po);
    }

 exit:
    NPY_END_THREADS;
    if (zeros) {
        for(jj = 0; jj < rank; jj++)
            free(zeros[jj]);
        free(zeros);
    }
    if (offsets) {
        for(jj = 0; jj < rank; jj++)
            free(offsets[jj]);
        free(offsets);
    }
    if (splvals) {
        for(jj = 0; jj < rank; jj++) {
            if (splvals[jj]) {
                for(hh = 0; hh < odimensions[jj]; hh++)
                    free(splvals[jj][hh]);
                free(splvals[jj]);
            }
        }
        free(splvals);
    }
    if (edge_offsets) {
        for(jj = 0; jj < rank; jj++) {
            if (edge_offsets[jj]) {
                for(hh = 0; hh < odimensions[jj]; hh++)
                    free(edge_offsets[jj][hh]);
                free(edge_offsets[jj]);
            }
        }
        free(edge_offsets);
    }
    if (edge_grid_const) {
        for(jj = 0; jj < rank; jj++) {
            if (edge_grid_const[jj]) {
                for(hh = 0; hh < odimensions[jj]; hh++)
                    free(edge_grid_const[jj][hh]);
                free(edge_grid_const[jj]);
            }
        }
        free(edge_grid_const);
    }
    free(foffsets);
    free(fcoordinates);
    return PyErr_Occurred() ? 0 : 1;
}
