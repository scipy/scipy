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
#include <string.h>
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


static void
NI_CopyToDoubleLine(char *input_data, int input_type, npy_intp start,
                    npy_intp stride, npy_intp len, double *line)
{
    npy_intp ii;

    if (input_type == NPY_FLOAT) {
        const npy_float *input = (const npy_float *)input_data;
        for (ii = 0; ii < len; ii++) {
            line[ii] = input[start + ii * stride];
        }
    } else {
        const npy_double *input = (const npy_double *)input_data;
        for (ii = 0; ii < len; ii++) {
            line[ii] = input[start + ii * stride];
        }
    }
}

static void
NI_CopyDoubleLineToOutput(double *output, npy_intp start, npy_intp stride,
                          npy_intp len, const double *line)
{
    npy_intp ii;

    for (ii = 0; ii < len; ii++) {
        output[start + ii * stride] = line[ii];
    }
}

static void
NI_FilterContiguousDoubleLine(char *input_data, int input_type, double *output,
                              npy_intp start, npy_intp len,
                              const double *poles, int npoles,
                              NI_ExtendMode mode)
{
    if (input_type == NPY_FLOAT) {
        const npy_float *input = (const npy_float *)input_data;
        npy_intp ii;
        for (ii = 0; ii < len; ii++) {
            output[start + ii] = input[start + ii];
        }
    } else {
        const npy_double *input = (const npy_double *)input_data;
        if (input != output) {
            memcpy(output + start, input + start, len * sizeof(double));
        }
    }

    if (len > 1) {
        apply_filter(output + start, len, poles, npoles, mode);
    }
}

static void
NI_GetCContiguousLineGeometry(npy_intp *dims, int rank, int axis,
                              npy_intp line_index, npy_intp *start,
                              npy_intp *stride)
{
    if (rank == 1) {
        *start = 0;
        *stride = 1;
    } else if (rank == 2) {
        if (axis == 0) {
            *start = line_index;
            *stride = dims[1];
        } else {
            *start = line_index * dims[1];
            *stride = 1;
        }
    } else {
        if (axis == 0) {
            *start = line_index;
            *stride = dims[1] * dims[2];
        } else if (axis == 1) {
            const npy_intp ii = line_index / dims[2];
            const npy_intp kk = line_index - ii * dims[2];
            *start = ii * dims[1] * dims[2] + kk;
            *stride = dims[2];
        } else {
            *start = line_index * dims[2];
            *stride = 1;
        }
    }
}

static int
NI_SplineFilter1DFastPath(PyArrayObject *input, int order, int axis,
                          NI_ExtendMode mode, PyArrayObject *output,
                          const double *poles, int npoles, int *handled)
{
    double *line = NULL;
    char *input_data = NULL;
    double *output_data = NULL;
    npy_intp dims[NPY_MAXDIMS], ii, len, line_count, line_index, size;
    int rank, input_type, jj;
    NPY_BEGIN_THREADS_DEF;

    *handled = 0;

    if (order != 3 || (mode != NI_EXTEND_MIRROR &&
            mode != NI_EXTEND_REFLECT)) {
        return 1;
    }

    rank = PyArray_NDIM(input);
    if (rank < 1 || rank > 3 || PyArray_NDIM(output) != rank) {
        return 1;
    }
    if (!PyArray_IS_C_CONTIGUOUS(input) ||
            !PyArray_IS_C_CONTIGUOUS(output)) {
        return 1;
    }

    input_type = PyArray_TYPE(input);
    if ((input_type != NPY_FLOAT && input_type != NPY_DOUBLE) ||
            PyArray_TYPE(output) != NPY_DOUBLE) {
        return 1;
    }

    for (jj = 0; jj < rank; jj++) {
        dims[jj] = PyArray_DIM(input, jj);
        if (dims[jj] != PyArray_DIM(output, jj)) {
            return 1;
        }
    }

    len = dims[axis];
    size = PyArray_SIZE(input);
    if (size == 0) {
        *handled = 1;
        return 1;
    }
    if (len <= 1) {
        input_data = (char *)PyArray_DATA(input);
        output_data = (double *)PyArray_DATA(output);
        if (input_type == NPY_FLOAT) {
            const npy_float *input_float = (const npy_float *)input_data;
            for (ii = 0; ii < size; ii++) {
                output_data[ii] = input_float[ii];
            }
        } else {
            const npy_double *input_double = (const npy_double *)input_data;
            if (input_double != output_data) {
                memcpy(output_data, input_double, size * sizeof(double));
            }
        }
        *handled = 1;
        return 1;
    }

    line_count = size / len;
    input_data = (char *)PyArray_DATA(input);
    output_data = (double *)PyArray_DATA(output);

    if (axis != rank - 1) {
        line = malloc(len * sizeof(double));
        if (NPY_UNLIKELY(!line)) {
            PyErr_NoMemory();
            return 0;
        }
    }

    NPY_BEGIN_THREADS;
    for (line_index = 0; line_index < line_count; line_index++) {
        npy_intp start, stride;
        NI_GetCContiguousLineGeometry(dims, rank, axis, line_index,
                                      &start, &stride);
        if (stride == 1) {
            NI_FilterContiguousDoubleLine(input_data, input_type, output_data,
                                          start, len, poles, npoles, mode);
        } else {
            NI_CopyToDoubleLine(input_data, input_type, start, stride, len,
                                line);
            apply_filter(line, len, poles, npoles, mode);
            NI_CopyDoubleLineToOutput(output_data, start, stride, len, line);
        }
    }
    NPY_END_THREADS;

    free(line);
    *handled = 1;
    return 1;
}

/* one-dimensional spline filter: */
int NI_SplineFilter1D(PyArrayObject *input, int order, int axis,
                      NI_ExtendMode mode, PyArrayObject *output)
{
    int npoles = 0, more, fast_path_handled = 0;
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

    if (!NI_SplineFilter1DFastPath(input, order, axis, mode, output,
                                   poles, npoles, &fast_path_handled)) {
        goto exit;
    }
    if (fast_path_handled) {
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

#define CASE_ZOOM_LINEAR_FAST(_TYPE, _type)                                \
case _TYPE:                                                                \
{                                                                          \
    const _type *in = (const _type *)PyArray_DATA(input);                   \
    _type *out = (_type *)PyArray_DATA(output);                             \
    npy_intp ii, jj, kk;                                                    \
                                                                           \
    if (rank == 1) {                                                        \
        for (ii = 0; ii < odimensions[0]; ii++) {                           \
            double t = 0.0, coeff;                                          \
            coeff = in[offset0[0][ii]];                                     \
            coeff *= weight0[0][ii];                                        \
            t += coeff;                                                     \
            coeff = in[offset1[0][ii]];                                     \
            coeff *= weight1[0][ii];                                        \
            t += coeff;                                                     \
            out[ii] = (_type)t;                                             \
        }                                                                  \
    } else if (rank == 2) {                                                 \
        npy_intp oo = 0;                                                    \
        for (ii = 0; ii < odimensions[0]; ii++) {                           \
            const npy_intp off00 = offset0[0][ii];                          \
            const npy_intp off01 = offset1[0][ii];                          \
            const double w00 = weight0[0][ii];                              \
            const double w01 = weight1[0][ii];                              \
            for (jj = 0; jj < odimensions[1]; jj++, oo++) {                 \
                double t = 0.0, coeff;                                      \
                coeff = in[off00 + offset0[1][jj]];                         \
                coeff *= w00;                                               \
                coeff *= weight0[1][jj];                                    \
                t += coeff;                                                 \
                coeff = in[off00 + offset1[1][jj]];                         \
                coeff *= w00;                                               \
                coeff *= weight1[1][jj];                                    \
                t += coeff;                                                 \
                coeff = in[off01 + offset0[1][jj]];                         \
                coeff *= w01;                                               \
                coeff *= weight0[1][jj];                                    \
                t += coeff;                                                 \
                coeff = in[off01 + offset1[1][jj]];                         \
                coeff *= w01;                                               \
                coeff *= weight1[1][jj];                                    \
                t += coeff;                                                 \
                out[oo] = (_type)t;                                         \
            }                                                              \
        }                                                                  \
    } else {                                                               \
        npy_intp oo = 0;                                                    \
        for (ii = 0; ii < odimensions[0]; ii++) {                           \
            const npy_intp off00 = offset0[0][ii];                          \
            const npy_intp off01 = offset1[0][ii];                          \
            const double w00 = weight0[0][ii];                              \
            const double w01 = weight1[0][ii];                              \
            for (jj = 0; jj < odimensions[1]; jj++) {                       \
                const npy_intp off10 = offset0[1][jj];                      \
                const npy_intp off11 = offset1[1][jj];                      \
                const double w10 = weight0[1][jj];                          \
                const double w11 = weight1[1][jj];                          \
                for (kk = 0; kk < odimensions[2]; kk++, oo++) {             \
                    double t = 0.0, coeff;                                  \
                    coeff = in[off00 + off10 + offset0[2][kk]];             \
                    coeff *= w00;                                           \
                    coeff *= w10;                                           \
                    coeff *= weight0[2][kk];                                \
                    t += coeff;                                             \
                    coeff = in[off00 + off10 + offset1[2][kk]];             \
                    coeff *= w00;                                           \
                    coeff *= w10;                                           \
                    coeff *= weight1[2][kk];                                \
                    t += coeff;                                             \
                    coeff = in[off00 + off11 + offset0[2][kk]];             \
                    coeff *= w00;                                           \
                    coeff *= w11;                                           \
                    coeff *= weight0[2][kk];                                \
                    t += coeff;                                             \
                    coeff = in[off00 + off11 + offset1[2][kk]];             \
                    coeff *= w00;                                           \
                    coeff *= w11;                                           \
                    coeff *= weight1[2][kk];                                \
                    t += coeff;                                             \
                    coeff = in[off01 + off10 + offset0[2][kk]];             \
                    coeff *= w01;                                           \
                    coeff *= w10;                                           \
                    coeff *= weight0[2][kk];                                \
                    t += coeff;                                             \
                    coeff = in[off01 + off10 + offset1[2][kk]];             \
                    coeff *= w01;                                           \
                    coeff *= w10;                                           \
                    coeff *= weight1[2][kk];                                \
                    t += coeff;                                             \
                    coeff = in[off01 + off11 + offset0[2][kk]];             \
                    coeff *= w01;                                           \
                    coeff *= w11;                                           \
                    coeff *= weight0[2][kk];                                \
                    t += coeff;                                             \
                    coeff = in[off01 + off11 + offset1[2][kk]];             \
                    coeff *= w01;                                           \
                    coeff *= w11;                                           \
                    coeff *= weight1[2][kk];                                \
                    t += coeff;                                             \
                    out[oo] = (_type)t;                                     \
                }                                                          \
            }                                                              \
        }                                                                  \
    }                                                                      \
}                                                                          \
break

#define ZOOM_CUBIC_FAST(_in_type, _out_type)                               \
do {                                                                       \
    const _in_type *in = (const _in_type *)PyArray_DATA(input);             \
    _out_type *out = (_out_type *)PyArray_DATA(output);                     \
    npy_intp ii, jj, kk, ll, mm, nn;                                        \
                                                                           \
    if (rank == 1) {                                                        \
        for (ii = 0; ii < odimensions[0]; ii++) {                           \
            const npy_intp base0 = 4 * ii;                                  \
            double t = 0.0;                                                 \
            for (ll = 0; ll < 4; ll++) {                                    \
                double coeff = in[offsets[0][base0 + ll]];                  \
                coeff *= weights[0][base0 + ll];                            \
                t += coeff;                                                 \
            }                                                              \
            out[ii] = (_out_type)t;                                         \
        }                                                                  \
    } else if (rank == 2) {                                                 \
        npy_intp oo = 0;                                                    \
        for (ii = 0; ii < odimensions[0]; ii++) {                           \
            const npy_intp base0 = 4 * ii;                                  \
            for (jj = 0; jj < odimensions[1]; jj++, oo++) {                 \
                const npy_intp base1 = 4 * jj;                              \
                double t = 0.0;                                             \
                for (ll = 0; ll < 4; ll++) {                                \
                    const npy_intp off0 = offsets[0][base0 + ll];           \
                    const double w0 = weights[0][base0 + ll];               \
                    double row = 0.0;                                       \
                    for (mm = 0; mm < 4; mm++) {                            \
                        double coeff = in[off0 + offsets[1][base1 + mm]];   \
                        coeff *= weights[1][base1 + mm];                    \
                        row += coeff;                                       \
                    }                                                      \
                    t += row * w0;                                          \
                }                                                          \
                out[oo] = (_out_type)t;                                     \
            }                                                              \
        }                                                                  \
    } else {                                                               \
        npy_intp oo = 0;                                                    \
        for (ii = 0; ii < odimensions[0]; ii++) {                           \
            const npy_intp base0 = 4 * ii;                                  \
            for (jj = 0; jj < odimensions[1]; jj++) {                       \
                const npy_intp base1 = 4 * jj;                              \
                for (kk = 0; kk < odimensions[2]; kk++, oo++) {             \
                    const npy_intp base2 = 4 * kk;                          \
                    double t = 0.0;                                         \
                    for (ll = 0; ll < 4; ll++) {                            \
                        const npy_intp off0 = offsets[0][base0 + ll];       \
                        const double w0 = weights[0][base0 + ll];           \
                        double plane = 0.0;                                 \
                        for (mm = 0; mm < 4; mm++) {                        \
                            const npy_intp off1 = offsets[1][base1 + mm];   \
                            const double w1 = weights[1][base1 + mm];       \
                            double row = 0.0;                               \
                            for (nn = 0; nn < 4; nn++) {                    \
                                double coeff = in[off0 + off1 +             \
                                                  offsets[2][base2 + nn]];  \
                                coeff *= weights[2][base2 + nn];            \
                                row += coeff;                               \
                            }                                              \
                            plane += row * w1;                              \
                        }                                                  \
                        t += plane * w0;                                    \
                    }                                                      \
                    out[oo] = (_out_type)t;                                 \
                }                                                          \
            }                                                              \
        }                                                                  \
    }                                                                      \
} while (0)

static void
NI_FreeZoomLinearFastData(npy_intp **offset0, npy_intp **offset1,
                          double **weight0, double **weight1, int rank)
{
    int ii;
    for (ii = 0; ii < rank; ii++) {
        free(offset0[ii]);
        free(offset1[ii]);
        free(weight0[ii]);
        free(weight1[ii]);
    }
}

static int
NI_InitZoomLinearFastAxis(npy_intp odim, npy_intp idim, npy_intp stride,
                          double zoom, int mode, npy_intp *offset0,
                          npy_intp *offset1, double *weight0,
                          double *weight1)
{
    const int spline_mode = _get_spline_boundary_mode(mode);
    npy_intp kk;

    for (kk = 0; kk < odim; kk++) {
        double cc = (double)kk * zoom;
        double cc_floor;
        npy_intp start, idx0, idx1;

        cc = map_coordinate(cc, idim, mode);
        cc_floor = floor(cc);
        start = (npy_intp)cc_floor;
        idx0 = start;
        idx1 = start + 1;

        if (start < 0 || start + 1 >= idim) {
            idx0 = (npy_intp)map_coordinate(idx0, idim, spline_mode);
            idx1 = (npy_intp)map_coordinate(idx1, idim, spline_mode);
        }

        offset0[kk] = stride * idx0;
        offset1[kk] = stride * idx1;
        weight1[kk] = cc - cc_floor;
        weight0[kk] = 1.0 - weight1[kk];
    }

    return 1;
}

static int
NI_ZoomShiftLinearFastPath(PyArrayObject *input, PyArrayObject *zoom_ar,
                           PyArrayObject *shift_ar, PyArrayObject *output,
                           int order, int mode, int nprepad,
                           int grid_mode, int *handled)
{
    npy_intp *offset0[NPY_MAXDIMS] = {NULL};
    npy_intp *offset1[NPY_MAXDIMS] = {NULL};
    double *weight0[NPY_MAXDIMS] = {NULL};
    double *weight1[NPY_MAXDIMS] = {NULL};
    npy_intp idimensions[NPY_MAXDIMS], odimensions[NPY_MAXDIMS];
    npy_intp istrides[NPY_MAXDIMS];
    npy_double *zooms = NULL;
    int rank, type_num, jj;
    NPY_BEGIN_THREADS_DEF;

    *handled = 0;

    if (order != 1 || zoom_ar == NULL || shift_ar != NULL || nprepad != 0 ||
            grid_mode) {
        return 1;
    }
    if (mode != NI_EXTEND_MIRROR && mode != NI_EXTEND_REFLECT) {
        return 1;
    }

    rank = PyArray_NDIM(input);
    if (rank < 1 || rank > 3 || PyArray_NDIM(output) != rank) {
        return 1;
    }

    type_num = PyArray_TYPE(input);
    if (type_num != PyArray_TYPE(output) ||
            (type_num != NPY_FLOAT && type_num != NPY_DOUBLE)) {
        return 1;
    }
    if (!PyArray_IS_C_CONTIGUOUS(input) ||
            !PyArray_IS_C_CONTIGUOUS(output)) {
        return 1;
    }
    if (PyArray_TYPE(zoom_ar) != NPY_DOUBLE ||
            PyArray_SIZE(zoom_ar) != rank) {
        return 1;
    }

    if (PyArray_SIZE(output) == 0) {
        *handled = 1;
        return 1;
    }

    for (jj = 0; jj < rank; jj++) {
        idimensions[jj] = PyArray_DIM(input, jj);
        odimensions[jj] = PyArray_DIM(output, jj);
        istrides[jj] = PyArray_STRIDE(input, jj) / PyArray_ITEMSIZE(input);
        if (idimensions[jj] < 1) {
            return 1;
        }
    }

    zooms = (npy_double *)PyArray_DATA(zoom_ar);
    for (jj = 0; jj < rank; jj++) {
        if (odimensions[jj] == 0) {
            continue;
        }
        offset0[jj] = malloc(odimensions[jj] * sizeof(npy_intp));
        offset1[jj] = malloc(odimensions[jj] * sizeof(npy_intp));
        weight0[jj] = malloc(odimensions[jj] * sizeof(double));
        weight1[jj] = malloc(odimensions[jj] * sizeof(double));
        if (NPY_UNLIKELY(!offset0[jj] || !offset1[jj] ||
                !weight0[jj] || !weight1[jj])) {
            NI_FreeZoomLinearFastData(offset0, offset1, weight0, weight1,
                                      rank);
            PyErr_NoMemory();
            return 0;
        }
        NI_InitZoomLinearFastAxis(odimensions[jj], idimensions[jj],
                                  istrides[jj], zooms[jj], mode,
                                  offset0[jj], offset1[jj],
                                  weight0[jj], weight1[jj]);
    }

    NPY_BEGIN_THREADS;
    switch (type_num) {
        CASE_ZOOM_LINEAR_FAST(NPY_FLOAT, npy_float);
        CASE_ZOOM_LINEAR_FAST(NPY_DOUBLE, npy_double);
    default:
        break;
    }
    NPY_END_THREADS;

    NI_FreeZoomLinearFastData(offset0, offset1, weight0, weight1, rank);
    *handled = 1;
    return 1;
}

static void
NI_FreeZoomCubicFastData(npy_intp **offsets, double **weights, int rank)
{
    int ii;
    for (ii = 0; ii < rank; ii++) {
        free(offsets[ii]);
        free(weights[ii]);
    }
}

static int
NI_InitZoomCubicFastAxis(npy_intp odim, npy_intp idim, npy_intp stride,
                         double zoom, int mode, npy_intp *offsets,
                         double *weights)
{
    const int spline_mode = _get_spline_boundary_mode(mode);
    npy_intp kk, hh;

    for (kk = 0; kk < odim; kk++) {
        double cc = (double)kk * zoom;
        npy_intp start;
        npy_intp base = 4 * kk;

        cc = map_coordinate(cc, idim, mode);
        start = (npy_intp)floor(cc) - 1;
        get_spline_interpolation_weights(cc, 3, weights + base);

        if (start < 0 || start + 3 >= idim) {
            for (hh = 0; hh < 4; hh++) {
                npy_intp idx = start + hh;
                idx = (npy_intp)map_coordinate(idx, idim, spline_mode);
                offsets[base + hh] = stride * idx;
            }
        } else {
            for (hh = 0; hh < 4; hh++) {
                offsets[base + hh] = stride * (start + hh);
            }
        }
    }

    return 1;
}

static int
NI_ZoomShiftCubicFastPath(PyArrayObject *input, PyArrayObject *zoom_ar,
                          PyArrayObject *shift_ar, PyArrayObject *output,
                          int order, int mode, int nprepad,
                          int grid_mode, int *handled)
{
    npy_intp *offsets[NPY_MAXDIMS] = {NULL};
    double *weights[NPY_MAXDIMS] = {NULL};
    npy_intp idimensions[NPY_MAXDIMS], odimensions[NPY_MAXDIMS];
    npy_intp istrides[NPY_MAXDIMS];
    npy_double *zooms = NULL;
    int rank, input_type, output_type, jj;
    NPY_BEGIN_THREADS_DEF;

    *handled = 0;

    if (order != 3 || zoom_ar == NULL || shift_ar != NULL || nprepad != 0 ||
            grid_mode) {
        return 1;
    }
    if (mode != NI_EXTEND_MIRROR && mode != NI_EXTEND_REFLECT) {
        return 1;
    }

    rank = PyArray_NDIM(input);
    if (rank < 1 || rank > 3 || PyArray_NDIM(output) != rank) {
        return 1;
    }

    input_type = PyArray_TYPE(input);
    output_type = PyArray_TYPE(output);
    if ((input_type != NPY_FLOAT && input_type != NPY_DOUBLE) ||
            (output_type != NPY_FLOAT && output_type != NPY_DOUBLE)) {
        return 1;
    }
    if (!PyArray_IS_C_CONTIGUOUS(input) ||
            !PyArray_IS_C_CONTIGUOUS(output)) {
        return 1;
    }
    if (PyArray_TYPE(zoom_ar) != NPY_DOUBLE ||
            PyArray_SIZE(zoom_ar) != rank) {
        return 1;
    }

    if (PyArray_SIZE(output) == 0) {
        *handled = 1;
        return 1;
    }

    for (jj = 0; jj < rank; jj++) {
        idimensions[jj] = PyArray_DIM(input, jj);
        odimensions[jj] = PyArray_DIM(output, jj);
        istrides[jj] = PyArray_STRIDE(input, jj) / PyArray_ITEMSIZE(input);
        if (idimensions[jj] < 1) {
            return 1;
        }
    }

    zooms = (npy_double *)PyArray_DATA(zoom_ar);
    for (jj = 0; jj < rank; jj++) {
        if (odimensions[jj] == 0) {
            continue;
        }
        offsets[jj] = malloc(4 * odimensions[jj] * sizeof(npy_intp));
        weights[jj] = malloc(4 * odimensions[jj] * sizeof(double));
        if (NPY_UNLIKELY(!offsets[jj] || !weights[jj])) {
            NI_FreeZoomCubicFastData(offsets, weights, rank);
            PyErr_NoMemory();
            return 0;
        }
        NI_InitZoomCubicFastAxis(odimensions[jj], idimensions[jj],
                                 istrides[jj], zooms[jj], mode,
                                 offsets[jj], weights[jj]);
    }

    NPY_BEGIN_THREADS;
    if (input_type == NPY_DOUBLE && output_type == NPY_DOUBLE) {
        ZOOM_CUBIC_FAST(npy_double, npy_double);
    } else if (input_type == NPY_DOUBLE && output_type == NPY_FLOAT) {
        ZOOM_CUBIC_FAST(npy_double, npy_float);
    } else if (input_type == NPY_FLOAT && output_type == NPY_FLOAT) {
        ZOOM_CUBIC_FAST(npy_float, npy_float);
    } else if (input_type == NPY_FLOAT && output_type == NPY_DOUBLE) {
        ZOOM_CUBIC_FAST(npy_float, npy_double);
    }
    NPY_END_THREADS;

    NI_FreeZoomCubicFastData(offsets, weights, rank);
    *handled = 1;
    return 1;
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
    double **splvals = NULL, icoor[NPY_MAXDIMS] = {0}, tmp;
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
    if (NPY_UNLIKELY(!edge_offsets)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }
    data_offsets = malloc(irank * sizeof(npy_intp*));
    if (NPY_UNLIKELY(!data_offsets)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }
    for(jj = 0; jj < irank; jj++)
        data_offsets[jj] = NULL;

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
            /* mapping is from a coordinates array: */
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
        } else {
            NPY_END_THREADS;
            PyErr_SetString(PyExc_RuntimeError,
                            "One of `map`, `matrix` or `coordinates` must be provided");
            goto exit;
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
    int rank = 0, fast_path_handled = 0;
    NPY_BEGIN_THREADS_DEF;

    if (!NI_ZoomShiftLinearFastPath(input, zoom_ar, shift_ar, output,
                                    order, mode, nprepad, grid_mode,
                                    &fast_path_handled)) {
        return 0;
    }
    if (fast_path_handled) {
        return 1;
    }
    if (!NI_ZoomShiftCubicFastPath(input, zoom_ar, shift_ar, output,
                                   order, mode, nprepad, grid_mode,
                                   &fast_path_handled)) {
        return 0;
    }
    if (fast_path_handled) {
        return 1;
    }

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
    if (NPY_UNLIKELY(!offsets)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }
    for(jj = 0; jj < rank; jj++)
        offsets[jj] = NULL;

    /* store spline coefficients, along each axis: */
    splvals = malloc(rank * sizeof(double**));
    if (NPY_UNLIKELY(!splvals)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }
    for(jj = 0; jj < rank; jj++)
        splvals[jj] = NULL;

    /* store offsets at all edges: */
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
