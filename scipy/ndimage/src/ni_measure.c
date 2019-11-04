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
#include "ni_measure.h"
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <assert.h>


#define CASE_FIND_OBJECT_POINT(_TYPE, _type, _pi, _regions, _array,  \
                               _max_label, _ii)                      \
case _TYPE:                                                          \
{                                                                    \
    int _kk;                                                         \
    npy_intp _rank = PyArray_NDIM(_array);                           \
    npy_intp _sindex = *(_type *)_pi - 1;                            \
    if (_sindex >= 0 && _sindex < _max_label) {                      \
        if (_rank > 0) {                                             \
            _sindex *= 2 * _rank;                                    \
            if (_regions[_sindex] < 0) {                             \
                for (_kk = 0; _kk < _rank; _kk++) {                  \
                    npy_intp _cc = _ii.coordinates[_kk];             \
                    _regions[_sindex + _kk] = _cc;                   \
                    _regions[_sindex + _kk + _rank] = _cc + 1;       \
                }                                                    \
            }                                                        \
            else {                                                   \
                for(_kk = 0; _kk < _rank; _kk++) {                   \
                    npy_intp _cc = _ii.coordinates[_kk];             \
                    if (_cc < _regions[_sindex + _kk]) {             \
                        _regions[_sindex + _kk] = _cc;               \
                    }                                                \
                    if (_cc + 1 > _regions[_sindex + _kk + _rank]) { \
                        _regions[_sindex + _kk + _rank] = _cc + 1;   \
                    }                                                \
                }                                                    \
            }                                                        \
        }                                                            \
        else {                                                       \
            _regions[_sindex] = 1;                                   \
        }                                                            \
    }                                                                \
}                                                                    \
break

int NI_FindObjects(PyArrayObject* input, npy_intp max_label,
                                     npy_intp* regions)
{
    npy_intp size, jj;
    NI_Iterator ii;
    char *pi;
    NPY_BEGIN_THREADS_DEF;

    NPY_BEGIN_THREADS;

    /* get input data, size and iterator: */
    pi = (void *)PyArray_DATA(input);
    size = PyArray_SIZE(input);
    if (!NI_InitPointIterator(input, &ii))
        goto exit;
    if (PyArray_NDIM(input) > 0) {
        for (jj = 0; jj < 2 * PyArray_NDIM(input) * max_label; jj++) {
            regions[jj] = -1;
        }
    } else {
        for(jj = 0; jj < max_label; jj++)
            regions[jj] = -1;
    }
    /* iterate over all points: */
    for(jj = 0 ; jj < size; jj++) {
        switch (PyArray_TYPE(input)) {
            CASE_FIND_OBJECT_POINT(NPY_BOOL, npy_bool,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_UBYTE, npy_ubyte,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_USHORT, npy_ushort,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_UINT, npy_uint,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_ULONG, npy_ulong,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_ULONGLONG, npy_ulonglong,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_BYTE, npy_byte,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_SHORT, npy_short,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_INT, npy_int,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_LONG, npy_long,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_LONGLONG, npy_longlong,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_FLOAT, npy_float,
                                   pi, regions, input, max_label, ii);
            CASE_FIND_OBJECT_POINT(NPY_DOUBLE, npy_double,
                                   pi, regions, input, max_label, ii);
        default:
            NPY_END_THREADS;
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        NI_ITERATOR_NEXT(ii, pi);
    }
 exit:
    NPY_END_THREADS;
    return PyErr_Occurred() ? 0 : 1;
}


#define WS_GET_INDEX(_TYPE, _type, _index, _c_strides,      \
                     _b_strides, _rank, _out,  _contiguous) \
do {                                                        \
    if (_contiguous) {                                      \
        _out = _index * sizeof(_type);                      \
    }                                                       \
    else {                                                  \
        int _qq;                                            \
        npy_intp _cc;                                       \
        npy_intp _idx = _index;                             \
        _out = 0;                                           \
        for (_qq = 0; _qq < _rank; _qq++) {                 \
            _cc = _idx / _c_strides[_qq];                   \
            _idx -= _cc * _c_strides[_qq];                  \
            _out += _b_strides[_qq] * _cc;                  \
        }                                                   \
    }                                                       \
} while(0)

#define CASE_GET_INPUT(_TYPE, _type, _ival, _pi) \
case _TYPE:                                      \
    _ival = *(_type *)_pi;                     \
    break

#define CASE_GET_LABEL(_TYPE, _type, _label, _pm) \
case _TYPE:                                       \
    _label = *(_type *)_pm;                       \
    break

#define CASE_PUT_LABEL(_TYPE, _type, _label, _pl) \
case _TYPE:                                       \
    *(_type *)_pl = _label;                       \
    break

#define CASE_WINDEX1(_TYPE, _type, _v_index, _p_index, _strides, _istrides, \
                     _irank, _icont, _p_idx, _v_idx, _pi, _vval, _pval)     \
case _TYPE:                                                                 \
    WS_GET_INDEX(_TYPE, _type, _v_index, _strides, _istrides, _irank,       \
                 _p_idx, _icont);                                           \
    WS_GET_INDEX(_TYPE, _type, _p_index, _strides, _istrides, _irank,       \
                 _v_idx, _icont);                                           \
    _vval = *(_type *)(_pi + _v_idx);                                       \
    _pval = *(_type *)(_pi + _p_idx);                                       \
    break

#define CASE_WINDEX2(_TYPE, _type, _v_index, _strides, _ostrides, \
                     _irank, _idx, _ocont, _label, _pl)           \
case _TYPE:                                                       \
    WS_GET_INDEX(_TYPE, _type, _v_index, _strides, _ostrides,     \
                 _irank, _idx, _ocont);                           \
    _label = *(_type *)(_pl + _idx);                              \
    break

#define CASE_WINDEX3(_TYPE, _type, _p_index, _strides, _ostrides, \
                     _irank, _idx, _ocont, _label, _pl)           \
case _TYPE:                                                       \
    WS_GET_INDEX(_TYPE, _type, _p_index, _strides, _ostrides,     \
                 _irank, _idx, _ocont);                           \
    *(_type *)(_pl + _idx) = _label;                              \
break

#define WS_MAXDIM 7

typedef struct {
    npy_intp index;
    void *next, *prev;
    npy_uint16 done;
    npy_uint8 cost;
} NI_WatershedElement;

int NI_WatershedIFT(PyArrayObject* input, PyArrayObject* markers,
                                        PyArrayObject* strct, PyArrayObject* output)
{
    char *pl, *pm, *pi;
    int ll;
    npy_intp size, jj, hh, kk, maxval;
    npy_intp strides[WS_MAXDIM], coordinates[WS_MAXDIM];
    npy_intp *nstrides = NULL, nneigh, ssize;
    int i_contiguous, o_contiguous;
    NI_WatershedElement *temp = NULL, **first = NULL, **last = NULL;
    npy_bool *ps = NULL;
    NI_Iterator mi, ii, li;
    NPY_BEGIN_THREADS_DEF;

    i_contiguous = PyArray_ISCONTIGUOUS(input);
    o_contiguous = PyArray_ISCONTIGUOUS(output);
    ssize = PyArray_SIZE(strct);
    if (PyArray_NDIM(input) > WS_MAXDIM) {
        PyErr_SetString(PyExc_RuntimeError, "too many dimensions");
        goto exit;
    }
    size = PyArray_SIZE(input);
    /* Storage for the temporary queue data. */
    temp = malloc(size * sizeof(NI_WatershedElement));
    if (!temp) {
        PyErr_NoMemory();
        goto exit;
    }

    NPY_BEGIN_THREADS;

    pi = (void *)PyArray_DATA(input);
    if (!NI_InitPointIterator(input, &ii))
        goto exit;
    /* Initialization and find the maximum of the input. */
    maxval = 0;
    for(jj = 0; jj < size; jj++) {
        npy_intp ival = 0;
        switch (PyArray_TYPE(input)) {
            CASE_GET_INPUT(NPY_UINT8, npy_uint8, ival, pi);
            CASE_GET_INPUT(NPY_UINT16, npy_uint16, ival, pi);
        default:
            NPY_END_THREADS;
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        temp[jj].index = jj;
        temp[jj].done = 0;
        if (ival > maxval)
            maxval = ival;
        NI_ITERATOR_NEXT(ii, pi);
    }
    pi = (void *)PyArray_DATA(input);
    /* Allocate and initialize the storage for the queue. */
    first = malloc((maxval + 1) * sizeof(NI_WatershedElement*));
    last = malloc((maxval + 1) * sizeof(NI_WatershedElement*));
    if (NPY_UNLIKELY(!first || !last)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }
    for(hh = 0; hh <= maxval; hh++) {
        first[hh] = NULL;
        last[hh] = NULL;
    }
    if (!NI_InitPointIterator(markers, &mi))
        goto exit;
    if (!NI_InitPointIterator(output, &li))
        goto exit;
    pm = (void *)PyArray_DATA(markers);
    pl = (void *)PyArray_DATA(output);
    /* initialize all nodes */
    for (ll = 0; ll < PyArray_NDIM(input); ll++) {
        coordinates[ll] = 0;
    }
    for(jj = 0; jj < size; jj++) {
        /* get marker */
        int label = 0;
        switch (PyArray_TYPE(markers)) {
            CASE_GET_LABEL(NPY_UBYTE, npy_ubyte, label, pm);
            CASE_GET_LABEL(NPY_USHORT, npy_ushort, label, pm);
            CASE_GET_LABEL(NPY_UINT, npy_uint, label, pm);
            CASE_GET_LABEL(NPY_ULONG, npy_ulong, label, pm);
            CASE_GET_LABEL(NPY_ULONGLONG, npy_ulonglong, label, pm);
            CASE_GET_LABEL(NPY_BYTE, npy_byte, label, pm);
            CASE_GET_LABEL(NPY_SHORT, npy_short, label, pm);
            CASE_GET_LABEL(NPY_INT, npy_int, label, pm);
            CASE_GET_LABEL(NPY_LONG, npy_long, label, pm);
            CASE_GET_LABEL(NPY_LONGLONG, npy_longlong, label, pm);
        default:
            NPY_END_THREADS;
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        switch (PyArray_TYPE(output)) {
            CASE_PUT_LABEL(NPY_UBYTE, npy_ubyte, label, pl);
            CASE_PUT_LABEL(NPY_USHORT, npy_ushort, label, pl);
            CASE_PUT_LABEL(NPY_UINT, npy_uint, label, pl);
            CASE_PUT_LABEL(NPY_ULONG, npy_ulong, label, pl);
            CASE_PUT_LABEL(NPY_ULONGLONG, npy_ulonglong, label, pl);
            CASE_PUT_LABEL(NPY_BYTE, npy_byte, label, pl);
            CASE_PUT_LABEL(NPY_SHORT, npy_short, label, pl);
            CASE_PUT_LABEL(NPY_INT, npy_int, label, pl);
            CASE_PUT_LABEL(NPY_LONG, npy_long, label, pl);
            CASE_PUT_LABEL(NPY_LONGLONG, npy_longlong, label, pl);
        default:
            NPY_END_THREADS;
            PyErr_SetString(PyExc_RuntimeError, "data type not supported");
            goto exit;
        }
        NI_ITERATOR_NEXT2(mi, li, pm, pl);
        if (label != 0) {
            /* This node is a marker */
            temp[jj].cost = 0;
            if (!first[0]) {
                first[0] = &(temp[jj]);
                first[0]->next = NULL;
                first[0]->prev = NULL;
                last[0] = first[0];
            } else {
                if (label > 0) {
                    /* object markers are enqueued at the beginning, so they
                       are processed first. */
                    temp[jj].next = first[0];
                    temp[jj].prev = NULL;
                    first[0]->prev = &(temp[jj]);
                    first[0] = &(temp[jj]);
                } else {
                    /* background markers are enqueued at the end, so they are
                         processed after the object markers. */
                    temp[jj].next = NULL;
                    temp[jj].prev = last[0];
                    last[0]->next = &(temp[jj]);
                    last[0] = &(temp[jj]);
                }
            }
        } else {
            /* This node is not a marker */
            temp[jj].cost = maxval + 1;
            temp[jj].next = NULL;
            temp[jj].prev = NULL;
        }
        for (ll = PyArray_NDIM(input) - 1; ll >= 0; ll--) {
            if (coordinates[ll] < PyArray_DIMS(input)[ll] - 1) {
                coordinates[ll]++;
                break;
            } else {
                coordinates[ll] = 0;
            }
        }
    }

    pl = (void *)PyArray_DATA(output);
    ps = (npy_bool*)PyArray_DATA(strct);
    nneigh = 0;
    for (kk = 0; kk < ssize; kk++)
        if (ps[kk] && kk != (ssize / 2))
            ++nneigh;
    nstrides = malloc(nneigh * sizeof(npy_intp));
    if (NPY_UNLIKELY(!nstrides)) {
        NPY_END_THREADS;
        PyErr_NoMemory();
        goto exit;
    }
    strides[PyArray_NDIM(input) - 1] = 1;
    for (ll = PyArray_NDIM(input) - 2; ll >= 0; ll--) {
        strides[ll] = PyArray_DIM(input, ll + 1) * strides[ll + 1];
    }
    for (ll = 0; ll < PyArray_NDIM(input); ll++) {
        coordinates[ll] = -1;
    }
    for(kk = 0; kk < nneigh; kk++)
        nstrides[kk] = 0;
    jj = 0;
    for(kk = 0; kk < ssize; kk++) {
        if (ps[kk]) {
            int offset = 0;
            for (ll = 0; ll < PyArray_NDIM(input); ll++) {
                offset += coordinates[ll] * strides[ll];
            }
            if (offset != 0)
                nstrides[jj++] += offset;
        }
        for (ll = PyArray_NDIM(input) - 1; ll >= 0; ll--) {
            if (coordinates[ll] < 1) {
                coordinates[ll]++;
                break;
            } else {
                coordinates[ll] = -1;
            }
        }
    }
    /* Propagation phase: */
    for(jj = 0; jj <= maxval; jj++) {
        while (first[jj]) {
            /* dequeue first element: */
            NI_WatershedElement *v = first[jj];
            first[jj] = first[jj]->next;
            if (first[jj])
                first[jj]->prev = NULL;
            v->prev = NULL;
            v->next = NULL;
            /* Mark element as done: */
            v->done = 1;
            /* Iterate over the neighbors of the element: */
            for(hh = 0; hh < nneigh; hh++) {
                npy_intp v_index = v->index, p_index = v->index, idx, cc;
                int qq, outside = 0;
                p_index += nstrides[hh];
                /* check if the neighbor is within the extent of the array: */
                idx = p_index;
                for (qq = 0; qq < PyArray_NDIM(input); qq++) {
                    cc = idx / strides[qq];
                    if (cc < 0 || cc >= PyArray_DIM(input, qq)) {
                        outside = 1;
                        break;
                    }
                    idx -= cc * strides[qq];
                }
                if (!outside) {
                    NI_WatershedElement *p = &(temp[p_index]);
                    if (!(p->done)) {
                        /* If the neighbor was not processed yet: */
                        int max, pval, vval, wvp, pcost, label, p_idx, v_idx;
                        switch (PyArray_TYPE(input)) {
                            CASE_WINDEX1(NPY_UBYTE, npy_ubyte,
                                v_index, p_index, strides,
                                PyArray_STRIDES(input), PyArray_NDIM(input),
                                i_contiguous, p_idx, v_idx, pi, vval, pval);
                            CASE_WINDEX1(NPY_USHORT, npy_ushort,
                                v_index, p_index, strides,
                                PyArray_STRIDES(input), PyArray_NDIM(input),
                                i_contiguous, p_idx, v_idx, pi, vval, pval);
                        default:
                            NPY_END_THREADS;
                            PyErr_SetString(PyExc_RuntimeError,
                                            "data type not supported");
                            goto exit;
                        }
                        /* Calculate cost: */
                        wvp = pval - vval;
                        if (wvp < 0)
                            wvp = -wvp;
                        /* Find the maximum of this cost and the current
                             element cost: */
                        pcost = p->cost;
                        max = v->cost > wvp ? v->cost : wvp;
                        if (max < pcost) {
                            /* If this maximum is less than the neighbors cost,
                                 adapt the cost and the label of the neighbor: */
                            int idx;
                            p->cost = max;
                            switch (PyArray_TYPE(output)) {
                                CASE_WINDEX2(NPY_UBYTE, npy_ubyte,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX2(NPY_USHORT, npy_ushort,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX2(NPY_UINT, npy_uint,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX2(NPY_ULONG, npy_ulong,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX2(NPY_ULONGLONG, npy_ulonglong,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX2(NPY_BYTE, npy_byte,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX2(NPY_SHORT, npy_short,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX2(NPY_INT, npy_int,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX2(NPY_LONG, npy_long,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX2(NPY_LONGLONG, npy_longlong,
                                             v_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                            default:
                                NPY_END_THREADS;
                                PyErr_SetString(PyExc_RuntimeError,
                                                "data type not supported");
                                goto exit;
                            }
                            switch (PyArray_TYPE(output)) {
                                CASE_WINDEX3(NPY_UBYTE, npy_ubyte,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX3(NPY_USHORT, npy_ushort,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX3(NPY_UINT, npy_uint,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX3(NPY_ULONG, npy_ulong,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX3(NPY_ULONGLONG, npy_ulonglong,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX3(NPY_BYTE, npy_byte,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX3(NPY_SHORT, npy_short,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX3(NPY_INT, npy_int,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX3(NPY_LONG, npy_long,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                                CASE_WINDEX3(NPY_LONGLONG, npy_longlong,
                                             p_index, strides,
                                             PyArray_STRIDES(output),
                                             PyArray_NDIM(input),
                                             idx, o_contiguous, label, pl);
                            default:
                                NPY_END_THREADS;
                                PyErr_SetString(PyExc_RuntimeError,
                                                "data type not supported");
                                goto exit;
                            }
                            /* If the neighbor is in a queue, remove it: */
                            if (p->next || p->prev) {
                                NI_WatershedElement *prev = p->prev, *next = p->next;
                                if (first[pcost] == p)
                                    first[pcost] = next;
                                if (last[pcost] == p)
                                    last[pcost] = prev;
                                if (prev)
                                    prev->next = next;
                                if (next)
                                    next->prev = prev;
                            }
                            /* Insert the neighbor in the appropriate queue: */
                            if (label < 0) {
                                p->prev = last[max];
                                p->next = NULL;
                                if (last[max])
                                    last[max]->next = p;
                                last[max] = p;
                                if (!first[max])
                                    first[max] = p;
                            } else {
                                p->next = first[max];
                                p->prev = NULL;
                                if (first[max])
                                    first[max]->prev = p;
                                first[max] = p;
                                if (!last[max])
                                    last[max] = p;
                            }
                        }
                    }
                }
            }
        }
    }
 exit:
    NPY_END_THREADS;
    free(temp);
    free(first);
    free(last);
    free(nstrides);
    return PyErr_Occurred() ? 0 : 1;
}
