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

/* initialize iterations over single array elements: */
int NI_InitPointIterator(PyArrayObject *array, NI_Iterator *iterator)
{
    int ii;

    iterator->rank_m1 = array->nd - 1;
    for(ii = 0; ii < array->nd; ii++) {
        /* adapt dimensions for use in the macros: */
        iterator->dimensions[ii] = array->dimensions[ii] - 1;
        /* initialize coordinates: */
        iterator->coordinates[ii] = 0;
        /* initialize strides: */
        iterator->strides[ii] = array->strides[ii];
        /* calculate the strides to move back at the end of an axis: */
        iterator->backstrides[ii] =
                array->strides[ii] * iterator->dimensions[ii];
    }
    return 1;
}


/* initialize iteration over a lower sub-space: */
int NI_SubspaceIterator(NI_Iterator *iterator, UInt32 axes)
{
    int ii, last = 0;

    for(ii = 0; ii <= iterator->rank_m1; ii++) {
        if (axes & (((UInt32)1) << ii)) {
            if (last != ii) {
                iterator->dimensions[last] = iterator->dimensions[ii];
                iterator->strides[last] = iterator->strides[ii];
                iterator->backstrides[last] = iterator->backstrides[ii];
            }
            ++last;
        }
    }
    iterator->rank_m1 = last - 1;
    return 1;
}

/* initialize iteration over array lines: */
int NI_LineIterator(NI_Iterator *iterator, int axis)
{
    UInt32 axes = ((UInt32)1) << axis;
    return NI_SubspaceIterator(iterator, ~axes);
}


/******************************************************************/
/* Line buffers */
/******************************************************************/

/* Allocate line buffer data */
int NI_AllocateLineBuffer(PyArrayObject* array, int axis, npy_intp size1,
        npy_intp size2, npy_intp *lines, npy_intp max_size, double **buffer)
{
    npy_intp line_size, max_lines;
    int ii;

    /* the number of lines of the array is an upper limit for the
         number of lines in the buffer: */
    max_lines = 1;
    for(ii = 0; ii < array->nd; ii++)
        max_lines *= array->dimensions[ii];
    if (array->nd > 0 && array->dimensions[axis] > 0)
        max_lines /= array->dimensions[axis];
    /* calculate the space needed for one line, including space to
         support the boundary conditions: */
    line_size = sizeof(double) * (array->dimensions[axis] + size1 + size2);
    /* if *lines < 1, no number of lines is proposed, so we calculate it
         from the maximum size allowed: */
    if (*lines < 1) {
        *lines = line_size > 0 ? max_size / line_size : 0;
        if (*lines < 1)
            *lines = 1;
    }
    /* no need to allocate too many lines: */
    if (*lines > max_lines)
        *lines = max_lines;
    /* allocate data for the buffer: */
    *buffer = malloc(*lines * line_size);
    if (!*buffer) {
        PyErr_NoMemory();
        return 0;
    }
    return 1;
}

/* Some NumPy types are ambiguous */
int NI_CanonicalType(int type_num)
{
    switch (type_num) {
        case NPY_INT:
            return NPY_INT32;

        case NPY_LONG:
#if NPY_SIZEOF_LONG == 4
            return NPY_INT32;
#else
            return NPY_INT64;
#endif

        case NPY_LONGLONG:
            return NPY_INT64;

        case NPY_UINT:
            return NPY_UINT32;

        case NPY_ULONG:
#if NPY_SIZEOF_LONG == 4
            return NPY_UINT32;
#else
            return NPY_UINT64;
#endif

        case NPY_ULONGLONG:
            return NPY_UINT64;

        default:
            return type_num;
    }
}

/* Initialize a line buffer */
int NI_InitLineBuffer(PyArrayObject *array, int axis, npy_intp size1,
        npy_intp size2, npy_intp buffer_lines, double *buffer_data,
        NI_ExtendMode extend_mode, double extend_value, NI_LineBuffer *buffer)
{
    npy_intp line_length = 0, array_lines = 0, size;
    int ii;

    size = 1;
    for(ii = 0; ii < array->nd; ii++)
        size *= array->dimensions[ii];
    /* check if the buffer is big enough: */
    if (size > 0 && buffer_lines < 1) {
        PyErr_SetString(PyExc_RuntimeError, "buffer too small");
        return 0;
    }
    /* Initialize a line iterator to move over the array: */
    if (!NI_InitPointIterator(array, &(buffer->iterator)))
        return 0;
    if (!NI_LineIterator(&(buffer->iterator), axis))
        return 0;
    line_length = array->nd > 0 ? array->dimensions[axis] : 1;
    if (line_length > 0)
        array_lines = line_length > 0 ? size / line_length : 1;
    /* initialize the buffer structure: */
    buffer->array_data = (void *)PyArray_DATA(array);
    buffer->buffer_data = buffer_data;
    buffer->buffer_lines = buffer_lines;
    buffer->array_type = NI_CanonicalType(PyArray_DESCR(array)->type_num);
    buffer->array_lines = array_lines;
    buffer->next_line = 0;
    buffer->size1 = size1;
    buffer->size2 = size2;
    buffer->line_length = line_length;
    buffer->line_stride = array->nd > 0 ? array->strides[axis] : 0;
    buffer->extend_mode = extend_mode;
    buffer->extend_value = extend_value;
    return 1;
}

/* Extend a line in memory to implement boundary conditions: */
int NI_ExtendLine(double *line, npy_intp length, npy_intp size1,
                  npy_intp size2, NI_ExtendMode mode, double constant_value, char * errmsg)
{
    npy_intp ii, jj, length1, nextend, rextend;
    double *l1, *l2, *l3, val;

    switch (mode) {
    case NI_EXTEND_WRAP:
        /* deal with situation where data is shorter than needed
           for filling the line */
        nextend = size1 / length;
        rextend = size1 - nextend * length;
        l1 = line + size1 + length - rextend;
        l2 = line;
        for(ii = 0; ii < rextend; ii++)
            *l2++ = *l1++;
        for(ii = 0; ii < nextend; ii++) {
            l1 = line + size1;
            for(jj = 0; jj < length; jj++)
                *l2++ = *l1++;
        }
        nextend = size2 / length;
        rextend = size2 - nextend * length;
        l1 = line + size1;
        l2 = line + size1 + length;
        for(ii = 0; ii < nextend; ii++) {
            l3 = l1;
            for(jj = 0; jj < length; jj++)
                *l2++ = *l3++;
        }
        for(ii = 0; ii < rextend; ii++)
            *l2++ = *l1++;
        break;
    case NI_EXTEND_MIRROR:
        if (length == 1) {
            l1 = line;
            val = line[size1];
            for(ii = 0; ii < size1; ii++)
                *l1++ = val;
            l1 = line + size1 + length;
            val = line[size1 + length - 1];
            for(ii = 0; ii < size2; ii++)
                *l1++ = val;
        } else {
            length1 = length - 1;
            nextend = size1 / length1;
            rextend = size1 - nextend * length1;
            l1 = line + size1 + 1;
            l2 = l1 - 2;
            for(ii = 0; ii < nextend; ii++) {
                l3 = l1;
                for(jj = 0; jj < length1; jj++)
                    *l2-- = *l3++;
                l1 -= length1;
            }
            for(ii = 0; ii < rextend; ii++)
                *l2-- = *l1++;
            nextend = size2 / length1;
            rextend = size2 - nextend * length1;
            l1 = line + size1 + length1 - 1;
            l2 = l1 + 2;
            for(ii = 0; ii < nextend; ii++) {
                l3 = l1;
                for(jj = 0; jj < length1; jj++)
                    *l2++ = *l3--;
                l1 += length1;
            }
            for(ii = 0; ii < rextend; ii++)
                *l2++ = *l1--;
        }
        break;
    case NI_EXTEND_REFLECT:
        nextend = size1 / length;
        rextend = size1 - nextend * length;
        l1 = line + size1;
        l2 = l1 - 1;
        for(ii = 0; ii < nextend; ii++) {
            l3 = l1;
            for(jj = 0; jj < length; jj++)
                *l2-- = *l3++;
            l1 -= length;
        }
        l3 = l1;
        for(ii = 0; ii < rextend; ii++)
            *l2-- = *l3++;
        nextend = size2 / length;
        rextend = size2 - nextend * length;
        l1 = line + size1 + length - 1;
        l2 = l1 + 1;
        for(ii = 0; ii < nextend; ii++) {
            l3 = l1;
            for(jj = 0; jj < length; jj++)
                *l2++ = *l3--;
            l1 += length;
        }
        for(ii = 0; ii < rextend; ii++)
            *l2++ = *l1--;
        break;
    case NI_EXTEND_NEAREST:
        l1 = line;
        val = line[size1];
        for(ii = 0; ii < size1; ii++)
            *l1++ = val;
        l1 = line + size1 + length;
        val = line[size1 + length - 1];
        for(ii = 0; ii < size2; ii++)
            *l1++ = val;
        break;
    case NI_EXTEND_CONSTANT:
        l1 = line;
        for(ii = 0; ii < size1; ii++)
            *l1++ = constant_value;
        l1 = line + size1 + length;
        for(ii = 0; ii < size2; ii++)
            *l1++ = constant_value;
        break;
    default:
        PyOS_snprintf(errmsg, NI_MAX_ERR_MSG, "mode not supported");
        return 0;
    }
    return 1;
}


#define CASE_COPY_DATA_TO_LINE(_pi, _po, _length, _stride, _type) \
case t ## _type:                                                  \
{                                                                 \
    npy_intp _ii;                                                  \
    for(_ii = 0; _ii < _length; _ii++) {                            \
        _po[_ii] = (double)*(_type*)_pi;                              \
        _pi += _stride;                                               \
    }                                                               \
}                                                                 \
break


/* Copy a line from an array to a buffer: */
int NI_ArrayToLineBuffer(NI_LineBuffer *buffer,
                         npy_intp *number_of_lines, int *more,
                         char * errmsg)
{
    double *pb = buffer->buffer_data;
    char *pa;
    npy_intp length = buffer->line_length;

    pb += buffer->size1;
    *number_of_lines = 0;
    /* fill until all lines in the array have been processed, or until
         the buffer is full: */
    while (buffer->next_line < buffer->array_lines &&
                 *number_of_lines < buffer->buffer_lines) {
        pa = buffer->array_data;
        /* copy the data from the array to the buffer: */
        switch (buffer->array_type) {
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, Bool);
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, UInt8);
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, UInt16);
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, UInt32);
#if HAS_UINT64
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, UInt64);
#endif
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, Int8);
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, Int16);
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, Int32);
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, Int64);
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, Float32);
            CASE_COPY_DATA_TO_LINE(pa, pb, length, buffer->line_stride, Float64);
        default:
            PyOS_snprintf(errmsg, NI_MAX_ERR_MSG, "array type %d not supported",
                          buffer->array_type);
            return 0;
        }
        /* goto next line in the array: */
        NI_ITERATOR_NEXT(buffer->iterator, buffer->array_data);
        /* implement boundary conditions to the line: */
        if (buffer->size1 + buffer->size2 > 0)
            if (!NI_ExtendLine(pb - buffer->size1, length, buffer->size1,
                                                 buffer->size2, buffer->extend_mode,
                                                 buffer->extend_value,
                                                 errmsg))
                return 0;
        /* The number of the array lines copied: */
        ++(buffer->next_line);
        /* keep track of (and return) the number of lines in the buffer: */
        ++(*number_of_lines);
        pb += buffer->line_length + buffer->size1 + buffer->size2;
    }
    /* if not all array lines were processed, *more is set true: */
    *more = buffer->next_line < buffer->array_lines;
    return 1;
}

#define CASE_COPY_LINE_TO_DATA(_pi, _po, _length, _stride, _type) \
case t ## _type:                                                  \
{                                                                 \
    npy_intp _ii;                                                  \
    for(_ii = 0; _ii < _length; _ii++) {                            \
        *(_type*)_po = (_type)_pi[_ii];                               \
        _po += _stride;                                               \
    }                                                               \
}                                                                 \
break

/* Copy a line from a buffer to an array: */
int NI_LineBufferToArray(NI_LineBuffer *buffer, char * errmsg)
{
    double *pb = buffer->buffer_data;
    char *pa;
    npy_intp jj, length = buffer->line_length;

    pb += buffer->size1;
    for(jj = 0; jj < buffer->buffer_lines; jj++) {
        /* if all array lines are copied return: */
        if (buffer->next_line == buffer->array_lines)
            break;
        pa = buffer->array_data;
        /* copy data from the buffer to the array: */
        switch (buffer->array_type) {
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, Bool);
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, UInt8);
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, UInt16);
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, UInt32);
#if HAS_UINT64
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, UInt64);
#endif
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, Int8);
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, Int16);
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, Int32);
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, Int64);
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, Float32);
            CASE_COPY_LINE_TO_DATA(pb, pa, length, buffer->line_stride, Float64);
        default:
            PyOS_snprintf(errmsg, NI_MAX_ERR_MSG, "array type not supported");
            return 0;
        }
        /* move to the next line in the array: */
        NI_ITERATOR_NEXT(buffer->iterator, buffer->array_data);
        /* number of lines copied: */
        ++(buffer->next_line);
        /* move the buffer data pointer to the next line: */
        pb += buffer->line_length + buffer->size1 + buffer->size2;
    }
    return 1;
}

/******************************************************************/
/* Multi-dimensional filter support functions */
/******************************************************************/

/* Initialize a filter iterator: */
int
NI_InitFilterIterator(int rank, npy_intp *filter_shape,
                    npy_intp filter_size, npy_intp *array_shape,
                    npy_intp *origins, NI_FilterIterator *iterator)
{
    int ii;
    npy_intp fshape[MAXDIM], forigins[MAXDIM];

    for(ii = 0; ii < rank; ii++) {
        fshape[ii] = *filter_shape++;
        forigins[ii] = origins ? *origins++ : 0;
    }
    /* calculate the strides, used to move the offsets pointer through
         the offsets table: */
    if (rank > 0) {
        iterator->strides[rank - 1] = filter_size;
        for(ii = rank - 2; ii >= 0; ii--) {
            npy_intp step = array_shape[ii + 1] < fshape[ii + 1] ?
                                                                         array_shape[ii + 1] : fshape[ii + 1];
            iterator->strides[ii] =  iterator->strides[ii + 1] * step;
        }
    }
    for(ii = 0; ii < rank; ii++) {
        npy_intp step = array_shape[ii] < fshape[ii] ?
                                                                                         array_shape[ii] : fshape[ii];
        npy_intp orgn = fshape[ii] / 2 + forigins[ii];
        /* stride for stepping back to previous offsets: */
        iterator->backstrides[ii] = (step - 1) * iterator->strides[ii];
        /* initialize boundary extension sizes: */
        iterator->bound1[ii] = orgn;
        iterator->bound2[ii] = array_shape[ii] - fshape[ii] + orgn;
    }
    return 1;
}

/* Calculate the offsets to the filter points, for all border regions and
     the interior of the array: */
int NI_InitFilterOffsets(PyArrayObject *array, Bool *footprint,
         npy_intp *filter_shape, npy_intp* origins,
         NI_ExtendMode mode, npy_intp **offsets, npy_intp *border_flag_value,
         npy_intp **coordinate_offsets)
{
    int rank, ii;
    npy_intp kk, ll, filter_size = 1, offsets_size = 1, max_size = 0;
    npy_intp max_stride = 0, *ashape = NULL, *astrides = NULL;
    npy_intp footprint_size = 0, coordinates[MAXDIM], position[MAXDIM];
    npy_intp fshape[MAXDIM], forigins[MAXDIM], *po, *pc = NULL;

    rank = array->nd;
    ashape = array->dimensions;
    astrides = array->strides;
    for(ii = 0; ii < rank; ii++) {
        fshape[ii] = *filter_shape++;
        forigins[ii] = origins ? *origins++ : 0;
    }
    /* the size of the footprint array: */
    for(ii = 0; ii < rank; ii++)
        filter_size *= fshape[ii];
    /* calculate the number of non-zero elements in the footprint: */
    if (footprint) {
        for(kk = 0; kk < filter_size; kk++)
            if (footprint[kk])
                ++footprint_size;
    } else {
        footprint_size = filter_size;
    }
    /* calculate how many sets of offsets must be stored: */
    for(ii = 0; ii < rank; ii++)
        offsets_size *= (ashape[ii] < fshape[ii] ? ashape[ii] : fshape[ii]);
    /* allocate offsets data: */
    *offsets = malloc(offsets_size * footprint_size * sizeof(npy_intp));
    if (!*offsets) {
        PyErr_NoMemory();
        goto exit;
    }
    if (coordinate_offsets) {
        *coordinate_offsets = malloc(offsets_size * rank
                                     * footprint_size * sizeof(npy_intp));
        if (!*coordinate_offsets) {
            PyErr_NoMemory();
            goto exit;
        }
    }
    for(ii = 0; ii < rank; ii++) {
        npy_intp stride;
        /* find maximum axis size: */
        if (ashape[ii] > max_size)
            max_size = ashape[ii];
        /* find maximum stride: */
        stride = astrides[ii] < 0 ? -astrides[ii] : astrides[ii];
        if (stride > max_stride)
            max_stride = stride;
        /* coordinates for iterating over the kernel elements: */
        coordinates[ii] = 0;
        /* keep track of the kernel position: */
        position[ii] = 0;
    }
    /* the flag to indicate that we are outside the border must have a
         value that is larger than any possible offset: */
    *border_flag_value = max_size * max_stride + 1;
    /* calculate all possible offsets to elements in the filter kernel,
         for all regions in the array (interior and border regions): */
    po = *offsets;
    if (coordinate_offsets) {
        pc = *coordinate_offsets;
    }
    /* iterate over all regions: */
    for(ll = 0; ll < offsets_size; ll++) {
        /* iterate over the elements in the footprint array: */
        for(kk = 0; kk < filter_size; kk++) {
            npy_intp offset = 0;
            /* only calculate an offset if the footprint is 1: */
            if (!footprint || footprint[kk]) {
                /* find offsets along all axes: */
                for(ii = 0; ii < rank; ii++) {
                    npy_intp orgn = fshape[ii] / 2 + forigins[ii];
                    npy_intp cc = coordinates[ii] - orgn + position[ii];
                    npy_intp len = ashape[ii];
                    /* apply boundary conditions, if necessary: */
                    switch (mode) {
                    case NI_EXTEND_MIRROR:
                        if (cc < 0) {
                            if (len <= 1) {
                                cc = 0;
                            } else {
                                int sz2 = 2 * len - 2;
                                cc = sz2 * (int)(-cc / sz2) + cc;
                                cc = cc <= 1 - len ? cc + sz2 : -cc;
                            }
                        } else if (cc >= len) {
                            if (len <= 1) {
                                cc = 0;
                            } else {
                                int sz2 = 2 * len - 2;
                                cc -= sz2 * (int)(cc / sz2);
                                if (cc >= len)
                                    cc = sz2 - cc;
                            }
                        }
                        break;
                    case NI_EXTEND_REFLECT:
                        if (cc < 0) {
                            if (len <= 1) {
                                cc = 0;
                            } else {
                                int sz2 = 2 * len;
                                if (cc < -sz2)
                                    cc = sz2 * (int)(-cc / sz2) + cc;
                                cc = cc < -len ? cc + sz2 : -cc - 1;
                            }
                        } else if (cc >= len) {
                            if (len <= 1) {cc = 0;
                            } else {
                                int sz2 = 2 * len;
                                cc -= sz2 * (int)(cc / sz2);
                                if (cc >= len)
                                    cc = sz2 - cc - 1;
                            }
                        }
                        break;
                    case NI_EXTEND_WRAP:
                        if (cc < 0) {
                            if (len <= 1) {
                                cc = 0;
                            } else {
                                int sz = len;
                                cc += sz * (int)(-cc / sz);
                                if (cc < 0)
                                    cc += sz;
                            }
                        } else if (cc >= len) {
                            if (len <= 1) {
                                cc = 0;
                            } else {
                                int sz = len;
                                cc -= sz * (int)(cc / sz);
                            }
                        }
                        break;
                    case NI_EXTEND_NEAREST:
                        if (cc < 0) {
                            cc = 0;
                        } else if (cc >= len) {
                            cc = len - 1;
                        }
                        break;
                    case NI_EXTEND_CONSTANT:
                        if (cc < 0 || cc >= len)
                            cc = *border_flag_value;
                        break;
                    default:
                    PyErr_SetString(PyExc_RuntimeError,
                                                                                    "boundary mode not supported");
                        goto exit;
                    }

                    /* calculate offset along current axis: */
                    if (cc == *border_flag_value) {
                        /* just flag that we are outside the border */
                        offset = *border_flag_value;
                        if (coordinate_offsets)
                            pc[ii] = 0;
                        break;
                    } else {
                        /* use an offset that is possibly mapped from outside
                           the border: */
                        cc = cc - position[ii];
                        offset += astrides[ii] * cc;
                        if (coordinate_offsets)
                            pc[ii] = cc;
                    }
                }
                /* store the offset */
                *po++ = offset;
                if (coordinate_offsets)
                    pc += rank;
            }
            /* next point in the filter: */
            for(ii = rank - 1; ii >= 0; ii--) {
                if (coordinates[ii] < fshape[ii] - 1) {
                    coordinates[ii]++;
                    break;
                } else {
                    coordinates[ii] = 0;
                }
            }
        }

        /* move to the next array region: */
        for(ii = rank - 1; ii >= 0; ii--) {
            int orgn = fshape[ii] / 2 + forigins[ii];
            if (position[ii] == orgn) {
                position[ii] += ashape[ii] - fshape[ii] + 1;
                if (position[ii] <= orgn)
                    position[ii] = orgn + 1;
            } else {
                position[ii]++;
            }
            if (position[ii] < ashape[ii]) {
                break;
            } else {
                position[ii] = 0;
            }
        }
    }

 exit:
    if (PyErr_Occurred()) {
        free(*offsets);
        if (coordinate_offsets) {
            free(*coordinate_offsets);
        }
        return 0;
    } else {
        return 1;
    }
}

NI_CoordinateList* NI_InitCoordinateList(int size, int rank)
{
    NI_CoordinateList *list = malloc(sizeof(NI_CoordinateList));
    if (!list) {
        return NULL;
    }
    list->block_size = size;
    list->rank = rank;
    list->blocks = NULL;
    return list;
}

int NI_CoordinateListStealBlocks(NI_CoordinateList *list1,
                                                                 NI_CoordinateList *list2)
{
    if (list1->block_size != list2->block_size ||
            list1->rank != list2->rank) {
        PyErr_SetString(PyExc_RuntimeError, "coordinate lists not compatible");
        return 1;
    }
    if (list1->blocks) {
        PyErr_SetString(PyExc_RuntimeError, "first is list not empty");
        return 1;
    }
    list1->blocks = list2->blocks;
    list2->blocks = NULL;
    return 0;
}

NI_CoordinateBlock* NI_CoordinateListAddBlock(NI_CoordinateList *list)
{
    NI_CoordinateBlock* block = NULL;
    block = malloc(sizeof(NI_CoordinateBlock));
    if (!block) {
        return NULL;
    }
    block->coordinates = malloc(list->block_size * list->rank
                                * sizeof(npy_intp));
    if (!block->coordinates) {
        free(block);
        return NULL;
    }
    block->next = list->blocks;
    list->blocks = block;
    block->size = 0;

    return block;
}

NI_CoordinateBlock* NI_CoordinateListDeleteBlock(NI_CoordinateList *list)
{
    NI_CoordinateBlock* block = list->blocks;
    if (block) {
        list->blocks = block->next;
        free(block->coordinates);
        free(block);
    }
    return list->blocks;
}

void NI_FreeCoordinateList(NI_CoordinateList *list)
{
    if (list) {
        NI_CoordinateBlock *block = list->blocks;
        while (block) {
            NI_CoordinateBlock *tmp = block;
            block = block->next;
            free(tmp->coordinates);
            free(tmp);
        }
        list->blocks = NULL;
        free(list);
    }
}
