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

#ifndef NI_SUPPORT_H
#define NI_SUPPORT_H

#include "nd_image.h"
#include <stdlib.h>
#include <float.h>
#include <limits.h>
#include <assert.h>

/* The different boundary conditions. The mirror condition is not used
     by the python code, but C code is kept around in case we might wish
     to add it. */
typedef enum {
    NI_EXTEND_FIRST = 0,
    NI_EXTEND_NEAREST = 0,
    NI_EXTEND_WRAP = 1,
    NI_EXTEND_REFLECT = 2,
    NI_EXTEND_MIRROR = 3,
    NI_EXTEND_CONSTANT = 4,
    NI_EXTEND_LAST = NI_EXTEND_CONSTANT,
    NI_EXTEND_DEFAULT = NI_EXTEND_MIRROR
} NI_ExtendMode;

/******************************************************************/
/* Iterators */
/******************************************************************/

/******************************************************************/
/* Iterators */
/******************************************************************/

/* the iterator structure: */
typedef struct {
    int rank_m1;
    maybelong dimensions[MAXDIM];
    maybelong coordinates[MAXDIM];
    maybelong strides[MAXDIM];
    maybelong backstrides[MAXDIM];
} NI_Iterator;

/* initialize iterations over single array elements: */
int NI_InitPointIterator(PyArrayObject*, NI_Iterator*);

/* initialize iterations over an arbritrary sub-space: */
int NI_SubspaceIterator(NI_Iterator*, UInt32);

/* initialize iteration over array lines: */
int NI_LineIterator(NI_Iterator*, int);

/* reset an iterator */
#define NI_ITERATOR_RESET(iterator)              \
{                                                \
    int _ii;                                       \
    for(_ii = 0; _ii <= (iterator).rank_m1; _ii++) \
        (iterator).coordinates[_ii] = 0;             \
}

/* go to the next point in a single array */
#define NI_ITERATOR_NEXT(iterator, pointer)                         \
{                                                                   \
    int _ii;                                                          \
    for(_ii = (iterator).rank_m1; _ii >= 0; _ii--)                    \
        if ((iterator).coordinates[_ii] < (iterator).dimensions[_ii]) { \
            (iterator).coordinates[_ii]++;                                \
            pointer += (iterator).strides[_ii];                           \
            break;                                                        \
        } else {                                                        \
            (iterator).coordinates[_ii] = 0;                              \
            pointer -= (iterator).backstrides[_ii];                       \
        }                                                               \
}

/* go to the next point in two arrays of the same size */
#define NI_ITERATOR_NEXT2(iterator1, iterator2,  pointer1, pointer2)  \
{                                                                     \
    int _ii;                                                            \
    for(_ii = (iterator1).rank_m1; _ii >= 0; _ii--)                     \
        if ((iterator1).coordinates[_ii] < (iterator1).dimensions[_ii]) { \
            (iterator1).coordinates[_ii]++;                                 \
            pointer1 += (iterator1).strides[_ii];                           \
            pointer2 += (iterator2).strides[_ii];                           \
            break;                                                          \
        } else {                                                          \
            (iterator1).coordinates[_ii] = 0;                               \
            pointer1 -= (iterator1).backstrides[_ii];                       \
            pointer2 -= (iterator2).backstrides[_ii];                       \
        }                                                                 \
}

/* go to the next point in three arrays of the same size */
#define NI_ITERATOR_NEXT3(iterator1, iterator2,  iterator3,           \
                                                    pointer1, pointer2, pointer3)               \
{                                                                     \
    int _ii;                                                            \
    for(_ii = (iterator1).rank_m1; _ii >= 0; _ii--)                     \
        if ((iterator1).coordinates[_ii] < (iterator1).dimensions[_ii]) { \
            (iterator1).coordinates[_ii]++;                                 \
            pointer1 += (iterator1).strides[_ii];                           \
            pointer2 += (iterator2).strides[_ii];                           \
            pointer3 += (iterator3).strides[_ii];                           \
            break;                                                          \
        } else {                                                          \
            (iterator1).coordinates[_ii] = 0;                               \
            pointer1 -= (iterator1).backstrides[_ii];                       \
            pointer2 -= (iterator2).backstrides[_ii];                       \
            pointer3 -= (iterator3).backstrides[_ii];                       \
        }                                                                 \
}

/* go to an arbitrary point in a single array */
#define NI_ITERATOR_GOTO(iterator, destination, base, pointer) \
{                                                              \
    int _ii;                                                     \
    pointer = base;                                              \
    for(_ii = (iterator).rank_m1; _ii >= 0; _ii--) {             \
        pointer += destination[_ii] * (iterator).strides[_ii];     \
        (iterator).coordinates[_ii] = destination[_ii];            \
    }                                                            \
}

/******************************************************************/
/* Line buffers */
/******************************************************************/

/* the linebuffer structure: */
typedef struct {
    double *buffer_data;
    maybelong buffer_lines, line_length, line_stride;
    maybelong size1, size2, array_lines, next_line;
    NI_Iterator iterator;
    char* array_data;
    NumarrayType array_type;
    NI_ExtendMode extend_mode;
    double extend_value;
} NI_LineBuffer;

/* Get the next line being processed: */
#define NI_GET_LINE(_buffer, _line)                                      \
    ((_buffer).buffer_data + (_line) * ((_buffer).line_length +            \
                                                                            (_buffer).size1 + (_buffer).size2))
/* Allocate line buffer data */
int NI_AllocateLineBuffer(PyArrayObject*, int, maybelong, maybelong,
                                                    maybelong*, maybelong, double**);

/* Initialize a line buffer */
int NI_InitLineBuffer(PyArrayObject*, int, maybelong, maybelong, maybelong,
                                            double*, NI_ExtendMode, double, NI_LineBuffer*);

/* Extend a line in memory to implement boundary conditions: */
int NI_ExtendLine(double*, maybelong, maybelong, maybelong, NI_ExtendMode, double);

/* Copy a line from an array to a buffer: */
int NI_ArrayToLineBuffer(NI_LineBuffer*, maybelong*, int*);

/* Copy a line from a buffer to an array: */
int NI_LineBufferToArray(NI_LineBuffer*);

/******************************************************************/
/* Multi-dimensional filter support functions */
/******************************************************************/

/* the filter iterator structure: */
typedef struct {
    maybelong strides[MAXDIM], backstrides[MAXDIM];
    maybelong bound1[MAXDIM], bound2[MAXDIM];
} NI_FilterIterator;

/* Initialize a filter iterator: */
int NI_InitFilterIterator(int, maybelong*, maybelong, maybelong*,
                                                                                    maybelong*, NI_FilterIterator*);

/* Calculate the offsets to the filter points, for all border regions and
     the interior of the array: */
int NI_InitFilterOffsets(PyArrayObject*, Bool*, maybelong*,
                    maybelong*, NI_ExtendMode, maybelong**, maybelong*, maybelong**);

/* Move to the next point in an array, possible changing the filter
     offsets, to adapt to boundary conditions: */
#define NI_FILTER_NEXT(iteratorf, iterator1, pointerf, pointer1)  \
{                                                                 \
    int _ii;                                                        \
    for(_ii = (iterator1).rank_m1; _ii >= 0; _ii--) {               \
        maybelong _pp = (iterator1).coordinates[_ii];                 \
        if (_pp < (iterator1).dimensions[_ii]) {                      \
            if (_pp < (iteratorf).bound1[_ii] ||                        \
                                                                    _pp >= (iteratorf).bound2[_ii]) \
                pointerf += (iteratorf).strides[_ii];                     \
            (iterator1).coordinates[_ii]++;                             \
            pointer1 += (iterator1).strides[_ii];                       \
            break;                                                      \
        } else {                                                      \
            (iterator1).coordinates[_ii] = 0;                           \
            pointer1 -= (iterator1).backstrides[_ii];                   \
            pointerf -= (iteratorf).backstrides[_ii];                   \
        }                                                             \
    }                                                               \
}

/* Move to the next point in two arrays, possible changing the pointer
     to the filter offsets when moving into a different region in the
     array: */
#define NI_FILTER_NEXT2(iteratorf, iterator1, iterator2,    \
                                                pointerf, pointer1, pointer2)       \
{                                                           \
    int _ii;                                                  \
    for(_ii = (iterator1).rank_m1; _ii >= 0; _ii--) {         \
        maybelong _pp = (iterator1).coordinates[_ii];           \
        if (_pp < (iterator1).dimensions[_ii]) {                \
            if (_pp < (iteratorf).bound1[_ii] ||                  \
                                                        _pp >= (iteratorf).bound2[_ii]) \
                pointerf += (iteratorf).strides[_ii];               \
            (iterator1).coordinates[_ii]++;                       \
            pointer1 += (iterator1).strides[_ii];                 \
            pointer2 += (iterator2).strides[_ii];                 \
            break;                                                \
        } else {                                                \
            (iterator1).coordinates[_ii] = 0;                     \
            pointer1 -= (iterator1).backstrides[_ii];             \
            pointer2 -= (iterator2).backstrides[_ii];             \
            pointerf -= (iteratorf).backstrides[_ii];             \
        }                                                       \
    }                                                         \
}

/* Move to the next point in three arrays, possible changing the pointer
     to the filter offsets when moving into a different region in the
     array: */
#define NI_FILTER_NEXT3(iteratorf, iterator1, iterator2, iterator3,  \
                                                pointerf, pointer1, pointer2, pointer3)      \
{                                                                    \
    int _ii;                                                           \
    for(_ii = (iterator1).rank_m1; _ii >= 0; _ii--) {                  \
        maybelong _pp = (iterator1).coordinates[_ii];                    \
        if (_pp < (iterator1).dimensions[_ii]) {                         \
            if (_pp < (iteratorf).bound1[_ii] ||                           \
                                                                         _pp >= (iteratorf).bound2[_ii]) \
                pointerf += (iteratorf).strides[_ii];                        \
            (iterator1).coordinates[_ii]++;                                \
            pointer1 += (iterator1).strides[_ii];                          \
            pointer2 += (iterator2).strides[_ii];                          \
            pointer3 += (iterator3).strides[_ii];                          \
            break;                                                         \
        } else {                                                         \
            (iterator1).coordinates[_ii] = 0;                              \
            pointer1 -= (iterator1).backstrides[_ii];                      \
            pointer2 -= (iterator2).backstrides[_ii];                      \
            pointer3 -= (iterator3).backstrides[_ii];                      \
            pointerf -= (iteratorf).backstrides[_ii];                      \
        }                                                                \
    }                                                                  \
}

/* Move the pointer to the filter offsets according to the given
    coordinates: */
#define NI_FILTER_GOTO(iteratorf, iterator, fbase, pointerf) \
{                                                            \
    int _ii;                                                   \
    maybelong _jj;                                             \
    pointerf = fbase;                                          \
    for(_ii = iterator.rank_m1; _ii >= 0; _ii--) {             \
        maybelong _pp = iterator.coordinates[_ii];               \
        maybelong b1 = (iteratorf).bound1[_ii];                  \
        maybelong b2 = (iteratorf).bound2[_ii];                  \
        if (_pp < b1) {                                          \
                _jj = _pp;                                           \
        } else if (_pp > b2 && b2 >= b1) {                       \
                _jj = _pp + b1 - b2;                                 \
        } else {                                                 \
                _jj = b1;                                            \
        }                                                        \
        pointerf += (iteratorf).strides[_ii] * _jj;              \
    }                                                          \
}

typedef struct {
        maybelong *coordinates;
        int size;
        void *next;
} NI_CoordinateBlock;

typedef struct {
        int block_size, rank;
        void *blocks;
} NI_CoordinateList;

NI_CoordinateList* NI_InitCoordinateList(int, int);
int NI_CoordinateListStealBlocks(NI_CoordinateList*, NI_CoordinateList*);
NI_CoordinateBlock* NI_CoordinateListAddBlock(NI_CoordinateList*);
NI_CoordinateBlock* NI_CoordinateListDeleteBlock(NI_CoordinateList*);
void NI_FreeCoordinateList(NI_CoordinateList*);

#endif
