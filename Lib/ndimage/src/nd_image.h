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

#ifndef ND_IMAGE_H
#define ND_IMAGE_H

#include "Python.h"
#include <numpy/noprefix.h>

#define NI_MAXDIM NPY_MAXDIMS

typedef npy_intp maybelong;
#define MAXDIM NPY_MAXDIMS

typedef enum
{ 
     tAny=-1,
     tBool=PyArray_BOOL,
     tInt8=PyArray_INT8,
     tUInt8=PyArray_UINT8,
     tInt16=PyArray_INT16,
     tUInt16=PyArray_UINT16,
     tInt32=PyArray_INT32,
     tUInt32=PyArray_UINT32,
     tInt64=PyArray_INT64,
     tUInt64=PyArray_UINT64,
     tFloat32=PyArray_FLOAT32,
     tFloat64=PyArray_FLOAT64,
     tComplex32=PyArray_COMPLEX64,
     tComplex64=PyArray_COMPLEX128,
     tObject=PyArray_OBJECT,        /* placeholder... does nothing */
     tMaxType=PyArray_NTYPES,
     tDefault = tFloat64,
#if NPY_BITSOF_LONG == 64
     tLong = tInt64,
#else
     tLong = tInt32,
#endif
} NumarrayType;

/* int NI_GetArrayRank(PyArrayObject*);
NumarrayType NI_GetArrayType(PyArrayObject*);
void NI_GetArrayDimensions(PyArrayObject*, int*);
void NI_GetArrayStrides(PyArrayObject*, int*);
char* NI_GetArrayData(PyArrayObject*);
int NI_ShapeEqual(PyArrayObject*, PyArrayObject*);
int NI_CheckArray(PyArrayObject*, NumarrayType, int, int*); */

#endif
