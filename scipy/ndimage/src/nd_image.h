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

#define PY_ARRAY_UNIQUE_SYMBOL _scipy_ndimage_ARRAY_API
#include <numpy/noprefix.h>

#include "numpy/npy_3kcompat.h"

/* Eventually get rid of everything below this line */

typedef enum
{
         tAny=-1,
         tBool=NPY_BOOL,
         tInt8=NPY_INT8,
         tUInt8=NPY_UINT8,
         tInt16=NPY_INT16,
         tUInt16=NPY_UINT16,
         tInt32=NPY_INT32,
         tUInt32=NPY_UINT32,
         tInt64=NPY_INT64,
         tUInt64=NPY_UINT64,
         tFloat32=NPY_FLOAT32,
         tFloat64=NPY_FLOAT64,
         tComplex64=NPY_COMPLEX64,
         tComplex128=NPY_COMPLEX128,
         tObject=NPY_OBJECT,        /* placeholder... does nothing */
         tMaxType=NPY_NTYPES,
         tDefault=NPY_FLOAT64
} NumarrayType;

#define HAS_UINT64 1

#endif /* ND_IMAGE_H */
