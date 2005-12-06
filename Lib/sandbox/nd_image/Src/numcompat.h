
#ifndef NUMCOMPAT_H
#define NUMCOMPAT_H

#define PY_ARRAY_UNIQUE_SYMBOL NUMCOMPAT
#include "scipy/arrayobject.h"

#define SZ_BUF  79
#define MAXDIM MAX_DIMS

#define maybelong intp

typedef enum 
{
  tAny,
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
  tDefault = tFloat64,
#if BITSOF_LONG == 64
  tLong = tInt64,
#else
  tLong = tInt32,
#endif
  tMaxType  
} NumarrayType;

#define nNumarrayType 16

typedef enum 
{
	NUM_LITTLE_ENDIAN=0,
	NUM_BIG_ENDIAN = 1
} NumarrayByteOrder;

#define IS_CARRAY (CONTIGUOUS | ALIGNED | NOTSWAPPED)

typedef enum
{
        NUM_CONTIGUOUS=CONTIGUOUS,
        NUM_NOTSWAPPED=NOTSWAPPED,
        NUM_ALIGNED=ALIGNED,
        NUM_WRITABLE=WRITEABLE,
        NUM_COPY=ENSURECOPY,

        NUM_C_ARRAY  = (NUM_CONTIGUOUS | NUM_ALIGNED | NUM_NOTSWAPPED),
        NUM_UNCONVERTED = 0
} NumRequirements;


#define UNCONVERTED 0
#define C_ARRAY     (NUM_CONTIGUOUS | NUM_NOTSWAPPED | NUM_ALIGNED)

#define import_libnumarray() import_array()
#define Complex64 Complex64_

typedef struct { Float32 r, i; } Complex32;
typedef struct { Float64 r, i; } Complex64;

#define _NAtype_toDescr(type) (((type)==tAny) ? NULL :	\
			       PyArray_DescrFromType(type))

#define NA_InputArray(obj, type, flags) \
	(PyArrayObject *)\
		PyArray_FromAny(obj, _NAtype_toDescr(type), 0, 0, flags)

#define NA_OutputArray(obj, type, flags)			\
	(PyArrayObject *) \
		((PyArray_Check(obj) && PyArray_CHKFLAGS(obj, flags)) ? \
		 PyArray_Empty(PyArray_NDIM(obj),			\
			       PyArray_DIMS(obj),			\
			       ((type == tAny) && \
				(Py_INCREF(PyArray_DESCR(obj)) || 1)) ? \
			       PyArray_DESCR(obj) :			\
			       PyArray_DescrFromType(type), 0) :	\
		 PyArray_FromAny(obj, _NAtype_toDescr(type), 0, 0, flags))

#define NA_IoArray(obj, type, flags) \
	(PyArrayObject *)						\
		PyArray_FromAny(obj, _NAtype_toDescr(type), 0, 0, flags | \
				UPDATEIFCOPY)

#define NA_vNewArray(buffer, type, ndim, shape) \
	(PyArrayObject *)						\
		PyArray_NewFromDescr(&PyArray_Type,			\
				     PyArray_DescrFromType(type),	\
				     ndim, shape, NULL, buffer, CARRAY_FLAGS, \
				     NULL)

#define NA_elements PyArray_SIZE

#define NA_OFFSETDATA(a) ((void *) PyArray_DATA(a))

PyArrayObject *NA_NewArray(void *buffer, NumarrayType type, int ndim, ...);

#endif
