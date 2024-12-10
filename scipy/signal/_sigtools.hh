#ifndef _SCIPY_PRIVATE_SIGNAL__SIGTOOLS_H_
#define _SCIPY_PRIVATE_SIGNAL__SIGTOOLS_H_

#include "Python.h"

#define BOUNDARY_MASK 12
#define OUTSIZE_MASK 3
#define FLIP_MASK  16
#define TYPE_MASK  (32+64+128+256+512)
#define TYPE_SHIFT 5

#define FULL  2
#define SAME  1
#define VALID 0

#define CIRCULAR 8
#define REFLECT  4
#define PAD      0

#define MAXTYPES 21


/* Generally useful structures for passing data into and out of
   subroutines.  Used in the generic routines instead of the
   Python Specific structures so that the routines can be easily
   grabbed and used in another scripting language */

typedef struct {
  char *data;
  int elsize;
} Generic_ptr;

typedef struct {
  char *data;
  npy_intp numels;
  int elsize;
  char *zero;        /* Pointer to Representation of zero */
} Generic_Vector;

typedef struct {
  char *data;
  int  nd;
  npy_intp  *dimensions;
  int  elsize;
  npy_intp  *strides;
  char *zero;         /* Pointer to Representation of zero */
} Generic_Array;

typedef void (MultAddFunction) (char *, npy_intp, char *, npy_intp, char *,
                                npy_intp *, npy_intp *, int, npy_intp, int,
                                npy_intp *, npy_intp *, npy_uintp *);

PyObject*
scipy_signal__sigtools_linear_filter(PyObject * NPY_UNUSED(dummy), PyObject * args);

PyObject*
scipy_signal__sigtools_correlateND(PyObject *NPY_UNUSED(dummy), PyObject *args);

void
scipy_signal__sigtools_linear_filter_module_init(void);

/*
static int index_out_of_bounds(int *, int *, int );
static long compute_offsets (unsigned long *, long *, int *, int *, int *, int *, int);
static int increment(int *, int, int *);
static void convolveND(Generic_Array *, Generic_Array *, Generic_Array *, MultAddFunction *, int);
static void RawFilter(Generic_Vector, Generic_Vector, Generic_Array, Generic_Array, Generic_Array *, Generic_Array *, BasicFilterFunction *, int);
*/


template<int TYPECODE> struct traits {};

// integer types
template<> struct traits<NPY_UBYTE> {
    typedef npy_ubyte real_type;
};
template<> struct traits<NPY_BYTE> {
    typedef npy_byte real_type;
};
template<> struct traits<NPY_USHORT> {
    typedef npy_ushort real_type;
};
template<> struct traits<NPY_SHORT> {
    typedef short real_type;
};
template<> struct traits<NPY_UINT> {
    typedef npy_uint real_type;
};
template<> struct traits<NPY_INT> {
    typedef int real_type;
};
template<> struct traits<NPY_ULONG> {
    typedef npy_ulong real_type;
};
template<> struct traits<NPY_LONG> {
    typedef npy_long real_type;
};
template<> struct traits<NPY_ULONGLONG> {
    typedef npy_ulonglong real_type;
};
template<> struct traits<NPY_LONGLONG> {
    typedef npy_longlong real_type;
};

// real types
template<> struct traits<NPY_FLOAT> {
    typedef float real_type;
};
template<> struct traits<NPY_DOUBLE> {
    typedef double real_type;
};
template<> struct traits<NPY_LONGDOUBLE> {
    typedef npy_longdouble real_type;
};

// complex
template<> struct traits<NPY_CFLOAT> {
    typedef float real_type;
};
template<> struct traits<NPY_CDOUBLE> {
    typedef double real_type;
};
template<> struct traits<NPY_CLONGDOUBLE> {
    typedef npy_longdouble real_type;
};

#endif
