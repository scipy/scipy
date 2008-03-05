/* -*- C -*- */
#ifdef SWIGPYTHON

%module _umfpack

/*
  See umfpack.py for more information.

  Created by: Robert Cimrman
*/

%{
#include <umfpack.h>
#include "numpy/arrayobject.h"
%}

%feature("autodoc", "1");

#include <umfpack.h>

%init %{
    import_array();
%}

%{
/*!
  Appends @a what to @a where. On input, @a where need not to be a tuple, but on
  return it always is.

  @par Revision history:
  - 17.02.2005, c
*/
PyObject *helper_appendToTuple( PyObject *where, PyObject *what ) {
  PyObject *o2, *o3;

  if ((!where) || (where == Py_None)) {
    where = what;
  } else {
    if (!PyTuple_Check( where )) {
      o2 = where;
      where = PyTuple_New( 1 );
      PyTuple_SetItem( where, 0, o2 );
    }
    o3 = PyTuple_New( 1 );
    PyTuple_SetItem( o3, 0, what );
    o2 = where;
    where = PySequence_Concat( o2, o3 );
    Py_DECREF( o2 );
    Py_DECREF( o3 );
  }
  return where;
}

/*!
  Gets PyArrayObject from a PyObject.

  @par Revision history:
  - 22.02.2005, c
  - 03.03.2005
  - 25.11.2005
  - 30.11.2005
  - 01.12.2005
*/
PyArrayObject *helper_getCArrayObject( PyObject *input, int type,
				       int minDim, int maxDim ) {
  PyArrayObject *obj;

  if (PyArray_Check( input )) {
    obj = (PyArrayObject *) input;
    if (!PyArray_ISCARRAY( obj )) {
      PyErr_SetString( PyExc_TypeError, "not a C array" );
      return NULL;
    }
    obj = (PyArrayObject *)
      PyArray_ContiguousFromAny( input, type, minDim, maxDim );
    if (!obj) return NULL;
  } else {
    PyErr_SetString( PyExc_TypeError, "not an array" );
    return NULL;
  }
  return obj;
}
%}

/*!
  Use for arrays as input arguments. Could be also used for changing an array
  in place.

  @a rtype ... return this C data type
  @a ctype ... C data type of the C function
  @a atype ... PyArray_* suffix

  @par Revision history:
  - 30.11.2005, c
*/
#define ARRAY_IN( rtype, ctype, atype ) \
%typemap( in ) (ctype *array) { \
  PyArrayObject *obj; \
  obj = helper_getCArrayObject( $input, PyArray_##atype, 1, 1 ); \
  if (!obj) return NULL; \
  $1 = (rtype *) obj->data; \
  Py_DECREF( obj ); \
};

/*!
  @par Revision history:
  - 30.11.2005, c
*/
#define CONF_IN( arSize ) \
%typemap( in ) (double conf [arSize]) { \
  PyArrayObject *obj; \
  obj = helper_getCArrayObject( $input, PyArray_DOUBLE, 1, 1 ); \
  if (!obj) return NULL; \
  if ((obj->nd != 1) || (obj->dimensions[0] != arSize)) { \
    PyErr_SetString( PyExc_ValueError, "wrong Control/Info array size" ); \
    Py_DECREF( obj ); \
    return NULL; \
  } \
  $1 = (double *) obj->data; \
  Py_DECREF( obj ); \
};

/*!
  @par Revision history:
  - 01.12.2005, c
  - 02.12.2005
*/
#define OPAQUE_ARGOUT( ttype ) \
%typemap( in, numinputs=0 ) ttype* opaque_argout( ttype tmp ) { \
  $1 = &tmp; \
}; \
%typemap( argout ) ttype* opaque_argout { \
  PyObject *obj; \
  obj = SWIG_NewPointerObj( (ttype) (*$1), $*1_descriptor, 0 ); \
  $result = helper_appendToTuple( $result, obj ); \
};

/*!
  @par Revision history:
  - 02.12.2005, c
*/
#define OPAQUE_ARGINOUT( ttype ) \
%typemap( in ) ttype* opaque_arginout( ttype tmp ) { \
  if ((SWIG_ConvertPtr( $input,(void **) &tmp, $*1_descriptor, \
			SWIG_POINTER_EXCEPTION)) == -1) return NULL; \
  $1 = &tmp; \
}; \
%typemap( argout ) ttype* opaque_arginout { \
  PyObject *obj; \
  obj = SWIG_NewPointerObj( (ttype) (*$1), $*1_descriptor, 0 ); \
  $result = helper_appendToTuple( $result, obj ); \
};

ARRAY_IN( int, const int, INT )
%apply const int *array {
    const int Ap [ ],
    const int Ai [ ]
};

ARRAY_IN( long, const long, LONG )
%apply const long *array {
    const long Ap [ ],
    const long Ai [ ]
};

ARRAY_IN( double, const double, DOUBLE )
%apply const double *array {
    const double Ax [ ],
    const double Az [ ],
    const double B [ ],
    const double Bx [ ],
    const double Bz [ ]
};

ARRAY_IN( double, double, DOUBLE )
%apply double *array {
    double X [ ],
    double Xx [ ],
    double Xz [ ]
};

CONF_IN( UMFPACK_CONTROL )
%apply (double conf [UMFPACK_CONTROL]) {
    double Control [ANY]
};

CONF_IN( UMFPACK_INFO )
%apply double conf [UMFPACK_INFO] {
    double Info [ANY]
};

%include <umfpack.h>
%include <umfpack_solve.h>
%include <umfpack_defaults.h>
%include <umfpack_triplet_to_col.h>
%include <umfpack_col_to_triplet.h>
%include <umfpack_transpose.h>
%include <umfpack_scale.h>

%include <umfpack_report_symbolic.h>
%include <umfpack_report_numeric.h>
%include <umfpack_report_info.h>
%include <umfpack_report_control.h>

/*
  The order is important below!
*/

OPAQUE_ARGOUT( void * )
%apply  void ** opaque_argout {
    void **Symbolic,
    void **Numeric
}

%include <umfpack_symbolic.h>
%include <umfpack_numeric.h>


OPAQUE_ARGINOUT( void * )
%apply  void ** opaque_arginout {
    void **Symbolic,
    void **Numeric
}

%include <umfpack_free_symbolic.h>
%include <umfpack_free_numeric.h>



/*
 * wnbell - attempt to get L,U,P,Q out
 */
%include "typemaps.i"
%apply int  *OUTPUT {
    int *lnz,
    int *unz,
    int *n_row,
    int *n_col,
    int *nz_udiag
};
%apply long *OUTPUT {
    long *lnz,
    long *unz,
    long *n_row,
    long *n_col,
    long *nz_udiag
};
%include <umfpack_get_lunz.h>


ARRAY_IN( double, double, DOUBLE )
%apply double *array {
    double Lx [ ],
    double Lz [ ],
    double Ux [ ],
    double Uz [ ],
    double Dx [ ],
    double Dz [ ],
    double Rs [ ]
};

ARRAY_IN( int, int, INT )
%apply int *array {
    int Lp [ ],
    int Lj [ ],
    int Up [ ],
    int Ui [ ],
    int P [ ],
    int Q [ ]
};
%apply int  *OUTPUT { int *do_recip};
%include <umfpack_get_numeric.h>

#endif
