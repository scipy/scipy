/* -*- C -*- */
%{
#ifndef SWIG_FILE_WITH_INIT
#  define NO_IMPORT_ARRAY
#endif
#include "stdio.h"
#include <numpy/arrayobject.h>
#include "complex_ops.h"
#include "bool_ops.h"


/* The following code originally appeared in
 * enthought/kiva/agg/src/numeric.i written by Eric Jones.  It was
 * translated from C++ to C by John Hunter.  Bill Spotz has modified
 * it slightly to fix some minor bugs, upgrade to numpy (all
 * versions), add some comments and some functionality.
 */

/* Macros to extract array attributes.
 */
#define is_array(a)            ((a) && PyArray_Check((PyArrayObject *)a))
#define array_type(a)          (int)(PyArray_TYPE(a))
#define array_numdims(a)       (((PyArrayObject *)a)->nd)
#define array_dimensions(a)    (((PyArrayObject *)a)->dimensions)
#define array_size(a,i)        (((PyArrayObject *)a)->dimensions[i])
#define array_data(a)          (((PyArrayObject *)a)->data)
#define array_is_contiguous(a) (PyArray_ISCONTIGUOUS(a))
#define array_is_native(a)     (PyArray_ISNOTSWAPPED(a))

/* Given a PyObject, return a string describing its type.
 */
static const char* pytype_string(PyObject* py_obj) {
  if (py_obj == NULL          ) return "C NULL value";
  if (py_obj == Py_None       ) return "Python None" ;
  if (PyCallable_Check(py_obj)) return "callable"    ;
  if (PyString_Check(  py_obj)) return "string"      ;
  if (PyInt_Check(     py_obj)) return "int"         ;
  if (PyFloat_Check(   py_obj)) return "float"       ;
  if (PyDict_Check(    py_obj)) return "dict"        ;
  if (PyList_Check(    py_obj)) return "list"        ;
  if (PyTuple_Check(   py_obj)) return "tuple"       ;
  if (PyFile_Check(    py_obj)) return "file"        ;
  if (PyModule_Check(  py_obj)) return "module"      ;
  if (PyInstance_Check(py_obj)) return "instance"    ;

  return "unknown type";
}

/* Given a NumPy typecode, return a string describing the type.
 */
static const char* typecode_string(int typecode) {
  static const char* type_names[25] = {"bool", "byte", "unsigned byte",
				 "short", "unsigned short", "int",
				 "unsigned int", "long", "unsigned long",
				 "long long", "unsigned long long",
				 "float", "double", "long double",
				 "complex float", "complex double",
				 "complex long double", "object",
				 "string", "unicode", "void", "ntypes",
				 "notype", "char", "unknown"};
  return typecode < 24 ? type_names[typecode] : type_names[24];
}

/* Make sure input has correct numpy type.  Allow character and byte
 * to match.  Also allow int and long to match.  This is deprecated.
 * You should use PyArray_EquivTypenums() instead.
 */
static int type_match(int actual_type, int desired_type) {
  return PyArray_EquivTypenums(actual_type, desired_type);
}

/* Given a PyObject pointer, cast it to a PyArrayObject pointer if
 * legal.  If not, set the python error string appropriately and
 * return NULL.
 */
static PyArrayObject* obj_to_array_no_conversion(PyObject* input, int typecode) {
  PyArrayObject* ary = NULL;
  if (is_array(input) && (typecode == NPY_NOTYPE ||
			  PyArray_EquivTypenums(array_type(input), typecode))) {
    ary = (PyArrayObject*) input;
  }
  else if is_array(input) {
    const char* desired_type = typecode_string(typecode);
    const char* actual_type  = typecode_string(array_type(input));
    PyErr_Format(PyExc_TypeError, 
		 "Array of type '%s' required.  Array of type '%s' given", 
		 desired_type, actual_type);
    ary = NULL;
  }
  else {
    const char * desired_type = typecode_string(typecode);
    const char * actual_type  = pytype_string(input);
    PyErr_Format(PyExc_TypeError, 
		 "Array of type '%s' required.  A '%s' was given", 
		 desired_type, actual_type);
    ary = NULL;
  }
  return ary;
}

/* Convert the given PyObject to a NumPy array with the given
 * typecode.  On success, return a valid PyArrayObject* with the
 * correct type.  On failure, the python error string will be set and
 * the routine returns NULL.
 */
static PyArrayObject* obj_to_array_allow_conversion(PyObject* input, int typecode,
                                             int* is_new_object) {
  PyArrayObject* ary = NULL;
  PyObject* py_obj;
  if (is_array(input) && (typecode == NPY_NOTYPE ||
			  PyArray_EquivTypenums(array_type(input),typecode))) {
    ary = (PyArrayObject*) input;
    *is_new_object = 0;
  }
  else {
    py_obj = PyArray_FromObject(input, typecode, 0, 0);
    /* If NULL, PyArray_FromObject will have set python error value.*/
    ary = (PyArrayObject*) py_obj;
    *is_new_object = 1;
  }
  return ary;
}

/* Given a PyArrayObject, check to see if it is contiguous.  If so,
 * return the input pointer and flag it as not a new object.  If it is
 * not contiguous, create a new PyArrayObject using the original data,
 * flag it as a new object and return the pointer.
 */
static PyArrayObject* make_contiguous(PyArrayObject* ary, int* is_new_object,
                               int min_dims, int max_dims) {
  PyArrayObject* result;
  if (array_is_contiguous(ary)) {
    result = ary;
    *is_new_object = 0;
  }
  else {
    result = (PyArrayObject*) PyArray_ContiguousFromObject((PyObject*)ary, 
							   array_type(ary), 
							   min_dims,
							   max_dims);
    *is_new_object = 1;
  }
  return result;
}

/* Convert a given PyObject to a contiguous PyArrayObject of the
 * specified type.  If the input object is not a contiguous
 * PyArrayObject, a new one will be created and the new object flag
 * will be set.
 */
static PyArrayObject* obj_to_array_contiguous_allow_conversion(PyObject* input,
                                                        int typecode,
                                                        int* is_new_object) {
  int is_new1 = 0;
  int is_new2 = 0;
  PyArrayObject* ary2;
  PyArrayObject* ary1 = obj_to_array_allow_conversion(input, typecode, &is_new1);
  if (ary1) {
    ary2 = make_contiguous(ary1, &is_new2, 0, 0);
    if ( is_new1 && is_new2) {
      Py_DECREF(ary1);
    }
    ary1 = ary2;    
  }
  *is_new_object = is_new1 || is_new2;
  return ary1;
}

/* Test whether a python object is contiguous.  If array is
 * contiguous, return 1.  Otherwise, set the python error string and
 * return 0.
 */
static int require_contiguous(PyArrayObject* ary) {
  int contiguous = 1;
  if (!array_is_contiguous(ary)) {
    PyErr_SetString(PyExc_TypeError,
		    "Array must be contiguous.  A non-contiguous array was given");
    contiguous = 0;
  }
  return contiguous;
}

/* Require that a numpy array is not byte-swapped.  If the array is
 * not byte-swapped, return 1.  Otherwise, set the python error string
 * and return 0.
 */
static int require_native(PyArrayObject* ary) {
  int native = 1;
  if (!array_is_native(ary)) {
    PyErr_SetString(PyExc_TypeError,
		    "Array must have native byteorder.  A byte-swapped array was given");
    native = 0;
  }
  return native;
}

/* Require the given PyArrayObject to have a specified number of
 * dimensions.  If the array has the specified number of dimensions,
 * return 1.  Otherwise, set the python error string and return 0.
 */
static int require_dimensions(PyArrayObject* ary, int exact_dimensions) {
  int success = 1;
  if (array_numdims(ary) != exact_dimensions) {
    PyErr_Format(PyExc_TypeError, 
		 "Array must have %d dimensions.  Given array has %d dimensions", 
		 exact_dimensions, array_numdims(ary));
    success = 0;
  }
  return success;
}

/* Require the given PyArrayObject to have one of a list of specified
 * number of dimensions.  If the array has one of the specified number
 * of dimensions, return 1.  Otherwise, set the python error string
 * and return 0.
 */
static int require_dimensions_n(PyArrayObject* ary, int* exact_dimensions, int n) {
  int success = 0;
  int i;
  char dims_str[255] = "";
  char s[255];
  for (i = 0; i < n && !success; i++) {
    if (array_numdims(ary) == exact_dimensions[i]) {
      success = 1;
    }
  }
  if (!success) {
    for (i = 0; i < n-1; i++) {
      sprintf(s, "%d, ", exact_dimensions[i]);                
      strcat(dims_str,s);
    }
    sprintf(s, " or %d", exact_dimensions[n-1]);            
    strcat(dims_str,s);
    PyErr_Format(PyExc_TypeError, 
		 "Array must be have %s dimensions.  Given array has %d dimensions",
		 dims_str, array_numdims(ary));
  }
  return success;
}    

/* Require the given PyArrayObject to have a specified shape.  If the
 * array has the specified shape, return 1.  Otherwise, set the python
 * error string and return 0.
 */
static int require_size(PyArrayObject* ary, npy_intp* size, int n) {
  int i;
  int success = 1;
  int len;
  char desired_dims[255] = "[";
  char s[255];
  char actual_dims[255] = "[";
  for(i=0; i < n;i++) {
    if (size[i] != -1 &&  size[i] != array_size(ary,i)) {
      success = 0;    
    }
  }
  if (!success) {
    for (i = 0; i < n; i++) {
      if (size[i] == -1) {
	sprintf(s, "*,");                
      }
      else
      {
	sprintf(s,"%" NPY_INTP_FMT ",", size[i]);                
      }    
      strcat(desired_dims,s);
    }
    len = strlen(desired_dims);
    desired_dims[len-1] = ']';
    for (i = 0; i < n; i++) {
      sprintf(s,"%" NPY_INTP_FMT ",", array_size(ary,i));                            
      strcat(actual_dims,s);
    }
    len = strlen(actual_dims);
    actual_dims[len-1] = ']';
    PyErr_Format(PyExc_TypeError, 
		 "Array must be have shape of %s.  Given array has shape of %s",
		 desired_dims, actual_dims);
  }
  return success;
}
/* End John Hunter translation (with modifications by Bill Spotz) */





/*!
  Appends @a what to @a where. On input, @a where need not to be a tuple, but on
  return it always is.

  @par Revision history:
  - 17.02.2005, c
*/
static PyObject *helper_appendToTuple( PyObject *where, PyObject *what ) {
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




%}

/* TYPEMAP_IN macros
 *
 * This family of typemaps allows pure input C arguments of the form
 *
 *     (type* IN_ARRAY1, int DIM1)
 *     (type* IN_ARRAY2, int DIM1, int DIM2)
 *
 * where "type" is any type supported by the Numeric module, to be
 * called in python with an argument list of a single array (or any
 * python object that can be passed to the Numeric.array constructor
 * to produce an arrayof te specified shape).  This can be applied to
 * a existing functions using the %apply directive:
 *
 *     %apply (double* IN_ARRAY1, int DIM1) {double* series, int length}
 *     %apply (double* IN_ARRAY2, int DIM1, int DIM2) {double* mx, int rows, int cols}
 *     double sum(double* series, int length);
 *     double max(double* mx, int rows, int cols);
 *
 * or with
 *
 *     double sum(double* IN_ARRAY1, int DIM1);
 *     double max(double* IN_ARRAY2, int DIM1, int DIM2);
 */

/* One dimensional input arrays */
%define TYPEMAP_IN1(type,typecode)
  %typemap(in) type* IN_ARRAY1 (PyArrayObject* array=NULL, int is_new_object) {
  npy_intp size[1] = {-1};
  array = obj_to_array_contiguous_allow_conversion($input, typecode, &is_new_object);
  if (!array || !require_dimensions(array,1) || !require_size(array,size,1)
             || !require_contiguous(array)   || !require_native(array)) SWIG_fail;

  $1 = (type*) array->data;
}
%typemap(freearg) type*  IN_ARRAY1 {
  if (is_new_object$argnum && array$argnum) { Py_DECREF(array$argnum); }
}
%enddef




 /* Two dimensional input arrays */
%define TYPEMAP_IN2(type,typecode)
  %typemap(in) (type* IN_ARRAY2)
               (PyArrayObject* array=NULL, int is_new_object) {
  npy_intp size[2] = {-1,-1};
  array = obj_to_array_contiguous_allow_conversion($input, typecode, &is_new_object);
  if (!array || !require_dimensions(array,2) || !require_size(array,size,1) 
             || !require_contiguous(array)   || !require_native(array)) SWIG_fail;
  $1 = (type*) array->data;
}
%typemap(freearg) (type* IN_ARRAY2) {
  if (is_new_object$argnum && array$argnum) { Py_DECREF(array$argnum); }
}
%enddef


/* TYPEMAP_INPLACE macros
 *
 * This family of typemaps allows input/output C arguments of the form
 *
 *     (type* INPLACE_ARRAY1, int DIM1)
 *     (type* INPLACE_ARRAY2, int DIM1, int DIM2)
 *
 * where "type" is any type supported by the Numeric module, to be
 * called in python with an argument list of a single contiguous
 * Numeric array.  This can be applied to an existing function using
 * the %apply directive:
 *
 *     %apply (double* INPLACE_ARRAY1, int DIM1) {double* series, int length}
 *     %apply (double* INPLACE_ARRAY2, int DIM1, int DIM2) {double* mx, int rows, int cols}
 *     void negate(double* series, int length);
 *     void normalize(double* mx, int rows, int cols);
 *     
 *
 * or with
 *
 *     void sum(double* INPLACE_ARRAY1, int DIM1);
 *     void sum(double* INPLACE_ARRAY2, int DIM1, int DIM2);
 */

 /* One dimensional input/output arrays */
%define TYPEMAP_INPLACE1(type,typecode)
%typemap(in) (type* INPLACE_ARRAY) (PyArrayObject* temp=NULL) {
  temp = obj_to_array_no_conversion($input,typecode);
  if (!temp  || !require_contiguous(temp) || !require_native(temp)) SWIG_fail;
  $1 = (type*) array_data(temp);
}
%enddef


 /* Two dimensional input/output arrays */
%define TYPEMAP_INPLACE2(type,typecode)
  %typemap(in) (type* INPLACE_ARRAY2) (PyArrayObject* temp=NULL) {
  temp = obj_to_array_no_conversion($input,typecode);
  if (!temp || !require_contiguous(temp) || !require_native(temp)) SWIG_fail;
  $1 = (type*) array_data(temp);
}
%enddef




/* TYPEMAP_ARGOUT macros
 *
 * This family of typemaps allows output C arguments of the form
 *
 *     (type* ARGOUT_ARRAY[ANY])
 *     (type* ARGOUT_ARRAY[ANY][ANY])
 *
 * where "type" is any type supported by the Numeric module, to be
 * called in python with an argument list of a single contiguous
 * Numeric array.  This can be applied to an existing function using
 * the %apply directive:
 *
 *     %apply (double* ARGOUT_ARRAY[ANY] {double series, int length}
 *     %apply (double* ARGOUT_ARRAY[ANY][ANY]) {double* mx, int rows, int cols}
 *     void negate(double* series, int length);
 *     void normalize(double* mx, int rows, int cols);
 *     
 *
 * or with
 *
 *     void sum(double* ARGOUT_ARRAY[ANY]);
 *     void sum(double* ARGOUT_ARRAY[ANY][ANY]);
 */

 /* One dimensional input/output arrays */
/*%define TYPEMAP_ARGOUT1(type,typecode)
%typemap(in,numinputs=0) type ARGOUT_ARRAY[ANY] {
  $1 = (type*) malloc($1_dim0*sizeof(type));
  if (!$1) {
    PyErr_SetString(PyExc_RuntimeError, "Failed to allocate memory");
    SWIG_fail;
  }
}
%typemap(argout) ARGOUT_ARRAY[ANY] {
  int dimensions[1] = {$1_dim0};
  PyObject* outArray = PyArray_FromDimsAndData(1, dimensions, typecode, (char*)$1);
}
%enddef
*/

 /* Two dimensional input/output arrays */
/*%define TYPEMAP_ARGOUT2(type,typecode)
  %typemap(in) (type* ARGOUT_ARRAY2, int DIM1, int DIM2) (PyArrayObject* temp=NULL) {
  temp = obj_to_array_no_conversion($input,typecode);
  if (!temp || !require_contiguous(temp) || !require_native(temp)) SWIG_fail;
  $1 = (type*) array(temp);
  $2 = temp->dimensions[0];
  $3 = temp->dimensions[1];
}
%enddef
*/





/*
 * WNBELL additions
 */


/*
 * Use STL vectors for ARGOUTs
 */
%define VEC_ARRAY_ARGOUT( ctype, atype  ) 
%typemap( in, numinputs=0 ) std::vector<ctype>* array_argout( std::vector<ctype>* tmp ) {
  tmp = new std::vector<ctype>(); 
  $1 = tmp; 
}; 
%typemap( argout ) std::vector<ctype>* array_argout { 
  npy_intp length = ($1)->size(); 
  PyObject *obj = PyArray_SimpleNew(1, &length, ##atype); 
  if (length > 0) {
    memcpy(PyArray_DATA(obj), &((*($1))[0]), sizeof(ctype)*length);
  }
  delete $1; 
  $result = helper_appendToTuple( $result, (PyObject *)obj ); 
}; 
%enddef



/*
  * make typechecks - used for overloading
  */
%include "typemaps.i"

%define NPY_TYPECHECK( ctype, atype )
%typemap(typecheck) ctype *, const ctype *, ctype [], const ctype []
{
    $1 = (is_array($input) && PyArray_CanCastSafely(PyArray_TYPE($input), ##atype)) ? 1 : 0;
};
%enddef



%define INSTANTIATE_TYPEMAPS(type,typecode)
TYPEMAP_IN1(      type,typecode)
TYPEMAP_IN1(const type,typecode)
TYPEMAP_IN2(      type,typecode)
TYPEMAP_IN2(const type,typecode)
TYPEMAP_INPLACE1(type,typecode)
TYPEMAP_INPLACE2(type,typecode)
/*TYPEMAP_ARGOUT1(type, typecode)
TYPEMAP_ARGOUT2(type, typecode)*/
VEC_ARRAY_ARGOUT(type, typecode)
NPY_TYPECHECK(type, typecode)
%enddef


INSTANTIATE_TYPEMAPS(npy_bool_wrapper,        NPY_BOOL       )
INSTANTIATE_TYPEMAPS(char,                    NPY_CHAR       )
INSTANTIATE_TYPEMAPS(unsigned char,           NPY_UBYTE      )
INSTANTIATE_TYPEMAPS(signed char,             NPY_BYTE       )
INSTANTIATE_TYPEMAPS(short,                   NPY_SHORT      )
INSTANTIATE_TYPEMAPS(unsigned short,          NPY_USHORT     )
INSTANTIATE_TYPEMAPS(int,                     NPY_INT        )
INSTANTIATE_TYPEMAPS(unsigned int,            NPY_UINT       )
INSTANTIATE_TYPEMAPS(long,                    NPY_LONG       )
INSTANTIATE_TYPEMAPS(long long,               NPY_LONGLONG   )
INSTANTIATE_TYPEMAPS(unsigned long long,      NPY_ULONGLONG  )
INSTANTIATE_TYPEMAPS(float,                   NPY_FLOAT      )
INSTANTIATE_TYPEMAPS(double,                  NPY_DOUBLE     )
INSTANTIATE_TYPEMAPS(long double,             NPY_LONGDOUBLE )
INSTANTIATE_TYPEMAPS(npy_cfloat_wrapper,      NPY_CFLOAT     )
INSTANTIATE_TYPEMAPS(npy_cdouble_wrapper,     NPY_CDOUBLE    )
INSTANTIATE_TYPEMAPS(npy_clongdouble_wrapper, NPY_CLONGDOUBLE)
INSTANTIATE_TYPEMAPS(PyObject,                NPY_OBJECT     )



#undef TYPEMAP_IN1
#undef TYPEMAP_IN2
#undef TYPEMAP_INPLACE1
#undef TYPEMAP_INPLACE2
#undef TYPEMAP_ARGOUT1
#undef TYPEMAP_ARGOUT2
#under NPY_TYPECHECK

