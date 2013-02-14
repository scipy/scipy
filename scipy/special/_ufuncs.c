

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#ifndef Py_PYTHON_H
    #error Python headers needed to compile C extensions, please install development version of Python.
#elif PY_VERSION_HEX < 0x02040000
    #error Cython requires Python 2.4+.
#else
#include <stddef.h> 
#ifndef offsetof
#define offsetof(type, member) ( (size_t) & ((type*)0) -> member )
#endif
#if !defined(WIN32) && !defined(MS_WINDOWS)
  #ifndef __stdcall
    #define __stdcall
  #endif
  #ifndef __cdecl
    #define __cdecl
  #endif
  #ifndef __fastcall
    #define __fastcall
  #endif
#endif
#ifndef DL_IMPORT
  #define DL_IMPORT(t) t
#endif
#ifndef DL_EXPORT
  #define DL_EXPORT(t) t
#endif
#ifndef PY_LONG_LONG
  #define PY_LONG_LONG LONG_LONG
#endif
#ifndef Py_HUGE_VAL
  #define Py_HUGE_VAL HUGE_VAL
#endif
#ifdef PYPY_VERSION
#define CYTHON_COMPILING_IN_PYPY 1
#define CYTHON_COMPILING_IN_CPYTHON 0
#else
#define CYTHON_COMPILING_IN_PYPY 0
#define CYTHON_COMPILING_IN_CPYTHON 1
#endif
#if PY_VERSION_HEX < 0x02050000
  typedef int Py_ssize_t;
  #define PY_SSIZE_T_MAX INT_MAX
  #define PY_SSIZE_T_MIN INT_MIN
  #define PY_FORMAT_SIZE_T ""
  #define CYTHON_FORMAT_SSIZE_T ""
  #define PyInt_FromSsize_t(z) PyInt_FromLong(z)
  #define PyInt_AsSsize_t(o)   __Pyx_PyInt_AsInt(o)
  #define PyNumber_Index(o)    ((PyNumber_Check(o) && !PyFloat_Check(o)) ? PyNumber_Int(o) : \
                                (PyErr_Format(PyExc_TypeError, \
                                              "expected index value, got %.200s", Py_TYPE(o)->tp_name), \
                                 (PyObject*)0))
  #define __Pyx_PyIndex_Check(o) (PyNumber_Check(o) && !PyFloat_Check(o) && \
                                  !PyComplex_Check(o))
  #define PyIndex_Check __Pyx_PyIndex_Check
  #define PyErr_WarnEx(category, message, stacklevel) PyErr_Warn(category, message)
  #define __PYX_BUILD_PY_SSIZE_T "i"
#else
  #define __PYX_BUILD_PY_SSIZE_T "n"
  #define CYTHON_FORMAT_SSIZE_T "z"
  #define __Pyx_PyIndex_Check PyIndex_Check
#endif
#if PY_VERSION_HEX < 0x02060000
  #define Py_REFCNT(ob) (((PyObject*)(ob))->ob_refcnt)
  #define Py_TYPE(ob)   (((PyObject*)(ob))->ob_type)
  #define Py_SIZE(ob)   (((PyVarObject*)(ob))->ob_size)
  #define PyVarObject_HEAD_INIT(type, size) \
          PyObject_HEAD_INIT(type) size,
  #define PyType_Modified(t)
  typedef struct {
     void *buf;
     PyObject *obj;
     Py_ssize_t len;
     Py_ssize_t itemsize;
     int readonly;
     int ndim;
     char *format;
     Py_ssize_t *shape;
     Py_ssize_t *strides;
     Py_ssize_t *suboffsets;
     void *internal;
  } Py_buffer;
  #define PyBUF_SIMPLE 0
  #define PyBUF_WRITABLE 0x0001
  #define PyBUF_FORMAT 0x0004
  #define PyBUF_ND 0x0008
  #define PyBUF_STRIDES (0x0010 | PyBUF_ND)
  #define PyBUF_C_CONTIGUOUS (0x0020 | PyBUF_STRIDES)
  #define PyBUF_F_CONTIGUOUS (0x0040 | PyBUF_STRIDES)
  #define PyBUF_ANY_CONTIGUOUS (0x0080 | PyBUF_STRIDES)
  #define PyBUF_INDIRECT (0x0100 | PyBUF_STRIDES)
  #define PyBUF_RECORDS (PyBUF_STRIDES | PyBUF_FORMAT | PyBUF_WRITABLE)
  #define PyBUF_FULL (PyBUF_INDIRECT | PyBUF_FORMAT | PyBUF_WRITABLE)
  typedef int (*getbufferproc)(PyObject *, Py_buffer *, int);
  typedef void (*releasebufferproc)(PyObject *, Py_buffer *);
#endif
#if PY_MAJOR_VERSION < 3
  #define __Pyx_BUILTIN_MODULE_NAME "__builtin__"
  #define __Pyx_PyCode_New(a, k, l, s, f, code, c, n, v, fv, cell, fn, name, fline, lnos) \
          PyCode_New(a, l, s, f, code, c, n, v, fv, cell, fn, name, fline, lnos)
#else
  #define __Pyx_BUILTIN_MODULE_NAME "builtins"
  #define __Pyx_PyCode_New(a, k, l, s, f, code, c, n, v, fv, cell, fn, name, fline, lnos) \
          PyCode_New(a, k, l, s, f, code, c, n, v, fv, cell, fn, name, fline, lnos)
#endif
#if PY_MAJOR_VERSION < 3 && PY_MINOR_VERSION < 6
  #define PyUnicode_FromString(s) PyUnicode_Decode(s, strlen(s), "UTF-8", "strict")
#endif
#if PY_MAJOR_VERSION >= 3
  #define Py_TPFLAGS_CHECKTYPES 0
  #define Py_TPFLAGS_HAVE_INDEX 0
#endif
#if (PY_VERSION_HEX < 0x02060000) || (PY_MAJOR_VERSION >= 3)
  #define Py_TPFLAGS_HAVE_NEWBUFFER 0
#endif
#if PY_VERSION_HEX > 0x03030000 && defined(PyUnicode_KIND)
  #define CYTHON_PEP393_ENABLED 1
  #define __Pyx_PyUnicode_READY(op)       (likely(PyUnicode_IS_READY(op)) ? \
                                              0 : _PyUnicode_Ready((PyObject *)(op)))
  #define __Pyx_PyUnicode_GET_LENGTH(u)   PyUnicode_GET_LENGTH(u)
  #define __Pyx_PyUnicode_READ_CHAR(u, i) PyUnicode_READ_CHAR(u, i)
  #define __Pyx_PyUnicode_READ(k, d, i)   PyUnicode_READ(k, d, i)
#else
  #define CYTHON_PEP393_ENABLED 0
  #define __Pyx_PyUnicode_READY(op)       (0)
  #define __Pyx_PyUnicode_GET_LENGTH(u)   PyUnicode_GET_SIZE(u)
  #define __Pyx_PyUnicode_READ_CHAR(u, i) ((Py_UCS4)(PyUnicode_AS_UNICODE(u)[i]))
  #define __Pyx_PyUnicode_READ(k, d, i)   ((k=k), (Py_UCS4)(((Py_UNICODE*)d)[i]))
#endif
#if PY_MAJOR_VERSION >= 3
  #define PyBaseString_Type            PyUnicode_Type
  #define PyStringObject               PyUnicodeObject
  #define PyString_Type                PyUnicode_Type
  #define PyString_Check               PyUnicode_Check
  #define PyString_CheckExact          PyUnicode_CheckExact
#endif
#if PY_VERSION_HEX < 0x02060000
  #define PyBytesObject                PyStringObject
  #define PyBytes_Type                 PyString_Type
  #define PyBytes_Check                PyString_Check
  #define PyBytes_CheckExact           PyString_CheckExact
  #define PyBytes_FromString           PyString_FromString
  #define PyBytes_FromStringAndSize    PyString_FromStringAndSize
  #define PyBytes_FromFormat           PyString_FromFormat
  #define PyBytes_DecodeEscape         PyString_DecodeEscape
  #define PyBytes_AsString             PyString_AsString
  #define PyBytes_AsStringAndSize      PyString_AsStringAndSize
  #define PyBytes_Size                 PyString_Size
  #define PyBytes_AS_STRING            PyString_AS_STRING
  #define PyBytes_GET_SIZE             PyString_GET_SIZE
  #define PyBytes_Repr                 PyString_Repr
  #define PyBytes_Concat               PyString_Concat
  #define PyBytes_ConcatAndDel         PyString_ConcatAndDel
#endif
#if PY_VERSION_HEX < 0x02060000
  #define PySet_Check(obj)             PyObject_TypeCheck(obj, &PySet_Type)
  #define PyFrozenSet_Check(obj)       PyObject_TypeCheck(obj, &PyFrozenSet_Type)
#endif
#ifndef PySet_CheckExact
  #define PySet_CheckExact(obj)        (Py_TYPE(obj) == &PySet_Type)
#endif
#define __Pyx_TypeCheck(obj, type) PyObject_TypeCheck(obj, (PyTypeObject *)type)
#if PY_MAJOR_VERSION >= 3
  #define PyIntObject                  PyLongObject
  #define PyInt_Type                   PyLong_Type
  #define PyInt_Check(op)              PyLong_Check(op)
  #define PyInt_CheckExact(op)         PyLong_CheckExact(op)
  #define PyInt_FromString             PyLong_FromString
  #define PyInt_FromUnicode            PyLong_FromUnicode
  #define PyInt_FromLong               PyLong_FromLong
  #define PyInt_FromSize_t             PyLong_FromSize_t
  #define PyInt_FromSsize_t            PyLong_FromSsize_t
  #define PyInt_AsLong                 PyLong_AsLong
  #define PyInt_AS_LONG                PyLong_AS_LONG
  #define PyInt_AsSsize_t              PyLong_AsSsize_t
  #define PyInt_AsUnsignedLongMask     PyLong_AsUnsignedLongMask
  #define PyInt_AsUnsignedLongLongMask PyLong_AsUnsignedLongLongMask
#endif
#if PY_MAJOR_VERSION >= 3
  #define PyBoolObject                 PyLongObject
#endif
#if PY_VERSION_HEX < 0x03020000
  typedef long Py_hash_t;
  #define __Pyx_PyInt_FromHash_t PyInt_FromLong
  #define __Pyx_PyInt_AsHash_t   PyInt_AsLong
#else
  #define __Pyx_PyInt_FromHash_t PyInt_FromSsize_t
  #define __Pyx_PyInt_AsHash_t   PyInt_AsSsize_t
#endif
#if (PY_MAJOR_VERSION < 3) || (PY_VERSION_HEX >= 0x03010300)
  #define __Pyx_PySequence_GetSlice(obj, a, b) PySequence_GetSlice(obj, a, b)
  #define __Pyx_PySequence_SetSlice(obj, a, b, value) PySequence_SetSlice(obj, a, b, value)
  #define __Pyx_PySequence_DelSlice(obj, a, b) PySequence_DelSlice(obj, a, b)
#else
  #define __Pyx_PySequence_GetSlice(obj, a, b) (unlikely(!(obj)) ? \
        (PyErr_SetString(PyExc_SystemError, "null argument to internal routine"), (PyObject*)0) : \
        (likely((obj)->ob_type->tp_as_mapping) ? (PySequence_GetSlice(obj, a, b)) : \
            (PyErr_Format(PyExc_TypeError, "'%.200s' object is unsliceable", (obj)->ob_type->tp_name), (PyObject*)0)))
  #define __Pyx_PySequence_SetSlice(obj, a, b, value) (unlikely(!(obj)) ? \
        (PyErr_SetString(PyExc_SystemError, "null argument to internal routine"), -1) : \
        (likely((obj)->ob_type->tp_as_mapping) ? (PySequence_SetSlice(obj, a, b, value)) : \
            (PyErr_Format(PyExc_TypeError, "'%.200s' object doesn't support slice assignment", (obj)->ob_type->tp_name), -1)))
  #define __Pyx_PySequence_DelSlice(obj, a, b) (unlikely(!(obj)) ? \
        (PyErr_SetString(PyExc_SystemError, "null argument to internal routine"), -1) : \
        (likely((obj)->ob_type->tp_as_mapping) ? (PySequence_DelSlice(obj, a, b)) : \
            (PyErr_Format(PyExc_TypeError, "'%.200s' object doesn't support slice deletion", (obj)->ob_type->tp_name), -1)))
#endif
#if PY_MAJOR_VERSION >= 3
  #define PyMethod_New(func, self, klass) ((self) ? PyMethod_New(func, self) : PyInstanceMethod_New(func))
#endif
#if PY_VERSION_HEX < 0x02050000
  #define __Pyx_GetAttrString(o,n)   PyObject_GetAttrString((o),((char *)(n)))
  #define __Pyx_SetAttrString(o,n,a) PyObject_SetAttrString((o),((char *)(n)),(a))
  #define __Pyx_DelAttrString(o,n)   PyObject_DelAttrString((o),((char *)(n)))
#else
  #define __Pyx_GetAttrString(o,n)   PyObject_GetAttrString((o),(n))
  #define __Pyx_SetAttrString(o,n,a) PyObject_SetAttrString((o),(n),(a))
  #define __Pyx_DelAttrString(o,n)   PyObject_DelAttrString((o),(n))
#endif
#if PY_VERSION_HEX < 0x02050000
  #define __Pyx_NAMESTR(n) ((char *)(n))
  #define __Pyx_DOCSTR(n)  ((char *)(n))
#else
  #define __Pyx_NAMESTR(n) (n)
  #define __Pyx_DOCSTR(n)  (n)
#endif


#if PY_MAJOR_VERSION >= 3
  #define __Pyx_PyNumber_Divide(x,y)         PyNumber_TrueDivide(x,y)
  #define __Pyx_PyNumber_InPlaceDivide(x,y)  PyNumber_InPlaceTrueDivide(x,y)
#else
  #define __Pyx_PyNumber_Divide(x,y)         PyNumber_Divide(x,y)
  #define __Pyx_PyNumber_InPlaceDivide(x,y)  PyNumber_InPlaceDivide(x,y)
#endif

#ifndef __PYX_EXTERN_C
  #ifdef __cplusplus
    #define __PYX_EXTERN_C extern "C"
  #else
    #define __PYX_EXTERN_C extern
  #endif
#endif

#if defined(WIN32) || defined(MS_WINDOWS)
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#define __PYX_HAVE__scipy__special___ufuncs
#define __PYX_HAVE_API__scipy__special___ufuncs
#include "_complexstuff.h"
#include "numpy/npy_math.h"
#include "stdio.h"
#include "stdlib.h"
#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"
#include "sf_error.h"
#include "math.h"
#include "_ufuncs_defs.h"
#include "cephes.h"
#include "specfun_wrappers.h"
#include "c_misc/misc.h"
#ifdef _OPENMP
#include <omp.h>
#endif 

#ifdef PYREX_WITHOUT_ASSERTIONS
#define CYTHON_WITHOUT_ASSERTIONS
#endif



#ifndef CYTHON_INLINE
  #if defined(__GNUC__)
    #define CYTHON_INLINE __inline__
  #elif defined(_MSC_VER)
    #define CYTHON_INLINE __inline
  #elif defined (__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
    #define CYTHON_INLINE inline
  #else
    #define CYTHON_INLINE
  #endif
#endif


#ifndef CYTHON_UNUSED
# if defined(__GNUC__)
#   if !(defined(__cplusplus)) || (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 4))
#     define CYTHON_UNUSED __attribute__ ((__unused__))
#   else
#     define CYTHON_UNUSED
#   endif
# elif defined(__ICC) || (defined(__INTEL_COMPILER) && !defined(_MSC_VER))
#   define CYTHON_UNUSED __attribute__ ((__unused__))
# else
#   define CYTHON_UNUSED
# endif
#endif

typedef struct {PyObject **p; char *s; const long n; const char* encoding; const char is_unicode; const char is_str; const char intern; } __Pyx_StringTabEntry; 




#define __Pyx_PyBytes_FromUString(s) PyBytes_FromString((char*)s)
#define __Pyx_PyBytes_AsUString(s)   ((unsigned char*) PyBytes_AsString(s))

#define __Pyx_Owned_Py_None(b) (Py_INCREF(Py_None), Py_None)
#define __Pyx_PyBool_FromLong(b) ((b) ? (Py_INCREF(Py_True), Py_True) : (Py_INCREF(Py_False), Py_False))
static CYTHON_INLINE int __Pyx_PyObject_IsTrue(PyObject*);
static CYTHON_INLINE PyObject* __Pyx_PyNumber_Int(PyObject* x);

static CYTHON_INLINE Py_ssize_t __Pyx_PyIndex_AsSsize_t(PyObject*);
static CYTHON_INLINE PyObject * __Pyx_PyInt_FromSize_t(size_t);
static CYTHON_INLINE size_t __Pyx_PyInt_AsSize_t(PyObject*);

#if CYTHON_COMPILING_IN_CPYTHON
#define __pyx_PyFloat_AsDouble(x) (PyFloat_CheckExact(x) ? PyFloat_AS_DOUBLE(x) : PyFloat_AsDouble(x))
#else
#define __pyx_PyFloat_AsDouble(x) PyFloat_AsDouble(x)
#endif
#define __pyx_PyFloat_AsFloat(x) ((float) __pyx_PyFloat_AsDouble(x))

#ifdef __GNUC__
  
  #if __GNUC__ > 2 || (__GNUC__ == 2 && (__GNUC_MINOR__ > 95))
    #define likely(x)   __builtin_expect(!!(x), 1)
    #define unlikely(x) __builtin_expect(!!(x), 0)
  #else 
    #define likely(x)   (x)
    #define unlikely(x) (x)
  #endif 
#else 
  #define likely(x)   (x)
  #define unlikely(x) (x)
#endif 
    
static PyObject *__pyx_m;
static PyObject *__pyx_b;
static PyObject *__pyx_empty_tuple;
static PyObject *__pyx_empty_bytes;
static int __pyx_lineno;
static int __pyx_clineno = 0;
static const char * __pyx_cfilenm= __FILE__;
static const char *__pyx_filename;

#if !defined(CYTHON_CCOMPLEX)
  #if defined(__cplusplus)
    #define CYTHON_CCOMPLEX 1
  #elif defined(_Complex_I)
    #define CYTHON_CCOMPLEX 1
  #else
    #define CYTHON_CCOMPLEX 0
  #endif
#endif
#if CYTHON_CCOMPLEX
  #ifdef __cplusplus
    #include <complex>
  #else
    #include <complex.h>
  #endif
#endif
#if CYTHON_CCOMPLEX && !defined(__cplusplus) && defined(__sun__) && defined(__GNUC__)
  #undef _Complex_I
  #define _Complex_I 1.0fj
#endif


static const char *__pyx_f[] = {
  "_ufuncs.pyx",
  "_legacy.pxd",
  "numpy.pxd",
  "type.pxd",
};


typedef npy_int8 __pyx_t_5numpy_int8_t;


typedef npy_int16 __pyx_t_5numpy_int16_t;


typedef npy_int32 __pyx_t_5numpy_int32_t;


typedef npy_int64 __pyx_t_5numpy_int64_t;


typedef npy_uint8 __pyx_t_5numpy_uint8_t;


typedef npy_uint16 __pyx_t_5numpy_uint16_t;


typedef npy_uint32 __pyx_t_5numpy_uint32_t;


typedef npy_uint64 __pyx_t_5numpy_uint64_t;


typedef npy_float32 __pyx_t_5numpy_float32_t;


typedef npy_float64 __pyx_t_5numpy_float64_t;


typedef npy_long __pyx_t_5numpy_int_t;


typedef npy_longlong __pyx_t_5numpy_long_t;


typedef npy_longlong __pyx_t_5numpy_longlong_t;


typedef npy_ulong __pyx_t_5numpy_uint_t;


typedef npy_ulonglong __pyx_t_5numpy_ulong_t;


typedef npy_ulonglong __pyx_t_5numpy_ulonglong_t;


typedef npy_intp __pyx_t_5numpy_intp_t;


typedef npy_uintp __pyx_t_5numpy_uintp_t;


typedef npy_double __pyx_t_5numpy_float_t;


typedef npy_double __pyx_t_5numpy_double_t;


typedef npy_longdouble __pyx_t_5numpy_longdouble_t;
#if CYTHON_CCOMPLEX
  #ifdef __cplusplus
    typedef ::std::complex< double > __pyx_t_double_complex;
  #else
    typedef double _Complex __pyx_t_double_complex;
  #endif
#else
    typedef struct { double real, imag; } __pyx_t_double_complex;
#endif

#if CYTHON_CCOMPLEX
  #ifdef __cplusplus
    typedef ::std::complex< float > __pyx_t_float_complex;
  #else
    typedef float _Complex __pyx_t_float_complex;
  #endif
#else
    typedef struct { float real, imag; } __pyx_t_float_complex;
#endif





typedef npy_cfloat __pyx_t_5numpy_cfloat_t;


typedef npy_cdouble __pyx_t_5numpy_cdouble_t;


typedef npy_clongdouble __pyx_t_5numpy_clongdouble_t;


typedef npy_cdouble __pyx_t_5numpy_complex_t;


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_lambertw_scalar_t(__pyx_t_double_complex, long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_bdtr_unsafe_t(double, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_bdtrc_unsafe_t(double, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_bdtri_unsafe_t(double, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_binom_t(double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyc_double__t(double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyc_double_complex__t(double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyc_l_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebys_double__t(double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebys_double_complex__t(double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebys_l_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyt_double__t(double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyt_double_complex__t(double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyt_l_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyu_double__t(double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyu_double_complex__t(double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyu_l_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_double__t(double, double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_double_complex__t(double, double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_l_t(long, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_double__t(double, double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_double_complex__t(double, double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_l_t(long, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_hermite_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_hermitenorm_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_jacobi_double__t(double, double, double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_jacobi_double_complex__t(double, double, double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_jacobi_l_t(long, double, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_laguerre_double__t(double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_laguerre_double_complex__t(double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_laguerre_l_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_legendre_double__t(double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_legendre_double_complex__t(double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_legendre_l_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_double__t(double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_double_complex__t(double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_l_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_double__t(double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_double_complex__t(double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_l_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_double__t(double, double, double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_double_complex__t(double, double, double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_l_t(long, double, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_double__t(double, double);


typedef __pyx_t_double_complex __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_double_complex__t(double, __pyx_t_double_complex);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_l_t(long, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_expn_unsafe_t(double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_hyp2f0_unsafe_t(double, double, double, double, double *);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_kn_unsafe_t(double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_nbdtr_unsafe_t(double, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_nbdtrc_unsafe_t(double, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_nbdtri_unsafe_t(double, double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_pdtr_unsafe_t(double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_pdtrc_unsafe_t(double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_pdtri_unsafe_t(double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_smirnov_unsafe_t(double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_smirnovi_unsafe_t(double, double);


typedef double __pyx_t_5scipy_7special_7_ufuncs__proto_yn_unsafe_t(double, double);
#ifndef CYTHON_REFNANNY
  #define CYTHON_REFNANNY 0
#endif
#if CYTHON_REFNANNY
  typedef struct {
    void (*INCREF)(void*, PyObject*, int);
    void (*DECREF)(void*, PyObject*, int);
    void (*GOTREF)(void*, PyObject*, int);
    void (*GIVEREF)(void*, PyObject*, int);
    void* (*SetupContext)(const char*, int, const char*);
    void (*FinishContext)(void**);
  } __Pyx_RefNannyAPIStruct;
  static __Pyx_RefNannyAPIStruct *__Pyx_RefNanny = NULL;
  static __Pyx_RefNannyAPIStruct *__Pyx_RefNannyImportAPI(const char *modname); 
  #define __Pyx_RefNannyDeclarations void *__pyx_refnanny = NULL;
#ifdef WITH_THREAD
  #define __Pyx_RefNannySetupContext(name, acquire_gil) \
          if (acquire_gil) { \
              PyGILState_STATE __pyx_gilstate_save = PyGILState_Ensure(); \
              __pyx_refnanny = __Pyx_RefNanny->SetupContext((name), __LINE__, __FILE__); \
              PyGILState_Release(__pyx_gilstate_save); \
          } else { \
              __pyx_refnanny = __Pyx_RefNanny->SetupContext((name), __LINE__, __FILE__); \
          }
#else
  #define __Pyx_RefNannySetupContext(name, acquire_gil) \
          __pyx_refnanny = __Pyx_RefNanny->SetupContext((name), __LINE__, __FILE__)
#endif
  #define __Pyx_RefNannyFinishContext() \
          __Pyx_RefNanny->FinishContext(&__pyx_refnanny)
  #define __Pyx_INCREF(r)  __Pyx_RefNanny->INCREF(__pyx_refnanny, (PyObject *)(r), __LINE__)
  #define __Pyx_DECREF(r)  __Pyx_RefNanny->DECREF(__pyx_refnanny, (PyObject *)(r), __LINE__)
  #define __Pyx_GOTREF(r)  __Pyx_RefNanny->GOTREF(__pyx_refnanny, (PyObject *)(r), __LINE__)
  #define __Pyx_GIVEREF(r) __Pyx_RefNanny->GIVEREF(__pyx_refnanny, (PyObject *)(r), __LINE__)
  #define __Pyx_XINCREF(r)  do { if((r) != NULL) {__Pyx_INCREF(r); }} while(0)
  #define __Pyx_XDECREF(r)  do { if((r) != NULL) {__Pyx_DECREF(r); }} while(0)
  #define __Pyx_XGOTREF(r)  do { if((r) != NULL) {__Pyx_GOTREF(r); }} while(0)
  #define __Pyx_XGIVEREF(r) do { if((r) != NULL) {__Pyx_GIVEREF(r);}} while(0)
#else
  #define __Pyx_RefNannyDeclarations
  #define __Pyx_RefNannySetupContext(name, acquire_gil)
  #define __Pyx_RefNannyFinishContext()
  #define __Pyx_INCREF(r) Py_INCREF(r)
  #define __Pyx_DECREF(r) Py_DECREF(r)
  #define __Pyx_GOTREF(r)
  #define __Pyx_GIVEREF(r)
  #define __Pyx_XINCREF(r) Py_XINCREF(r)
  #define __Pyx_XDECREF(r) Py_XDECREF(r)
  #define __Pyx_XGOTREF(r)
  #define __Pyx_XGIVEREF(r)
#endif 
#define __Pyx_CLEAR(r)    do { PyObject* tmp = ((PyObject*)(r)); r = NULL; __Pyx_DECREF(tmp);} while(0)
#define __Pyx_XCLEAR(r)   do { if((r) != NULL) {PyObject* tmp = ((PyObject*)(r)); r = NULL; __Pyx_DECREF(tmp);}} while(0)

static PyObject *__Pyx_GetName(PyObject *dict, PyObject *name); 

static void __Pyx_RaiseDoubleKeywordsError(const char* func_name, PyObject* kw_name); 

static int __Pyx_ParseOptionalKeywords(PyObject *kwds, PyObject **argnames[], \
    PyObject *kwds2, PyObject *values[], Py_ssize_t num_pos_args, \
    const char* function_name); 

static void __Pyx_RaiseArgtupleInvalid(const char* func_name, int exact,
    Py_ssize_t num_min, Py_ssize_t num_max, Py_ssize_t num_found); 

#ifndef __PYX_FORCE_INIT_THREADS
  #define __PYX_FORCE_INIT_THREADS 0
#endif

static CYTHON_INLINE long __Pyx_mod_long(long, long); 

static CYTHON_INLINE long __Pyx_div_long(long, long); 

static CYTHON_INLINE void __Pyx_ErrRestore(PyObject *type, PyObject *value, PyObject *tb); 
static CYTHON_INLINE void __Pyx_ErrFetch(PyObject **type, PyObject **value, PyObject **tb); 

static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb, PyObject *cause); 

static CYTHON_INLINE void __Pyx_RaiseTooManyValuesError(Py_ssize_t expected);

static CYTHON_INLINE void __Pyx_RaiseNeedMoreValuesError(Py_ssize_t index);

static CYTHON_INLINE void __Pyx_RaiseNoneNotIterableError(void);

static CYTHON_INLINE int __Pyx_IterFinish(void); 

static int __Pyx_IternextUnpackEndCheck(PyObject *retval, Py_ssize_t expected); 

static CYTHON_INLINE int __Pyx_TypeTest(PyObject *obj, PyTypeObject *type); 

#if CYTHON_CCOMPLEX
  #ifdef __cplusplus
    #define __Pyx_CREAL(z) ((z).real())
    #define __Pyx_CIMAG(z) ((z).imag())
  #else
    #define __Pyx_CREAL(z) (__real__(z))
    #define __Pyx_CIMAG(z) (__imag__(z))
  #endif
#else
    #define __Pyx_CREAL(z) ((z).real)
    #define __Pyx_CIMAG(z) ((z).imag)
#endif
#if defined(_WIN32) && defined(__cplusplus) && CYTHON_CCOMPLEX
    #define __Pyx_SET_CREAL(z,x) ((z).real(x))
    #define __Pyx_SET_CIMAG(z,y) ((z).imag(y))
#else
    #define __Pyx_SET_CREAL(z,x) __Pyx_CREAL(z) = (x)
    #define __Pyx_SET_CIMAG(z,y) __Pyx_CIMAG(z) = (y)
#endif

static CYTHON_INLINE __pyx_t_double_complex __pyx_t_double_complex_from_parts(double, double);

#if CYTHON_CCOMPLEX
    #define __Pyx_c_eq(a, b)   ((a)==(b))
    #define __Pyx_c_sum(a, b)  ((a)+(b))
    #define __Pyx_c_diff(a, b) ((a)-(b))
    #define __Pyx_c_prod(a, b) ((a)*(b))
    #define __Pyx_c_quot(a, b) ((a)/(b))
    #define __Pyx_c_neg(a)     (-(a))
  #ifdef __cplusplus
    #define __Pyx_c_is_zero(z) ((z)==(double)0)
    #define __Pyx_c_conj(z)    (::std::conj(z))
    #if 1
        #define __Pyx_c_abs(z)     (::std::abs(z))
        #define __Pyx_c_pow(a, b)  (::std::pow(a, b))
    #endif
  #else
    #define __Pyx_c_is_zero(z) ((z)==0)
    #define __Pyx_c_conj(z)    (conj(z))
    #if 1
        #define __Pyx_c_abs(z)     (cabs(z))
        #define __Pyx_c_pow(a, b)  (cpow(a, b))
    #endif
 #endif
#else
    static CYTHON_INLINE int __Pyx_c_eq(__pyx_t_double_complex, __pyx_t_double_complex);
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_sum(__pyx_t_double_complex, __pyx_t_double_complex);
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_diff(__pyx_t_double_complex, __pyx_t_double_complex);
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_prod(__pyx_t_double_complex, __pyx_t_double_complex);
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_quot(__pyx_t_double_complex, __pyx_t_double_complex);
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_neg(__pyx_t_double_complex);
    static CYTHON_INLINE int __Pyx_c_is_zero(__pyx_t_double_complex);
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_conj(__pyx_t_double_complex);
    #if 1
        static CYTHON_INLINE double __Pyx_c_abs(__pyx_t_double_complex);
        static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_pow(__pyx_t_double_complex, __pyx_t_double_complex);
    #endif
#endif

static CYTHON_INLINE PyObject *__Pyx_PyInt_to_py_Py_intptr_t(Py_intptr_t);

static CYTHON_INLINE Py_intptr_t __Pyx_PyInt_from_py_Py_intptr_t(PyObject *);

static CYTHON_INLINE __pyx_t_float_complex __pyx_t_float_complex_from_parts(float, float);

#if CYTHON_CCOMPLEX
    #define __Pyx_c_eqf(a, b)   ((a)==(b))
    #define __Pyx_c_sumf(a, b)  ((a)+(b))
    #define __Pyx_c_difff(a, b) ((a)-(b))
    #define __Pyx_c_prodf(a, b) ((a)*(b))
    #define __Pyx_c_quotf(a, b) ((a)/(b))
    #define __Pyx_c_negf(a)     (-(a))
  #ifdef __cplusplus
    #define __Pyx_c_is_zerof(z) ((z)==(float)0)
    #define __Pyx_c_conjf(z)    (::std::conj(z))
    #if 1
        #define __Pyx_c_absf(z)     (::std::abs(z))
        #define __Pyx_c_powf(a, b)  (::std::pow(a, b))
    #endif
  #else
    #define __Pyx_c_is_zerof(z) ((z)==0)
    #define __Pyx_c_conjf(z)    (conjf(z))
    #if 1
        #define __Pyx_c_absf(z)     (cabsf(z))
        #define __Pyx_c_powf(a, b)  (cpowf(a, b))
    #endif
 #endif
#else
    static CYTHON_INLINE int __Pyx_c_eqf(__pyx_t_float_complex, __pyx_t_float_complex);
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_sumf(__pyx_t_float_complex, __pyx_t_float_complex);
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_difff(__pyx_t_float_complex, __pyx_t_float_complex);
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_prodf(__pyx_t_float_complex, __pyx_t_float_complex);
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_quotf(__pyx_t_float_complex, __pyx_t_float_complex);
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_negf(__pyx_t_float_complex);
    static CYTHON_INLINE int __Pyx_c_is_zerof(__pyx_t_float_complex);
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_conjf(__pyx_t_float_complex);
    #if 1
        static CYTHON_INLINE float __Pyx_c_absf(__pyx_t_float_complex);
        static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_powf(__pyx_t_float_complex, __pyx_t_float_complex);
    #endif
#endif

static CYTHON_INLINE long __Pyx_pow_long(long, long); 

static CYTHON_INLINE unsigned char __Pyx_PyInt_AsUnsignedChar(PyObject *);

static CYTHON_INLINE unsigned short __Pyx_PyInt_AsUnsignedShort(PyObject *);

static CYTHON_INLINE unsigned int __Pyx_PyInt_AsUnsignedInt(PyObject *);

static CYTHON_INLINE char __Pyx_PyInt_AsChar(PyObject *);

static CYTHON_INLINE short __Pyx_PyInt_AsShort(PyObject *);

static CYTHON_INLINE int __Pyx_PyInt_AsInt(PyObject *);

static CYTHON_INLINE signed char __Pyx_PyInt_AsSignedChar(PyObject *);

static CYTHON_INLINE signed short __Pyx_PyInt_AsSignedShort(PyObject *);

static CYTHON_INLINE signed int __Pyx_PyInt_AsSignedInt(PyObject *);

static CYTHON_INLINE int __Pyx_PyInt_AsLongDouble(PyObject *);

static CYTHON_INLINE unsigned long __Pyx_PyInt_AsUnsignedLong(PyObject *);

static CYTHON_INLINE unsigned PY_LONG_LONG __Pyx_PyInt_AsUnsignedLongLong(PyObject *);

static CYTHON_INLINE long __Pyx_PyInt_AsLong(PyObject *);

static CYTHON_INLINE PY_LONG_LONG __Pyx_PyInt_AsLongLong(PyObject *);

static CYTHON_INLINE signed long __Pyx_PyInt_AsSignedLong(PyObject *);

static CYTHON_INLINE signed PY_LONG_LONG __Pyx_PyInt_AsSignedLongLong(PyObject *);

static void __Pyx_WriteUnraisable(const char *name, int clineno,
                                  int lineno, const char *filename); 

static int __Pyx_check_binary_version(void);

#if !defined(__Pyx_PyIdentifier_FromString)
#if PY_MAJOR_VERSION < 3
  #define __Pyx_PyIdentifier_FromString(s) PyString_FromString(s)
#else
  #define __Pyx_PyIdentifier_FromString(s) PyUnicode_FromString(s)
#endif
#endif

static PyObject *__Pyx_ImportModule(const char *name); 

static PyTypeObject *__Pyx_ImportType(const char *module_name, const char *class_name, size_t size, int strict);  

typedef struct {
    int code_line;
    PyCodeObject* code_object;
} __Pyx_CodeObjectCacheEntry;
struct __Pyx_CodeObjectCache {
    int count;
    int max_count;
    __Pyx_CodeObjectCacheEntry* entries;
};
static struct __Pyx_CodeObjectCache __pyx_code_cache = {0,0,NULL};
static int __pyx_bisect_code_objects(__Pyx_CodeObjectCacheEntry* entries, int count, int code_line);
static PyCodeObject *__pyx_find_code_object(int code_line);
static void __pyx_insert_code_object(int code_line, PyCodeObject* code_object);

static void __Pyx_AddTraceback(const char *funcname, int c_line,
                               int py_line, const char *filename); 

static int __Pyx_InitStrings(__Pyx_StringTabEntry *t); 













static PyTypeObject *__pyx_ptype_7cpython_4type_type = 0;






static PyTypeObject *__pyx_ptype_5numpy_dtype = 0;
static PyTypeObject *__pyx_ptype_5numpy_flatiter = 0;
static PyTypeObject *__pyx_ptype_5numpy_broadcast = 0;
static PyTypeObject *__pyx_ptype_5numpy_ndarray = 0;
static PyTypeObject *__pyx_ptype_5numpy_ufunc = 0;
static CYTHON_INLINE char *__pyx_f_5numpy__util_dtypestring(PyArray_Descr *, char *, char *, int *); 








static CYTHON_INLINE int __pyx_f_5scipy_7special_13_complexstuff_zisnan(__pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_13_complexstuff_zabs(__pyx_t_double_complex); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_f_5scipy_7special_13_complexstuff_zlog(__pyx_t_double_complex); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_f_5scipy_7special_13_complexstuff_zexp(__pyx_t_double_complex); 


static CYTHON_INLINE __pyx_t_double_complex __pyx_f_5scipy_7special_8lambertw_lambertw_scalar(__pyx_t_double_complex, long, double); 


static CYTHON_INLINE void __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(char *, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_bdtrc_unsafe(double, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_bdtr_unsafe(double, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_bdtri_unsafe(double, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_expn_unsafe(double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_hyp2f0_unsafe(double, double, double, double, double *); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_nbdtrc_unsafe(double, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_nbdtr_unsafe(double, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_nbdtri_unsafe(double, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_pdtrc_unsafe(double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_pdtr_unsafe(double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_pdtri_unsafe(double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_kn_unsafe(double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_yn_unsafe(double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_smirnov_unsafe(double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_smirnovi_unsafe(double, double); 




static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_binom(double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi_l(long, double, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi_l(long, double, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer_l(long, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt_l(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu_l(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys_l(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc_l(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt_l(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu_l(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre_l(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre_l(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre_l(long, double, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre_l(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_hermite(long, double); 
static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_hermitenorm(long, double); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(double, double, double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(double, double, double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp1f1(double, double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp1f1(double, double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi(double, double, double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi(double, double, double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi(double, double, double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi(double, double, double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer(double, double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer(double, double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt(double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt(double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu(double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu(double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys(double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys(double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc(double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc(double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt(double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt(double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu(double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu(double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre(double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre(double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre(double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre(double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre(double, double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre(double, double, __pyx_t_double_complex); 
static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre(double, double); 
static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre(double, __pyx_t_double_complex); 


static __pyx_t_5scipy_7special_7_ufuncs__proto_lambertw_scalar_t *__pyx_v_5scipy_7special_7_ufuncs__proto_lambertw_scalar_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_bdtr_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_bdtr_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_bdtrc_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_bdtrc_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_bdtri_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_bdtri_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_binom_t *__pyx_v_5scipy_7special_7_ufuncs__proto_binom_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyc_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyc_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyc_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyc_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyc_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyc_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebys_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebys_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebys_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebys_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebys_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebys_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyt_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyt_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyt_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyt_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyt_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyt_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyu_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyu_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyu_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyu_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_chebyu_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyu_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_hermite_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_hermite_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_hermitenorm_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_hermitenorm_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_jacobi_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_jacobi_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_jacobi_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_jacobi_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_jacobi_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_jacobi_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_laguerre_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_laguerre_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_laguerre_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_laguerre_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_laguerre_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_laguerre_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_legendre_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_legendre_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_legendre_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_legendre_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_legendre_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_legendre_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_double__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_double__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_double_complex__t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_double_complex__t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_l_t *__pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_l_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_expn_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_expn_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_hyp2f0_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_hyp2f0_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_kn_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_kn_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_nbdtr_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_nbdtr_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_nbdtrc_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_nbdtrc_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_nbdtri_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_nbdtri_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_pdtr_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_pdtr_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_pdtrc_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_pdtrc_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_pdtri_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_pdtri_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_smirnov_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_smirnov_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_smirnovi_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_smirnovi_unsafe_t_var;
static __pyx_t_5scipy_7special_7_ufuncs__proto_yn_unsafe_t *__pyx_v_5scipy_7special_7_ufuncs__proto_yn_unsafe_t_var;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[20];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[20];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[16];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[16];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[16];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[18];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[18];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[18];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[18];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[24];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[24];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[30];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[18];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[18];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[18];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[18];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[30];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_loops[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[12];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_data[6];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[18];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_loops[3];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_data[3];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[16];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[24];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[20];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_loops[3];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[6];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_data[3];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[16];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[16];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[16];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[10];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[14];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[14];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[14];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[14];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[14];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[14];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types[8];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_round_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_round_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_loops[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[8];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_data[4];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[12];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_types[6];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_doc;
static PyUFuncGenericFunction __pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_loops[2];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_ptr[4];
static void *__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_data[2];
static char __pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_types[4];
static char *__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_doc;
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_f_f__As_f_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_ddddd_dd(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_d_dd(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_D_D(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddi_d_As_fffl_ff(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DD_As_F_FF(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DDDD_As_f_FFFF(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_f_ffff(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddi_d_As_dddl_dd(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ddD_D(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dddd_As_ff_ffff(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_dd_dd(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_fff_ff(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ldd__As_ldd_d(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_dddD_D(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_d_dddd(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_fffff_ff(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DDDD_As_d_DDDD(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_ddd_dd(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_Dld__As_Dld_D(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_ff_ff(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_fffF_F(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_lld_d(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DD_As_d_DD(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_F_F(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_llf_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_lddd__As_lddd_d(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ffF_F(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DDDD_As_F_FFFF(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_lddd__As_lfff_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DD_As_D_DD(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_Dld__As_Flf_F(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dddd_As_dd_dddd(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_f_ff(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DD_As_f_FF(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DDDD_As_D_DDDD(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_g_g__As_g_g(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ldd__As_lff_f(char **, npy_intp *, npy_intp *, void *); 
static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d(char **, npy_intp *, npy_intp *, void *); 
__PYX_EXTERN_C DL_EXPORT(int) wrap_PyUFunc_getfperr(void); 
#define __Pyx_MODULE_NAME "scipy.special._ufuncs"
int __pyx_module_is_main_scipy__special___ufuncs = 0;


static PyObject *__pyx_builtin_range;
static PyObject *__pyx_builtin_FutureWarning;
static PyObject *__pyx_builtin_ValueError;
static PyObject *__pyx_builtin_RuntimeError;
static PyObject *__pyx_pf_5scipy_7special_7_ufuncs__errprint(CYTHON_UNUSED PyObject *__pyx_self, PyObject *__pyx_v_inflag); 
static int __pyx_pf_5numpy_7ndarray___getbuffer__(PyArrayObject *__pyx_v_self, Py_buffer *__pyx_v_info, int __pyx_v_flags); 
static void __pyx_pf_5numpy_7ndarray_2__releasebuffer__(PyArrayObject *__pyx_v_self, Py_buffer *__pyx_v_info); 
static char __pyx_k_1[] = "invalid input argument";
static char __pyx_k_2[] = "scipy.special.";
static char __pyx_k_3[] = ": floating point number truncated to an integer ";
static char __pyx_k_4[] = "(this behavior changes in a future Scipy release)";
static char __pyx_k_5[] = "iteration failed to converge: %g + %gj";
static char __pyx_k_6[] = "ndarray is not C contiguous";
static char __pyx_k_8[] = "ndarray is not Fortran contiguous";
static char __pyx_k_10[] = "Non-native byte order not supported";
static char __pyx_k_12[] = "unknown dtype code in numpy.pxd (%d)";
static char __pyx_k_13[] = "Format string allocated too short, see comment in numpy.pxd";
static char __pyx_k_16[] = "Format string allocated too short.";
static char __pyx_k_18[] = "Internal function, use `lambertw` instead.";
static char __pyx_k_19[] = "(Ai,Aip,Bi,Bip)=airy(z) calculates the Airy functions and their derivatives\nevaluated at real or complex number z.  The Airy functions Ai and Bi\nare two independent solutions of y''(x)=xy.  Aip and Bip are the first derivatives\nevaluated at x of Ai and Bi respectively.";
static char __pyx_k_20[] = "(Aie,Aipe,Bie,Bipe)=airye(z) calculates the exponentially scaled Airy functions and\ntheir derivatives evaluated at real or complex number z.\nairye(z)[0:1] = airy(z)[0:1] * exp(2.0/3.0*z*sqrt(z))\nairye(z)[2:3] = airy(z)[2:3] * exp(-abs((2.0/3.0*z*sqrt(z)).real))";
static char __pyx_k_21[] = "y=bdtr(k,n,p) returns the sum of the terms 0 through k of the\nBinomial probability density:  sum(nCj p**j (1-p)**(n-j),j=0..k)";
static char __pyx_k_22[] = "y=bdtrc(k,n,p) returns the sum of the terms k+1 through n of the\nBinomial probability density: sum(nCj p**j (1-p)**(n-j), j=k+1..n)";
static char __pyx_k_23[] = "p=bdtri(k,n,y) finds the probability p such that the sum of the\nterms 0 through k of the Binomial probability density is equal to the\ngiven cumulative probability y.";
static char __pyx_k_24[] = "";
static char __pyx_k_25[] = "y=bei(x) returns the Kelvin function bei x";
static char __pyx_k_26[] = "y=beip(x) returns the derivative of the Kelvin function bei x";
static char __pyx_k_27[] = "y=ber(x) returns the Kelvin function ber x";
static char __pyx_k_28[] = "y=berp(x) returns the derivative of the Kelvin function ber x";
static char __pyx_k_29[] = "y=besselpoly(a,lam,nu) returns the value of the integral:\nintegral(x**lam * jv(nu,2*a*x),x=0..1).";
static char __pyx_k_30[] = "y=beta(a,b) returns gamma(a) * gamma(b) / gamma(a+b)";
static char __pyx_k_31[] = "betainc(a, b, x)\n\nCompute the incomplete beta integral of the arguments, evaluated\nfrom zero to x::\n\n    gamma(a+b) / (gamma(a)*gamma(b)) * integral(t**(a-1) (1-t)**(b-1), t=0..x).\n\nNotes\n-----\nThe incomplete beta is also sometimes defined without the terms\nin gamma, in which case the above definition is the so-called regularized\nincomplete beta. Under this definition, you can get the incomplete beta by\nmultiplying the result of the scipy function by beta(a, b).";
static char __pyx_k_32[] = "betaincinv(a,b,y)\n\nCompute x such that betainc(a,b,x) = y.";
static char __pyx_k_33[] = "y=betaln(a,b) returns the natural logarithm of the absolute value of\nbeta: ln(abs(beta(x))).";
static char __pyx_k_34[] = "binom(n, k)\n\nBinomial coefficient";
static char __pyx_k_35[] = "y=btdtr(a,b,x) returns the area from zero to x under the beta\ndensity function: gamma(a+b)/(gamma(a)*gamma(b)))*integral(t**(a-1)\n(1-t)**(b-1), t=0..x).  SEE ALSO betainc";
static char __pyx_k_36[] = "x=btdtri(a,b,p) returns the pth quantile of the beta distribution.  It is\neffectively the inverse of btdtr returning the value of x for which\nbtdtr(a,b,x) = p.   SEE ALSO betaincinv";
static char __pyx_k_37[] = "y=cbrt(x) returns the real cube root of x.";
static char __pyx_k_38[] = "p=chdtr(v,x) Returns the area under the left hand tail (from 0 to x) of the Chi\nsquare probability density function with v degrees of freedom:\n1/(2**(v/2) * gamma(v/2)) * integral(t**(v/2-1) * exp(-t/2), t=0..x)";
static char __pyx_k_39[] = "p=chdtrc(v,x) returns the area under the right hand tail (from x to\ninfinity) of the Chi square probability density function with v\ndegrees of freedom:\n1/(2**(v/2) * gamma(v/2)) * integral(t**(v/2-1) * exp(-t/2), t=x..inf)";
static char __pyx_k_40[] = "x=chdtri(v,p) returns the argument x such that chdtrc(v,x) is equal\nto p.";
static char __pyx_k_41[] = "y=cosdg(x) calculates the cosine of the angle x given in degrees.";
static char __pyx_k_42[] = "y=calculates cos(x) - 1 for use when x is near zero.";
static char __pyx_k_43[] = "y=cotdg(x) calculates the cotangent of the angle x given in degrees.";
static char __pyx_k_44[] = "y=ellipe(m) returns the complete integral of the second kind:\nintegral(sqrt(1-m*sin(t)**2),t=0..pi/2)";
static char __pyx_k_45[] = "y=ellipeinc(phi,m) returns the incomplete elliptic integral of the\nsecond kind: integral(sqrt(1-m*sin(t)**2),t=0..phi)";
static char __pyx_k_46[] = "(sn,cn,dn,ph)=ellipj(u,m) calculates the Jacobian elliptic functions of\nparameter m between 0 and 1, and real u.  The returned functions are\noften written sn(u|m), cn(u|m), and dn(u|m).  The value of ph is such\nthat if u = ellik(ph,m), then sn(u|m) = sin(ph) and cn(u|m) = cos(ph).";
static char __pyx_k_47[] = "y=ellipkinc(phi,m) returns the incomplete elliptic integral of the first\nkind: integral(1/sqrt(1-m*sin(t)**2),t=0..phi)";
static char __pyx_k_48[] = "y=ellipkm1(1 - m) returns the complete integral of the first kind:\nintegral(1/sqrt(1-m*sin(t)**2),t=0..pi/2)";
static char __pyx_k_49[] = "eval_chebyc(n, x, out=None)\n\nEvaluate Chebyshev C polynomial at a point.";
static char __pyx_k_50[] = "eval_chebys(n, x, out=None)\n\nEvaluate Chebyshev S polynomial at a point.";
static char __pyx_k_51[] = "eval_chebyt(n, x, out=None)\n\nEvaluate Chebyshev T polynomial at a point.\n\nThis routine is numerically stable for `x` in ``[-1, 1]`` at least\nup to order ``10000``.";
static char __pyx_k_52[] = "eval_chebyu(n, x, out=None)\n\nEvaluate Chebyshev U polynomial at a point.";
static char __pyx_k_53[] = "eval_gegenbauer(n, alpha, x, out=None)\n\nEvaluate Gegenbauer polynomial at a point.";
static char __pyx_k_54[] = "eval_genlaguerre(n, alpha, x, out=None)\n\nEvaluate generalized Laguerre polynomial at a point.";
static char __pyx_k_55[] = "eval_hermite(n, x, out=None)\n\nEvaluate Hermite polynomial at a point.";
static char __pyx_k_56[] = "eval_hermitenorm(n, x, out=None)\n\nEvaluate normalized Hermite polynomial at a point.";
static char __pyx_k_57[] = "eval_jacobi(n, alpha, beta, x, out=None)\n\nEvaluate Jacobi polynomial at a point.";
static char __pyx_k_58[] = "eval_laguerre(n, x, out=None)\n\nEvaluate Laguerre polynomial at a point.";
static char __pyx_k_59[] = "eval_legendre(n, x, out=None)\n\nEvaluate Legendre polynomial at a point.";
static char __pyx_k_60[] = "eval_sh_chebyt(n, x, out=None)\n\nEvaluate shifted Chebyshev T polynomial at a point.";
static char __pyx_k_61[] = "eval_sh_chebyu(n, x, out=None)\n\nEvaluate shifted Chebyshev U polynomial at a point.";
static char __pyx_k_62[] = "eval_sh_jacobi(n, p, q, x, out=None)\n\nEvaluate shifted Jacobi polynomial at a point.";
static char __pyx_k_63[] = "eval_sh_legendre(n, x, out=None)\n\nEvaluate shifted Legendre polynomial at a point.";
static char __pyx_k_64[] = "y=exp1(z) returns the exponential integral (n=1) of complex argument\nz: integral(exp(-z*t)/t,t=1..inf).";
static char __pyx_k_65[] = "y=exp10(x) returns 10 raised to the x power.";
static char __pyx_k_66[] = "y=exp2(x) returns 2 raised to the x power.";
static char __pyx_k_67[] = "y=expi(x) returns an exponential integral of argument x defined as\nintegral(exp(t)/t,t=-inf..x).  See expn for a different exponential\nintegral.";
static char __pyx_k_68[] = "Expit ufunc for ndarrays.\n\nThe expit function is defined as expit(x) = 1/(1+exp(-x)).\nNote that expit is the inverse logit function.\n\nParameters\n----------\nx : ndarray\n    The ndarray to apply expit to element-wise.\n\nReturns\n-------\nout : ndarray\n    An ndarray of the same shape as x. Its entries\n    are expit of the corresponding entry of x.\n\nNotes\n-----\nAs a ufunc logit takes a number of optional\nkeywork arguments. For more information\nsee `ufuncs <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_";
static char __pyx_k_69[] = "y=expm1(x) calculates exp(x) - 1 for use when x is near zero.";
static char __pyx_k_70[] = "y=expn(n,x) returns the exponential integral for integer n and\nnon-negative x and n: integral(exp(-x*t) / t**n, t=1..inf).";
static char __pyx_k_71[] = "y=fdtr(dfn,dfd,x) returns the area from zero to x under the F density\nfunction (also known as Snedcor's density or the variance ratio\ndensity).  This is the density of X = (unum/dfn)/(uden/dfd), where unum and\nuden are random variables having Chi square distributions with dfn and\ndfd degrees of freedom, respectively.";
static char __pyx_k_72[] = "y=fdtrc(dfn,dfd,x) returns the complemented F distribution function.";
static char __pyx_k_73[] = "x=fdtri(dfn,dfd,p) finds the F density argument x such that\nfdtr(dfn,dfd,x)=p.";
static char __pyx_k_74[] = "x=fdtridfd(dfn,p,x) finds the F density argument dfd such that\nfdtr(dfn,dfd,x)=p.";
static char __pyx_k_75[] = "(ssa,cca)=fresnel(z) returns the Fresnel sin and cos integrals: integral(sin(pi/2\n* t**2),t=0..z) and integral(cos(pi/2 * t**2),t=0..z) for real or\ncomplex z.";
static char __pyx_k_76[] = "y=gamma(z) returns the gamma function of the argument.  The gamma\nfunction is often referred to as the generalized factorial since\nz*gamma(z) = gamma(z+1) and gamma(n+1) = n! for natural number n.";
static char __pyx_k_77[] = "y=gammainc(a,x) returns the incomplete gamma integral defined as\n1 / gamma(a) * integral(exp(-t) * t**(a-1), t=0..x).  a must be\npositive and x must be >= 0.";
static char __pyx_k_78[] = "y=gammaincc(a,x) returns the complemented incomplete gamma integral\ndefined as 1 / gamma(a) * integral(exp(-t) * t**(a-1), t=x..inf) = 1 -\ngammainc(a,x).  a must be positive and x must be >= 0.";
static char __pyx_k_79[] = "x=gammainccinv(a,y) returns x such that gammaincc(a,x) = y.";
static char __pyx_k_80[] = "gammaincinv(a, y) returns x such that gammainc(a, x) = y.";
static char __pyx_k_81[] = "y=gammaln(z) returns the base e logarithm of the absolute value of the\ngamma function of z: ln(abs(gamma(z)))\n\nSee Also\n--------\ngammasgn";
static char __pyx_k_82[] = "y=gammasgn(x) returns the sign of the gamma function.\n\nSee Also\n--------\ngammaln";
static char __pyx_k_83[] = "y=gdtr(a,b,x) returns the integral from zero to x of the gamma\nprobability density function: a**b / gamma(b) * integral(t**(b-1) exp(-at),t=0..x).\nThe arguments a and b are used differently here than in other definitions.";
static char __pyx_k_84[] = "y=gdtrc(a,b,x) returns the integral from x to infinity of the gamma\nprobability density function.  SEE gdtr, gdtri";
static char __pyx_k_85[] = "y=hankel1(v,z) returns the Hankel function of the first kind for real order v and complex argument z.";
static char __pyx_k_86[] = "y=hankel1e(v,z) returns the exponentially scaled Hankel function of the first\nkind for real order v and complex argument z:\nhankel1e(v,z) = hankel1(v,z) * exp(-1j * z)";
static char __pyx_k_87[] = "y=hankel2(v,z) returns the Hankel function of the second kind for real order v and complex argument z.";
static char __pyx_k_88[] = "y=hankel2e(v,z) returns the exponentially scaled Hankel function of the second\nkind for real order v and complex argument z:\nhankel1e(v,z) = hankel1(v,z) * exp(1j * z)";
static char __pyx_k_89[] = "y=hyp1f1(a,b,x) returns the confluent hypergeometeric function\n( 1F1(a,b;x) ) evaluated at the values a, b, and x.";
static char __pyx_k_90[] = "(y,err)=hyp1f2(a,b,c,x) returns (y,err) with the hypergeometric function 1F2 in y and an error estimate in err.";
static char __pyx_k_91[] = "(y,err)=hyp2f0(a,b,x,type) returns (y,err) with the hypergeometric function 2F0 in y and an error estimate in err.  The input type determines a convergence factor and\ncan be either 1 or 2.";
static char __pyx_k_92[] = "y=hyp2f1(a,b,c,z) returns the Gauss hypergeometric function\n( 2F1(a,b;c;z) ).";
static char __pyx_k_93[] = "(y,err)=hyp3f0(a,b,c,x) returns (y,err) with the hypergeometric function 3F0 in y and an error estimate in err.";
static char __pyx_k_94[] = "y=hyperu(a,b,x) returns the confluent hypergeometric function of the\nsecond kind U(a,b,x).";
static char __pyx_k_95[] = "y=i0(x) returns the modified Bessel function of order 0 at x.";
static char __pyx_k_96[] = "y=i0e(x) returns the exponentially scaled modified Bessel function\nof order 0 at x.  i0e(x) = exp(-abs(x)) * i0(x).";
static char __pyx_k_97[] = "y=i1(x) returns the modified Bessel function of order 1 at x.";
static char __pyx_k_98[] = "y=i1e(x) returns the exponentially scaled modified Bessel function\nof order 0 at x.  i1e(x) = exp(-abs(x)) * i1(x).";
static char __pyx_k_99[] = "(ii0,ik0)=it2i0k0(x) returns the integrals int((i0(t)-1)/t,t=0..x) and\nint(k0(t)/t,t=x..infinitity).";
static char __pyx_k__B[] = "B";
static char __pyx_k__H[] = "H";
static char __pyx_k__I[] = "I";
static char __pyx_k__L[] = "L";
static char __pyx_k__O[] = "O";
static char __pyx_k__Q[] = "Q";
static char __pyx_k__b[] = "b";
static char __pyx_k__d[] = "d";
static char __pyx_k__f[] = "f";
static char __pyx_k__g[] = "g";
static char __pyx_k__h[] = "h";
static char __pyx_k__i[] = "i";
static char __pyx_k__l[] = "l";
static char __pyx_k__q[] = "q";
static char __pyx_k_100[] = "(ij0,iy0)=it2j0y0(x) returns the integrals int((1-j0(t))/t,t=0..x) and\nint(y0(t)/t,t=x..infinitity).";
static char __pyx_k_101[] = "y=it2struve0(x) returns the integral of the Struve function of order 0\ndivided by t from x to infinity:  integral(H0(t)/t, t=x..inf).";
static char __pyx_k_102[] = "(Apt,Bpt,Ant,Bnt)=itairy(x) calculates the integral of Airy functions from 0 to x\nfor positive (Apt, Bpt) and negative (Ant, Bnt) arguments.";
static char __pyx_k_103[] = "(ii0,ik0)=iti0k0(x) returns simple integrals from 0 to x of the zeroth order\nmodified Bessel functions i0 and k0.";
static char __pyx_k_104[] = "(ij0,iy0)=itj0y0(x) returns simple integrals from 0 to x of the zeroth order\nBessel functions j0 and y0.";
static char __pyx_k_105[] = "y=itmodstruve0(x) returns the integral of the modified Struve function\nof order 0 from 0 to x:  integral(L0(t), t=0..x).";
static char __pyx_k_106[] = "y=itstruve0(x) returns the integral of the Struve function of order 0\nfrom 0 to x:  integral(H0(t), t=0..x).";
static char __pyx_k_107[] = "y=iv(v,z) returns the modified Bessel function of real order v of\nz.  If z is of real type and negative, v must be integer valued.";
static char __pyx_k_108[] = "y=ive(v,z) returns the exponentially scaled modified Bessel function of\nreal order v and complex z: ive(v,z) = iv(v,z) * exp(-abs(z.real))";
static char __pyx_k_109[] = "y=j0(x) returns the Bessel function of order 0 at x.";
static char __pyx_k_110[] = "y=j1(x) returns the Bessel function of order 1 at x.";
static char __pyx_k_111[] = "y=jv(v,z) returns the Bessel function of real order v at complex z.";
static char __pyx_k_112[] = "y=jve(v,z) returns the exponentially scaled Bessel function of real order\nv at complex z: jve(v,z) = jv(v,z) * exp(-abs(z.imag))";
static char __pyx_k_113[] = "y=k0(x) returns the modified Bessel function of the second kind (sometimes called the third kind) of\norder 0 at x.";
static char __pyx_k_114[] = "y=k0e(x) returns the exponentially scaled modified Bessel function\nof the second kind (sometimes called the third kind) of order 0 at x.  k0e(x) = exp(x) * k0(x).";
static char __pyx_k_115[] = "y=i1(x) returns the modified Bessel function of the second kind (sometimes called the third kind) of\norder 1 at x.";
static char __pyx_k_116[] = "y=k1e(x) returns the exponentially scaled modified Bessel function\nof the second kind (sometimes called the third kind) of order 1 at x.  k1e(x) = exp(x) * k1(x)";
static char __pyx_k_117[] = "y=kei(x) returns the Kelvin function ker x";
static char __pyx_k_118[] = "y=keip(x) returns the derivative of the Kelvin function kei x";
static char __pyx_k_119[] = "(Be, Ke, Bep, Kep)=kelvin(x) returns the tuple (Be, Ke, Bep, Kep) which contains\ncomplex numbers representing the real and imaginary Kelvin functions\nand their derivatives evaluated at x.  For example,\nkelvin(x)[0].real = ber x and kelvin(x)[0].imag = bei x with similar\nrelationships for ker and kei.";
static char __pyx_k_120[] = "y=ker(x) returns the Kelvin function ker x";
static char __pyx_k_121[] = "y=kerp(x) returns the derivative of the Kelvin function ker x";
static char __pyx_k_122[] = "y=kn(n,x) returns the modified Bessel function of the second kind (sometimes called the third kind) for\ninteger order n at x.";
static char __pyx_k_123[] = "y=kolmogi(p) returns y such that kolmogorov(y) = p";
static char __pyx_k_124[] = "p=kolmogorov(y) returns the complementary cumulative distribution\nfunction of Kolmogorov's limiting distribution (Kn* for large n)\nof a two-sided test for equality between an empirical and a theoretical\ndistribution. It is equal to the (limit as n->infinity of the) probability\nthat sqrt(n) * max absolute deviation > y.";
static char __pyx_k_125[] = "y=kv(v,z) returns the modified Bessel function of the second kind (sometimes called the third kind) for\nreal order v at complex z.";
static char __pyx_k_126[] = "y=kve(v,z) returns the exponentially scaled, modified Bessel function\nof the second kind (sometimes called the third kind) for real order v at complex z: kve(v,z) = kv(v,z) * exp(z)";
static char __pyx_k_127[] = "y=log1p(x) calculates log(1+x) for use when x is near zero.";
static char __pyx_k_128[] = "y=log_ndtr(x) returns the log of the area under the standard Gaussian probability\ndensity function, integrated from minus infinity to x:\n1/sqrt(2*pi) * integral(exp(-t**2 / 2),t=-inf..x)";
static char __pyx_k_129[] = "Logit ufunc for ndarrays.\n\nThe logit function is defined as logit(p) = log(p/(1-p)).\nNote that logit(0) = -inf, logit(1) = inf, and logit(p)\nfor p<0 or p>1 yields nan.\n\nParameters\n----------\nx : ndarray\n    The ndarray to apply logit to element-wise.\n\nReturns\n-------\nout : ndarray\n    An ndarray of the same shape as x. Its entries\n    are logit of the corresponding entry of x.\n\nNotes\n-----\nAs a ufunc logit takes a number of optional\nkeywork arguments. For more information\nsee `ufuncs <http://docs.scipy.org/doc/numpy/reference/ufuncs.html>`_";
static char __pyx_k_130[] = "y=lpmv(m,v,x) returns the associated legendre function of integer order\nm and real degree v (s.t. v>-m-1 or v<m): ``|x| <= 1``.";
static char __pyx_k_131[] = "lmbda=mathieu_a(m,q) returns the characteristic value for the even solution,\nce_m(z,q), of Mathieu's equation";
static char __pyx_k_132[] = "lmbda=mathieu_b(m,q) returns the characteristic value for the odd solution,\nse_m(z,q), of Mathieu's equation";
static char __pyx_k_133[] = "(y,yp)=mathieu_cem(m,q,x) returns the even Mathieu function, ce_m(x,q),\nof order m and parameter q evaluated at x (given in degrees).\nAlso returns the derivative with respect to x of ce_m(x,q)";
static char __pyx_k_134[] = "(y,yp)=mathieu_modcem1(m,q,x) evaluates the even modified Mathieu function\nof the first kind, Mc1m(x,q), and its derivative at x for order m and\nparameter q.";
static char __pyx_k_135[] = "(y,yp)=mathieu_modcem2(m,q,x) evaluates the even modified Mathieu function\nof the second kind, Mc2m(x,q), and its derivative at x (given in degrees)\nfor order m and parameter q.";
static char __pyx_k_136[] = "(y,yp)=mathieu_modsem1(m,q,x) evaluates the odd modified Mathieu function\nof the first kind, Ms1m(x,q), and its derivative at x (given in degrees)\nfor order m and parameter q.";
static char __pyx_k_137[] = "(y,yp)=mathieu_modsem2(m,q,x) evaluates the odd modified Mathieu function\nof the second kind, Ms2m(x,q), and its derivative at x (given in degrees)\nfor order m and parameter q.";
static char __pyx_k_138[] = "(y,yp)=mathieu_sem(m,q,x) returns the odd Mathieu function, se_m(x,q),\nof order m and parameter q evaluated at x (given in degrees).\nAlso returns the derivative with respect to x of se_m(x,q).";
static char __pyx_k_139[] = "(fm,km)=modfresnelp(x) returns the modified Fresnel integrals ``F_-(x)`` and ``K_-(x)``\nas ``fp=integral(exp(-1j*t*t),t=x..inf)`` and ``kp=1/sqrt(pi)*exp(1j*(x*x+pi/4))*fp``";
static char __pyx_k_140[] = "(fp,kp)=modfresnelp(x) returns the modified Fresnel integrals F_+(x) and K_+(x)\nas fp=integral(exp(1j*t*t),t=x..inf) and kp=1/sqrt(pi)*exp(-1j*(x*x+pi/4))*fp";
static char __pyx_k_141[] = "y=modstruve(v,x) returns the modified Struve function Lv(x) of order\nv at x, x must be positive unless v is an integer and it is recommended\nthat ``|v| <= 20``.";
static char __pyx_k_142[] = "y=nbdtr(k,n,p) returns the sum of the terms 0 through k of the\nnegative binomial distribution: sum((n+j-1)Cj p**n (1-p)**j,j=0..k).\nIn a sequence of Bernoulli trials this is the probability that k or\nfewer failures precede the nth success.";
static char __pyx_k_143[] = "y=nbdtrc(k,n,p) returns the sum of the terms k+1 to infinity of the\nnegative binomial distribution.";
static char __pyx_k_144[] = "p=nbdtri(k,n,y) finds the argument p such that nbdtr(k,n,p)=y.";
static char __pyx_k_145[] = "k=nbdtrik(y,n,p) finds the argument k such that nbdtr(k,n,p)=y.";
static char __pyx_k_146[] = "n=nbdtrin(k,y,p) finds the argument n such that nbdtr(k,n,p)=y.";
static char __pyx_k_147[] = "y=ndtr(x) returns the area under the standard Gaussian probability\ndensity function, integrated from minus infinity to x:\n1/sqrt(2*pi) * integral(exp(-t**2 / 2),t=-inf..x)";
static char __pyx_k_148[] = "x=ndtri(y) returns the argument x for which the area udnder the\nGaussian probability density function (integrated from minus infinity\nto x) is equal to y.";
static char __pyx_k_149[] = "(s,sp)=obl_ang1(m,n,c,x) computes the oblate sheroidal angular function\nof the first kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.";
static char __pyx_k_150[] = "(s,sp)=obl_ang1_cv(m,n,c,cv,x) computes the oblate sheroidal angular function\nof the first kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed\ncharacteristic value.";
static char __pyx_k_151[] = "cv=obl_cv(m,n,c) computes the characteristic value of oblate spheroidal\nwave functions of order m,n (n>=m) and spheroidal parameter c.";
static char __pyx_k_152[] = "(s,sp)=obl_rad1(m,n,c,x) computes the oblate sheroidal radial function\nof the first kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.";
static char __pyx_k_153[] = "(s,sp)=obl_rad1_cv(m,n,c,cv,x) computes the oblate sheroidal radial function\nof the first kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed\ncharacteristic value.";
static char __pyx_k_154[] = "(s,sp)=obl_rad2(m,n,c,x) computes the oblate sheroidal radial function\nof the second kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.";
static char __pyx_k_155[] = "(s,sp)=obl_rad2_cv(m,n,c,cv,x) computes the oblate sheroidal radial function\nof the second kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed\ncharacteristic value.";
static char __pyx_k_156[] = "(d,dp)=pbdv(v,x) returns (d,dp) with the parabolic cylinder function Dv(x) in\nd and the derivative, Dv'(x) in dp.";
static char __pyx_k_157[] = "(v,vp)=pbvv(v,x) returns (v,vp) with the parabolic cylinder function Vv(x) in\nv and the derivative, Vv'(x) in vp.";
static char __pyx_k_158[] = "(w,wp)=pbwa(a,x) returns (w,wp) with the parabolic cylinder function W(a,x) in\nw and the derivative, W'(a,x) in wp.  May not be accurate for large (>5)\narguments in a and/or x.";
static char __pyx_k_159[] = "y=pdtr(k,m) returns the sum of the first k terms of the Poisson\ndistribution: sum(exp(-m) * m**j / j!, j=0..k) = gammaincc( k+1, m).\nArguments must both be positive and k an integer.";
static char __pyx_k_160[] = "y=pdtrc(k,m) returns the sum of the terms from k+1 to infinity of the\nPoisson distribution: sum(exp(-m) * m**j / j!, j=k+1..inf) = gammainc( k+1, m).\nArguments must both be positive and k an integer.";
static char __pyx_k_161[] = "m=pdtri(k,y) returns the Poisson variable m such that the sum\nfrom 0 to k of the Poisson density is equal to the given probability\ny:  calculated by gammaincinv( k+1, y).  k must be a nonnegative integer and\ny between 0 and 1.";
static char __pyx_k_162[] = "k=pdtrik(p,m) returns the quantile k such that pdtr(k,m)=p";
static char __pyx_k_163[] = "(s,sp)=pro_ang1(m,n,c,x) computes the prolate sheroidal angular function\nof the first kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.";
static char __pyx_k_164[] = "(s,sp)=pro_ang1_cv(m,n,c,cv,x) computes the prolate sheroidal angular function\nof the first kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed\ncharacteristic value.";
static char __pyx_k_165[] = "cv=pro_cv(m,n,c) computes the characteristic value of prolate spheroidal\nwave functions of order m,n (n>=m) and spheroidal parameter c.";
static char __pyx_k_166[] = "(s,sp)=pro_rad1(m,n,c,x) computes the prolate sheroidal radial function\nof the first kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``.";
static char __pyx_k_167[] = "(s,sp)=pro_rad1_cv(m,n,c,cv,x) computes the prolate sheroidal radial function\nof the first kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed\ncharacteristic value.";
static char __pyx_k_168[] = "(s,sp)=pro_rad2(m,n,c,x) computes the prolate sheroidal radial function\nof the second kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and |x|<1.0.";
static char __pyx_k_169[] = "(s,sp)=pro_rad2_cv(m,n,c,cv,x) computes the prolate sheroidal radial function\nof the second kind and its derivative (with respect to x) for mode parameters\nm>=0 and n>=m, spheroidal parameter c and ``|x| < 1.0``. Requires pre-computed\ncharacteristic value.";
static char __pyx_k_170[] = "y=psi(z) is the derivative of the logarithm of the gamma function\nevaluated at z (also called the digamma function).";
static char __pyx_k_171[] = "y=radian(d,m,s) returns the angle given in (d)egrees, (m)inutes, and\n(s)econds in radians.";
static char __pyx_k_172[] = "y=rgamma(z) returns one divided by the gamma function of x.";
static char __pyx_k_173[] = "y=Returns the nearest integer to x as a double precision\nfloating point result.  If x ends in 0.5 exactly, the\nnearest even integer is chosen.";
static char __pyx_k_174[] = "(shi,chi)=shichi(x) returns the hyperbolic sine and cosine integrals:\nintegral(sinh(t)/t,t=0..x) and eul + ln x +\nintegral((cosh(t)-1)/t,t=0..x) where eul is Euler's Constant.";
static char __pyx_k_175[] = "(si,ci)=sici(x) returns in si the integral of the sinc function from 0 to x:\nintegral(sin(t)/t,t=0..x).  It returns in ci the cosine integral: eul + ln x +\nintegral((cos(t) - 1)/t,t=0..x).";
static char __pyx_k_176[] = "y=sindg(x) calculates the sine of the angle x given in degrees.";
static char __pyx_k_177[] = "y=smirnov(n,e) returns the exact Kolmogorov-Smirnov complementary\ncumulative distribution function (Dn+ or Dn-) for a one-sided test of\nequality between an empirical and a theoretical distribution. It is equal\nto the probability that the maximum difference between a theoretical\ndistribution and an empirical one based on n samples is greater than e.";
static char __pyx_k_178[] = "e=smirnovi(n,y) returns e such that smirnov(n,e) = y.";
static char __pyx_k_179[] = "y=spence(x) returns the dilogarithm integral: -integral(log t /\n(t-1),t=1..x)";
static char __pyx_k_180[] = "p=stdtr(df,t) returns the integral from minus infinity to t of the Student t\ndistribution with df > 0 degrees of freedom:\ngamma((df+1)/2)/(sqrt(df*pi)*gamma(df/2)) * integral((1+x**2/df)**(-df/2-1/2),\nx=-inf..t)";
static char __pyx_k_181[] = "t=stdtridf(p,t) returns the argument df such that stdtr(df,t) is equal to p.";
static char __pyx_k_182[] = "t=stdtrit(df,p) returns the argument t such that stdtr(df,t) is equal to p.";
static char __pyx_k_183[] = "y=struve(v,x) returns the Struve function Hv(x) of order v at x, x\nmust be positive unless v is an integer.";
static char __pyx_k_184[] = "y=tandg(x) calculates the tangent of the angle x given in degrees.";
static char __pyx_k_185[] = "y=y0(x) returns the Bessel function of the second kind of order 0 at x.";
static char __pyx_k_186[] = "y=y1(x) returns the Bessel function of the second kind of order 1 at x.";
static char __pyx_k_187[] = "y=yn(n,x) returns the Bessel function of the second kind of integer\norder n at x.";
static char __pyx_k_188[] = "y=yv(v,z) returns the Bessel function of the second kind of real\norder v at complex z.";
static char __pyx_k_189[] = "y=yve(v,z) returns the exponentially scaled Bessel function of the second\nkind of real order v at complex z: yve(v,z) = yv(v,z) * exp(-abs(z.imag))";
static char __pyx_k_190[] = "y=zeta(x,q) returns the Riemann zeta function of two arguments:\nsum((k+q)**(-x),k=0..inf)";
static char __pyx_k_191[] = "y=zetac(x) returns 1.0 - the Riemann zeta function: sum(k**(-x), k=2..inf)";
static char __pyx_k_194[] = "/home/pauli/prj/scipy/scipy/scipy/special/_ufuncs.pyx";
static char __pyx_k_195[] = "scipy.special._ufuncs";
static char __pyx_k__Zd[] = "Zd";
static char __pyx_k__Zf[] = "Zf";
static char __pyx_k__Zg[] = "Zg";
static char __pyx_k__i0[] = "i0";
static char __pyx_k__i1[] = "i1";
static char __pyx_k__iv[] = "iv";
static char __pyx_k__j0[] = "j0";
static char __pyx_k__j1[] = "j1";
static char __pyx_k__jn[] = "jn";
static char __pyx_k__jv[] = "jv";
static char __pyx_k__k0[] = "k0";
static char __pyx_k__k1[] = "k1";
static char __pyx_k__kn[] = "kn";
static char __pyx_k__kv[] = "kv";
static char __pyx_k__y0[] = "y0";
static char __pyx_k__y1[] = "y1";
static char __pyx_k__yn[] = "yn";
static char __pyx_k__yv[] = "yv";
static char __pyx_k__bei[] = "bei";
static char __pyx_k__ber[] = "ber";
static char __pyx_k__i0e[] = "i0e";
static char __pyx_k__i1e[] = "i1e";
static char __pyx_k__ive[] = "ive";
static char __pyx_k__jve[] = "jve";
static char __pyx_k__k0e[] = "k0e";
static char __pyx_k__k1e[] = "k1e";
static char __pyx_k__kei[] = "kei";
static char __pyx_k__ker[] = "ker";
static char __pyx_k__kve[] = "kve";
static char __pyx_k__psi[] = "psi";
static char __pyx_k__yve[] = "yve";
static char __pyx_k__airy[] = "airy";
static char __pyx_k__bdtr[] = "bdtr";
static char __pyx_k__beip[] = "beip";
static char __pyx_k__berp[] = "berp";
static char __pyx_k__beta[] = "beta";
static char __pyx_k__cbrt[] = "cbrt";
static char __pyx_k__exp1[] = "exp1";
static char __pyx_k__exp2[] = "exp2";
static char __pyx_k__expi[] = "expi";
static char __pyx_k__expn[] = "expn";
static char __pyx_k__fdtr[] = "fdtr";
static char __pyx_k__gdtr[] = "gdtr";
static char __pyx_k__keip[] = "keip";
static char __pyx_k__kerp[] = "kerp";
static char __pyx_k__lpmv[] = "lpmv";
static char __pyx_k__ndtr[] = "ndtr";
static char __pyx_k__pbdv[] = "pbdv";
static char __pyx_k__pbvv[] = "pbvv";
static char __pyx_k__pbwa[] = "pbwa";
static char __pyx_k__pdtr[] = "pdtr";
static char __pyx_k__sici[] = "sici";
static char __pyx_k__zeta[] = "zeta";
static char __pyx_k__airye[] = "airye";
static char __pyx_k__bdtrc[] = "bdtrc";
static char __pyx_k__bdtri[] = "bdtri";
static char __pyx_k__binom[] = "binom";
static char __pyx_k__btdtr[] = "btdtr";
static char __pyx_k__chdtr[] = "chdtr";
static char __pyx_k__cosdg[] = "cosdg";
static char __pyx_k__cosm1[] = "cosm1";
static char __pyx_k__cotdg[] = "cotdg";
static char __pyx_k__exp10[] = "exp10";
static char __pyx_k__expit[] = "expit";
static char __pyx_k__expm1[] = "expm1";
static char __pyx_k__fdtrc[] = "fdtrc";
static char __pyx_k__fdtri[] = "fdtri";
static char __pyx_k__gamma[] = "gamma";
static char __pyx_k__gdtrc[] = "gdtrc";
static char __pyx_k__log1p[] = "log1p";
static char __pyx_k__logit[] = "logit";
static char __pyx_k__nbdtr[] = "nbdtr";
static char __pyx_k__ndtri[] = "ndtri";
static char __pyx_k__pdtrc[] = "pdtrc";
static char __pyx_k__pdtri[] = "pdtri";
static char __pyx_k__range[] = "range";
static char __pyx_k__round[] = "round";
static char __pyx_k__sindg[] = "sindg";
static char __pyx_k__stdtr[] = "stdtr";
static char __pyx_k__tandg[] = "tandg";
static char __pyx_k__zetac[] = "zetac";
static char __pyx_k__bdtrik[] = "bdtrik";
static char __pyx_k__bdtrin[] = "bdtrin";
static char __pyx_k__betaln[] = "betaln";
static char __pyx_k__btdtri[] = "btdtri";
static char __pyx_k__chdtrc[] = "chdtrc";
static char __pyx_k__chdtri[] = "chdtri";
static char __pyx_k__chndtr[] = "chndtr";
static char __pyx_k__ellipe[] = "ellipe";
static char __pyx_k__ellipj[] = "ellipj";
static char __pyx_k__gdtria[] = "gdtria";
static char __pyx_k__gdtrib[] = "gdtrib";
static char __pyx_k__gdtrix[] = "gdtrix";
static char __pyx_k__hyp1f1[] = "hyp1f1";
static char __pyx_k__hyp1f2[] = "hyp1f2";
static char __pyx_k__hyp2f0[] = "hyp2f0";
static char __pyx_k__hyp2f1[] = "hyp2f1";
static char __pyx_k__hyp3f0[] = "hyp3f0";
static char __pyx_k__hyperu[] = "hyperu";
static char __pyx_k__inflag[] = "inflag";
static char __pyx_k__itairy[] = "itairy";
static char __pyx_k__iti0k0[] = "iti0k0";
static char __pyx_k__itj0y0[] = "itj0y0";
static char __pyx_k__kelvin[] = "kelvin";
static char __pyx_k__nbdtrc[] = "nbdtrc";
static char __pyx_k__nbdtri[] = "nbdtri";
static char __pyx_k__ncfdtr[] = "ncfdtr";
static char __pyx_k__nctdtr[] = "nctdtr";
static char __pyx_k__obl_cv[] = "obl_cv";
static char __pyx_k__pdtrik[] = "pdtrik";
static char __pyx_k__pro_cv[] = "pro_cv";
static char __pyx_k__radian[] = "radian";
static char __pyx_k__rgamma[] = "rgamma";
static char __pyx_k__shichi[] = "shichi";
static char __pyx_k__spence[] = "spence";
static char __pyx_k__struve[] = "struve";
static char __pyx_k__betainc[] = "betainc";
static char __pyx_k__btdtria[] = "btdtria";
static char __pyx_k__btdtrib[] = "btdtrib";
static char __pyx_k__chdtriv[] = "chdtriv";
static char __pyx_k__fresnel[] = "fresnel";
static char __pyx_k__gammaln[] = "gammaln";
static char __pyx_k__hankel1[] = "hankel1";
static char __pyx_k__hankel2[] = "hankel2";
static char __pyx_k__it2i0k0[] = "it2i0k0";
static char __pyx_k__it2j0y0[] = "it2j0y0";
static char __pyx_k__kolmogi[] = "kolmogi";
static char __pyx_k__nbdtrik[] = "nbdtrik";
static char __pyx_k__nbdtrin[] = "nbdtrin";
static char __pyx_k__ncfdtri[] = "ncfdtri";
static char __pyx_k__smirnov[] = "smirnov";
static char __pyx_k__stdtrit[] = "stdtrit";
static char __pyx_k__tklmbda[] = "tklmbda";
static char __pyx_k____main__[] = "__main__";
static char __pyx_k____test__[] = "__test__";
static char __pyx_k__chndtrix[] = "chndtrix";
static char __pyx_k__ellipkm1[] = "ellipkm1";
static char __pyx_k__fdtridfd[] = "fdtridfd";
static char __pyx_k__gammainc[] = "gammainc";
static char __pyx_k__gammasgn[] = "gammasgn";
static char __pyx_k__hankel1e[] = "hankel1e";
static char __pyx_k__hankel2e[] = "hankel2e";
static char __pyx_k__lambertw[] = "lambertw";
static char __pyx_k__log_ndtr[] = "log_ndtr";
static char __pyx_k__nctdtrit[] = "nctdtrit";
static char __pyx_k__nrdtrimn[] = "nrdtrimn";
static char __pyx_k__nrdtrisd[] = "nrdtrisd";
static char __pyx_k__obl_ang1[] = "obl_ang1";
static char __pyx_k__obl_rad1[] = "obl_rad1";
static char __pyx_k__obl_rad2[] = "obl_rad2";
static char __pyx_k__pro_ang1[] = "pro_ang1";
static char __pyx_k__pro_rad1[] = "pro_rad1";
static char __pyx_k__pro_rad2[] = "pro_rad2";
static char __pyx_k__smirnovi[] = "smirnovi";
static char __pyx_k__stdtridf[] = "stdtridf";
static char __pyx_k___errprint[] = "_errprint";
static char __pyx_k___lambertw[] = "_lambertw";
static char __pyx_k__chndtridf[] = "chndtridf";
static char __pyx_k__chndtrinc[] = "chndtrinc";
static char __pyx_k__ellipeinc[] = "ellipeinc";
static char __pyx_k__ellipkinc[] = "ellipkinc";
static char __pyx_k__gammaincc[] = "gammaincc";
static char __pyx_k__itstruve0[] = "itstruve0";
static char __pyx_k__mathieu_a[] = "mathieu_a";
static char __pyx_k__mathieu_b[] = "mathieu_b";
static char __pyx_k__modstruve[] = "modstruve";
static char __pyx_k__ncfdtrinc[] = "ncfdtrinc";
static char __pyx_k__nctdtridf[] = "nctdtridf";
static char __pyx_k__nctdtrinc[] = "nctdtrinc";
static char __pyx_k__ValueError[] = "ValueError";
static char __pyx_k__besselpoly[] = "besselpoly";
static char __pyx_k__betaincinv[] = "betaincinv";
static char __pyx_k__it2struve0[] = "it2struve0";
static char __pyx_k__kolmogorov[] = "kolmogorov";
static char __pyx_k__ncfdtridfd[] = "ncfdtridfd";
static char __pyx_k__ncfdtridfn[] = "ncfdtridfn";
static char __pyx_k__eval_chebyc[] = "eval_chebyc";
static char __pyx_k__eval_chebys[] = "eval_chebys";
static char __pyx_k__eval_chebyt[] = "eval_chebyt";
static char __pyx_k__eval_chebyu[] = "eval_chebyu";
static char __pyx_k__eval_jacobi[] = "eval_jacobi";
static char __pyx_k__gammaincinv[] = "gammaincinv";
static char __pyx_k__mathieu_cem[] = "mathieu_cem";
static char __pyx_k__mathieu_sem[] = "mathieu_sem";
static char __pyx_k__modfresnelm[] = "modfresnelm";
static char __pyx_k__modfresnelp[] = "modfresnelp";
static char __pyx_k__obl_ang1_cv[] = "obl_ang1_cv";
static char __pyx_k__obl_rad1_cv[] = "obl_rad1_cv";
static char __pyx_k__obl_rad2_cv[] = "obl_rad2_cv";
static char __pyx_k__pro_ang1_cv[] = "pro_ang1_cv";
static char __pyx_k__pro_rad1_cv[] = "pro_rad1_cv";
static char __pyx_k__pro_rad2_cv[] = "pro_rad2_cv";
static char __pyx_k__RuntimeError[] = "RuntimeError";
static char __pyx_k__eval_hermite[] = "eval_hermite";
static char __pyx_k__gammainccinv[] = "gammainccinv";
static char __pyx_k__itmodstruve0[] = "itmodstruve0";
static char __pyx_k__FutureWarning[] = "FutureWarning";
static char __pyx_k__eval_laguerre[] = "eval_laguerre";
static char __pyx_k__eval_legendre[] = "eval_legendre";
static char __pyx_k__eval_sh_chebyt[] = "eval_sh_chebyt";
static char __pyx_k__eval_sh_chebyu[] = "eval_sh_chebyu";
static char __pyx_k__eval_sh_jacobi[] = "eval_sh_jacobi";
static char __pyx_k__eval_gegenbauer[] = "eval_gegenbauer";
static char __pyx_k__mathieu_modcem1[] = "mathieu_modcem1";
static char __pyx_k__mathieu_modcem2[] = "mathieu_modcem2";
static char __pyx_k__mathieu_modsem1[] = "mathieu_modsem1";
static char __pyx_k__mathieu_modsem2[] = "mathieu_modsem2";
static char __pyx_k__eval_genlaguerre[] = "eval_genlaguerre";
static char __pyx_k__eval_hermitenorm[] = "eval_hermitenorm";
static char __pyx_k__eval_sh_legendre[] = "eval_sh_legendre";
static PyObject *__pyx_kp_u_10;
static PyObject *__pyx_kp_u_12;
static PyObject *__pyx_kp_u_13;
static PyObject *__pyx_kp_u_16;
static PyObject *__pyx_kp_s_194;
static PyObject *__pyx_n_s_195;
static PyObject *__pyx_kp_b_2;
static PyObject *__pyx_kp_b_3;
static PyObject *__pyx_kp_b_4;
static PyObject *__pyx_kp_u_6;
static PyObject *__pyx_kp_u_8;
static PyObject *__pyx_n_s__FutureWarning;
static PyObject *__pyx_n_s__RuntimeError;
static PyObject *__pyx_n_s__ValueError;
static PyObject *__pyx_n_s____main__;
static PyObject *__pyx_n_s____test__;
static PyObject *__pyx_n_s___errprint;
static PyObject *__pyx_n_s___lambertw;
static PyObject *__pyx_n_s__airy;
static PyObject *__pyx_n_s__airye;
static PyObject *__pyx_n_s__bdtr;
static PyObject *__pyx_n_s__bdtrc;
static PyObject *__pyx_n_s__bdtri;
static PyObject *__pyx_n_s__bdtrik;
static PyObject *__pyx_n_s__bdtrin;
static PyObject *__pyx_n_s__bei;
static PyObject *__pyx_n_s__beip;
static PyObject *__pyx_n_s__ber;
static PyObject *__pyx_n_s__berp;
static PyObject *__pyx_n_s__besselpoly;
static PyObject *__pyx_n_s__beta;
static PyObject *__pyx_n_s__betainc;
static PyObject *__pyx_n_s__betaincinv;
static PyObject *__pyx_n_s__betaln;
static PyObject *__pyx_n_s__binom;
static PyObject *__pyx_n_s__btdtr;
static PyObject *__pyx_n_s__btdtri;
static PyObject *__pyx_n_s__btdtria;
static PyObject *__pyx_n_s__btdtrib;
static PyObject *__pyx_n_s__cbrt;
static PyObject *__pyx_n_s__chdtr;
static PyObject *__pyx_n_s__chdtrc;
static PyObject *__pyx_n_s__chdtri;
static PyObject *__pyx_n_s__chdtriv;
static PyObject *__pyx_n_s__chndtr;
static PyObject *__pyx_n_s__chndtridf;
static PyObject *__pyx_n_s__chndtrinc;
static PyObject *__pyx_n_s__chndtrix;
static PyObject *__pyx_n_s__cosdg;
static PyObject *__pyx_n_s__cosm1;
static PyObject *__pyx_n_s__cotdg;
static PyObject *__pyx_n_s__ellipe;
static PyObject *__pyx_n_s__ellipeinc;
static PyObject *__pyx_n_s__ellipj;
static PyObject *__pyx_n_s__ellipkinc;
static PyObject *__pyx_n_s__ellipkm1;
static PyObject *__pyx_n_s__eval_chebyc;
static PyObject *__pyx_n_s__eval_chebys;
static PyObject *__pyx_n_s__eval_chebyt;
static PyObject *__pyx_n_s__eval_chebyu;
static PyObject *__pyx_n_s__eval_gegenbauer;
static PyObject *__pyx_n_s__eval_genlaguerre;
static PyObject *__pyx_n_s__eval_hermite;
static PyObject *__pyx_n_s__eval_hermitenorm;
static PyObject *__pyx_n_s__eval_jacobi;
static PyObject *__pyx_n_s__eval_laguerre;
static PyObject *__pyx_n_s__eval_legendre;
static PyObject *__pyx_n_s__eval_sh_chebyt;
static PyObject *__pyx_n_s__eval_sh_chebyu;
static PyObject *__pyx_n_s__eval_sh_jacobi;
static PyObject *__pyx_n_s__eval_sh_legendre;
static PyObject *__pyx_n_s__exp1;
static PyObject *__pyx_n_s__exp10;
static PyObject *__pyx_n_s__exp2;
static PyObject *__pyx_n_s__expi;
static PyObject *__pyx_n_s__expit;
static PyObject *__pyx_n_s__expm1;
static PyObject *__pyx_n_s__expn;
static PyObject *__pyx_n_s__fdtr;
static PyObject *__pyx_n_s__fdtrc;
static PyObject *__pyx_n_s__fdtri;
static PyObject *__pyx_n_s__fdtridfd;
static PyObject *__pyx_n_s__fresnel;
static PyObject *__pyx_n_s__gamma;
static PyObject *__pyx_n_s__gammainc;
static PyObject *__pyx_n_s__gammaincc;
static PyObject *__pyx_n_s__gammainccinv;
static PyObject *__pyx_n_s__gammaincinv;
static PyObject *__pyx_n_s__gammaln;
static PyObject *__pyx_n_s__gammasgn;
static PyObject *__pyx_n_s__gdtr;
static PyObject *__pyx_n_s__gdtrc;
static PyObject *__pyx_n_s__gdtria;
static PyObject *__pyx_n_s__gdtrib;
static PyObject *__pyx_n_s__gdtrix;
static PyObject *__pyx_n_s__hankel1;
static PyObject *__pyx_n_s__hankel1e;
static PyObject *__pyx_n_s__hankel2;
static PyObject *__pyx_n_s__hankel2e;
static PyObject *__pyx_n_s__hyp1f1;
static PyObject *__pyx_n_s__hyp1f2;
static PyObject *__pyx_n_s__hyp2f0;
static PyObject *__pyx_n_s__hyp2f1;
static PyObject *__pyx_n_s__hyp3f0;
static PyObject *__pyx_n_s__hyperu;
static PyObject *__pyx_n_s__i0;
static PyObject *__pyx_n_s__i0e;
static PyObject *__pyx_n_s__i1;
static PyObject *__pyx_n_s__i1e;
static PyObject *__pyx_n_s__inflag;
static PyObject *__pyx_n_s__it2i0k0;
static PyObject *__pyx_n_s__it2j0y0;
static PyObject *__pyx_n_s__it2struve0;
static PyObject *__pyx_n_s__itairy;
static PyObject *__pyx_n_s__iti0k0;
static PyObject *__pyx_n_s__itj0y0;
static PyObject *__pyx_n_s__itmodstruve0;
static PyObject *__pyx_n_s__itstruve0;
static PyObject *__pyx_n_s__iv;
static PyObject *__pyx_n_s__ive;
static PyObject *__pyx_n_s__j0;
static PyObject *__pyx_n_s__j1;
static PyObject *__pyx_n_s__jn;
static PyObject *__pyx_n_s__jv;
static PyObject *__pyx_n_s__jve;
static PyObject *__pyx_n_s__k0;
static PyObject *__pyx_n_s__k0e;
static PyObject *__pyx_n_s__k1;
static PyObject *__pyx_n_s__k1e;
static PyObject *__pyx_n_s__kei;
static PyObject *__pyx_n_s__keip;
static PyObject *__pyx_n_s__kelvin;
static PyObject *__pyx_n_s__ker;
static PyObject *__pyx_n_s__kerp;
static PyObject *__pyx_n_s__kn;
static PyObject *__pyx_n_s__kolmogi;
static PyObject *__pyx_n_s__kolmogorov;
static PyObject *__pyx_n_s__kv;
static PyObject *__pyx_n_s__kve;
static PyObject *__pyx_n_s__log1p;
static PyObject *__pyx_n_s__log_ndtr;
static PyObject *__pyx_n_s__logit;
static PyObject *__pyx_n_s__lpmv;
static PyObject *__pyx_n_s__mathieu_a;
static PyObject *__pyx_n_s__mathieu_b;
static PyObject *__pyx_n_s__mathieu_cem;
static PyObject *__pyx_n_s__mathieu_modcem1;
static PyObject *__pyx_n_s__mathieu_modcem2;
static PyObject *__pyx_n_s__mathieu_modsem1;
static PyObject *__pyx_n_s__mathieu_modsem2;
static PyObject *__pyx_n_s__mathieu_sem;
static PyObject *__pyx_n_s__modfresnelm;
static PyObject *__pyx_n_s__modfresnelp;
static PyObject *__pyx_n_s__modstruve;
static PyObject *__pyx_n_s__nbdtr;
static PyObject *__pyx_n_s__nbdtrc;
static PyObject *__pyx_n_s__nbdtri;
static PyObject *__pyx_n_s__nbdtrik;
static PyObject *__pyx_n_s__nbdtrin;
static PyObject *__pyx_n_s__ncfdtr;
static PyObject *__pyx_n_s__ncfdtri;
static PyObject *__pyx_n_s__ncfdtridfd;
static PyObject *__pyx_n_s__ncfdtridfn;
static PyObject *__pyx_n_s__ncfdtrinc;
static PyObject *__pyx_n_s__nctdtr;
static PyObject *__pyx_n_s__nctdtridf;
static PyObject *__pyx_n_s__nctdtrinc;
static PyObject *__pyx_n_s__nctdtrit;
static PyObject *__pyx_n_s__ndtr;
static PyObject *__pyx_n_s__ndtri;
static PyObject *__pyx_n_s__nrdtrimn;
static PyObject *__pyx_n_s__nrdtrisd;
static PyObject *__pyx_n_s__obl_ang1;
static PyObject *__pyx_n_s__obl_ang1_cv;
static PyObject *__pyx_n_s__obl_cv;
static PyObject *__pyx_n_s__obl_rad1;
static PyObject *__pyx_n_s__obl_rad1_cv;
static PyObject *__pyx_n_s__obl_rad2;
static PyObject *__pyx_n_s__obl_rad2_cv;
static PyObject *__pyx_n_s__pbdv;
static PyObject *__pyx_n_s__pbvv;
static PyObject *__pyx_n_s__pbwa;
static PyObject *__pyx_n_s__pdtr;
static PyObject *__pyx_n_s__pdtrc;
static PyObject *__pyx_n_s__pdtri;
static PyObject *__pyx_n_s__pdtrik;
static PyObject *__pyx_n_s__pro_ang1;
static PyObject *__pyx_n_s__pro_ang1_cv;
static PyObject *__pyx_n_s__pro_cv;
static PyObject *__pyx_n_s__pro_rad1;
static PyObject *__pyx_n_s__pro_rad1_cv;
static PyObject *__pyx_n_s__pro_rad2;
static PyObject *__pyx_n_s__pro_rad2_cv;
static PyObject *__pyx_n_s__psi;
static PyObject *__pyx_n_s__radian;
static PyObject *__pyx_n_s__range;
static PyObject *__pyx_n_s__rgamma;
static PyObject *__pyx_n_s__round;
static PyObject *__pyx_n_s__shichi;
static PyObject *__pyx_n_s__sici;
static PyObject *__pyx_n_s__sindg;
static PyObject *__pyx_n_s__smirnov;
static PyObject *__pyx_n_s__smirnovi;
static PyObject *__pyx_n_s__spence;
static PyObject *__pyx_n_s__stdtr;
static PyObject *__pyx_n_s__stdtridf;
static PyObject *__pyx_n_s__stdtrit;
static PyObject *__pyx_n_s__struve;
static PyObject *__pyx_n_s__tandg;
static PyObject *__pyx_n_s__tklmbda;
static PyObject *__pyx_n_s__y0;
static PyObject *__pyx_n_s__y1;
static PyObject *__pyx_n_s__yn;
static PyObject *__pyx_n_s__yv;
static PyObject *__pyx_n_s__yve;
static PyObject *__pyx_n_s__zeta;
static PyObject *__pyx_n_s__zetac;
static PyObject *__pyx_int_15;
static PyObject *__pyx_k_tuple_7;
static PyObject *__pyx_k_tuple_9;
static PyObject *__pyx_k_tuple_11;
static PyObject *__pyx_k_tuple_14;
static PyObject *__pyx_k_tuple_15;
static PyObject *__pyx_k_tuple_17;
static PyObject *__pyx_k_tuple_192;
static PyObject *__pyx_k_codeobj_193;



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double, double, double, double))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])), ((double)(((float *)__pyx_v_ip3)[0])));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_f_f__As_f_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  float __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((float (*)(float))__pyx_v_func)(((float)(((float *)__pyx_v_ip0)[0])));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_ddddd_dd(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_ip4;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_ip4 = (__pyx_v_args[4]);
  __pyx_v_op0 = (__pyx_v_args[5]);
  __pyx_v_op1 = (__pyx_v_args[6]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double, double, double, double, double *, double *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])), ((double)(((double *)__pyx_v_ip3)[0])), ((double)(((double *)__pyx_v_ip4)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    (((double *)__pyx_v_op1)[0]) = ((double)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_ip4 = (__pyx_v_ip4 + (__pyx_v_steps[4]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[5]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[6]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(double, __pyx_t_double_complex))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_double_complex *)__pyx_v_ip1)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_double_complex *)__pyx_v_ip1)[0])))));

    
    (((__pyx_t_double_complex *)__pyx_v_op0)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov0)), ((double)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_d_dd(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double *, double *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    (((double *)__pyx_v_op1)[0]) = ((double)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_D_D(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(__pyx_t_double_complex))__pyx_v_func)(__pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_double_complex *)__pyx_v_ip0)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_double_complex *)__pyx_v_ip0)[0])))));

    
    (((__pyx_t_double_complex *)__pyx_v_op0)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov0)), ((double)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddi_d_As_fffl_ff(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;
  int __pyx_t_3;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);
  __pyx_v_op1 = (__pyx_v_args[5]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_t_3 = (((int)(((long *)__pyx_v_ip3)[0])) == (((long *)__pyx_v_ip3)[0]));
    if (__pyx_t_3) {

      
      __pyx_v_ov0 = ((double (*)(double, double, double, int, double *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])), ((int)(((long *)__pyx_v_ip3)[0])), (&__pyx_v_ov1));
      goto __pyx_L5;
    }
     {

      
      sf_error(__pyx_v_func_name, SF_ERROR_DOMAIN, __pyx_k_1);

      
      __pyx_v_ov0 = ((double)NPY_NAN);

      
      __pyx_v_ov1 = ((double)NPY_NAN);
    }
    __pyx_L5:;

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    (((float *)__pyx_v_op1)[0]) = ((float)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[5]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DD_As_F_FF(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  __pyx_t_double_complex __pyx_v_ov0;
  __pyx_t_double_complex __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(__pyx_t_double_complex, __pyx_t_double_complex *, __pyx_t_double_complex *))__pyx_v_func)(__pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_float_complex *)__pyx_v_ip0)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_float_complex *)__pyx_v_ip0)[0])))), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((__pyx_t_float_complex *)__pyx_v_op0)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov0)), ((float)__Pyx_CIMAG(__pyx_v_ov0)));

    
    (((__pyx_t_float_complex *)__pyx_v_op1)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov1)), ((float)__Pyx_CIMAG(__pyx_v_ov1)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DDDD_As_f_FFFF(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  char *__pyx_v_op2;
  char *__pyx_v_op3;
  __pyx_t_double_complex __pyx_v_ov0;
  __pyx_t_double_complex __pyx_v_ov1;
  __pyx_t_double_complex __pyx_v_ov2;
  __pyx_t_double_complex __pyx_v_ov3;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);
  __pyx_v_op2 = (__pyx_v_args[3]);
  __pyx_v_op3 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, __pyx_t_double_complex *, __pyx_t_double_complex *, __pyx_t_double_complex *, __pyx_t_double_complex *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1), (&__pyx_v_ov2), (&__pyx_v_ov3));

    
    (((__pyx_t_float_complex *)__pyx_v_op0)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov0)), ((float)__Pyx_CIMAG(__pyx_v_ov0)));

    
    (((__pyx_t_float_complex *)__pyx_v_op1)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov1)), ((float)__Pyx_CIMAG(__pyx_v_ov1)));

    
    (((__pyx_t_float_complex *)__pyx_v_op2)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov2)), ((float)__Pyx_CIMAG(__pyx_v_ov2)));

    
    (((__pyx_t_float_complex *)__pyx_v_op3)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov3)), ((float)__Pyx_CIMAG(__pyx_v_ov3)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));

    
    __pyx_v_op2 = (__pyx_v_op2 + (__pyx_v_steps[3]));

    
    __pyx_v_op3 = (__pyx_v_op3 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_f_ffff(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  char *__pyx_v_op2;
  char *__pyx_v_op3;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  double __pyx_v_ov2;
  double __pyx_v_ov3;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);
  __pyx_v_op2 = (__pyx_v_args[3]);
  __pyx_v_op3 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double *, double *, double *, double *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1), (&__pyx_v_ov2), (&__pyx_v_ov3));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    (((float *)__pyx_v_op1)[0]) = ((float)__pyx_v_ov1);

    
    (((float *)__pyx_v_op2)[0]) = ((float)__pyx_v_ov2);

    
    (((float *)__pyx_v_op3)[0]) = ((float)__pyx_v_ov3);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));

    
    __pyx_v_op2 = (__pyx_v_op2 + (__pyx_v_steps[3]));

    
    __pyx_v_op3 = (__pyx_v_op3 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddi_d_As_dddl_dd(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;
  int __pyx_t_3;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);
  __pyx_v_op1 = (__pyx_v_args[5]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_t_3 = (((int)(((long *)__pyx_v_ip3)[0])) == (((long *)__pyx_v_ip3)[0]));
    if (__pyx_t_3) {

      
      __pyx_v_ov0 = ((double (*)(double, double, double, int, double *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])), ((int)(((long *)__pyx_v_ip3)[0])), (&__pyx_v_ov1));
      goto __pyx_L5;
    }
     {

      
      sf_error(__pyx_v_func_name, SF_ERROR_DOMAIN, __pyx_k_1);

      
      __pyx_v_ov0 = ((double)NPY_NAN);

      
      __pyx_v_ov1 = ((double)NPY_NAN);
    }
    __pyx_L5:;

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    (((double *)__pyx_v_op1)[0]) = ((double)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[5]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ddD_D(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(double, double, __pyx_t_double_complex))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_double_complex *)__pyx_v_ip2)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_double_complex *)__pyx_v_ip2)[0])))));

    
    (((__pyx_t_double_complex *)__pyx_v_op0)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov0)), ((double)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double, double, double))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(double, __pyx_t_double_complex))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_float_complex *)__pyx_v_ip1)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_float_complex *)__pyx_v_ip1)[0])))));

    
    (((__pyx_t_float_complex *)__pyx_v_op0)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov0)), ((float)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dddd_As_ff_ffff(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  char *__pyx_v_op2;
  char *__pyx_v_op3;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  double __pyx_v_ov2;
  double __pyx_v_ov3;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);
  __pyx_v_op1 = (__pyx_v_args[3]);
  __pyx_v_op2 = (__pyx_v_args[4]);
  __pyx_v_op3 = (__pyx_v_args[5]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double, double *, double *, double *, double *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1), (&__pyx_v_ov2), (&__pyx_v_ov3));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    (((float *)__pyx_v_op1)[0]) = ((float)__pyx_v_ov1);

    
    (((float *)__pyx_v_op2)[0]) = ((float)__pyx_v_ov2);

    
    (((float *)__pyx_v_op3)[0]) = ((float)__pyx_v_ov3);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[3]));

    
    __pyx_v_op2 = (__pyx_v_op2 + (__pyx_v_steps[4]));

    
    __pyx_v_op3 = (__pyx_v_op3 + (__pyx_v_steps[5]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_dd_dd(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);
  __pyx_v_op1 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double, double *, double *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    (((double *)__pyx_v_op1)[0]) = ((double)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_fff_ff(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);
  __pyx_v_op1 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double, double, double *, double *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    (((float *)__pyx_v_op1)[0]) = ((float)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ldd__As_ldd_d(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(long, double, double))__pyx_v_func)(((long)(((long *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_dddD_D(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(double, double, double, __pyx_t_double_complex))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])), __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_double_complex *)__pyx_v_ip3)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_double_complex *)__pyx_v_ip3)[0])))));

    
    (((__pyx_t_double_complex *)__pyx_v_op0)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov0)), ((double)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_d_dddd(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  char *__pyx_v_op2;
  char *__pyx_v_op3;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  double __pyx_v_ov2;
  double __pyx_v_ov3;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);
  __pyx_v_op2 = (__pyx_v_args[3]);
  __pyx_v_op3 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double *, double *, double *, double *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1), (&__pyx_v_ov2), (&__pyx_v_ov3));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    (((double *)__pyx_v_op1)[0]) = ((double)__pyx_v_ov1);

    
    (((double *)__pyx_v_op2)[0]) = ((double)__pyx_v_ov2);

    
    (((double *)__pyx_v_op3)[0]) = ((double)__pyx_v_ov3);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));

    
    __pyx_v_op2 = (__pyx_v_op2 + (__pyx_v_steps[3]));

    
    __pyx_v_op3 = (__pyx_v_op3 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_fffff_ff(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_ip4;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_ip4 = (__pyx_v_args[4]);
  __pyx_v_op0 = (__pyx_v_args[5]);
  __pyx_v_op1 = (__pyx_v_args[6]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double, double, double, double, double *, double *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])), ((double)(((float *)__pyx_v_ip3)[0])), ((double)(((float *)__pyx_v_ip4)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    (((float *)__pyx_v_op1)[0]) = ((float)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_ip4 = (__pyx_v_ip4 + (__pyx_v_steps[4]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[5]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[6]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DDDD_As_d_DDDD(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  char *__pyx_v_op2;
  char *__pyx_v_op3;
  __pyx_t_double_complex __pyx_v_ov0;
  __pyx_t_double_complex __pyx_v_ov1;
  __pyx_t_double_complex __pyx_v_ov2;
  __pyx_t_double_complex __pyx_v_ov3;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);
  __pyx_v_op2 = (__pyx_v_args[3]);
  __pyx_v_op3 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, __pyx_t_double_complex *, __pyx_t_double_complex *, __pyx_t_double_complex *, __pyx_t_double_complex *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1), (&__pyx_v_ov2), (&__pyx_v_ov3));

    
    (((__pyx_t_double_complex *)__pyx_v_op0)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov0)), ((double)__Pyx_CIMAG(__pyx_v_ov0)));

    
    (((__pyx_t_double_complex *)__pyx_v_op1)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov1)), ((double)__Pyx_CIMAG(__pyx_v_ov1)));

    
    (((__pyx_t_double_complex *)__pyx_v_op2)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov2)), ((double)__Pyx_CIMAG(__pyx_v_ov2)));

    
    (((__pyx_t_double_complex *)__pyx_v_op3)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov3)), ((double)__Pyx_CIMAG(__pyx_v_ov3)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));

    
    __pyx_v_op2 = (__pyx_v_op2 + (__pyx_v_steps[3]));

    
    __pyx_v_op3 = (__pyx_v_op3 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_ddd_dd(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);
  __pyx_v_op1 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double, double, double *, double *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    (((double *)__pyx_v_op1)[0]) = ((double)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_Dld__As_Dld_D(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(__pyx_t_double_complex, long, double))__pyx_v_func)(__pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_double_complex *)__pyx_v_ip0)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_double_complex *)__pyx_v_ip0)[0])))), ((long)(((long *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])));

    
    (((__pyx_t_double_complex *)__pyx_v_op0)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov0)), ((double)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_ff_ff(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);
  __pyx_v_op1 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double, double *, double *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    (((float *)__pyx_v_op1)[0]) = ((float)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_fffF_F(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(double, double, double, __pyx_t_double_complex))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])), __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_float_complex *)__pyx_v_ip3)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_float_complex *)__pyx_v_ip3)[0])))));

    
    (((__pyx_t_float_complex *)__pyx_v_op0)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov0)), ((float)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;
  int __pyx_t_3;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_t_3 = (((int)(((long *)__pyx_v_ip0)[0])) == (((long *)__pyx_v_ip0)[0]));
    if (__pyx_t_3) {

      
      __pyx_v_ov0 = ((double (*)(int, double))__pyx_v_func)(((int)(((long *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])));
      goto __pyx_L5;
    }
     {

      
      sf_error(__pyx_v_func_name, SF_ERROR_DOMAIN, __pyx_k_1);

      
      __pyx_v_ov0 = ((double)NPY_NAN);
    }
    __pyx_L5:;

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double, double, double, double))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])), ((double)(((double *)__pyx_v_ip3)[0])));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_lld_d(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;
  int __pyx_t_3;
  int __pyx_t_4;
  int __pyx_t_5;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_t_3 = (((int)(((long *)__pyx_v_ip0)[0])) == (((long *)__pyx_v_ip0)[0]));
    if (__pyx_t_3) {
      __pyx_t_4 = (((int)(((long *)__pyx_v_ip1)[0])) == (((long *)__pyx_v_ip1)[0]));
      __pyx_t_5 = __pyx_t_4;
    } else {
      __pyx_t_5 = __pyx_t_3;
    }
    if (__pyx_t_5) {

      
      __pyx_v_ov0 = ((double (*)(int, int, double))__pyx_v_func)(((int)(((long *)__pyx_v_ip0)[0])), ((int)(((long *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])));
      goto __pyx_L5;
    }
     {

      
      sf_error(__pyx_v_func_name, SF_ERROR_DOMAIN, __pyx_k_1);

      
      __pyx_v_ov0 = ((double)NPY_NAN);
    }
    __pyx_L5:;

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DD_As_d_DD(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  __pyx_t_double_complex __pyx_v_ov0;
  __pyx_t_double_complex __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, __pyx_t_double_complex *, __pyx_t_double_complex *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((__pyx_t_double_complex *)__pyx_v_op0)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov0)), ((double)__Pyx_CIMAG(__pyx_v_ov0)));

    
    (((__pyx_t_double_complex *)__pyx_v_op1)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov1)), ((double)__Pyx_CIMAG(__pyx_v_ov1)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_F_F(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(__pyx_t_double_complex))__pyx_v_func)(__pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_float_complex *)__pyx_v_ip0)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_float_complex *)__pyx_v_ip0)[0])))));

    
    (((__pyx_t_float_complex *)__pyx_v_op0)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov0)), ((float)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);
  __pyx_v_op1 = (__pyx_v_args[5]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double, double, double, double, double *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])), ((double)(((float *)__pyx_v_ip3)[0])), (&__pyx_v_ov1));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    (((float *)__pyx_v_op1)[0]) = ((float)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[5]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_llf_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;
  int __pyx_t_3;
  int __pyx_t_4;
  int __pyx_t_5;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_t_3 = (((int)(((long *)__pyx_v_ip0)[0])) == (((long *)__pyx_v_ip0)[0]));
    if (__pyx_t_3) {
      __pyx_t_4 = (((int)(((long *)__pyx_v_ip1)[0])) == (((long *)__pyx_v_ip1)[0]));
      __pyx_t_5 = __pyx_t_4;
    } else {
      __pyx_t_5 = __pyx_t_3;
    }
    if (__pyx_t_5) {

      
      __pyx_v_ov0 = ((double (*)(int, int, double))__pyx_v_func)(((int)(((long *)__pyx_v_ip0)[0])), ((int)(((long *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])));
      goto __pyx_L5;
    }
     {

      
      sf_error(__pyx_v_func_name, SF_ERROR_DOMAIN, __pyx_k_1);

      
      __pyx_v_ov0 = ((double)NPY_NAN);
    }
    __pyx_L5:;

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_lddd__As_lddd_d(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(long, double, double, double))__pyx_v_func)(((long)(((long *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])), ((double)(((double *)__pyx_v_ip3)[0])));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ffF_F(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(double, double, __pyx_t_double_complex))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_float_complex *)__pyx_v_ip2)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_float_complex *)__pyx_v_ip2)[0])))));

    
    (((__pyx_t_float_complex *)__pyx_v_op0)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov0)), ((float)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DDDD_As_F_FFFF(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  char *__pyx_v_op2;
  char *__pyx_v_op3;
  __pyx_t_double_complex __pyx_v_ov0;
  __pyx_t_double_complex __pyx_v_ov1;
  __pyx_t_double_complex __pyx_v_ov2;
  __pyx_t_double_complex __pyx_v_ov3;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);
  __pyx_v_op2 = (__pyx_v_args[3]);
  __pyx_v_op3 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(__pyx_t_double_complex, __pyx_t_double_complex *, __pyx_t_double_complex *, __pyx_t_double_complex *, __pyx_t_double_complex *))__pyx_v_func)(__pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_float_complex *)__pyx_v_ip0)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_float_complex *)__pyx_v_ip0)[0])))), (&__pyx_v_ov0), (&__pyx_v_ov1), (&__pyx_v_ov2), (&__pyx_v_ov3));

    
    (((__pyx_t_float_complex *)__pyx_v_op0)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov0)), ((float)__Pyx_CIMAG(__pyx_v_ov0)));

    
    (((__pyx_t_float_complex *)__pyx_v_op1)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov1)), ((float)__Pyx_CIMAG(__pyx_v_ov1)));

    
    (((__pyx_t_float_complex *)__pyx_v_op2)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov2)), ((float)__Pyx_CIMAG(__pyx_v_ov2)));

    
    (((__pyx_t_float_complex *)__pyx_v_op3)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov3)), ((float)__Pyx_CIMAG(__pyx_v_ov3)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));

    
    __pyx_v_op2 = (__pyx_v_op2 + (__pyx_v_steps[3]));

    
    __pyx_v_op3 = (__pyx_v_op3 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_lddd__As_lfff_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(long, double, double, double))__pyx_v_func)(((long)(((long *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])), ((double)(((float *)__pyx_v_ip3)[0])));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_ip3;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_ip3 = (__pyx_v_args[3]);
  __pyx_v_op0 = (__pyx_v_args[4]);
  __pyx_v_op1 = (__pyx_v_args[5]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double, double, double, double, double *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])), ((double)(((double *)__pyx_v_ip3)[0])), (&__pyx_v_ov1));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    (((double *)__pyx_v_op1)[0]) = ((double)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_ip3 = (__pyx_v_ip3 + (__pyx_v_steps[3]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[4]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[5]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DD_As_D_DD(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  __pyx_t_double_complex __pyx_v_ov0;
  __pyx_t_double_complex __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(__pyx_t_double_complex, __pyx_t_double_complex *, __pyx_t_double_complex *))__pyx_v_func)(__pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_double_complex *)__pyx_v_ip0)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_double_complex *)__pyx_v_ip0)[0])))), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((__pyx_t_double_complex *)__pyx_v_op0)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov0)), ((double)__Pyx_CIMAG(__pyx_v_ov0)));

    
    (((__pyx_t_double_complex *)__pyx_v_op1)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov1)), ((double)__Pyx_CIMAG(__pyx_v_ov1)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_D_Dld__As_Flf_F(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  __pyx_t_double_complex __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((__pyx_t_double_complex (*)(__pyx_t_double_complex, long, double))__pyx_v_func)(__pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_float_complex *)__pyx_v_ip0)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_float_complex *)__pyx_v_ip0)[0])))), ((long)(((long *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])));

    
    (((__pyx_t_float_complex *)__pyx_v_op0)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov0)), ((float)__Pyx_CIMAG(__pyx_v_ov0)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dddd_As_dd_dddd(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  char *__pyx_v_op2;
  char *__pyx_v_op3;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  double __pyx_v_ov2;
  double __pyx_v_ov3;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);
  __pyx_v_op1 = (__pyx_v_args[3]);
  __pyx_v_op2 = (__pyx_v_args[4]);
  __pyx_v_op3 = (__pyx_v_args[5]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double, double *, double *, double *, double *))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1), (&__pyx_v_ov2), (&__pyx_v_ov3));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    (((double *)__pyx_v_op1)[0]) = ((double)__pyx_v_ov1);

    
    (((double *)__pyx_v_op2)[0]) = ((double)__pyx_v_ov2);

    
    (((double *)__pyx_v_op3)[0]) = ((double)__pyx_v_ov3);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[3]));

    
    __pyx_v_op2 = (__pyx_v_op2 + (__pyx_v_steps[4]));

    
    __pyx_v_op3 = (__pyx_v_op3 + (__pyx_v_steps[5]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double, double))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_f_ff(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  double __pyx_v_ov0;
  double __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, double *, double *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    (((float *)__pyx_v_op1)[0]) = ((float)__pyx_v_ov1);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(long, double))__pyx_v_func)(((long)(((long *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;
  int __pyx_t_3;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_t_3 = (((int)(((long *)__pyx_v_ip0)[0])) == (((long *)__pyx_v_ip0)[0]));
    if (__pyx_t_3) {

      
      __pyx_v_ov0 = ((double (*)(int, double))__pyx_v_func)(((int)(((long *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])));
      goto __pyx_L5;
    }
     {

      
      sf_error(__pyx_v_func_name, SF_ERROR_DOMAIN, __pyx_k_1);

      
      __pyx_v_ov0 = ((double)NPY_NAN);
    }
    __pyx_L5:;

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double, double))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DD_As_f_FF(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  __pyx_t_double_complex __pyx_v_ov0;
  __pyx_t_double_complex __pyx_v_ov1;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(double, __pyx_t_double_complex *, __pyx_t_double_complex *))__pyx_v_func)(((double)(((float *)__pyx_v_ip0)[0])), (&__pyx_v_ov0), (&__pyx_v_ov1));

    
    (((__pyx_t_float_complex *)__pyx_v_op0)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov0)), ((float)__Pyx_CIMAG(__pyx_v_ov0)));

    
    (((__pyx_t_float_complex *)__pyx_v_op1)[0]) = __pyx_t_float_complex_from_parts(((float)__Pyx_CREAL(__pyx_v_ov1)), ((float)__Pyx_CIMAG(__pyx_v_ov1)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(double, double, double))__pyx_v_func)(((double)(((double *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])), ((double)(((double *)__pyx_v_ip2)[0])));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DDDD_As_D_DDDD(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  char *__pyx_v_op1;
  char *__pyx_v_op2;
  char *__pyx_v_op3;
  __pyx_t_double_complex __pyx_v_ov0;
  __pyx_t_double_complex __pyx_v_ov1;
  __pyx_t_double_complex __pyx_v_ov2;
  __pyx_t_double_complex __pyx_v_ov3;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);
  __pyx_v_op1 = (__pyx_v_args[2]);
  __pyx_v_op2 = (__pyx_v_args[3]);
  __pyx_v_op3 = (__pyx_v_args[4]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    ((int (*)(__pyx_t_double_complex, __pyx_t_double_complex *, __pyx_t_double_complex *, __pyx_t_double_complex *, __pyx_t_double_complex *))__pyx_v_func)(__pyx_t_double_complex_from_parts(((double)__Pyx_CREAL((((__pyx_t_double_complex *)__pyx_v_ip0)[0]))), ((double)__Pyx_CIMAG((((__pyx_t_double_complex *)__pyx_v_ip0)[0])))), (&__pyx_v_ov0), (&__pyx_v_ov1), (&__pyx_v_ov2), (&__pyx_v_ov3));

    
    (((__pyx_t_double_complex *)__pyx_v_op0)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov0)), ((double)__Pyx_CIMAG(__pyx_v_ov0)));

    
    (((__pyx_t_double_complex *)__pyx_v_op1)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov1)), ((double)__Pyx_CIMAG(__pyx_v_ov1)));

    
    (((__pyx_t_double_complex *)__pyx_v_op2)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov2)), ((double)__Pyx_CIMAG(__pyx_v_ov2)));

    
    (((__pyx_t_double_complex *)__pyx_v_op3)[0]) = __pyx_t_double_complex_from_parts(((double)__Pyx_CREAL(__pyx_v_ov3)), ((double)__Pyx_CIMAG(__pyx_v_ov3)));

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));

    
    __pyx_v_op1 = (__pyx_v_op1 + (__pyx_v_steps[2]));

    
    __pyx_v_op2 = (__pyx_v_op2 + (__pyx_v_steps[3]));

    
    __pyx_v_op3 = (__pyx_v_op3 + (__pyx_v_steps[4]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_g_g__As_g_g(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_op0;
  long double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_op0 = (__pyx_v_args[1]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((long double (*)(long double))__pyx_v_func)(((long double)(((long double *)__pyx_v_ip0)[0])));

    
    (((long double *)__pyx_v_op0)[0]) = ((long double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[1]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ldd__As_lff_f(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_ip2;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_ip2 = (__pyx_v_args[2]);
  __pyx_v_op0 = (__pyx_v_args[3]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(long, double, double))__pyx_v_func)(((long)(((long *)__pyx_v_ip0)[0])), ((double)(((float *)__pyx_v_ip1)[0])), ((double)(((float *)__pyx_v_ip2)[0])));

    
    (((float *)__pyx_v_op0)[0]) = ((float)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_ip2 = (__pyx_v_ip2 + (__pyx_v_steps[2]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[3]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



static void __pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d(char **__pyx_v_args, npy_intp *__pyx_v_dims, npy_intp *__pyx_v_steps, void *__pyx_v_data) {
  CYTHON_UNUSED npy_intp __pyx_v_i;
  npy_intp __pyx_v_n;
  void *__pyx_v_func;
  char *__pyx_v_func_name;
  char *__pyx_v_ip0;
  char *__pyx_v_ip1;
  char *__pyx_v_op0;
  double __pyx_v_ov0;
  npy_intp __pyx_t_1;
  npy_intp __pyx_t_2;

  
  __pyx_v_n = (__pyx_v_dims[0]);

  
  __pyx_v_func = (((void **)__pyx_v_data)[0]);

  
  __pyx_v_func_name = ((char *)(((void **)__pyx_v_data)[1]));

  
  __pyx_v_ip0 = (__pyx_v_args[0]);
  __pyx_v_ip1 = (__pyx_v_args[1]);
  __pyx_v_op0 = (__pyx_v_args[2]);

  
  __pyx_t_1 = __pyx_v_n;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_v_ov0 = ((double (*)(long, double))__pyx_v_func)(((long)(((long *)__pyx_v_ip0)[0])), ((double)(((double *)__pyx_v_ip1)[0])));

    
    (((double *)__pyx_v_op0)[0]) = ((double)__pyx_v_ov0);

    
    __pyx_v_ip0 = (__pyx_v_ip0 + (__pyx_v_steps[0]));

    
    __pyx_v_ip1 = (__pyx_v_ip1 + (__pyx_v_steps[1]));

    
    __pyx_v_op0 = (__pyx_v_op0 + (__pyx_v_steps[2]));
  }

  
  sf_error_check_fpe(__pyx_v_func_name);

}



int wrap_PyUFunc_getfperr(void) {
  int __pyx_r;

  
  __pyx_r = PyUFunc_getfperr();
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}


static PyObject *__pyx_pw_5scipy_7special_7_ufuncs_1_errprint(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); 
static char __pyx_doc_5scipy_7special_7_ufuncs__errprint[] = "\n    errprint(flag=None)\n\n    Sets or returns the error printing flag for special functions.\n\n    Parameters\n    ----------\n    flag : bool, optional\n        Whether warnings concerning evaluation of special functions in\n        scipy.special are shown. If omitted, no change is made to the\n        current setting.\n\n    Returns\n    -------\n    old_flag\n        Previous value of the error flag\n\n    ";
static PyMethodDef __pyx_mdef_5scipy_7special_7_ufuncs_1_errprint = {__Pyx_NAMESTR("_errprint"), (PyCFunction)__pyx_pw_5scipy_7special_7_ufuncs_1_errprint, METH_VARARGS|METH_KEYWORDS, __Pyx_DOCSTR(__pyx_doc_5scipy_7special_7_ufuncs__errprint)};
static PyObject *__pyx_pw_5scipy_7special_7_ufuncs_1_errprint(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyObject *__pyx_v_inflag = 0;
  PyObject *__pyx_r = 0;
  __Pyx_RefNannyDeclarations
  __Pyx_RefNannySetupContext("_errprint (wrapper)", 0);
  {
    static PyObject **__pyx_pyargnames[] = {&__pyx_n_s__inflag,0};
    PyObject* values[1] = {0};

    
    values[0] = ((PyObject *)Py_None);
    if (unlikely(__pyx_kwds)) {
      Py_ssize_t kw_args;
      const Py_ssize_t pos_args = PyTuple_GET_SIZE(__pyx_args);
      switch (pos_args) {
        case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
        case  0: break;
        default: goto __pyx_L5_argtuple_error;
      }
      kw_args = PyDict_Size(__pyx_kwds);
      switch (pos_args) {
        case  0:
        if (kw_args > 0) {
          PyObject* value = PyDict_GetItem(__pyx_kwds, __pyx_n_s__inflag);
          if (value) { values[0] = value; kw_args--; }
        }
      }
      if (unlikely(kw_args > 0)) {
        if (unlikely(__Pyx_ParseOptionalKeywords(__pyx_kwds, __pyx_pyargnames, 0, values, pos_args, "_errprint") < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7228; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
      }
    } else {
      switch (PyTuple_GET_SIZE(__pyx_args)) {
        case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
        case  0: break;
        default: goto __pyx_L5_argtuple_error;
      }
    }
    __pyx_v_inflag = values[0];
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L5_argtuple_error:;
  __Pyx_RaiseArgtupleInvalid("_errprint", 0, 0, 1, PyTuple_GET_SIZE(__pyx_args)); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7228; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  __pyx_L3_error:;
  __Pyx_AddTraceback("scipy.special._ufuncs._errprint", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __Pyx_RefNannyFinishContext();
  return NULL;
  __pyx_L4_argument_unpacking_done:;
  __pyx_r = __pyx_pf_5scipy_7special_7_ufuncs__errprint(__pyx_self, __pyx_v_inflag);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}

static PyObject *__pyx_pf_5scipy_7special_7_ufuncs__errprint(CYTHON_UNUSED PyObject *__pyx_self, PyObject *__pyx_v_inflag) {
  PyObject *__pyx_r = NULL;
  __Pyx_RefNannyDeclarations
  int __pyx_t_1;
  PyObject *__pyx_t_2 = NULL;
  int __pyx_lineno = 0;
  const char *__pyx_filename = NULL;
  int __pyx_clineno = 0;
  __Pyx_RefNannySetupContext("_errprint", 0);

  
  __pyx_t_1 = (__pyx_v_inflag != Py_None);
  if (__pyx_t_1) {

    
    __Pyx_XDECREF(__pyx_r);
    __pyx_t_1 = __Pyx_PyObject_IsTrue(__pyx_v_inflag); if (unlikely(__pyx_t_1 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7248; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __pyx_t_2 = PyInt_FromLong(sf_error_set_print(((int)(!(!__pyx_t_1))))); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7248; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_2);
    __pyx_r = __pyx_t_2;
    __pyx_t_2 = 0;
    goto __pyx_L0;
    goto __pyx_L3;
  }
   {

    
    __Pyx_XDECREF(__pyx_r);
    __pyx_t_2 = PyInt_FromLong(sf_error_get_print()); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7250; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_2);
    __pyx_r = __pyx_t_2;
    __pyx_t_2 = 0;
    goto __pyx_L0;
  }
  __pyx_L3:;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_AddTraceback("scipy.special._ufuncs._errprint", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = NULL;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE int __pyx_f_5scipy_7special_13_complexstuff_zisnan(__pyx_t_double_complex __pyx_v_x) {
  int __pyx_r;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;

  
  __pyx_t_1 = npy_isnan(__Pyx_CREAL(__pyx_v_x));
  if (!__pyx_t_1) {
    __pyx_t_2 = npy_isnan(__Pyx_CIMAG(__pyx_v_x));
    __pyx_t_3 = __pyx_t_2;
  } else {
    __pyx_t_3 = __pyx_t_1;
  }
  __pyx_r = __pyx_t_3;
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_13_complexstuff_zabs(__pyx_t_double_complex __pyx_v_x) {
  double __pyx_v_r;
  double __pyx_r;

  
  __pyx_v_r = npy_cabs((((npy_cdouble *)(&__pyx_v_x))[0]));

  
  __pyx_r = __pyx_v_r;
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_f_5scipy_7special_13_complexstuff_zlog(__pyx_t_double_complex __pyx_v_x) {
  npy_cdouble __pyx_v_r;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_r = npy_clog((((npy_cdouble *)(&__pyx_v_x))[0]));

  
  __pyx_r = (((__pyx_t_double_complex *)(&__pyx_v_r))[0]);
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_f_5scipy_7special_13_complexstuff_zexp(__pyx_t_double_complex __pyx_v_x) {
  npy_cdouble __pyx_v_r;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_r = npy_cexp((((npy_cdouble *)(&__pyx_v_x))[0]));

  
  __pyx_r = (((__pyx_t_double_complex *)(&__pyx_v_r))[0]);
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}

static CYTHON_INLINE void __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(char *__pyx_v_func_name, double __pyx_v_x, double __pyx_v_y) {
  PyObject *__pyx_v_msg = NULL;
  __Pyx_RefNannyDeclarations
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  char *__pyx_t_6;
  int __pyx_lineno = 0;
  const char *__pyx_filename = NULL;
  int __pyx_clineno = 0;
  #ifdef WITH_THREAD
  PyGILState_STATE __pyx_gilstate_save;
  #endif
  __Pyx_RefNannySetupContext("_legacy_cast_check", 1);
   {

    
    __pyx_t_1 = (((int)__pyx_v_x) != __pyx_v_x);
    if (!__pyx_t_1) {
      __pyx_t_2 = (((int)__pyx_v_y) != __pyx_v_y);
      __pyx_t_3 = __pyx_t_2;
    } else {
      __pyx_t_3 = __pyx_t_1;
    }
    if (__pyx_t_3) {

      
      {
          #ifdef WITH_THREAD
          PyGILState_STATE __pyx_gilstate_save = PyGILState_Ensure();
          #endif
           {

            
            __pyx_t_4 = PyBytes_FromString(__pyx_v_func_name); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 37; __pyx_clineno = __LINE__; goto __pyx_L8;}
            __Pyx_GOTREF(((PyObject *)__pyx_t_4));
            __pyx_t_5 = PyNumber_Add(((PyObject *)__pyx_kp_b_2), ((PyObject *)__pyx_t_4)); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 36; __pyx_clineno = __LINE__; goto __pyx_L8;}
            __Pyx_GOTREF(((PyObject *)__pyx_t_5));
            __Pyx_DECREF(((PyObject *)__pyx_t_4)); __pyx_t_4 = 0;
            __pyx_t_4 = PyNumber_Add(((PyObject *)__pyx_t_5), ((PyObject *)__pyx_kp_b_3)); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 37; __pyx_clineno = __LINE__; goto __pyx_L8;}
            __Pyx_GOTREF(((PyObject *)__pyx_t_4));
            __Pyx_DECREF(((PyObject *)__pyx_t_5)); __pyx_t_5 = 0;
            __pyx_t_5 = PyNumber_Add(((PyObject *)__pyx_t_4), ((PyObject *)__pyx_kp_b_4)); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 38; __pyx_clineno = __LINE__; goto __pyx_L8;}
            __Pyx_GOTREF(((PyObject *)__pyx_t_5));
            __Pyx_DECREF(((PyObject *)__pyx_t_4)); __pyx_t_4 = 0;
            __pyx_v_msg = __pyx_t_5;
            __pyx_t_5 = 0;

            
            __pyx_t_5 = __pyx_builtin_FutureWarning;
            __Pyx_INCREF(__pyx_t_5);
            __pyx_t_6 = PyBytes_AsString(((PyObject *)__pyx_v_msg)); if (unlikely((!__pyx_t_6) && PyErr_Occurred())) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 39; __pyx_clineno = __LINE__; goto __pyx_L8;}
            PyErr_WarnEx(__pyx_t_5, __pyx_t_6, 1);
            __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
          }

          
           {
            int __pyx_why;
            __pyx_why = 0; goto __pyx_L9;
            __pyx_L8: __pyx_why = 4; goto __pyx_L9;
            __pyx_L9:;
            #ifdef WITH_THREAD
            PyGILState_Release(__pyx_gilstate_save);
            #endif
            switch (__pyx_why) {
              case 4: goto __pyx_L4;
            }
          }
      }
      goto __pyx_L6;
    }
    __pyx_L6:;
  }

  
   {
    int __pyx_why;
    __pyx_why = 0; goto __pyx_L5;
    __pyx_L4: __pyx_why = 4; goto __pyx_L5;
    __pyx_L5:;
    #ifdef WITH_THREAD
    __pyx_gilstate_save = PyGILState_Ensure();
    #endif
    switch (__pyx_why) {
      case 4: goto __pyx_L1_error;
    }
  }

  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_WriteUnraisable("scipy.special._legacy._legacy_cast_check", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_L0:;
  __Pyx_XDECREF(__pyx_v_msg);
  #ifdef WITH_THREAD
  PyGILState_Release(__pyx_gilstate_save);
  #endif
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_bdtrc_unsafe(double __pyx_v_k, double __pyx_v_n, double __pyx_v_p) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__bdtrc, __pyx_v_k, __pyx_v_n);

  
  __pyx_r = bdtrc(((int)__pyx_v_k), ((int)__pyx_v_n), __pyx_v_p);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_bdtr_unsafe(double __pyx_v_k, double __pyx_v_n, double __pyx_v_p) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__bdtr, __pyx_v_k, __pyx_v_n);

  
  __pyx_r = bdtr(((int)__pyx_v_k), ((int)__pyx_v_n), __pyx_v_p);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_bdtri_unsafe(double __pyx_v_k, double __pyx_v_n, double __pyx_v_y) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__bdtri, __pyx_v_k, __pyx_v_n);

  
  __pyx_r = bdtri(((int)__pyx_v_k), ((int)__pyx_v_n), __pyx_v_y);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_expn_unsafe(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__expn, __pyx_v_n, 0.0);

  
  __pyx_r = expn(((int)__pyx_v_n), __pyx_v_x);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_hyp2f0_unsafe(double __pyx_v_a, double __pyx_v_b, double __pyx_v_x, double __pyx_v_type, double *__pyx_v_err) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__hyp2f0, __pyx_v_type, 0.0);

  
  __pyx_r = hyp2f0(__pyx_v_a, __pyx_v_b, __pyx_v_x, ((int)__pyx_v_type), __pyx_v_err);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_nbdtrc_unsafe(double __pyx_v_k, double __pyx_v_n, double __pyx_v_p) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__nbdtrc, __pyx_v_k, __pyx_v_n);

  
  __pyx_r = nbdtrc(((int)__pyx_v_k), ((int)__pyx_v_n), __pyx_v_p);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_nbdtr_unsafe(double __pyx_v_k, double __pyx_v_n, double __pyx_v_p) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__nbdtr, __pyx_v_k, __pyx_v_n);

  
  __pyx_r = nbdtr(((int)__pyx_v_k), ((int)__pyx_v_n), __pyx_v_p);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_nbdtri_unsafe(double __pyx_v_k, double __pyx_v_n, double __pyx_v_p) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__nbdtri, __pyx_v_k, __pyx_v_n);

  
  __pyx_r = nbdtri(((int)__pyx_v_k), ((int)__pyx_v_n), __pyx_v_p);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_pdtrc_unsafe(double __pyx_v_k, double __pyx_v_m) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__pdtrc, __pyx_v_k, 0.0);

  
  __pyx_r = pdtrc(((int)__pyx_v_k), __pyx_v_m);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_pdtr_unsafe(double __pyx_v_k, double __pyx_v_m) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__pdtr, __pyx_v_k, 0.0);

  
  __pyx_r = pdtr(((int)__pyx_v_k), __pyx_v_m);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_pdtri_unsafe(double __pyx_v_k, double __pyx_v_y) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__pdtri, __pyx_v_k, 0.0);

  
  __pyx_r = pdtri(((int)__pyx_v_k), __pyx_v_y);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_kn_unsafe(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__kn, __pyx_v_n, 0.0);

  
  __pyx_r = kn(((int)__pyx_v_n), __pyx_v_x);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_yn_unsafe(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__yn, __pyx_v_n, 0.0);

  
  __pyx_r = yn(((int)__pyx_v_n), __pyx_v_x);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_smirnov_unsafe(double __pyx_v_n, double __pyx_v_e) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__smirnov, __pyx_v_n, 0.0);

  
  __pyx_r = smirnov(((int)__pyx_v_n), __pyx_v_e);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_7_legacy_smirnovi_unsafe(double __pyx_v_n, double __pyx_v_p) {
  double __pyx_r;

  
  __pyx_f_5scipy_7special_7_legacy__legacy_cast_check(__pyx_k__smirnovi, __pyx_v_n, 0.0);

  
  __pyx_r = smirnovi(((int)__pyx_v_n), __pyx_v_p);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(double __pyx_v_a, double __pyx_v_b, double __pyx_v_c, double __pyx_v_z) {
  double __pyx_r;

  
  __pyx_r = hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_z);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(double __pyx_v_a, double __pyx_v_b, double __pyx_v_c, __pyx_t_double_complex __pyx_v_z) {
  npy_cdouble __pyx_v_r;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_r = chyp2f1_wrap(__pyx_v_a, __pyx_v_b, __pyx_v_c, (((npy_cdouble *)(&__pyx_v_z))[0]));

  
  __pyx_r = (((__pyx_t_double_complex *)(&__pyx_v_r))[0]);
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp1f1(double __pyx_v_a, double __pyx_v_b, double __pyx_v_z) {
  double __pyx_r;

  
  __pyx_r = hyp1f1_wrap(__pyx_v_a, __pyx_v_b, __pyx_v_z);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp1f1(double __pyx_v_a, double __pyx_v_b, __pyx_t_double_complex __pyx_v_z) {
  npy_cdouble __pyx_v_r;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_r = chyp1f1_wrap(__pyx_v_a, __pyx_v_b, (((npy_cdouble *)(&__pyx_v_z))[0]));

  
  __pyx_r = (((__pyx_t_double_complex *)(&__pyx_v_r))[0]);
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_binom(double __pyx_v_n, double __pyx_v_k) {
  double __pyx_v_kx;
  double __pyx_v_nx;
  double __pyx_v_num;
  double __pyx_v_den;
  int __pyx_v_i;
  double __pyx_r;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  int __pyx_t_4;
  long __pyx_t_5;
  int __pyx_t_6;

  
  __pyx_t_1 = (__pyx_v_n < 0.0);
  if (__pyx_t_1) {

    
    __pyx_v_nx = floor(__pyx_v_n);

    
    __pyx_t_1 = (__pyx_v_n == __pyx_v_nx);
    if (__pyx_t_1) {

      
      __pyx_r = NPY_NAN;
      goto __pyx_L0;
      goto __pyx_L4;
    }
    __pyx_L4:;
    goto __pyx_L3;
  }
  __pyx_L3:;

  
  __pyx_v_kx = floor(__pyx_v_k);

  
  __pyx_t_1 = (__pyx_v_k == __pyx_v_kx);
  if (__pyx_t_1) {

    
    __pyx_v_nx = floor(__pyx_v_n);

    
    __pyx_t_1 = (__pyx_v_nx == __pyx_v_n);
    if (__pyx_t_1) {
      __pyx_t_2 = (__pyx_v_kx > (__pyx_v_nx / 2.0));
      if (__pyx_t_2) {
        __pyx_t_3 = (__pyx_v_nx > 0.0);
        __pyx_t_4 = __pyx_t_3;
      } else {
        __pyx_t_4 = __pyx_t_2;
      }
      __pyx_t_2 = __pyx_t_4;
    } else {
      __pyx_t_2 = __pyx_t_1;
    }
    if (__pyx_t_2) {

      
      __pyx_v_kx = (__pyx_v_nx - __pyx_v_kx);
      goto __pyx_L6;
    }
    __pyx_L6:;

    
    __pyx_t_2 = (__pyx_v_kx >= 1.0);
    if (__pyx_t_2) {
      __pyx_t_1 = (__pyx_v_kx < 20.0);
      __pyx_t_4 = __pyx_t_1;
    } else {
      __pyx_t_4 = __pyx_t_2;
    }
    if (__pyx_t_4) {

      
      __pyx_v_num = 1.0;

      
      __pyx_v_den = 1.0;

      
      __pyx_t_5 = (1 + ((int)__pyx_v_kx));
      for (__pyx_t_6 = 1; __pyx_t_6 < __pyx_t_5; __pyx_t_6+=1) {
        __pyx_v_i = __pyx_t_6;

        
        __pyx_v_num = (__pyx_v_num * ((__pyx_v_i + __pyx_v_n) - __pyx_v_kx));

        
        __pyx_v_den = (__pyx_v_den * __pyx_v_i);

        
        __pyx_t_4 = (fabs(__pyx_v_num) > 1e50);
        if (__pyx_t_4) {

          
          __pyx_v_num = (__pyx_v_num / __pyx_v_den);

          
          __pyx_v_den = 1.0;
          goto __pyx_L10;
        }
        __pyx_L10:;
      }

      
      __pyx_r = (__pyx_v_num / __pyx_v_den);
      goto __pyx_L0;
      goto __pyx_L7;
    }
    __pyx_L7:;
    goto __pyx_L5;
  }
  __pyx_L5:;

  
  __pyx_r = ((1.0 / beta(((1.0 + __pyx_v_n) - __pyx_v_k), (1.0 + __pyx_v_k))) / (__pyx_v_n + 1.0));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi(double __pyx_v_n, double __pyx_v_alpha, double __pyx_v_beta, double __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  double __pyx_v_d;
  double __pyx_v_g;
  double __pyx_r;

  
  __pyx_v_d = __pyx_f_5scipy_7special_15orthogonal_eval_binom((__pyx_v_n + __pyx_v_alpha), __pyx_v_n);

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (((__pyx_v_n + __pyx_v_alpha) + __pyx_v_beta) + 1.0);

  
  __pyx_v_c = (__pyx_v_alpha + 1.0);

  
  __pyx_v_g = (0.5 * (1.0 - __pyx_v_x));

  
  __pyx_r = (__pyx_v_d * __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi(double __pyx_v_n, double __pyx_v_alpha, double __pyx_v_beta, __pyx_t_double_complex __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  double __pyx_v_d;
  __pyx_t_double_complex __pyx_v_g;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_d = __pyx_f_5scipy_7special_15orthogonal_eval_binom((__pyx_v_n + __pyx_v_alpha), __pyx_v_n);

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (((__pyx_v_n + __pyx_v_alpha) + __pyx_v_beta) + 1.0);

  
  __pyx_v_c = (__pyx_v_alpha + 1.0);

  
  __pyx_v_g = __Pyx_c_prod(__pyx_t_double_complex_from_parts(0.5, 0), __Pyx_c_diff(__pyx_t_double_complex_from_parts(1, 0), __pyx_v_x));

  
  __pyx_r = __Pyx_c_prod(__pyx_t_double_complex_from_parts(__pyx_v_d, 0), __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi_l(long __pyx_v_n, double __pyx_v_alpha, double __pyx_v_beta, double __pyx_v_x) {
  long __pyx_v_kk;
  double __pyx_v_p;
  double __pyx_v_d;
  double __pyx_v_k;
  double __pyx_v_t;
  double __pyx_r;
  int __pyx_t_1;
  long __pyx_t_2;
  long __pyx_t_3;

  
  __pyx_t_1 = (__pyx_v_n < 0);
  if (__pyx_t_1) {

    
    __pyx_r = 0.0;
    goto __pyx_L0;
    goto __pyx_L3;
  }

  
  __pyx_t_1 = (__pyx_v_n == 0);
  if (__pyx_t_1) {

    
    __pyx_r = 1.0;
    goto __pyx_L0;
    goto __pyx_L3;
  }

  
  __pyx_t_1 = (__pyx_v_n == 1);
  if (__pyx_t_1) {

    
    __pyx_r = (0.5 * ((2.0 * (__pyx_v_alpha + 1.0)) + (((__pyx_v_alpha + __pyx_v_beta) + 2.0) * (__pyx_v_x - 1.0))));
    goto __pyx_L0;
    goto __pyx_L3;
  }
   {

    
    __pyx_v_d = ((((__pyx_v_alpha + __pyx_v_beta) + 2.0) * (__pyx_v_x - 1.0)) / (2.0 * (__pyx_v_alpha + 1.0)));

    
    __pyx_v_p = (__pyx_v_d + 1.0);

    
    __pyx_t_2 = (__pyx_v_n - 1);
    for (__pyx_t_3 = 0; __pyx_t_3 < __pyx_t_2; __pyx_t_3+=1) {
      __pyx_v_kk = __pyx_t_3;

      
      __pyx_v_k = (__pyx_v_kk + 1.0);

      
      __pyx_v_t = (((2.0 * __pyx_v_k) + __pyx_v_alpha) + __pyx_v_beta);

      
      __pyx_v_d = ((((((__pyx_v_t * (__pyx_v_t + 1.0)) * (__pyx_v_t + 2.0)) * (__pyx_v_x - 1.0)) * __pyx_v_p) + ((((2.0 * __pyx_v_k) * (__pyx_v_k + __pyx_v_beta)) * (__pyx_v_t + 2.0)) * __pyx_v_d)) / (((2.0 * ((__pyx_v_k + __pyx_v_alpha) + 1.0)) * (((__pyx_v_k + __pyx_v_alpha) + __pyx_v_beta) + 1.0)) * __pyx_v_t));

      
      __pyx_v_p = (__pyx_v_d + __pyx_v_p);
    }

    
    __pyx_r = (__pyx_f_5scipy_7special_15orthogonal_eval_binom((__pyx_v_n + __pyx_v_alpha), __pyx_v_n) * __pyx_v_p);
    goto __pyx_L0;
  }
  __pyx_L3:;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi(double __pyx_v_n, double __pyx_v_p, double __pyx_v_q, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = (__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi(__pyx_v_n, (__pyx_v_p - __pyx_v_q), (__pyx_v_q - 1.0), ((2.0 * __pyx_v_x) - 1.0)) / __pyx_f_5scipy_7special_15orthogonal_eval_binom((((2.0 * __pyx_v_n) + __pyx_v_p) - 1.0), __pyx_v_n));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi(double __pyx_v_n, double __pyx_v_p, double __pyx_v_q, __pyx_t_double_complex __pyx_v_x) {
  __pyx_t_double_complex __pyx_r;

  
  __pyx_r = __Pyx_c_quot(__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi(__pyx_v_n, (__pyx_v_p - __pyx_v_q), (__pyx_v_q - 1.0), __Pyx_c_diff(__Pyx_c_prod(__pyx_t_double_complex_from_parts(2, 0), __pyx_v_x), __pyx_t_double_complex_from_parts(1, 0))), __pyx_t_double_complex_from_parts(__pyx_f_5scipy_7special_15orthogonal_eval_binom((((2.0 * __pyx_v_n) + __pyx_v_p) - 1.0), __pyx_v_n), 0));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi_l(long __pyx_v_n, double __pyx_v_p, double __pyx_v_q, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = (__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi_l(__pyx_v_n, (__pyx_v_p - __pyx_v_q), (__pyx_v_q - 1.0), ((2.0 * __pyx_v_x) - 1.0)) / __pyx_f_5scipy_7special_15orthogonal_eval_binom((((2 * __pyx_v_n) + __pyx_v_p) - 1.0), __pyx_v_n));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer(double __pyx_v_n, double __pyx_v_alpha, double __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  double __pyx_v_d;
  double __pyx_v_g;
  double __pyx_r;

  
  __pyx_v_d = ((Gamma((__pyx_v_n + (2.0 * __pyx_v_alpha))) / Gamma((1.0 + __pyx_v_n))) / Gamma((2.0 * __pyx_v_alpha)));

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (__pyx_v_n + (2.0 * __pyx_v_alpha));

  
  __pyx_v_c = (__pyx_v_alpha + 0.5);

  
  __pyx_v_g = ((1.0 - __pyx_v_x) / 2.0);

  
  __pyx_r = (__pyx_v_d * __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer(double __pyx_v_n, double __pyx_v_alpha, __pyx_t_double_complex __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  double __pyx_v_d;
  __pyx_t_double_complex __pyx_v_g;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_d = ((Gamma((__pyx_v_n + (2.0 * __pyx_v_alpha))) / Gamma((1.0 + __pyx_v_n))) / Gamma((2.0 * __pyx_v_alpha)));

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (__pyx_v_n + (2.0 * __pyx_v_alpha));

  
  __pyx_v_c = (__pyx_v_alpha + 0.5);

  
  __pyx_v_g = __Pyx_c_quot(__Pyx_c_diff(__pyx_t_double_complex_from_parts(1, 0), __pyx_v_x), __pyx_t_double_complex_from_parts(2.0, 0));

  
  __pyx_r = __Pyx_c_prod(__pyx_t_double_complex_from_parts(__pyx_v_d, 0), __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer_l(long __pyx_v_n, double __pyx_v_alpha, double __pyx_v_x) {
  long __pyx_v_kk;
  double __pyx_v_p;
  double __pyx_v_d;
  double __pyx_v_k;
  double __pyx_r;
  int __pyx_t_1;
  long __pyx_t_2;
  long __pyx_t_3;

  
  __pyx_t_1 = (__pyx_v_n < 0);
  if (__pyx_t_1) {

    
    __pyx_r = 0.0;
    goto __pyx_L0;
    goto __pyx_L3;
  }

  
  __pyx_t_1 = (__pyx_v_n == 0);
  if (__pyx_t_1) {

    
    __pyx_r = 1.0;
    goto __pyx_L0;
    goto __pyx_L3;
  }

  
  __pyx_t_1 = (__pyx_v_n == 1);
  if (__pyx_t_1) {

    
    __pyx_r = ((2.0 * __pyx_v_alpha) * __pyx_v_x);
    goto __pyx_L0;
    goto __pyx_L3;
  }

  
  __pyx_t_1 = (__pyx_v_alpha == 0.0);
  if (__pyx_t_1) {

    
    __pyx_r = __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer(__pyx_v_n, __pyx_v_alpha, __pyx_v_x);
    goto __pyx_L0;
    goto __pyx_L3;
  }
   {

    
    __pyx_v_d = (__pyx_v_x - 1.0);

    
    __pyx_v_p = __pyx_v_x;

    
    __pyx_t_2 = (__pyx_v_n - 1);
    for (__pyx_t_3 = 0; __pyx_t_3 < __pyx_t_2; __pyx_t_3+=1) {
      __pyx_v_kk = __pyx_t_3;

      
      __pyx_v_k = (__pyx_v_kk + 1.0);

      
      __pyx_v_d = (((((2.0 * (__pyx_v_k + __pyx_v_alpha)) / (__pyx_v_k + (2.0 * __pyx_v_alpha))) * (__pyx_v_x - 1.0)) * __pyx_v_p) + ((__pyx_v_k / (__pyx_v_k + (2.0 * __pyx_v_alpha))) * __pyx_v_d));

      
      __pyx_v_p = (__pyx_v_d + __pyx_v_p);
    }

    
    __pyx_r = (__pyx_f_5scipy_7special_15orthogonal_eval_binom(((__pyx_v_n + (2.0 * __pyx_v_alpha)) - 1.0), __pyx_v_n) * __pyx_v_p);
    goto __pyx_L0;
  }
  __pyx_L3:;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  CYTHON_UNUSED double __pyx_v_d;
  double __pyx_v_g;
  double __pyx_r;

  
  __pyx_v_d = 1.0;

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = __pyx_v_n;

  
  __pyx_v_c = 0.5;

  
  __pyx_v_g = (0.5 * (1.0 - __pyx_v_x));

  
  __pyx_r = __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt(double __pyx_v_n, __pyx_t_double_complex __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  CYTHON_UNUSED double __pyx_v_d;
  __pyx_t_double_complex __pyx_v_g;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_d = 1.0;

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = __pyx_v_n;

  
  __pyx_v_c = 0.5;

  
  __pyx_v_g = __Pyx_c_prod(__pyx_t_double_complex_from_parts(0.5, 0), __Pyx_c_diff(__pyx_t_double_complex_from_parts(1, 0), __pyx_v_x));

  
  __pyx_r = __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g);
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt_l(long __pyx_v_k, double __pyx_v_x) {
  CYTHON_UNUSED long __pyx_v_m;
  double __pyx_v_b2;
  double __pyx_v_b1;
  double __pyx_v_b0;
  double __pyx_r;
  long __pyx_t_1;
  long __pyx_t_2;

  
  __pyx_v_b2 = 0.0;

  
  __pyx_v_b1 = -1.0;

  
  __pyx_v_b0 = 0.0;

  
  __pyx_v_x = (2.0 * __pyx_v_x);

  
  __pyx_t_1 = (__pyx_v_k + 1);
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_m = __pyx_t_2;

    
    __pyx_v_b2 = __pyx_v_b1;

    
    __pyx_v_b1 = __pyx_v_b0;

    
    __pyx_v_b0 = ((__pyx_v_x * __pyx_v_b1) - __pyx_v_b2);
  }

  
  __pyx_r = ((__pyx_v_b0 - __pyx_v_b2) / 2.0);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  double __pyx_v_d;
  double __pyx_v_g;
  double __pyx_r;

  
  __pyx_v_d = (__pyx_v_n + 1.0);

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (__pyx_v_n + 2.0);

  
  __pyx_v_c = 1.5;

  
  __pyx_v_g = (0.5 * (1.0 - __pyx_v_x));

  
  __pyx_r = (__pyx_v_d * __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu(double __pyx_v_n, __pyx_t_double_complex __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  double __pyx_v_d;
  __pyx_t_double_complex __pyx_v_g;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_d = (__pyx_v_n + 1.0);

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (__pyx_v_n + 2.0);

  
  __pyx_v_c = 1.5;

  
  __pyx_v_g = __Pyx_c_prod(__pyx_t_double_complex_from_parts(0.5, 0), __Pyx_c_diff(__pyx_t_double_complex_from_parts(1, 0), __pyx_v_x));

  
  __pyx_r = __Pyx_c_prod(__pyx_t_double_complex_from_parts(__pyx_v_d, 0), __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu_l(long __pyx_v_k, double __pyx_v_x) {
  CYTHON_UNUSED long __pyx_v_m;
  double __pyx_v_b2;
  double __pyx_v_b1;
  double __pyx_v_b0;
  double __pyx_r;
  long __pyx_t_1;
  long __pyx_t_2;

  
  __pyx_v_b2 = 0.0;

  
  __pyx_v_b1 = -1.0;

  
  __pyx_v_b0 = 0.0;

  
  __pyx_v_x = (2.0 * __pyx_v_x);

  
  __pyx_t_1 = (__pyx_v_k + 1);
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_m = __pyx_t_2;

    
    __pyx_v_b2 = __pyx_v_b1;

    
    __pyx_v_b1 = __pyx_v_b0;

    
    __pyx_v_b0 = ((__pyx_v_x * __pyx_v_b1) - __pyx_v_b2);
  }

  
  __pyx_r = __pyx_v_b0;
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu(__pyx_v_n, (0.5 * __pyx_v_x));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys(double __pyx_v_n, __pyx_t_double_complex __pyx_v_x) {
  __pyx_t_double_complex __pyx_r;

  
  __pyx_r = __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu(__pyx_v_n, __Pyx_c_prod(__pyx_t_double_complex_from_parts(0.5, 0), __pyx_v_x));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys_l(long __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu_l(__pyx_v_n, (0.5 * __pyx_v_x));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = (2.0 * __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt(__pyx_v_n, (0.5 * __pyx_v_x)));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc(double __pyx_v_n, __pyx_t_double_complex __pyx_v_x) {
  __pyx_t_double_complex __pyx_r;

  
  __pyx_r = __Pyx_c_prod(__pyx_t_double_complex_from_parts(2, 0), __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt(__pyx_v_n, __Pyx_c_prod(__pyx_t_double_complex_from_parts(0.5, 0), __pyx_v_x)));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc_l(long __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = (2.0 * __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt_l(__pyx_v_n, (0.5 * __pyx_v_x)));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt(__pyx_v_n, ((2.0 * __pyx_v_x) - 1.0));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt(double __pyx_v_n, __pyx_t_double_complex __pyx_v_x) {
  __pyx_t_double_complex __pyx_r;

  
  __pyx_r = __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt(__pyx_v_n, __Pyx_c_diff(__Pyx_c_prod(__pyx_t_double_complex_from_parts(2, 0), __pyx_v_x), __pyx_t_double_complex_from_parts(1, 0)));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt_l(long __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt_l(__pyx_v_n, ((2.0 * __pyx_v_x) - 1.0));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu(__pyx_v_n, ((2.0 * __pyx_v_x) - 1.0));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu(double __pyx_v_n, __pyx_t_double_complex __pyx_v_x) {
  __pyx_t_double_complex __pyx_r;

  
  __pyx_r = __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu(__pyx_v_n, __Pyx_c_diff(__Pyx_c_prod(__pyx_t_double_complex_from_parts(2, 0), __pyx_v_x), __pyx_t_double_complex_from_parts(1, 0)));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu_l(long __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu_l(__pyx_v_n, ((2.0 * __pyx_v_x) - 1.0));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  double __pyx_v_d;
  double __pyx_v_g;
  double __pyx_r;

  
  __pyx_v_d = 1.0;

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (__pyx_v_n + 1.0);

  
  __pyx_v_c = 1.0;

  
  __pyx_v_g = (0.5 * (1.0 - __pyx_v_x));

  
  __pyx_r = (__pyx_v_d * __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre(double __pyx_v_n, __pyx_t_double_complex __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_c;
  double __pyx_v_d;
  __pyx_t_double_complex __pyx_v_g;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_d = 1.0;

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (__pyx_v_n + 1.0);

  
  __pyx_v_c = 1.0;

  
  __pyx_v_g = __Pyx_c_prod(__pyx_t_double_complex_from_parts(0.5, 0), __Pyx_c_diff(__pyx_t_double_complex_from_parts(1, 0), __pyx_v_x));

  
  __pyx_r = __Pyx_c_prod(__pyx_t_double_complex_from_parts(__pyx_v_d, 0), __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp2f1(__pyx_v_a, __pyx_v_b, __pyx_v_c, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre_l(long __pyx_v_n, double __pyx_v_x) {
  long __pyx_v_kk;
  double __pyx_v_p;
  double __pyx_v_d;
  double __pyx_v_k;
  double __pyx_r;
  int __pyx_t_1;
  long __pyx_t_2;
  long __pyx_t_3;

  
  __pyx_t_1 = (__pyx_v_n < 0);
  if (__pyx_t_1) {

    
    __pyx_r = 0.0;
    goto __pyx_L0;
    goto __pyx_L3;
  }

  
  __pyx_t_1 = (__pyx_v_n == 0);
  if (__pyx_t_1) {

    
    __pyx_r = 1.0;
    goto __pyx_L0;
    goto __pyx_L3;
  }

  
  __pyx_t_1 = (__pyx_v_n == 1);
  if (__pyx_t_1) {

    
    __pyx_r = __pyx_v_x;
    goto __pyx_L0;
    goto __pyx_L3;
  }
   {

    
    __pyx_v_d = (__pyx_v_x - 1.0);

    
    __pyx_v_p = __pyx_v_x;

    
    __pyx_t_2 = (__pyx_v_n - 1);
    for (__pyx_t_3 = 0; __pyx_t_3 < __pyx_t_2; __pyx_t_3+=1) {
      __pyx_v_kk = __pyx_t_3;

      
      __pyx_v_k = (__pyx_v_kk + 1.0);

      
      __pyx_v_d = ((((((2.0 * __pyx_v_k) + 1.0) / (__pyx_v_k + 1.0)) * (__pyx_v_x - 1.0)) * __pyx_v_p) + ((__pyx_v_k / (__pyx_v_k + 1.0)) * __pyx_v_d));

      
      __pyx_v_p = (__pyx_v_d + __pyx_v_p);
    }

    
    __pyx_r = __pyx_v_p;
    goto __pyx_L0;
  }
  __pyx_L3:;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre(__pyx_v_n, ((2.0 * __pyx_v_x) - 1.0));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre(double __pyx_v_n, __pyx_t_double_complex __pyx_v_x) {
  __pyx_t_double_complex __pyx_r;

  
  __pyx_r = __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre(__pyx_v_n, __Pyx_c_diff(__Pyx_c_prod(__pyx_t_double_complex_from_parts(2, 0), __pyx_v_x), __pyx_t_double_complex_from_parts(1, 0)));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre_l(long __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre_l(__pyx_v_n, ((2.0 * __pyx_v_x) - 1.0));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre(double __pyx_v_n, double __pyx_v_alpha, double __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_d;
  double __pyx_v_g;
  double __pyx_r;

  
  __pyx_v_d = __pyx_f_5scipy_7special_15orthogonal_eval_binom((__pyx_v_n + __pyx_v_alpha), __pyx_v_n);

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (__pyx_v_alpha + 1.0);

  
  __pyx_v_g = __pyx_v_x;

  
  __pyx_r = (__pyx_v_d * __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_hyp1f1(__pyx_v_a, __pyx_v_b, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre(double __pyx_v_n, double __pyx_v_alpha, __pyx_t_double_complex __pyx_v_x) {
  double __pyx_v_a;
  double __pyx_v_b;
  double __pyx_v_d;
  __pyx_t_double_complex __pyx_v_g;
  __pyx_t_double_complex __pyx_r;

  
  __pyx_v_d = __pyx_f_5scipy_7special_15orthogonal_eval_binom((__pyx_v_n + __pyx_v_alpha), __pyx_v_n);

  
  __pyx_v_a = (-__pyx_v_n);

  
  __pyx_v_b = (__pyx_v_alpha + 1.0);

  
  __pyx_v_g = __pyx_v_x;

  
  __pyx_r = __Pyx_c_prod(__pyx_t_double_complex_from_parts(__pyx_v_d, 0), __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_hyp1f1(__pyx_v_a, __pyx_v_b, __pyx_v_g));
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre_l(long __pyx_v_n, double __pyx_v_alpha, double __pyx_v_x) {
  long __pyx_v_kk;
  double __pyx_v_p;
  double __pyx_v_d;
  double __pyx_v_k;
  double __pyx_r;
  int __pyx_t_1;
  long __pyx_t_2;
  long __pyx_t_3;

  
  __pyx_t_1 = (__pyx_v_n < 0);
  if (__pyx_t_1) {

    
    __pyx_r = 0.0;
    goto __pyx_L0;
    goto __pyx_L3;
  }

  
  __pyx_t_1 = (__pyx_v_n == 0);
  if (__pyx_t_1) {

    
    __pyx_r = 1.0;
    goto __pyx_L0;
    goto __pyx_L3;
  }

  
  __pyx_t_1 = (__pyx_v_n == 1);
  if (__pyx_t_1) {

    
    __pyx_r = (((-__pyx_v_x) + __pyx_v_alpha) + 1.0);
    goto __pyx_L0;
    goto __pyx_L3;
  }
   {

    
    __pyx_v_d = ((-__pyx_v_x) / (__pyx_v_alpha + 1.0));

    
    __pyx_v_p = (__pyx_v_d + 1.0);

    
    __pyx_t_2 = (__pyx_v_n - 1);
    for (__pyx_t_3 = 0; __pyx_t_3 < __pyx_t_2; __pyx_t_3+=1) {
      __pyx_v_kk = __pyx_t_3;

      
      __pyx_v_k = (__pyx_v_kk + 1.0);

      
      __pyx_v_d = ((((-__pyx_v_x) / ((__pyx_v_k + __pyx_v_alpha) + 1.0)) * __pyx_v_p) + ((__pyx_v_k / ((__pyx_v_k + __pyx_v_alpha) + 1.0)) * __pyx_v_d));

      
      __pyx_v_p = (__pyx_v_d + __pyx_v_p);
    }

    
    __pyx_r = (__pyx_f_5scipy_7special_15orthogonal_eval_binom((__pyx_v_n + __pyx_v_alpha), __pyx_v_n) * __pyx_v_p);
    goto __pyx_L0;
  }
  __pyx_L3:;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre(double __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre(__pyx_v_n, 0., __pyx_v_x);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre(double __pyx_v_n, __pyx_t_double_complex __pyx_v_x) {
  __pyx_t_double_complex __pyx_r;

  
  __pyx_r = __pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre(__pyx_v_n, 0., __pyx_v_x);
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre_l(long __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = __pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre_l(__pyx_v_n, 0., __pyx_v_x);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_hermite(long __pyx_v_n, double __pyx_v_x) {
  long __pyx_v_m;
  double __pyx_r;
  int __pyx_t_1;

  
  __pyx_t_1 = (__Pyx_mod_long(__pyx_v_n, 2) == 0);
  if (__pyx_t_1) {

    
    __pyx_v_m = __Pyx_div_long(__pyx_v_n, 2);

    
    __pyx_r = (((__Pyx_pow_long(-1, __pyx_v_m) * __Pyx_pow_long(2, (2 * __pyx_v_m))) * Gamma((1 + __pyx_v_m))) * __pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre_l(__pyx_v_m, -0.5, pow(__pyx_v_x, 2.0)));
    goto __pyx_L0;
    goto __pyx_L3;
  }
   {

    
    __pyx_v_m = __Pyx_div_long((__pyx_v_n - 1), 2);

    
    __pyx_r = ((((__Pyx_pow_long(-1, __pyx_v_m) * __Pyx_pow_long(2, ((2 * __pyx_v_m) + 1))) * Gamma((1 + __pyx_v_m))) * __pyx_v_x) * __pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre_l(__pyx_v_m, 0.5, pow(__pyx_v_x, 2.0)));
    goto __pyx_L0;
  }
  __pyx_L3:;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE double __pyx_f_5scipy_7special_15orthogonal_eval_eval_hermitenorm(long __pyx_v_n, double __pyx_v_x) {
  double __pyx_r;

  
  __pyx_r = (__pyx_f_5scipy_7special_15orthogonal_eval_eval_hermite(__pyx_v_n, (__pyx_v_x / sqrt(2.0))) * pow(2.0, ((-__pyx_v_n) / 2.0)));
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static CYTHON_INLINE __pyx_t_double_complex __pyx_f_5scipy_7special_8lambertw_lambertw_scalar(__pyx_t_double_complex __pyx_v_z, long __pyx_v_k, double __pyx_v_tol) {
  __pyx_t_double_complex __pyx_v_w;
  double __pyx_v_u;
  double __pyx_v_absz;
  __pyx_t_double_complex __pyx_v_ew;
  __pyx_t_double_complex __pyx_v_wew;
  __pyx_t_double_complex __pyx_v_wewz;
  __pyx_t_double_complex __pyx_v_wn;
  CYTHON_UNUSED int __pyx_v_i;
  __pyx_t_double_complex __pyx_r;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  int __pyx_t_4;
  int __pyx_t_5;

  
  __pyx_t_1 = __pyx_f_5scipy_7special_13_complexstuff_zisnan(__pyx_v_z);
  if (__pyx_t_1) {

    
    __pyx_r = __pyx_v_z;
    goto __pyx_L0;
    goto __pyx_L3;
  }
  __pyx_L3:;

  
  __pyx_v_u = exp(-1.0);

  
  __pyx_v_absz = __pyx_f_5scipy_7special_13_complexstuff_zabs(__pyx_v_z);

  
  __pyx_t_1 = (__pyx_v_absz <= __pyx_v_u);
  if (__pyx_t_1) {

    
    __pyx_t_1 = (__Pyx_c_eq(__pyx_v_z, __pyx_t_double_complex_from_parts(0, 0)));
    if (__pyx_t_1) {

      
      __pyx_t_1 = (__pyx_v_k == 0);
      if (__pyx_t_1) {

        
        __pyx_r = __pyx_v_z;
        goto __pyx_L0;
        goto __pyx_L6;
      }
      __pyx_L6:;

      
      sf_error(__pyx_k__lambertw, SF_ERROR_SINGULAR, NULL);

      
      __pyx_r = __pyx_t_double_complex_from_parts((-NPY_INFINITY), 0);
      goto __pyx_L0;
      goto __pyx_L5;
    }
    __pyx_L5:;

    
    __pyx_t_1 = (__pyx_v_k == 0);
    if (__pyx_t_1) {

      
      __pyx_v_w = __pyx_v_z;
      goto __pyx_L7;
    }

    
    __pyx_t_1 = (__pyx_v_k == -1);
    if (__pyx_t_1) {
      __pyx_t_2 = (__Pyx_CIMAG(__pyx_v_z) == 0.0);
      if (__pyx_t_2) {
        __pyx_t_3 = (__Pyx_CREAL(__pyx_v_z) < 0.0);
        __pyx_t_4 = __pyx_t_3;
      } else {
        __pyx_t_4 = __pyx_t_2;
      }
      __pyx_t_2 = __pyx_t_4;
    } else {
      __pyx_t_2 = __pyx_t_1;
    }
    if (__pyx_t_2) {

      
      __pyx_v_w = __pyx_t_double_complex_from_parts(log((-__Pyx_CREAL(__pyx_v_z))), 0);
      goto __pyx_L7;
    }
     {

      
      __pyx_v_w = __pyx_f_5scipy_7special_13_complexstuff_zlog(__pyx_v_z);

      
      if (__pyx_v_k) {
        __pyx_v_w = __Pyx_c_sum(__pyx_v_w, __Pyx_c_prod(__pyx_t_double_complex_from_parts(((__pyx_v_k * 2) * NPY_PI), 0), __pyx_t_double_complex_from_parts(0, 1.0)));
        goto __pyx_L8;
      }
      __pyx_L8:;
    }
    __pyx_L7:;
    goto __pyx_L4;
  }

  
  __pyx_t_2 = (__pyx_v_k == 0);
  if (__pyx_t_2) {
    if ((__Pyx_CIMAG(__pyx_v_z) != 0)) {
      __pyx_t_1 = (__pyx_f_5scipy_7special_13_complexstuff_zabs(__pyx_v_z) <= 0.7);
      __pyx_t_4 = __pyx_t_1;
    } else {
      __pyx_t_4 = (__Pyx_CIMAG(__pyx_v_z) != 0);
    }
    __pyx_t_1 = __pyx_t_4;
  } else {
    __pyx_t_1 = __pyx_t_2;
  }
  if (__pyx_t_1) {

    
    __pyx_t_1 = (__pyx_f_5scipy_7special_13_complexstuff_zabs(__Pyx_c_sum(__pyx_v_z, __pyx_t_double_complex_from_parts(0.5, 0))) < 0.1);
    if (__pyx_t_1) {

      
      __pyx_t_1 = (__Pyx_CIMAG(__pyx_v_z) > 0.0);
      if (__pyx_t_1) {

        
        __pyx_v_w = __Pyx_c_sum(__pyx_t_double_complex_from_parts(0.7, 0), __pyx_t_double_complex_from_parts(0, 0.7));
        goto __pyx_L10;
      }
       {

        
        __pyx_v_w = __Pyx_c_diff(__pyx_t_double_complex_from_parts(0.7, 0), __pyx_t_double_complex_from_parts(0, 0.7));
      }
      __pyx_L10:;
      goto __pyx_L9;
    }
     {

      
      __pyx_v_w = __pyx_v_z;
    }
    __pyx_L9:;
    goto __pyx_L4;
  }
   {

    
    __pyx_t_1 = (__Pyx_CREAL(__pyx_v_z) == NPY_INFINITY);
    if (__pyx_t_1) {

      
      __pyx_t_1 = (__pyx_v_k == 0);
      if (__pyx_t_1) {

        
        __pyx_r = __pyx_v_z;
        goto __pyx_L0;
        goto __pyx_L12;
      }
       {

        
        __pyx_r = __Pyx_c_sum(__pyx_v_z, __Pyx_c_prod(__pyx_t_double_complex_from_parts(((2 * __pyx_v_k) * NPY_PI), 0), __pyx_t_double_complex_from_parts(0, 1.0)));
        goto __pyx_L0;
      }
      __pyx_L12:;
      goto __pyx_L11;
    }
    __pyx_L11:;

    
    __pyx_t_1 = (__Pyx_CREAL(__pyx_v_z) == (-NPY_INFINITY));
    if (__pyx_t_1) {

      
      __pyx_r = __Pyx_c_sum(__Pyx_c_neg(__pyx_v_z), __Pyx_c_prod(__pyx_t_double_complex_from_parts((((2 * __pyx_v_k) + 1) * NPY_PI), 0), __pyx_t_double_complex_from_parts(0, 1.0)));
      goto __pyx_L0;
      goto __pyx_L13;
    }
    __pyx_L13:;

    
    __pyx_v_w = __pyx_f_5scipy_7special_13_complexstuff_zlog(__pyx_v_z);

    
    if (__pyx_v_k) {
      __pyx_v_w = __Pyx_c_sum(__pyx_v_w, __Pyx_c_prod(__pyx_t_double_complex_from_parts(((__pyx_v_k * 2) * NPY_PI), 0), __pyx_t_double_complex_from_parts(0, 1.0)));
      goto __pyx_L14;
    }
    __pyx_L14:;
  }
  __pyx_L4:;

  
  for (__pyx_t_5 = 0; __pyx_t_5 < 100; __pyx_t_5+=1) {
    __pyx_v_i = __pyx_t_5;

    
    __pyx_v_ew = __pyx_f_5scipy_7special_13_complexstuff_zexp(__pyx_v_w);

    
    __pyx_v_wew = __Pyx_c_prod(__pyx_v_w, __pyx_v_ew);

    
    __pyx_v_wewz = __Pyx_c_diff(__pyx_v_wew, __pyx_v_z);

    
    __pyx_v_wn = __Pyx_c_diff(__pyx_v_w, __Pyx_c_quot(__pyx_v_wewz, __Pyx_c_diff(__Pyx_c_sum(__pyx_v_wew, __pyx_v_ew), __Pyx_c_quot(__Pyx_c_prod(__Pyx_c_sum(__pyx_v_w, __pyx_t_double_complex_from_parts(2, 0)), __pyx_v_wewz), __Pyx_c_sum(__Pyx_c_prod(__pyx_t_double_complex_from_parts(2, 0), __pyx_v_w), __pyx_t_double_complex_from_parts(2, 0))))));

    
    __pyx_t_1 = (__pyx_f_5scipy_7special_13_complexstuff_zabs(__Pyx_c_diff(__pyx_v_wn, __pyx_v_w)) < (__pyx_v_tol * __pyx_f_5scipy_7special_13_complexstuff_zabs(__pyx_v_wn)));
    if (__pyx_t_1) {

      
      __pyx_r = __pyx_v_wn;
      goto __pyx_L0;
      goto __pyx_L17;
    }
     {

      
      __pyx_v_w = __pyx_v_wn;
    }
    __pyx_L17:;
  }

  
  sf_error(__pyx_k__lambertw, SF_ERROR_SLOW, __pyx_k_5, ((double)__Pyx_CREAL(__pyx_v_z)), ((double)__Pyx_CIMAG(__pyx_v_z)));

  
  __pyx_r = __pyx_t_double_complex_from_parts(NPY_NAN, 0);
  goto __pyx_L0;

  __pyx_r = __pyx_t_double_complex_from_parts(0, 0);
  __pyx_L0:;
  return __pyx_r;
}


static int __pyx_pw_5numpy_7ndarray_1__getbuffer__(PyObject *__pyx_v_self, Py_buffer *__pyx_v_info, int __pyx_v_flags); 
static int __pyx_pw_5numpy_7ndarray_1__getbuffer__(PyObject *__pyx_v_self, Py_buffer *__pyx_v_info, int __pyx_v_flags) {
  int __pyx_r;
  __Pyx_RefNannyDeclarations
  __Pyx_RefNannySetupContext("__getbuffer__ (wrapper)", 0);
  __pyx_r = __pyx_pf_5numpy_7ndarray___getbuffer__(((PyArrayObject *)__pyx_v_self), ((Py_buffer *)__pyx_v_info), ((int)__pyx_v_flags));
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static int __pyx_pf_5numpy_7ndarray___getbuffer__(PyArrayObject *__pyx_v_self, Py_buffer *__pyx_v_info, int __pyx_v_flags) {
  int __pyx_v_copy_shape;
  int __pyx_v_i;
  int __pyx_v_ndim;
  int __pyx_v_endian_detector;
  int __pyx_v_little_endian;
  int __pyx_v_t;
  char *__pyx_v_f;
  PyArray_Descr *__pyx_v_descr = 0;
  int __pyx_v_offset;
  int __pyx_v_hasfields;
  int __pyx_r;
  __Pyx_RefNannyDeclarations
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  PyObject *__pyx_t_4 = NULL;
  int __pyx_t_5;
  int __pyx_t_6;
  int __pyx_t_7;
  PyObject *__pyx_t_8 = NULL;
  char *__pyx_t_9;
  int __pyx_lineno = 0;
  const char *__pyx_filename = NULL;
  int __pyx_clineno = 0;
  __Pyx_RefNannySetupContext("__getbuffer__", 0);
  if (__pyx_v_info != NULL) {
    __pyx_v_info->obj = Py_None; __Pyx_INCREF(Py_None);
    __Pyx_GIVEREF(__pyx_v_info->obj);
  }

  
  __pyx_t_1 = (__pyx_v_info == NULL);
  if (__pyx_t_1) {
    __pyx_r = 0;
    goto __pyx_L0;
    goto __pyx_L3;
  }
  __pyx_L3:;

  
  __pyx_v_endian_detector = 1;

  
  __pyx_v_little_endian = ((((char *)(&__pyx_v_endian_detector))[0]) != 0);

  
  __pyx_v_ndim = PyArray_NDIM(__pyx_v_self);

  
  __pyx_t_1 = ((sizeof(npy_intp)) != (sizeof(Py_ssize_t)));
  if (__pyx_t_1) {

    
    __pyx_v_copy_shape = 1;
    goto __pyx_L4;
  }
   {

    
    __pyx_v_copy_shape = 0;
  }
  __pyx_L4:;

  
  __pyx_t_1 = ((__pyx_v_flags & PyBUF_C_CONTIGUOUS) == PyBUF_C_CONTIGUOUS);
  if (__pyx_t_1) {

    
    __pyx_t_2 = (!PyArray_CHKFLAGS(__pyx_v_self, NPY_C_CONTIGUOUS));
    __pyx_t_3 = __pyx_t_2;
  } else {
    __pyx_t_3 = __pyx_t_1;
  }
  if (__pyx_t_3) {

    
    __pyx_t_4 = PyObject_Call(__pyx_builtin_ValueError, ((PyObject *)__pyx_k_tuple_7), NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 215; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __Pyx_Raise(__pyx_t_4, 0, 0, 0);
    __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
    {__pyx_filename = __pyx_f[2]; __pyx_lineno = 215; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L5;
  }
  __pyx_L5:;

  
  __pyx_t_3 = ((__pyx_v_flags & PyBUF_F_CONTIGUOUS) == PyBUF_F_CONTIGUOUS);
  if (__pyx_t_3) {

    
    __pyx_t_1 = (!PyArray_CHKFLAGS(__pyx_v_self, NPY_F_CONTIGUOUS));
    __pyx_t_2 = __pyx_t_1;
  } else {
    __pyx_t_2 = __pyx_t_3;
  }
  if (__pyx_t_2) {

    
    __pyx_t_4 = PyObject_Call(__pyx_builtin_ValueError, ((PyObject *)__pyx_k_tuple_9), NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 219; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __Pyx_Raise(__pyx_t_4, 0, 0, 0);
    __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
    {__pyx_filename = __pyx_f[2]; __pyx_lineno = 219; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L6;
  }
  __pyx_L6:;

  
  __pyx_v_info->buf = PyArray_DATA(__pyx_v_self);

  
  __pyx_v_info->ndim = __pyx_v_ndim;

  
  if (__pyx_v_copy_shape) {

    
    __pyx_v_info->strides = ((Py_ssize_t *)malloc((((sizeof(Py_ssize_t)) * ((size_t)__pyx_v_ndim)) * 2)));

    
    __pyx_v_info->shape = (__pyx_v_info->strides + __pyx_v_ndim);

    
    __pyx_t_5 = __pyx_v_ndim;
    for (__pyx_t_6 = 0; __pyx_t_6 < __pyx_t_5; __pyx_t_6+=1) {
      __pyx_v_i = __pyx_t_6;

      
      (__pyx_v_info->strides[__pyx_v_i]) = (PyArray_STRIDES(__pyx_v_self)[__pyx_v_i]);

      
      (__pyx_v_info->shape[__pyx_v_i]) = (PyArray_DIMS(__pyx_v_self)[__pyx_v_i]);
    }
    goto __pyx_L7;
  }
   {

    
    __pyx_v_info->strides = ((Py_ssize_t *)PyArray_STRIDES(__pyx_v_self));

    
    __pyx_v_info->shape = ((Py_ssize_t *)PyArray_DIMS(__pyx_v_self));
  }
  __pyx_L7:;

  
  __pyx_v_info->suboffsets = NULL;

  
  __pyx_v_info->itemsize = PyArray_ITEMSIZE(__pyx_v_self);

  
  __pyx_v_info->readonly = (!PyArray_ISWRITEABLE(__pyx_v_self));

  
  __pyx_v_f = NULL;

  
  __pyx_t_4 = ((PyObject *)__pyx_v_self->descr);
  __Pyx_INCREF(__pyx_t_4);
  __pyx_v_descr = ((PyArray_Descr *)__pyx_t_4);
  __pyx_t_4 = 0;

  
  __pyx_v_hasfields = PyDataType_HASFIELDS(__pyx_v_descr);

  
  __pyx_t_2 = (!__pyx_v_hasfields);
  if (__pyx_t_2) {
    __pyx_t_3 = (!__pyx_v_copy_shape);
    __pyx_t_1 = __pyx_t_3;
  } else {
    __pyx_t_1 = __pyx_t_2;
  }
  if (__pyx_t_1) {

    
    __Pyx_INCREF(Py_None);
    __Pyx_GIVEREF(Py_None);
    __Pyx_GOTREF(__pyx_v_info->obj);
    __Pyx_DECREF(__pyx_v_info->obj);
    __pyx_v_info->obj = Py_None;
    goto __pyx_L10;
  }
   {

    
    __Pyx_INCREF(((PyObject *)__pyx_v_self));
    __Pyx_GIVEREF(((PyObject *)__pyx_v_self));
    __Pyx_GOTREF(__pyx_v_info->obj);
    __Pyx_DECREF(__pyx_v_info->obj);
    __pyx_v_info->obj = ((PyObject *)__pyx_v_self);
  }
  __pyx_L10:;

  
  __pyx_t_1 = (!__pyx_v_hasfields);
  if (__pyx_t_1) {

    
    __pyx_t_5 = __pyx_v_descr->type_num;
    __pyx_v_t = __pyx_t_5;

    
    __pyx_t_1 = (__pyx_v_descr->byteorder == '>');
    if (__pyx_t_1) {
      __pyx_t_2 = __pyx_v_little_endian;
    } else {
      __pyx_t_2 = __pyx_t_1;
    }
    if (!__pyx_t_2) {

      
      __pyx_t_1 = (__pyx_v_descr->byteorder == '<');
      if (__pyx_t_1) {
        __pyx_t_3 = (!__pyx_v_little_endian);
        __pyx_t_7 = __pyx_t_3;
      } else {
        __pyx_t_7 = __pyx_t_1;
      }
      __pyx_t_1 = __pyx_t_7;
    } else {
      __pyx_t_1 = __pyx_t_2;
    }
    if (__pyx_t_1) {

      
      __pyx_t_4 = PyObject_Call(__pyx_builtin_ValueError, ((PyObject *)__pyx_k_tuple_11), NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 257; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_4);
      __Pyx_Raise(__pyx_t_4, 0, 0, 0);
      __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
      {__pyx_filename = __pyx_f[2]; __pyx_lineno = 257; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      goto __pyx_L12;
    }
    __pyx_L12:;

    
    __pyx_t_1 = (__pyx_v_t == NPY_BYTE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__b;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_UBYTE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__B;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_SHORT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__h;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_USHORT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__H;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_INT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__i;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_UINT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__I;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_LONG);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__l;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_ULONG);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__L;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_LONGLONG);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__q;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_ULONGLONG);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__Q;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_FLOAT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__f;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_DOUBLE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__d;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_LONGDOUBLE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__g;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_CFLOAT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__Zf;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_CDOUBLE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__Zd;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_CLONGDOUBLE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__Zg;
      goto __pyx_L13;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_OBJECT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__O;
      goto __pyx_L13;
    }
     {

      
      __pyx_t_4 = PyInt_FromLong(__pyx_v_t); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 276; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_4);
      __pyx_t_8 = PyNumber_Remainder(((PyObject *)__pyx_kp_u_12), __pyx_t_4); if (unlikely(!__pyx_t_8)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 276; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(((PyObject *)__pyx_t_8));
      __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
      __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 276; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_4);
      PyTuple_SET_ITEM(__pyx_t_4, 0, ((PyObject *)__pyx_t_8));
      __Pyx_GIVEREF(((PyObject *)__pyx_t_8));
      __pyx_t_8 = 0;
      __pyx_t_8 = PyObject_Call(__pyx_builtin_ValueError, ((PyObject *)__pyx_t_4), NULL); if (unlikely(!__pyx_t_8)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 276; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_8);
      __Pyx_DECREF(((PyObject *)__pyx_t_4)); __pyx_t_4 = 0;
      __Pyx_Raise(__pyx_t_8, 0, 0, 0);
      __Pyx_DECREF(__pyx_t_8); __pyx_t_8 = 0;
      {__pyx_filename = __pyx_f[2]; __pyx_lineno = 276; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    }
    __pyx_L13:;

    
    __pyx_v_info->format = __pyx_v_f;

    
    __pyx_r = 0;
    goto __pyx_L0;
    goto __pyx_L11;
  }
   {

    
    __pyx_v_info->format = ((char *)malloc(255));

    
    (__pyx_v_info->format[0]) = '^';

    
    __pyx_v_offset = 0;

    
    __pyx_t_9 = __pyx_f_5numpy__util_dtypestring(__pyx_v_descr, (__pyx_v_info->format + 1), (__pyx_v_info->format + 255), (&__pyx_v_offset)); if (unlikely(__pyx_t_9 == NULL)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 283; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __pyx_v_f = __pyx_t_9;

    
    (__pyx_v_f[0]) = '\x00';
  }
  __pyx_L11:;

  __pyx_r = 0;
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_8);
  __Pyx_AddTraceback("numpy.ndarray.__getbuffer__", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = -1;
  if (__pyx_v_info != NULL && __pyx_v_info->obj != NULL) {
    __Pyx_GOTREF(__pyx_v_info->obj);
    __Pyx_DECREF(__pyx_v_info->obj); __pyx_v_info->obj = NULL;
  }
  goto __pyx_L2;
  __pyx_L0:;
  if (__pyx_v_info != NULL && __pyx_v_info->obj == Py_None) {
    __Pyx_GOTREF(Py_None);
    __Pyx_DECREF(Py_None); __pyx_v_info->obj = NULL;
  }
  __pyx_L2:;
  __Pyx_XDECREF((PyObject *)__pyx_v_descr);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}


static void __pyx_pw_5numpy_7ndarray_3__releasebuffer__(PyObject *__pyx_v_self, Py_buffer *__pyx_v_info); 
static void __pyx_pw_5numpy_7ndarray_3__releasebuffer__(PyObject *__pyx_v_self, Py_buffer *__pyx_v_info) {
  __Pyx_RefNannyDeclarations
  __Pyx_RefNannySetupContext("__releasebuffer__ (wrapper)", 0);
  __pyx_pf_5numpy_7ndarray_2__releasebuffer__(((PyArrayObject *)__pyx_v_self), ((Py_buffer *)__pyx_v_info));
  __Pyx_RefNannyFinishContext();
}



static void __pyx_pf_5numpy_7ndarray_2__releasebuffer__(PyArrayObject *__pyx_v_self, Py_buffer *__pyx_v_info) {
  __Pyx_RefNannyDeclarations
  int __pyx_t_1;
  __Pyx_RefNannySetupContext("__releasebuffer__", 0);

  
  __pyx_t_1 = PyArray_HASFIELDS(__pyx_v_self);
  if (__pyx_t_1) {

    
    free(__pyx_v_info->format);
    goto __pyx_L3;
  }
  __pyx_L3:;

  
  __pyx_t_1 = ((sizeof(npy_intp)) != (sizeof(Py_ssize_t)));
  if (__pyx_t_1) {

    
    free(__pyx_v_info->strides);
    goto __pyx_L4;
  }
  __pyx_L4:;

  __Pyx_RefNannyFinishContext();
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew1(PyObject *__pyx_v_a) {
  PyObject *__pyx_r = NULL;
  __Pyx_RefNannyDeclarations
  PyObject *__pyx_t_1 = NULL;
  int __pyx_lineno = 0;
  const char *__pyx_filename = NULL;
  int __pyx_clineno = 0;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew1", 0);

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(1, ((void *)__pyx_v_a)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 769; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew1", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew2(PyObject *__pyx_v_a, PyObject *__pyx_v_b) {
  PyObject *__pyx_r = NULL;
  __Pyx_RefNannyDeclarations
  PyObject *__pyx_t_1 = NULL;
  int __pyx_lineno = 0;
  const char *__pyx_filename = NULL;
  int __pyx_clineno = 0;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew2", 0);

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(2, ((void *)__pyx_v_a), ((void *)__pyx_v_b)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 772; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew2", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew3(PyObject *__pyx_v_a, PyObject *__pyx_v_b, PyObject *__pyx_v_c) {
  PyObject *__pyx_r = NULL;
  __Pyx_RefNannyDeclarations
  PyObject *__pyx_t_1 = NULL;
  int __pyx_lineno = 0;
  const char *__pyx_filename = NULL;
  int __pyx_clineno = 0;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew3", 0);

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(3, ((void *)__pyx_v_a), ((void *)__pyx_v_b), ((void *)__pyx_v_c)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 775; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew3", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew4(PyObject *__pyx_v_a, PyObject *__pyx_v_b, PyObject *__pyx_v_c, PyObject *__pyx_v_d) {
  PyObject *__pyx_r = NULL;
  __Pyx_RefNannyDeclarations
  PyObject *__pyx_t_1 = NULL;
  int __pyx_lineno = 0;
  const char *__pyx_filename = NULL;
  int __pyx_clineno = 0;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew4", 0);

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(4, ((void *)__pyx_v_a), ((void *)__pyx_v_b), ((void *)__pyx_v_c), ((void *)__pyx_v_d)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 778; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew4", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew5(PyObject *__pyx_v_a, PyObject *__pyx_v_b, PyObject *__pyx_v_c, PyObject *__pyx_v_d, PyObject *__pyx_v_e) {
  PyObject *__pyx_r = NULL;
  __Pyx_RefNannyDeclarations
  PyObject *__pyx_t_1 = NULL;
  int __pyx_lineno = 0;
  const char *__pyx_filename = NULL;
  int __pyx_clineno = 0;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew5", 0);

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(5, ((void *)__pyx_v_a), ((void *)__pyx_v_b), ((void *)__pyx_v_c), ((void *)__pyx_v_d), ((void *)__pyx_v_e)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 781; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew5", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE char *__pyx_f_5numpy__util_dtypestring(PyArray_Descr *__pyx_v_descr, char *__pyx_v_f, char *__pyx_v_end, int *__pyx_v_offset) {
  PyArray_Descr *__pyx_v_child = 0;
  int __pyx_v_endian_detector;
  int __pyx_v_little_endian;
  PyObject *__pyx_v_fields = 0;
  PyObject *__pyx_v_childname = NULL;
  PyObject *__pyx_v_new_offset = NULL;
  PyObject *__pyx_v_t = NULL;
  char *__pyx_r;
  __Pyx_RefNannyDeclarations
  PyObject *__pyx_t_1 = NULL;
  Py_ssize_t __pyx_t_2;
  PyObject *__pyx_t_3 = NULL;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  PyObject *(*__pyx_t_6)(PyObject *);
  int __pyx_t_7;
  int __pyx_t_8;
  int __pyx_t_9;
  int __pyx_t_10;
  long __pyx_t_11;
  char *__pyx_t_12;
  int __pyx_lineno = 0;
  const char *__pyx_filename = NULL;
  int __pyx_clineno = 0;
  __Pyx_RefNannySetupContext("_util_dtypestring", 0);

  
  __pyx_v_endian_detector = 1;

  
  __pyx_v_little_endian = ((((char *)(&__pyx_v_endian_detector))[0]) != 0);

  
  if (unlikely(((PyObject *)__pyx_v_descr->names) == Py_None)) {
    PyErr_SetString(PyExc_TypeError, "'NoneType' object is not iterable");
    {__pyx_filename = __pyx_f[2]; __pyx_lineno = 794; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_1 = ((PyObject *)__pyx_v_descr->names); __Pyx_INCREF(__pyx_t_1); __pyx_t_2 = 0;
  for (;;) {
    if (__pyx_t_2 >= PyTuple_GET_SIZE(__pyx_t_1)) break;
    #if CYTHON_COMPILING_IN_CPYTHON
    __pyx_t_3 = PyTuple_GET_ITEM(__pyx_t_1, __pyx_t_2); __Pyx_INCREF(__pyx_t_3); __pyx_t_2++; if (unlikely(0 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 794; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    #else
    __pyx_t_3 = PySequence_ITEM(__pyx_t_1, __pyx_t_2); __pyx_t_2++; if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 794; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    #endif
    __Pyx_XDECREF(__pyx_v_childname);
    __pyx_v_childname = __pyx_t_3;
    __pyx_t_3 = 0;

    
    __pyx_t_3 = PyObject_GetItem(__pyx_v_descr->fields, __pyx_v_childname); if (!__pyx_t_3) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 795; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    if (!(likely(PyTuple_CheckExact(__pyx_t_3))||((__pyx_t_3) == Py_None)||(PyErr_Format(PyExc_TypeError, "Expected tuple, got %.200s", Py_TYPE(__pyx_t_3)->tp_name), 0))) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 795; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_XDECREF(((PyObject *)__pyx_v_fields));
    __pyx_v_fields = ((PyObject*)__pyx_t_3);
    __pyx_t_3 = 0;

    
    if (likely(PyTuple_CheckExact(((PyObject *)__pyx_v_fields)))) {
      PyObject* sequence = ((PyObject *)__pyx_v_fields);
      #if CYTHON_COMPILING_IN_CPYTHON
      Py_ssize_t size = Py_SIZE(sequence);
      #else
      Py_ssize_t size = PySequence_Size(sequence);
      #endif
      if (unlikely(size != 2)) {
        if (size > 2) __Pyx_RaiseTooManyValuesError(2);
        else if (size >= 0) __Pyx_RaiseNeedMoreValuesError(size);
        {__pyx_filename = __pyx_f[2]; __pyx_lineno = 796; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      #if CYTHON_COMPILING_IN_CPYTHON
      __pyx_t_3 = PyTuple_GET_ITEM(sequence, 0); 
      __pyx_t_4 = PyTuple_GET_ITEM(sequence, 1); 
      __Pyx_INCREF(__pyx_t_3);
      __Pyx_INCREF(__pyx_t_4);
      #else
      __pyx_t_3 = PySequence_ITEM(sequence, 0); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 796; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __pyx_t_4 = PySequence_ITEM(sequence, 1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 796; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      #endif
    } else if (1) {
      __Pyx_RaiseNoneNotIterableError(); {__pyx_filename = __pyx_f[2]; __pyx_lineno = 796; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    } else
    {
      Py_ssize_t index = -1;
      __pyx_t_5 = PyObject_GetIter(((PyObject *)__pyx_v_fields)); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 796; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_6 = Py_TYPE(__pyx_t_5)->tp_iternext;
      index = 0; __pyx_t_3 = __pyx_t_6(__pyx_t_5); if (unlikely(!__pyx_t_3)) goto __pyx_L5_unpacking_failed;
      __Pyx_GOTREF(__pyx_t_3);
      index = 1; __pyx_t_4 = __pyx_t_6(__pyx_t_5); if (unlikely(!__pyx_t_4)) goto __pyx_L5_unpacking_failed;
      __Pyx_GOTREF(__pyx_t_4);
      if (__Pyx_IternextUnpackEndCheck(__pyx_t_6(__pyx_t_5), 2) < 0) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 796; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __pyx_t_6 = NULL;
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      goto __pyx_L6_unpacking_done;
      __pyx_L5_unpacking_failed:;
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = NULL;
      if (__Pyx_IterFinish() == 0) __Pyx_RaiseNeedMoreValuesError(index);
      {__pyx_filename = __pyx_f[2]; __pyx_lineno = 796; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __pyx_L6_unpacking_done:;
    }
    if (!(likely(((__pyx_t_3) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_3, __pyx_ptype_5numpy_dtype))))) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 796; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_XDECREF(((PyObject *)__pyx_v_child));
    __pyx_v_child = ((PyArray_Descr *)__pyx_t_3);
    __pyx_t_3 = 0;
    __Pyx_XDECREF(__pyx_v_new_offset);
    __pyx_v_new_offset = __pyx_t_4;
    __pyx_t_4 = 0;

    
    __pyx_t_4 = PyInt_FromLong((__pyx_v_end - __pyx_v_f)); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 798; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __pyx_t_3 = PyInt_FromLong((__pyx_v_offset[0])); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 798; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __pyx_t_5 = PyNumber_Subtract(__pyx_v_new_offset, __pyx_t_3); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 798; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_5);
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    __pyx_t_3 = PyNumber_Subtract(__pyx_t_4, __pyx_t_5); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 798; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
    __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
    __pyx_t_5 = PyObject_RichCompare(__pyx_t_3, __pyx_int_15, Py_LT); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 798; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 798; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
    if (__pyx_t_7) {

      
      __pyx_t_5 = PyObject_Call(__pyx_builtin_RuntimeError, ((PyObject *)__pyx_k_tuple_14), NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 799; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_Raise(__pyx_t_5, 0, 0, 0);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      {__pyx_filename = __pyx_f[2]; __pyx_lineno = 799; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      goto __pyx_L7;
    }
    __pyx_L7:;

    
    __pyx_t_7 = (__pyx_v_child->byteorder == '>');
    if (__pyx_t_7) {
      __pyx_t_8 = __pyx_v_little_endian;
    } else {
      __pyx_t_8 = __pyx_t_7;
    }
    if (!__pyx_t_8) {

      
      __pyx_t_7 = (__pyx_v_child->byteorder == '<');
      if (__pyx_t_7) {
        __pyx_t_9 = (!__pyx_v_little_endian);
        __pyx_t_10 = __pyx_t_9;
      } else {
        __pyx_t_10 = __pyx_t_7;
      }
      __pyx_t_7 = __pyx_t_10;
    } else {
      __pyx_t_7 = __pyx_t_8;
    }
    if (__pyx_t_7) {

      
      __pyx_t_5 = PyObject_Call(__pyx_builtin_ValueError, ((PyObject *)__pyx_k_tuple_15), NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 803; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_Raise(__pyx_t_5, 0, 0, 0);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      {__pyx_filename = __pyx_f[2]; __pyx_lineno = 803; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      goto __pyx_L8;
    }
    __pyx_L8:;

    
    while (1) {
      __pyx_t_5 = PyInt_FromLong((__pyx_v_offset[0])); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 813; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_t_5, __pyx_v_new_offset, Py_LT); __Pyx_XGOTREF(__pyx_t_3); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 813; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 813; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (!__pyx_t_7) break;

      
      (__pyx_v_f[0]) = 120;

      
      __pyx_v_f = (__pyx_v_f + 1);

      
      __pyx_t_11 = 0;
      (__pyx_v_offset[__pyx_t_11]) = ((__pyx_v_offset[__pyx_t_11]) + 1);
    }

    
    __pyx_t_11 = 0;
    (__pyx_v_offset[__pyx_t_11]) = ((__pyx_v_offset[__pyx_t_11]) + __pyx_v_child->elsize);

    
    __pyx_t_7 = (!PyDataType_HASFIELDS(__pyx_v_child));
    if (__pyx_t_7) {

      
      __pyx_t_3 = PyInt_FromLong(__pyx_v_child->type_num); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 821; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_XDECREF(__pyx_v_t);
      __pyx_v_t = __pyx_t_3;
      __pyx_t_3 = 0;

      
      __pyx_t_7 = ((__pyx_v_end - __pyx_v_f) < 5);
      if (__pyx_t_7) {

        
        __pyx_t_3 = PyObject_Call(__pyx_builtin_RuntimeError, ((PyObject *)__pyx_k_tuple_17), NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 823; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        __Pyx_GOTREF(__pyx_t_3);
        __Pyx_Raise(__pyx_t_3, 0, 0, 0);
        __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
        {__pyx_filename = __pyx_f[2]; __pyx_lineno = 823; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        goto __pyx_L12;
      }
      __pyx_L12:;

      
      __pyx_t_3 = PyInt_FromLong(NPY_BYTE); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 826; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 826; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 826; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 98;
        goto __pyx_L13;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_UBYTE); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 827; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); __Pyx_XGOTREF(__pyx_t_3); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 827; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 827; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 66;
        goto __pyx_L13;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_SHORT); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 828; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 828; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 828; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 104;
        goto __pyx_L13;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_USHORT); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 829; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); __Pyx_XGOTREF(__pyx_t_3); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 829; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 829; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 72;
        goto __pyx_L13;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_INT); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 830; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 830; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 830; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 105;
        goto __pyx_L13;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_UINT); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 831; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); __Pyx_XGOTREF(__pyx_t_3); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 831; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 831; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 73;
        goto __pyx_L13;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_LONG); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 832; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 832; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 832; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 108;
        goto __pyx_L13;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_ULONG); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 833; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); __Pyx_XGOTREF(__pyx_t_3); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 833; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 833; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 76;
        goto __pyx_L13;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_LONGLONG); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 834; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 834; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 834; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 113;
        goto __pyx_L13;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_ULONGLONG); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 835; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); __Pyx_XGOTREF(__pyx_t_3); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 835; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 835; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 81;
        goto __pyx_L13;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_FLOAT); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 836; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 836; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 836; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 102;
        goto __pyx_L13;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_DOUBLE); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 837; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); __Pyx_XGOTREF(__pyx_t_3); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 837; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 837; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 100;
        goto __pyx_L13;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_LONGDOUBLE); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 838; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 838; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 838; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 103;
        goto __pyx_L13;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_CFLOAT); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 839; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); __Pyx_XGOTREF(__pyx_t_3); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 839; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 839; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 90;
        (__pyx_v_f[1]) = 102;
        __pyx_v_f = (__pyx_v_f + 1);
        goto __pyx_L13;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_CDOUBLE); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 840; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 840; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 840; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 90;
        (__pyx_v_f[1]) = 100;
        __pyx_v_f = (__pyx_v_f + 1);
        goto __pyx_L13;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_CLONGDOUBLE); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 841; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); __Pyx_XGOTREF(__pyx_t_3); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 841; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 841; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 90;
        (__pyx_v_f[1]) = 103;
        __pyx_v_f = (__pyx_v_f + 1);
        goto __pyx_L13;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_OBJECT); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 842; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); __Pyx_XGOTREF(__pyx_t_5); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 842; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_7 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 842; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_7) {
        (__pyx_v_f[0]) = 79;
        goto __pyx_L13;
      }
       {

        
        __pyx_t_5 = PyNumber_Remainder(((PyObject *)__pyx_kp_u_12), __pyx_v_t); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 844; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        __Pyx_GOTREF(((PyObject *)__pyx_t_5));
        __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 844; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        __Pyx_GOTREF(__pyx_t_3);
        PyTuple_SET_ITEM(__pyx_t_3, 0, ((PyObject *)__pyx_t_5));
        __Pyx_GIVEREF(((PyObject *)__pyx_t_5));
        __pyx_t_5 = 0;
        __pyx_t_5 = PyObject_Call(__pyx_builtin_ValueError, ((PyObject *)__pyx_t_3), NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 844; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        __Pyx_GOTREF(__pyx_t_5);
        __Pyx_DECREF(((PyObject *)__pyx_t_3)); __pyx_t_3 = 0;
        __Pyx_Raise(__pyx_t_5, 0, 0, 0);
        __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
        {__pyx_filename = __pyx_f[2]; __pyx_lineno = 844; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      __pyx_L13:;

      
      __pyx_v_f = (__pyx_v_f + 1);
      goto __pyx_L11;
    }
     {

      
      __pyx_t_12 = __pyx_f_5numpy__util_dtypestring(__pyx_v_child, __pyx_v_f, __pyx_v_end, __pyx_v_offset); if (unlikely(__pyx_t_12 == NULL)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 849; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __pyx_v_f = __pyx_t_12;
    }
    __pyx_L11:;
  }
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_r = __pyx_v_f;
  goto __pyx_L0;

  __pyx_r = 0;
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_AddTraceback("numpy._util_dtypestring", __pyx_clineno, __pyx_lineno, __pyx_filename);
  __pyx_r = NULL;
  __pyx_L0:;
  __Pyx_XDECREF((PyObject *)__pyx_v_child);
  __Pyx_XDECREF(__pyx_v_fields);
  __Pyx_XDECREF(__pyx_v_childname);
  __Pyx_XDECREF(__pyx_v_new_offset);
  __Pyx_XDECREF(__pyx_v_t);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE void __pyx_f_5numpy_set_array_base(PyArrayObject *__pyx_v_arr, PyObject *__pyx_v_base) {
  PyObject *__pyx_v_baseptr;
  __Pyx_RefNannyDeclarations
  int __pyx_t_1;
  __Pyx_RefNannySetupContext("set_array_base", 0);

  
  __pyx_t_1 = (__pyx_v_base == Py_None);
  if (__pyx_t_1) {

    
    __pyx_v_baseptr = NULL;
    goto __pyx_L3;
  }
   {

    
    Py_INCREF(__pyx_v_base);

    
    __pyx_v_baseptr = ((PyObject *)__pyx_v_base);
  }
  __pyx_L3:;

  
  Py_XDECREF(__pyx_v_arr->base);

  
  __pyx_v_arr->base = __pyx_v_baseptr;

  __Pyx_RefNannyFinishContext();
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_get_array_base(PyArrayObject *__pyx_v_arr) {
  PyObject *__pyx_r = NULL;
  __Pyx_RefNannyDeclarations
  int __pyx_t_1;
  __Pyx_RefNannySetupContext("get_array_base", 0);

  
  __pyx_t_1 = (__pyx_v_arr->base == NULL);
  if (__pyx_t_1) {

    
    __Pyx_XDECREF(__pyx_r);
    __Pyx_INCREF(Py_None);
    __pyx_r = Py_None;
    goto __pyx_L0;
    goto __pyx_L3;
  }
   {

    
    __Pyx_XDECREF(__pyx_r);
    __Pyx_INCREF(((PyObject *)__pyx_v_arr->base));
    __pyx_r = ((PyObject *)__pyx_v_arr->base);
    goto __pyx_L0;
  }
  __pyx_L3:;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}

static PyMethodDef __pyx_methods[] = {
  {0, 0, 0, 0}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef __pyx_moduledef = {
    PyModuleDef_HEAD_INIT,
    __Pyx_NAMESTR("_ufuncs"),
    0, 
    -1, 
    __pyx_methods ,
    NULL, 
    NULL, 
    NULL, 
    NULL 
};
#endif

static __Pyx_StringTabEntry __pyx_string_tab[] = {
  {&__pyx_kp_u_10, __pyx_k_10, sizeof(__pyx_k_10), 0, 1, 0, 0},
  {&__pyx_kp_u_12, __pyx_k_12, sizeof(__pyx_k_12), 0, 1, 0, 0},
  {&__pyx_kp_u_13, __pyx_k_13, sizeof(__pyx_k_13), 0, 1, 0, 0},
  {&__pyx_kp_u_16, __pyx_k_16, sizeof(__pyx_k_16), 0, 1, 0, 0},
  {&__pyx_kp_s_194, __pyx_k_194, sizeof(__pyx_k_194), 0, 0, 1, 0},
  {&__pyx_n_s_195, __pyx_k_195, sizeof(__pyx_k_195), 0, 0, 1, 1},
  {&__pyx_kp_b_2, __pyx_k_2, sizeof(__pyx_k_2), 0, 0, 0, 0},
  {&__pyx_kp_b_3, __pyx_k_3, sizeof(__pyx_k_3), 0, 0, 0, 0},
  {&__pyx_kp_b_4, __pyx_k_4, sizeof(__pyx_k_4), 0, 0, 0, 0},
  {&__pyx_kp_u_6, __pyx_k_6, sizeof(__pyx_k_6), 0, 1, 0, 0},
  {&__pyx_kp_u_8, __pyx_k_8, sizeof(__pyx_k_8), 0, 1, 0, 0},
  {&__pyx_n_s__FutureWarning, __pyx_k__FutureWarning, sizeof(__pyx_k__FutureWarning), 0, 0, 1, 1},
  {&__pyx_n_s__RuntimeError, __pyx_k__RuntimeError, sizeof(__pyx_k__RuntimeError), 0, 0, 1, 1},
  {&__pyx_n_s__ValueError, __pyx_k__ValueError, sizeof(__pyx_k__ValueError), 0, 0, 1, 1},
  {&__pyx_n_s____main__, __pyx_k____main__, sizeof(__pyx_k____main__), 0, 0, 1, 1},
  {&__pyx_n_s____test__, __pyx_k____test__, sizeof(__pyx_k____test__), 0, 0, 1, 1},
  {&__pyx_n_s___errprint, __pyx_k___errprint, sizeof(__pyx_k___errprint), 0, 0, 1, 1},
  {&__pyx_n_s___lambertw, __pyx_k___lambertw, sizeof(__pyx_k___lambertw), 0, 0, 1, 1},
  {&__pyx_n_s__airy, __pyx_k__airy, sizeof(__pyx_k__airy), 0, 0, 1, 1},
  {&__pyx_n_s__airye, __pyx_k__airye, sizeof(__pyx_k__airye), 0, 0, 1, 1},
  {&__pyx_n_s__bdtr, __pyx_k__bdtr, sizeof(__pyx_k__bdtr), 0, 0, 1, 1},
  {&__pyx_n_s__bdtrc, __pyx_k__bdtrc, sizeof(__pyx_k__bdtrc), 0, 0, 1, 1},
  {&__pyx_n_s__bdtri, __pyx_k__bdtri, sizeof(__pyx_k__bdtri), 0, 0, 1, 1},
  {&__pyx_n_s__bdtrik, __pyx_k__bdtrik, sizeof(__pyx_k__bdtrik), 0, 0, 1, 1},
  {&__pyx_n_s__bdtrin, __pyx_k__bdtrin, sizeof(__pyx_k__bdtrin), 0, 0, 1, 1},
  {&__pyx_n_s__bei, __pyx_k__bei, sizeof(__pyx_k__bei), 0, 0, 1, 1},
  {&__pyx_n_s__beip, __pyx_k__beip, sizeof(__pyx_k__beip), 0, 0, 1, 1},
  {&__pyx_n_s__ber, __pyx_k__ber, sizeof(__pyx_k__ber), 0, 0, 1, 1},
  {&__pyx_n_s__berp, __pyx_k__berp, sizeof(__pyx_k__berp), 0, 0, 1, 1},
  {&__pyx_n_s__besselpoly, __pyx_k__besselpoly, sizeof(__pyx_k__besselpoly), 0, 0, 1, 1},
  {&__pyx_n_s__beta, __pyx_k__beta, sizeof(__pyx_k__beta), 0, 0, 1, 1},
  {&__pyx_n_s__betainc, __pyx_k__betainc, sizeof(__pyx_k__betainc), 0, 0, 1, 1},
  {&__pyx_n_s__betaincinv, __pyx_k__betaincinv, sizeof(__pyx_k__betaincinv), 0, 0, 1, 1},
  {&__pyx_n_s__betaln, __pyx_k__betaln, sizeof(__pyx_k__betaln), 0, 0, 1, 1},
  {&__pyx_n_s__binom, __pyx_k__binom, sizeof(__pyx_k__binom), 0, 0, 1, 1},
  {&__pyx_n_s__btdtr, __pyx_k__btdtr, sizeof(__pyx_k__btdtr), 0, 0, 1, 1},
  {&__pyx_n_s__btdtri, __pyx_k__btdtri, sizeof(__pyx_k__btdtri), 0, 0, 1, 1},
  {&__pyx_n_s__btdtria, __pyx_k__btdtria, sizeof(__pyx_k__btdtria), 0, 0, 1, 1},
  {&__pyx_n_s__btdtrib, __pyx_k__btdtrib, sizeof(__pyx_k__btdtrib), 0, 0, 1, 1},
  {&__pyx_n_s__cbrt, __pyx_k__cbrt, sizeof(__pyx_k__cbrt), 0, 0, 1, 1},
  {&__pyx_n_s__chdtr, __pyx_k__chdtr, sizeof(__pyx_k__chdtr), 0, 0, 1, 1},
  {&__pyx_n_s__chdtrc, __pyx_k__chdtrc, sizeof(__pyx_k__chdtrc), 0, 0, 1, 1},
  {&__pyx_n_s__chdtri, __pyx_k__chdtri, sizeof(__pyx_k__chdtri), 0, 0, 1, 1},
  {&__pyx_n_s__chdtriv, __pyx_k__chdtriv, sizeof(__pyx_k__chdtriv), 0, 0, 1, 1},
  {&__pyx_n_s__chndtr, __pyx_k__chndtr, sizeof(__pyx_k__chndtr), 0, 0, 1, 1},
  {&__pyx_n_s__chndtridf, __pyx_k__chndtridf, sizeof(__pyx_k__chndtridf), 0, 0, 1, 1},
  {&__pyx_n_s__chndtrinc, __pyx_k__chndtrinc, sizeof(__pyx_k__chndtrinc), 0, 0, 1, 1},
  {&__pyx_n_s__chndtrix, __pyx_k__chndtrix, sizeof(__pyx_k__chndtrix), 0, 0, 1, 1},
  {&__pyx_n_s__cosdg, __pyx_k__cosdg, sizeof(__pyx_k__cosdg), 0, 0, 1, 1},
  {&__pyx_n_s__cosm1, __pyx_k__cosm1, sizeof(__pyx_k__cosm1), 0, 0, 1, 1},
  {&__pyx_n_s__cotdg, __pyx_k__cotdg, sizeof(__pyx_k__cotdg), 0, 0, 1, 1},
  {&__pyx_n_s__ellipe, __pyx_k__ellipe, sizeof(__pyx_k__ellipe), 0, 0, 1, 1},
  {&__pyx_n_s__ellipeinc, __pyx_k__ellipeinc, sizeof(__pyx_k__ellipeinc), 0, 0, 1, 1},
  {&__pyx_n_s__ellipj, __pyx_k__ellipj, sizeof(__pyx_k__ellipj), 0, 0, 1, 1},
  {&__pyx_n_s__ellipkinc, __pyx_k__ellipkinc, sizeof(__pyx_k__ellipkinc), 0, 0, 1, 1},
  {&__pyx_n_s__ellipkm1, __pyx_k__ellipkm1, sizeof(__pyx_k__ellipkm1), 0, 0, 1, 1},
  {&__pyx_n_s__eval_chebyc, __pyx_k__eval_chebyc, sizeof(__pyx_k__eval_chebyc), 0, 0, 1, 1},
  {&__pyx_n_s__eval_chebys, __pyx_k__eval_chebys, sizeof(__pyx_k__eval_chebys), 0, 0, 1, 1},
  {&__pyx_n_s__eval_chebyt, __pyx_k__eval_chebyt, sizeof(__pyx_k__eval_chebyt), 0, 0, 1, 1},
  {&__pyx_n_s__eval_chebyu, __pyx_k__eval_chebyu, sizeof(__pyx_k__eval_chebyu), 0, 0, 1, 1},
  {&__pyx_n_s__eval_gegenbauer, __pyx_k__eval_gegenbauer, sizeof(__pyx_k__eval_gegenbauer), 0, 0, 1, 1},
  {&__pyx_n_s__eval_genlaguerre, __pyx_k__eval_genlaguerre, sizeof(__pyx_k__eval_genlaguerre), 0, 0, 1, 1},
  {&__pyx_n_s__eval_hermite, __pyx_k__eval_hermite, sizeof(__pyx_k__eval_hermite), 0, 0, 1, 1},
  {&__pyx_n_s__eval_hermitenorm, __pyx_k__eval_hermitenorm, sizeof(__pyx_k__eval_hermitenorm), 0, 0, 1, 1},
  {&__pyx_n_s__eval_jacobi, __pyx_k__eval_jacobi, sizeof(__pyx_k__eval_jacobi), 0, 0, 1, 1},
  {&__pyx_n_s__eval_laguerre, __pyx_k__eval_laguerre, sizeof(__pyx_k__eval_laguerre), 0, 0, 1, 1},
  {&__pyx_n_s__eval_legendre, __pyx_k__eval_legendre, sizeof(__pyx_k__eval_legendre), 0, 0, 1, 1},
  {&__pyx_n_s__eval_sh_chebyt, __pyx_k__eval_sh_chebyt, sizeof(__pyx_k__eval_sh_chebyt), 0, 0, 1, 1},
  {&__pyx_n_s__eval_sh_chebyu, __pyx_k__eval_sh_chebyu, sizeof(__pyx_k__eval_sh_chebyu), 0, 0, 1, 1},
  {&__pyx_n_s__eval_sh_jacobi, __pyx_k__eval_sh_jacobi, sizeof(__pyx_k__eval_sh_jacobi), 0, 0, 1, 1},
  {&__pyx_n_s__eval_sh_legendre, __pyx_k__eval_sh_legendre, sizeof(__pyx_k__eval_sh_legendre), 0, 0, 1, 1},
  {&__pyx_n_s__exp1, __pyx_k__exp1, sizeof(__pyx_k__exp1), 0, 0, 1, 1},
  {&__pyx_n_s__exp10, __pyx_k__exp10, sizeof(__pyx_k__exp10), 0, 0, 1, 1},
  {&__pyx_n_s__exp2, __pyx_k__exp2, sizeof(__pyx_k__exp2), 0, 0, 1, 1},
  {&__pyx_n_s__expi, __pyx_k__expi, sizeof(__pyx_k__expi), 0, 0, 1, 1},
  {&__pyx_n_s__expit, __pyx_k__expit, sizeof(__pyx_k__expit), 0, 0, 1, 1},
  {&__pyx_n_s__expm1, __pyx_k__expm1, sizeof(__pyx_k__expm1), 0, 0, 1, 1},
  {&__pyx_n_s__expn, __pyx_k__expn, sizeof(__pyx_k__expn), 0, 0, 1, 1},
  {&__pyx_n_s__fdtr, __pyx_k__fdtr, sizeof(__pyx_k__fdtr), 0, 0, 1, 1},
  {&__pyx_n_s__fdtrc, __pyx_k__fdtrc, sizeof(__pyx_k__fdtrc), 0, 0, 1, 1},
  {&__pyx_n_s__fdtri, __pyx_k__fdtri, sizeof(__pyx_k__fdtri), 0, 0, 1, 1},
  {&__pyx_n_s__fdtridfd, __pyx_k__fdtridfd, sizeof(__pyx_k__fdtridfd), 0, 0, 1, 1},
  {&__pyx_n_s__fresnel, __pyx_k__fresnel, sizeof(__pyx_k__fresnel), 0, 0, 1, 1},
  {&__pyx_n_s__gamma, __pyx_k__gamma, sizeof(__pyx_k__gamma), 0, 0, 1, 1},
  {&__pyx_n_s__gammainc, __pyx_k__gammainc, sizeof(__pyx_k__gammainc), 0, 0, 1, 1},
  {&__pyx_n_s__gammaincc, __pyx_k__gammaincc, sizeof(__pyx_k__gammaincc), 0, 0, 1, 1},
  {&__pyx_n_s__gammainccinv, __pyx_k__gammainccinv, sizeof(__pyx_k__gammainccinv), 0, 0, 1, 1},
  {&__pyx_n_s__gammaincinv, __pyx_k__gammaincinv, sizeof(__pyx_k__gammaincinv), 0, 0, 1, 1},
  {&__pyx_n_s__gammaln, __pyx_k__gammaln, sizeof(__pyx_k__gammaln), 0, 0, 1, 1},
  {&__pyx_n_s__gammasgn, __pyx_k__gammasgn, sizeof(__pyx_k__gammasgn), 0, 0, 1, 1},
  {&__pyx_n_s__gdtr, __pyx_k__gdtr, sizeof(__pyx_k__gdtr), 0, 0, 1, 1},
  {&__pyx_n_s__gdtrc, __pyx_k__gdtrc, sizeof(__pyx_k__gdtrc), 0, 0, 1, 1},
  {&__pyx_n_s__gdtria, __pyx_k__gdtria, sizeof(__pyx_k__gdtria), 0, 0, 1, 1},
  {&__pyx_n_s__gdtrib, __pyx_k__gdtrib, sizeof(__pyx_k__gdtrib), 0, 0, 1, 1},
  {&__pyx_n_s__gdtrix, __pyx_k__gdtrix, sizeof(__pyx_k__gdtrix), 0, 0, 1, 1},
  {&__pyx_n_s__hankel1, __pyx_k__hankel1, sizeof(__pyx_k__hankel1), 0, 0, 1, 1},
  {&__pyx_n_s__hankel1e, __pyx_k__hankel1e, sizeof(__pyx_k__hankel1e), 0, 0, 1, 1},
  {&__pyx_n_s__hankel2, __pyx_k__hankel2, sizeof(__pyx_k__hankel2), 0, 0, 1, 1},
  {&__pyx_n_s__hankel2e, __pyx_k__hankel2e, sizeof(__pyx_k__hankel2e), 0, 0, 1, 1},
  {&__pyx_n_s__hyp1f1, __pyx_k__hyp1f1, sizeof(__pyx_k__hyp1f1), 0, 0, 1, 1},
  {&__pyx_n_s__hyp1f2, __pyx_k__hyp1f2, sizeof(__pyx_k__hyp1f2), 0, 0, 1, 1},
  {&__pyx_n_s__hyp2f0, __pyx_k__hyp2f0, sizeof(__pyx_k__hyp2f0), 0, 0, 1, 1},
  {&__pyx_n_s__hyp2f1, __pyx_k__hyp2f1, sizeof(__pyx_k__hyp2f1), 0, 0, 1, 1},
  {&__pyx_n_s__hyp3f0, __pyx_k__hyp3f0, sizeof(__pyx_k__hyp3f0), 0, 0, 1, 1},
  {&__pyx_n_s__hyperu, __pyx_k__hyperu, sizeof(__pyx_k__hyperu), 0, 0, 1, 1},
  {&__pyx_n_s__i0, __pyx_k__i0, sizeof(__pyx_k__i0), 0, 0, 1, 1},
  {&__pyx_n_s__i0e, __pyx_k__i0e, sizeof(__pyx_k__i0e), 0, 0, 1, 1},
  {&__pyx_n_s__i1, __pyx_k__i1, sizeof(__pyx_k__i1), 0, 0, 1, 1},
  {&__pyx_n_s__i1e, __pyx_k__i1e, sizeof(__pyx_k__i1e), 0, 0, 1, 1},
  {&__pyx_n_s__inflag, __pyx_k__inflag, sizeof(__pyx_k__inflag), 0, 0, 1, 1},
  {&__pyx_n_s__it2i0k0, __pyx_k__it2i0k0, sizeof(__pyx_k__it2i0k0), 0, 0, 1, 1},
  {&__pyx_n_s__it2j0y0, __pyx_k__it2j0y0, sizeof(__pyx_k__it2j0y0), 0, 0, 1, 1},
  {&__pyx_n_s__it2struve0, __pyx_k__it2struve0, sizeof(__pyx_k__it2struve0), 0, 0, 1, 1},
  {&__pyx_n_s__itairy, __pyx_k__itairy, sizeof(__pyx_k__itairy), 0, 0, 1, 1},
  {&__pyx_n_s__iti0k0, __pyx_k__iti0k0, sizeof(__pyx_k__iti0k0), 0, 0, 1, 1},
  {&__pyx_n_s__itj0y0, __pyx_k__itj0y0, sizeof(__pyx_k__itj0y0), 0, 0, 1, 1},
  {&__pyx_n_s__itmodstruve0, __pyx_k__itmodstruve0, sizeof(__pyx_k__itmodstruve0), 0, 0, 1, 1},
  {&__pyx_n_s__itstruve0, __pyx_k__itstruve0, sizeof(__pyx_k__itstruve0), 0, 0, 1, 1},
  {&__pyx_n_s__iv, __pyx_k__iv, sizeof(__pyx_k__iv), 0, 0, 1, 1},
  {&__pyx_n_s__ive, __pyx_k__ive, sizeof(__pyx_k__ive), 0, 0, 1, 1},
  {&__pyx_n_s__j0, __pyx_k__j0, sizeof(__pyx_k__j0), 0, 0, 1, 1},
  {&__pyx_n_s__j1, __pyx_k__j1, sizeof(__pyx_k__j1), 0, 0, 1, 1},
  {&__pyx_n_s__jn, __pyx_k__jn, sizeof(__pyx_k__jn), 0, 0, 1, 1},
  {&__pyx_n_s__jv, __pyx_k__jv, sizeof(__pyx_k__jv), 0, 0, 1, 1},
  {&__pyx_n_s__jve, __pyx_k__jve, sizeof(__pyx_k__jve), 0, 0, 1, 1},
  {&__pyx_n_s__k0, __pyx_k__k0, sizeof(__pyx_k__k0), 0, 0, 1, 1},
  {&__pyx_n_s__k0e, __pyx_k__k0e, sizeof(__pyx_k__k0e), 0, 0, 1, 1},
  {&__pyx_n_s__k1, __pyx_k__k1, sizeof(__pyx_k__k1), 0, 0, 1, 1},
  {&__pyx_n_s__k1e, __pyx_k__k1e, sizeof(__pyx_k__k1e), 0, 0, 1, 1},
  {&__pyx_n_s__kei, __pyx_k__kei, sizeof(__pyx_k__kei), 0, 0, 1, 1},
  {&__pyx_n_s__keip, __pyx_k__keip, sizeof(__pyx_k__keip), 0, 0, 1, 1},
  {&__pyx_n_s__kelvin, __pyx_k__kelvin, sizeof(__pyx_k__kelvin), 0, 0, 1, 1},
  {&__pyx_n_s__ker, __pyx_k__ker, sizeof(__pyx_k__ker), 0, 0, 1, 1},
  {&__pyx_n_s__kerp, __pyx_k__kerp, sizeof(__pyx_k__kerp), 0, 0, 1, 1},
  {&__pyx_n_s__kn, __pyx_k__kn, sizeof(__pyx_k__kn), 0, 0, 1, 1},
  {&__pyx_n_s__kolmogi, __pyx_k__kolmogi, sizeof(__pyx_k__kolmogi), 0, 0, 1, 1},
  {&__pyx_n_s__kolmogorov, __pyx_k__kolmogorov, sizeof(__pyx_k__kolmogorov), 0, 0, 1, 1},
  {&__pyx_n_s__kv, __pyx_k__kv, sizeof(__pyx_k__kv), 0, 0, 1, 1},
  {&__pyx_n_s__kve, __pyx_k__kve, sizeof(__pyx_k__kve), 0, 0, 1, 1},
  {&__pyx_n_s__log1p, __pyx_k__log1p, sizeof(__pyx_k__log1p), 0, 0, 1, 1},
  {&__pyx_n_s__log_ndtr, __pyx_k__log_ndtr, sizeof(__pyx_k__log_ndtr), 0, 0, 1, 1},
  {&__pyx_n_s__logit, __pyx_k__logit, sizeof(__pyx_k__logit), 0, 0, 1, 1},
  {&__pyx_n_s__lpmv, __pyx_k__lpmv, sizeof(__pyx_k__lpmv), 0, 0, 1, 1},
  {&__pyx_n_s__mathieu_a, __pyx_k__mathieu_a, sizeof(__pyx_k__mathieu_a), 0, 0, 1, 1},
  {&__pyx_n_s__mathieu_b, __pyx_k__mathieu_b, sizeof(__pyx_k__mathieu_b), 0, 0, 1, 1},
  {&__pyx_n_s__mathieu_cem, __pyx_k__mathieu_cem, sizeof(__pyx_k__mathieu_cem), 0, 0, 1, 1},
  {&__pyx_n_s__mathieu_modcem1, __pyx_k__mathieu_modcem1, sizeof(__pyx_k__mathieu_modcem1), 0, 0, 1, 1},
  {&__pyx_n_s__mathieu_modcem2, __pyx_k__mathieu_modcem2, sizeof(__pyx_k__mathieu_modcem2), 0, 0, 1, 1},
  {&__pyx_n_s__mathieu_modsem1, __pyx_k__mathieu_modsem1, sizeof(__pyx_k__mathieu_modsem1), 0, 0, 1, 1},
  {&__pyx_n_s__mathieu_modsem2, __pyx_k__mathieu_modsem2, sizeof(__pyx_k__mathieu_modsem2), 0, 0, 1, 1},
  {&__pyx_n_s__mathieu_sem, __pyx_k__mathieu_sem, sizeof(__pyx_k__mathieu_sem), 0, 0, 1, 1},
  {&__pyx_n_s__modfresnelm, __pyx_k__modfresnelm, sizeof(__pyx_k__modfresnelm), 0, 0, 1, 1},
  {&__pyx_n_s__modfresnelp, __pyx_k__modfresnelp, sizeof(__pyx_k__modfresnelp), 0, 0, 1, 1},
  {&__pyx_n_s__modstruve, __pyx_k__modstruve, sizeof(__pyx_k__modstruve), 0, 0, 1, 1},
  {&__pyx_n_s__nbdtr, __pyx_k__nbdtr, sizeof(__pyx_k__nbdtr), 0, 0, 1, 1},
  {&__pyx_n_s__nbdtrc, __pyx_k__nbdtrc, sizeof(__pyx_k__nbdtrc), 0, 0, 1, 1},
  {&__pyx_n_s__nbdtri, __pyx_k__nbdtri, sizeof(__pyx_k__nbdtri), 0, 0, 1, 1},
  {&__pyx_n_s__nbdtrik, __pyx_k__nbdtrik, sizeof(__pyx_k__nbdtrik), 0, 0, 1, 1},
  {&__pyx_n_s__nbdtrin, __pyx_k__nbdtrin, sizeof(__pyx_k__nbdtrin), 0, 0, 1, 1},
  {&__pyx_n_s__ncfdtr, __pyx_k__ncfdtr, sizeof(__pyx_k__ncfdtr), 0, 0, 1, 1},
  {&__pyx_n_s__ncfdtri, __pyx_k__ncfdtri, sizeof(__pyx_k__ncfdtri), 0, 0, 1, 1},
  {&__pyx_n_s__ncfdtridfd, __pyx_k__ncfdtridfd, sizeof(__pyx_k__ncfdtridfd), 0, 0, 1, 1},
  {&__pyx_n_s__ncfdtridfn, __pyx_k__ncfdtridfn, sizeof(__pyx_k__ncfdtridfn), 0, 0, 1, 1},
  {&__pyx_n_s__ncfdtrinc, __pyx_k__ncfdtrinc, sizeof(__pyx_k__ncfdtrinc), 0, 0, 1, 1},
  {&__pyx_n_s__nctdtr, __pyx_k__nctdtr, sizeof(__pyx_k__nctdtr), 0, 0, 1, 1},
  {&__pyx_n_s__nctdtridf, __pyx_k__nctdtridf, sizeof(__pyx_k__nctdtridf), 0, 0, 1, 1},
  {&__pyx_n_s__nctdtrinc, __pyx_k__nctdtrinc, sizeof(__pyx_k__nctdtrinc), 0, 0, 1, 1},
  {&__pyx_n_s__nctdtrit, __pyx_k__nctdtrit, sizeof(__pyx_k__nctdtrit), 0, 0, 1, 1},
  {&__pyx_n_s__ndtr, __pyx_k__ndtr, sizeof(__pyx_k__ndtr), 0, 0, 1, 1},
  {&__pyx_n_s__ndtri, __pyx_k__ndtri, sizeof(__pyx_k__ndtri), 0, 0, 1, 1},
  {&__pyx_n_s__nrdtrimn, __pyx_k__nrdtrimn, sizeof(__pyx_k__nrdtrimn), 0, 0, 1, 1},
  {&__pyx_n_s__nrdtrisd, __pyx_k__nrdtrisd, sizeof(__pyx_k__nrdtrisd), 0, 0, 1, 1},
  {&__pyx_n_s__obl_ang1, __pyx_k__obl_ang1, sizeof(__pyx_k__obl_ang1), 0, 0, 1, 1},
  {&__pyx_n_s__obl_ang1_cv, __pyx_k__obl_ang1_cv, sizeof(__pyx_k__obl_ang1_cv), 0, 0, 1, 1},
  {&__pyx_n_s__obl_cv, __pyx_k__obl_cv, sizeof(__pyx_k__obl_cv), 0, 0, 1, 1},
  {&__pyx_n_s__obl_rad1, __pyx_k__obl_rad1, sizeof(__pyx_k__obl_rad1), 0, 0, 1, 1},
  {&__pyx_n_s__obl_rad1_cv, __pyx_k__obl_rad1_cv, sizeof(__pyx_k__obl_rad1_cv), 0, 0, 1, 1},
  {&__pyx_n_s__obl_rad2, __pyx_k__obl_rad2, sizeof(__pyx_k__obl_rad2), 0, 0, 1, 1},
  {&__pyx_n_s__obl_rad2_cv, __pyx_k__obl_rad2_cv, sizeof(__pyx_k__obl_rad2_cv), 0, 0, 1, 1},
  {&__pyx_n_s__pbdv, __pyx_k__pbdv, sizeof(__pyx_k__pbdv), 0, 0, 1, 1},
  {&__pyx_n_s__pbvv, __pyx_k__pbvv, sizeof(__pyx_k__pbvv), 0, 0, 1, 1},
  {&__pyx_n_s__pbwa, __pyx_k__pbwa, sizeof(__pyx_k__pbwa), 0, 0, 1, 1},
  {&__pyx_n_s__pdtr, __pyx_k__pdtr, sizeof(__pyx_k__pdtr), 0, 0, 1, 1},
  {&__pyx_n_s__pdtrc, __pyx_k__pdtrc, sizeof(__pyx_k__pdtrc), 0, 0, 1, 1},
  {&__pyx_n_s__pdtri, __pyx_k__pdtri, sizeof(__pyx_k__pdtri), 0, 0, 1, 1},
  {&__pyx_n_s__pdtrik, __pyx_k__pdtrik, sizeof(__pyx_k__pdtrik), 0, 0, 1, 1},
  {&__pyx_n_s__pro_ang1, __pyx_k__pro_ang1, sizeof(__pyx_k__pro_ang1), 0, 0, 1, 1},
  {&__pyx_n_s__pro_ang1_cv, __pyx_k__pro_ang1_cv, sizeof(__pyx_k__pro_ang1_cv), 0, 0, 1, 1},
  {&__pyx_n_s__pro_cv, __pyx_k__pro_cv, sizeof(__pyx_k__pro_cv), 0, 0, 1, 1},
  {&__pyx_n_s__pro_rad1, __pyx_k__pro_rad1, sizeof(__pyx_k__pro_rad1), 0, 0, 1, 1},
  {&__pyx_n_s__pro_rad1_cv, __pyx_k__pro_rad1_cv, sizeof(__pyx_k__pro_rad1_cv), 0, 0, 1, 1},
  {&__pyx_n_s__pro_rad2, __pyx_k__pro_rad2, sizeof(__pyx_k__pro_rad2), 0, 0, 1, 1},
  {&__pyx_n_s__pro_rad2_cv, __pyx_k__pro_rad2_cv, sizeof(__pyx_k__pro_rad2_cv), 0, 0, 1, 1},
  {&__pyx_n_s__psi, __pyx_k__psi, sizeof(__pyx_k__psi), 0, 0, 1, 1},
  {&__pyx_n_s__radian, __pyx_k__radian, sizeof(__pyx_k__radian), 0, 0, 1, 1},
  {&__pyx_n_s__range, __pyx_k__range, sizeof(__pyx_k__range), 0, 0, 1, 1},
  {&__pyx_n_s__rgamma, __pyx_k__rgamma, sizeof(__pyx_k__rgamma), 0, 0, 1, 1},
  {&__pyx_n_s__round, __pyx_k__round, sizeof(__pyx_k__round), 0, 0, 1, 1},
  {&__pyx_n_s__shichi, __pyx_k__shichi, sizeof(__pyx_k__shichi), 0, 0, 1, 1},
  {&__pyx_n_s__sici, __pyx_k__sici, sizeof(__pyx_k__sici), 0, 0, 1, 1},
  {&__pyx_n_s__sindg, __pyx_k__sindg, sizeof(__pyx_k__sindg), 0, 0, 1, 1},
  {&__pyx_n_s__smirnov, __pyx_k__smirnov, sizeof(__pyx_k__smirnov), 0, 0, 1, 1},
  {&__pyx_n_s__smirnovi, __pyx_k__smirnovi, sizeof(__pyx_k__smirnovi), 0, 0, 1, 1},
  {&__pyx_n_s__spence, __pyx_k__spence, sizeof(__pyx_k__spence), 0, 0, 1, 1},
  {&__pyx_n_s__stdtr, __pyx_k__stdtr, sizeof(__pyx_k__stdtr), 0, 0, 1, 1},
  {&__pyx_n_s__stdtridf, __pyx_k__stdtridf, sizeof(__pyx_k__stdtridf), 0, 0, 1, 1},
  {&__pyx_n_s__stdtrit, __pyx_k__stdtrit, sizeof(__pyx_k__stdtrit), 0, 0, 1, 1},
  {&__pyx_n_s__struve, __pyx_k__struve, sizeof(__pyx_k__struve), 0, 0, 1, 1},
  {&__pyx_n_s__tandg, __pyx_k__tandg, sizeof(__pyx_k__tandg), 0, 0, 1, 1},
  {&__pyx_n_s__tklmbda, __pyx_k__tklmbda, sizeof(__pyx_k__tklmbda), 0, 0, 1, 1},
  {&__pyx_n_s__y0, __pyx_k__y0, sizeof(__pyx_k__y0), 0, 0, 1, 1},
  {&__pyx_n_s__y1, __pyx_k__y1, sizeof(__pyx_k__y1), 0, 0, 1, 1},
  {&__pyx_n_s__yn, __pyx_k__yn, sizeof(__pyx_k__yn), 0, 0, 1, 1},
  {&__pyx_n_s__yv, __pyx_k__yv, sizeof(__pyx_k__yv), 0, 0, 1, 1},
  {&__pyx_n_s__yve, __pyx_k__yve, sizeof(__pyx_k__yve), 0, 0, 1, 1},
  {&__pyx_n_s__zeta, __pyx_k__zeta, sizeof(__pyx_k__zeta), 0, 0, 1, 1},
  {&__pyx_n_s__zetac, __pyx_k__zetac, sizeof(__pyx_k__zetac), 0, 0, 1, 1},
  {0, 0, 0, 0, 0, 0, 0}
};
static int __Pyx_InitCachedBuiltins(void) {
  __pyx_builtin_range = __Pyx_GetName(__pyx_b, __pyx_n_s__range); if (!__pyx_builtin_range) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 34; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_builtin_FutureWarning = __Pyx_GetName(__pyx_b, __pyx_n_s__FutureWarning); if (!__pyx_builtin_FutureWarning) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 39; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_builtin_ValueError = __Pyx_GetName(__pyx_b, __pyx_n_s__ValueError); if (!__pyx_builtin_ValueError) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 215; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_builtin_RuntimeError = __Pyx_GetName(__pyx_b, __pyx_n_s__RuntimeError); if (!__pyx_builtin_RuntimeError) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 799; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  return 0;
  __pyx_L1_error:;
  return -1;
}

static int __Pyx_InitCachedConstants(void) {
  __Pyx_RefNannyDeclarations
  __Pyx_RefNannySetupContext("__Pyx_InitCachedConstants", 0);

  
  __pyx_k_tuple_7 = PyTuple_New(1); if (unlikely(!__pyx_k_tuple_7)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 215; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_k_tuple_7);
  __Pyx_INCREF(((PyObject *)__pyx_kp_u_6));
  PyTuple_SET_ITEM(__pyx_k_tuple_7, 0, ((PyObject *)__pyx_kp_u_6));
  __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_6));
  __Pyx_GIVEREF(((PyObject *)__pyx_k_tuple_7));

  
  __pyx_k_tuple_9 = PyTuple_New(1); if (unlikely(!__pyx_k_tuple_9)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 219; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_k_tuple_9);
  __Pyx_INCREF(((PyObject *)__pyx_kp_u_8));
  PyTuple_SET_ITEM(__pyx_k_tuple_9, 0, ((PyObject *)__pyx_kp_u_8));
  __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_8));
  __Pyx_GIVEREF(((PyObject *)__pyx_k_tuple_9));

  
  __pyx_k_tuple_11 = PyTuple_New(1); if (unlikely(!__pyx_k_tuple_11)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 257; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_k_tuple_11);
  __Pyx_INCREF(((PyObject *)__pyx_kp_u_10));
  PyTuple_SET_ITEM(__pyx_k_tuple_11, 0, ((PyObject *)__pyx_kp_u_10));
  __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_10));
  __Pyx_GIVEREF(((PyObject *)__pyx_k_tuple_11));

  
  __pyx_k_tuple_14 = PyTuple_New(1); if (unlikely(!__pyx_k_tuple_14)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 799; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_k_tuple_14);
  __Pyx_INCREF(((PyObject *)__pyx_kp_u_13));
  PyTuple_SET_ITEM(__pyx_k_tuple_14, 0, ((PyObject *)__pyx_kp_u_13));
  __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_13));
  __Pyx_GIVEREF(((PyObject *)__pyx_k_tuple_14));

  
  __pyx_k_tuple_15 = PyTuple_New(1); if (unlikely(!__pyx_k_tuple_15)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 803; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_k_tuple_15);
  __Pyx_INCREF(((PyObject *)__pyx_kp_u_10));
  PyTuple_SET_ITEM(__pyx_k_tuple_15, 0, ((PyObject *)__pyx_kp_u_10));
  __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_10));
  __Pyx_GIVEREF(((PyObject *)__pyx_k_tuple_15));

  
  __pyx_k_tuple_17 = PyTuple_New(1); if (unlikely(!__pyx_k_tuple_17)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 823; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_k_tuple_17);
  __Pyx_INCREF(((PyObject *)__pyx_kp_u_16));
  PyTuple_SET_ITEM(__pyx_k_tuple_17, 0, ((PyObject *)__pyx_kp_u_16));
  __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_16));
  __Pyx_GIVEREF(((PyObject *)__pyx_k_tuple_17));

  
  __pyx_k_tuple_192 = PyTuple_New(1); if (unlikely(!__pyx_k_tuple_192)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7228; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_k_tuple_192);
  __Pyx_INCREF(((PyObject *)__pyx_n_s__inflag));
  PyTuple_SET_ITEM(__pyx_k_tuple_192, 0, ((PyObject *)__pyx_n_s__inflag));
  __Pyx_GIVEREF(((PyObject *)__pyx_n_s__inflag));
  __Pyx_GIVEREF(((PyObject *)__pyx_k_tuple_192));
  __pyx_k_codeobj_193 = (PyObject*)__Pyx_PyCode_New(1, 0, 1, 0, 0, __pyx_empty_bytes, __pyx_empty_tuple, __pyx_empty_tuple, __pyx_k_tuple_192, __pyx_empty_tuple, __pyx_empty_tuple, __pyx_kp_s_194, __pyx_n_s___errprint, 7228, __pyx_empty_bytes); if (unlikely(!__pyx_k_codeobj_193)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7228; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_RefNannyFinishContext();
  return 0;
  __pyx_L1_error:;
  __Pyx_RefNannyFinishContext();
  return -1;
}

static int __Pyx_InitGlobals(void) {
  if (__Pyx_InitStrings(__pyx_string_tab) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  __pyx_int_15 = PyInt_FromLong(15); if (unlikely(!__pyx_int_15)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  return 0;
  __pyx_L1_error:;
  return -1;
}

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC init_ufuncs(void); 
PyMODINIT_FUNC init_ufuncs(void)
#else
PyMODINIT_FUNC PyInit__ufuncs(void); 
PyMODINIT_FUNC PyInit__ufuncs(void)
#endif
{
  PyObject *__pyx_t_1 = NULL;
  __Pyx_RefNannyDeclarations
  #if CYTHON_REFNANNY
  __Pyx_RefNanny = __Pyx_RefNannyImportAPI("refnanny");
  if (!__Pyx_RefNanny) {
      PyErr_Clear();
      __Pyx_RefNanny = __Pyx_RefNannyImportAPI("Cython.Runtime.refnanny");
      if (!__Pyx_RefNanny)
          Py_FatalError("failed to import 'refnanny' module");
  }
  #endif
  __Pyx_RefNannySetupContext("PyMODINIT_FUNC PyInit__ufuncs(void)", 0);
  if ( __Pyx_check_binary_version() < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_empty_tuple = PyTuple_New(0); if (unlikely(!__pyx_empty_tuple)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_empty_bytes = PyBytes_FromStringAndSize("", 0); if (unlikely(!__pyx_empty_bytes)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #ifdef __Pyx_CyFunction_USED
  if (__Pyx_CyFunction_init() < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #endif
  #ifdef __Pyx_FusedFunction_USED
  if (__pyx_FusedFunction_init() < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #endif
  #ifdef __Pyx_Generator_USED
  if (__pyx_Generator_init() < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #endif
  
  
  #if defined(__PYX_FORCE_INIT_THREADS) && __PYX_FORCE_INIT_THREADS
  #ifdef WITH_THREAD 
  PyEval_InitThreads();
  #endif
  #endif
  
  #if PY_MAJOR_VERSION < 3
  __pyx_m = Py_InitModule4(__Pyx_NAMESTR("_ufuncs"), __pyx_methods, 0, 0, PYTHON_API_VERSION); Py_XINCREF(__pyx_m);
  #else
  __pyx_m = PyModule_Create(&__pyx_moduledef);
  #endif
  if (unlikely(!__pyx_m)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #if PY_MAJOR_VERSION >= 3
  {
    PyObject *modules = PyImport_GetModuleDict(); if (unlikely(!modules)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    if (!PyDict_GetItemString(modules, "scipy.special._ufuncs")) {
      if (unlikely(PyDict_SetItemString(modules, "scipy.special._ufuncs", __pyx_m) < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    }
  }
  #endif
  __pyx_b = PyImport_AddModule(__Pyx_NAMESTR(__Pyx_BUILTIN_MODULE_NAME)); if (unlikely(!__pyx_b)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #if CYTHON_COMPILING_IN_PYPY
  Py_INCREF(__pyx_b);
  #endif
  if (__Pyx_SetAttrString(__pyx_m, "__builtins__", __pyx_b) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  
  if (unlikely(__Pyx_InitGlobals() < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__pyx_module_is_main_scipy__special___ufuncs) {
    if (__Pyx_SetAttrString(__pyx_m, "__name__", __pyx_n_s____main__) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  }
  
  if (unlikely(__Pyx_InitCachedBuiltins() < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  
  if (unlikely(__Pyx_InitCachedConstants() < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  
  
  
  
  
  __pyx_ptype_7cpython_4type_type = __Pyx_ImportType(__Pyx_BUILTIN_MODULE_NAME, "type", 
  #if CYTHON_COMPILING_IN_PYPY
  sizeof(PyTypeObject),
  #else
  sizeof(PyHeapTypeObject),
  #endif
  0); if (unlikely(!__pyx_ptype_7cpython_4type_type)) {__pyx_filename = __pyx_f[3]; __pyx_lineno = 9; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5numpy_dtype = __Pyx_ImportType("numpy", "dtype", sizeof(PyArray_Descr), 0); if (unlikely(!__pyx_ptype_5numpy_dtype)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 155; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5numpy_flatiter = __Pyx_ImportType("numpy", "flatiter", sizeof(PyArrayIterObject), 0); if (unlikely(!__pyx_ptype_5numpy_flatiter)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 165; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5numpy_broadcast = __Pyx_ImportType("numpy", "broadcast", sizeof(PyArrayMultiIterObject), 0); if (unlikely(!__pyx_ptype_5numpy_broadcast)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 169; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5numpy_ndarray = __Pyx_ImportType("numpy", "ndarray", sizeof(PyArrayObject), 0); if (unlikely(!__pyx_ptype_5numpy_ndarray)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 178; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5numpy_ufunc = __Pyx_ImportType("numpy", "ufunc", sizeof(PyUFuncObject), 0); if (unlikely(!__pyx_ptype_5numpy_ufunc)) {__pyx_filename = __pyx_f[2]; __pyx_lineno = 861; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  
  
  

  
  import_array();

  
  import_ufunc();

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_lambertw_scalar_t_var = (&__pyx_f_5scipy_7special_8lambertw_lambertw_scalar);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_bdtr_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_bdtr_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_bdtrc_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_bdtrc_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_bdtri_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_bdtri_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_binom_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_binom);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyc_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyc_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyc_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebys_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebys_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebys_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyt_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyt_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyt_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyu_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyu_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_chebyu_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_gegenbauer_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_genlaguerre_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_hermite_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_hermite);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_hermitenorm_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_hermitenorm);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_jacobi_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_jacobi_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_jacobi_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_laguerre_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_laguerre_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_laguerre_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_legendre_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_legendre_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_legendre_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyt_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_chebyu_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_jacobi_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_double__t_var = (&__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_double_complex__t_var = (&__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_eval_sh_legendre_l_t_var = (&__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre_l);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_expn_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_expn_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_hyp2f0_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_hyp2f0_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_kn_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_kn_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_nbdtr_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_nbdtr_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_nbdtrc_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_nbdtrc_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_nbdtri_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_nbdtri_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_pdtr_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_pdtr_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_pdtrc_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_pdtrc_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_pdtri_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_pdtri_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_smirnov_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_smirnov_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_smirnovi_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_smirnovi_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs__proto_yn_unsafe_t_var = (&__pyx_f_5scipy_7special_7_legacy_yn_unsafe);

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_doc = __pyx_k_18;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_Dld__As_Flf_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_Dld__As_Dld_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types[0]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types[1]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types[3]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types[4]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types[5]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types[7]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_ptr[0]) = ((void *)__pyx_f_5scipy_7special_8lambertw_lambertw_scalar);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_ptr[1]) = ((void *)((char *)__pyx_k___lambertw));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_ptr[2]) = ((void *)__pyx_f_5scipy_7special_8lambertw_lambertw_scalar);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_ptr[3]) = ((void *)((char *)__pyx_k___lambertw));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_types, 2, 3, 1, 0, __pyx_k___lambertw, __pyx_v_5scipy_7special_7_ufuncs_ufunc__lambertw_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1553; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s___lambertw, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1553; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_doc = __pyx_k_19;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_f_ffff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_d_dddd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DDDD_As_F_FFFF);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DDDD_As_D_DDDD);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[12]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[13]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[14]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[15]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[18]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types[19]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[0]) = ((void *)airy);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[1]) = ((void *)((char *)__pyx_k__airy));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[2]) = ((void *)airy);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[3]) = ((void *)((char *)__pyx_k__airy));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[4]) = ((void *)cairy_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[5]) = ((void *)((char *)__pyx_k__airy));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[6]) = ((void *)cairy_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[7]) = ((void *)((char *)__pyx_k__airy));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_types, 4, 1, 4, 0, __pyx_k__airy, __pyx_v_5scipy_7special_7_ufuncs_ufunc_airy_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1600; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__airy, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1600; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_doc = __pyx_k_20;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_f_ffff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_d_dddd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DDDD_As_F_FFFF);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DDDD_As_D_DDDD);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[12]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[13]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[14]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[15]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[18]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types[19]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[0]) = ((void *)cairy_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[1]) = ((void *)((char *)__pyx_k__airye));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[2]) = ((void *)cairy_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[3]) = ((void *)((char *)__pyx_k__airye));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[4]) = ((void *)cairy_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[5]) = ((void *)((char *)__pyx_k__airye));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[6]) = ((void *)cairy_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[7]) = ((void *)((char *)__pyx_k__airye));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_types, 4, 1, 4, 0, __pyx_k__airye, __pyx_v_5scipy_7special_7_ufuncs_ufunc_airye_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1647; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__airye, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1647; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_doc = __pyx_k_21;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_llf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_lld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[1]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[4]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[5]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[0]) = ((void *)bdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[1]) = ((void *)((char *)__pyx_k__bdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[2]) = ((void *)bdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[3]) = ((void *)((char *)__pyx_k__bdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_bdtr_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[5]) = ((void *)((char *)__pyx_k__bdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_bdtr_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[7]) = ((void *)((char *)__pyx_k__bdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_types, 4, 3, 1, 0, __pyx_k__bdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1688; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__bdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1688; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_doc = __pyx_k_22;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_llf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_lld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[1]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[4]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[5]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[0]) = ((void *)bdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[1]) = ((void *)((char *)__pyx_k__bdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[2]) = ((void *)bdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[3]) = ((void *)((char *)__pyx_k__bdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_bdtrc_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[5]) = ((void *)((char *)__pyx_k__bdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_bdtrc_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[7]) = ((void *)((char *)__pyx_k__bdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_types, 4, 3, 1, 0, __pyx_k__bdtrc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1729; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__bdtrc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1729; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_doc = __pyx_k_23;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_llf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_lld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[1]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[4]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[5]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[0]) = ((void *)bdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[1]) = ((void *)((char *)__pyx_k__bdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[2]) = ((void *)bdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[3]) = ((void *)((char *)__pyx_k__bdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_bdtri_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[5]) = ((void *)((char *)__pyx_k__bdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_bdtri_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[7]) = ((void *)((char *)__pyx_k__bdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_types, 4, 3, 1, 0, __pyx_k__bdtri, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtri_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1771; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__bdtri, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1771; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_ptr[0]) = ((void *)cdfbin2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_ptr[1]) = ((void *)((char *)__pyx_k__bdtrik));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_ptr[2]) = ((void *)cdfbin2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_ptr[3]) = ((void *)((char *)__pyx_k__bdtrik));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_types, 2, 3, 1, 0, __pyx_k__bdtrik, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrik_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1795; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__bdtrik, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1795; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_ptr[0]) = ((void *)cdfbin3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_ptr[1]) = ((void *)((char *)__pyx_k__bdtrin));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_ptr[2]) = ((void *)cdfbin3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_ptr[3]) = ((void *)((char *)__pyx_k__bdtrin));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_types, 2, 3, 1, 0, __pyx_k__bdtrin, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bdtrin_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1819; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__bdtrin, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1819; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_doc = __pyx_k_25;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_ptr[0]) = ((void *)bei_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_ptr[1]) = ((void *)((char *)__pyx_k__bei));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_ptr[2]) = ((void *)bei_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_ptr[3]) = ((void *)((char *)__pyx_k__bei));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_types, 2, 1, 1, 0, __pyx_k__bei, __pyx_v_5scipy_7special_7_ufuncs_ufunc_bei_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1839; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__bei, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1839; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_doc = __pyx_k_26;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_ptr[0]) = ((void *)beip_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_ptr[1]) = ((void *)((char *)__pyx_k__beip));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_ptr[2]) = ((void *)beip_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_ptr[3]) = ((void *)((char *)__pyx_k__beip));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_types, 2, 1, 1, 0, __pyx_k__beip, __pyx_v_5scipy_7special_7_ufuncs_ufunc_beip_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1859; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__beip, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1859; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_doc = __pyx_k_27;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_ptr[0]) = ((void *)ber_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_ptr[1]) = ((void *)((char *)__pyx_k__ber));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_ptr[2]) = ((void *)ber_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_ptr[3]) = ((void *)((char *)__pyx_k__ber));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_types, 2, 1, 1, 0, __pyx_k__ber, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ber_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1879; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ber, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1879; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_doc = __pyx_k_28;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_ptr[0]) = ((void *)berp_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_ptr[1]) = ((void *)((char *)__pyx_k__berp));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_ptr[2]) = ((void *)berp_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_ptr[3]) = ((void *)((char *)__pyx_k__berp));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_types, 2, 1, 1, 0, __pyx_k__berp, __pyx_v_5scipy_7special_7_ufuncs_ufunc_berp_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1899; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__berp, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1899; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_doc = __pyx_k_29;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_ptr[0]) = ((void *)besselpoly);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_ptr[1]) = ((void *)((char *)__pyx_k__besselpoly));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_ptr[2]) = ((void *)besselpoly);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_ptr[3]) = ((void *)((char *)__pyx_k__besselpoly));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_types, 2, 3, 1, 0, __pyx_k__besselpoly, __pyx_v_5scipy_7special_7_ufuncs_ufunc_besselpoly_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1924; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__besselpoly, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1924; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_doc = __pyx_k_30;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_ptr[0]) = ((void *)beta);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_ptr[1]) = ((void *)((char *)__pyx_k__beta));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_ptr[2]) = ((void *)beta);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_ptr[3]) = ((void *)((char *)__pyx_k__beta));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_types, 2, 2, 1, 0, __pyx_k__beta, __pyx_v_5scipy_7special_7_ufuncs_ufunc_beta_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1946; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__beta, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1946; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_doc = __pyx_k_31;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_ptr[0]) = ((void *)incbet);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_ptr[1]) = ((void *)((char *)__pyx_k__betainc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_ptr[2]) = ((void *)incbet);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_ptr[3]) = ((void *)((char *)__pyx_k__betainc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_types, 2, 3, 1, 0, __pyx_k__betainc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_betainc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1982; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__betainc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1982; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_doc = __pyx_k_32;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_ptr[0]) = ((void *)incbi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_ptr[1]) = ((void *)((char *)__pyx_k__betaincinv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_ptr[2]) = ((void *)incbi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_ptr[3]) = ((void *)((char *)__pyx_k__betaincinv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_types, 2, 3, 1, 0, __pyx_k__betaincinv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaincinv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2008; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__betaincinv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2008; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_doc = __pyx_k_33;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_ptr[0]) = ((void *)lbeta);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_ptr[1]) = ((void *)((char *)__pyx_k__betaln));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_ptr[2]) = ((void *)lbeta);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_ptr[3]) = ((void *)((char *)__pyx_k__betaln));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_types, 2, 2, 1, 0, __pyx_k__betaln, __pyx_v_5scipy_7special_7_ufuncs_ufunc_betaln_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2031; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__betaln, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2031; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_doc = __pyx_k_34;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_binom);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_ptr[1]) = ((void *)((char *)__pyx_k__binom));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_binom);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_ptr[3]) = ((void *)((char *)__pyx_k__binom));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_types, 2, 2, 1, 0, __pyx_k__binom, __pyx_v_5scipy_7special_7_ufuncs_ufunc_binom_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2055; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__binom, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2055; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_doc = __pyx_k_35;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_ptr[0]) = ((void *)btdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_ptr[1]) = ((void *)((char *)__pyx_k__btdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_ptr[2]) = ((void *)btdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_ptr[3]) = ((void *)((char *)__pyx_k__btdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_types, 2, 3, 1, 0, __pyx_k__btdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2081; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__btdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2081; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_doc = __pyx_k_36;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_ptr[0]) = ((void *)incbi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_ptr[1]) = ((void *)((char *)__pyx_k__btdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_ptr[2]) = ((void *)incbi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_ptr[3]) = ((void *)((char *)__pyx_k__btdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_types, 2, 3, 1, 0, __pyx_k__btdtri, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtri_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2107; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__btdtri, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2107; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_ptr[0]) = ((void *)cdfbet3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_ptr[1]) = ((void *)((char *)__pyx_k__btdtria));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_ptr[2]) = ((void *)cdfbet3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_ptr[3]) = ((void *)((char *)__pyx_k__btdtria));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_types, 2, 3, 1, 0, __pyx_k__btdtria, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtria_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2131; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__btdtria, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2131; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_ptr[0]) = ((void *)cdfbet4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_ptr[1]) = ((void *)((char *)__pyx_k__btdtrib));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_ptr[2]) = ((void *)cdfbet4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_ptr[3]) = ((void *)((char *)__pyx_k__btdtrib));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_types, 2, 3, 1, 0, __pyx_k__btdtrib, __pyx_v_5scipy_7special_7_ufuncs_ufunc_btdtrib_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2155; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__btdtrib, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2155; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_doc = __pyx_k_37;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_ptr[0]) = ((void *)cbrt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_ptr[1]) = ((void *)((char *)__pyx_k__cbrt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_ptr[2]) = ((void *)cbrt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_ptr[3]) = ((void *)((char *)__pyx_k__cbrt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_types, 2, 1, 1, 0, __pyx_k__cbrt, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cbrt_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2175; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__cbrt, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2175; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_doc = __pyx_k_38;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_ptr[0]) = ((void *)chdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_ptr[1]) = ((void *)((char *)__pyx_k__chdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_ptr[2]) = ((void *)chdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_ptr[3]) = ((void *)((char *)__pyx_k__chdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_types, 2, 2, 1, 0, __pyx_k__chdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2199; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__chdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2199; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_doc = __pyx_k_39;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_ptr[0]) = ((void *)chdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_ptr[1]) = ((void *)((char *)__pyx_k__chdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_ptr[2]) = ((void *)chdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_ptr[3]) = ((void *)((char *)__pyx_k__chdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_types, 2, 2, 1, 0, __pyx_k__chdtrc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtrc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2224; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__chdtrc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2224; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_doc = __pyx_k_40;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_ptr[0]) = ((void *)chdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_ptr[1]) = ((void *)((char *)__pyx_k__chdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_ptr[2]) = ((void *)chdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_ptr[3]) = ((void *)((char *)__pyx_k__chdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_types, 2, 2, 1, 0, __pyx_k__chdtri, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtri_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2247; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__chdtri, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2247; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_ptr[0]) = ((void *)cdfchi3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_ptr[1]) = ((void *)((char *)__pyx_k__chdtriv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_ptr[2]) = ((void *)cdfchi3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_ptr[3]) = ((void *)((char *)__pyx_k__chdtriv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_types, 2, 2, 1, 0, __pyx_k__chdtriv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chdtriv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2269; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__chdtriv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2269; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_ptr[0]) = ((void *)cdfchn1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_ptr[1]) = ((void *)((char *)__pyx_k__chndtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_ptr[2]) = ((void *)cdfchn1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_ptr[3]) = ((void *)((char *)__pyx_k__chndtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_types, 2, 3, 1, 0, __pyx_k__chndtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2293; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__chndtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2293; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_ptr[0]) = ((void *)cdfchn3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_ptr[1]) = ((void *)((char *)__pyx_k__chndtridf));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_ptr[2]) = ((void *)cdfchn3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_ptr[3]) = ((void *)((char *)__pyx_k__chndtridf));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_types, 2, 3, 1, 0, __pyx_k__chndtridf, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtridf_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2317; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__chndtridf, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2317; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_ptr[0]) = ((void *)cdfchn4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_ptr[1]) = ((void *)((char *)__pyx_k__chndtrinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_ptr[2]) = ((void *)cdfchn4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_ptr[3]) = ((void *)((char *)__pyx_k__chndtrinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_types, 2, 3, 1, 0, __pyx_k__chndtrinc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrinc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2341; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__chndtrinc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2341; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_ptr[0]) = ((void *)cdfchn2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_ptr[1]) = ((void *)((char *)__pyx_k__chndtrix));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_ptr[2]) = ((void *)cdfchn2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_ptr[3]) = ((void *)((char *)__pyx_k__chndtrix));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_types, 2, 3, 1, 0, __pyx_k__chndtrix, __pyx_v_5scipy_7special_7_ufuncs_ufunc_chndtrix_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2365; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__chndtrix, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2365; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_doc = __pyx_k_41;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_ptr[0]) = ((void *)cosdg);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_ptr[1]) = ((void *)((char *)__pyx_k__cosdg));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_ptr[2]) = ((void *)cosdg);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_ptr[3]) = ((void *)((char *)__pyx_k__cosdg));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_types, 2, 1, 1, 0, __pyx_k__cosdg, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosdg_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2385; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__cosdg, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2385; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_doc = __pyx_k_42;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_ptr[0]) = ((void *)cosm1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_ptr[1]) = ((void *)((char *)__pyx_k__cosm1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_ptr[2]) = ((void *)cosm1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_ptr[3]) = ((void *)((char *)__pyx_k__cosm1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_types, 2, 1, 1, 0, __pyx_k__cosm1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cosm1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2405; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__cosm1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2405; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_doc = __pyx_k_43;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_ptr[0]) = ((void *)cotdg);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_ptr[1]) = ((void *)((char *)__pyx_k__cotdg));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_ptr[2]) = ((void *)cotdg);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_ptr[3]) = ((void *)((char *)__pyx_k__cotdg));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_types, 2, 1, 1, 0, __pyx_k__cotdg, __pyx_v_5scipy_7special_7_ufuncs_ufunc_cotdg_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2425; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__cotdg, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2425; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_doc = __pyx_k_44;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_ptr[0]) = ((void *)ellpe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_ptr[1]) = ((void *)((char *)__pyx_k__ellipe));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_ptr[2]) = ((void *)ellpe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_ptr[3]) = ((void *)((char *)__pyx_k__ellipe));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_types, 2, 1, 1, 0, __pyx_k__ellipe, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipe_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2446; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ellipe, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2446; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_doc = __pyx_k_45;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_ptr[0]) = ((void *)ellie);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_ptr[1]) = ((void *)((char *)__pyx_k__ellipeinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_ptr[2]) = ((void *)ellie);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_ptr[3]) = ((void *)((char *)__pyx_k__ellipeinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_types, 2, 2, 1, 0, __pyx_k__ellipeinc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipeinc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2469; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ellipeinc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2469; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_doc = __pyx_k_46;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dddd_As_ff_ffff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dddd_As_dd_dddd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_ptr[0]) = ((void *)ellpj);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_ptr[1]) = ((void *)((char *)__pyx_k__ellipj));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_ptr[2]) = ((void *)ellpj);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_ptr[3]) = ((void *)((char *)__pyx_k__ellipj));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_types, 2, 2, 4, 0, __pyx_k__ellipj, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipj_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2500; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ellipj, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2500; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_doc = __pyx_k_47;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_ptr[0]) = ((void *)ellik);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_ptr[1]) = ((void *)((char *)__pyx_k__ellipkinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_ptr[2]) = ((void *)ellik);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_ptr[3]) = ((void *)((char *)__pyx_k__ellipkinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_types, 2, 2, 1, 0, __pyx_k__ellipkinc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkinc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2523; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ellipkinc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2523; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_doc = __pyx_k_48;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_ptr[0]) = ((void *)ellpk);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_ptr[1]) = ((void *)((char *)__pyx_k__ellipkm1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_ptr[2]) = ((void *)ellpk);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_ptr[3]) = ((void *)((char *)__pyx_k__ellipkm1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_types, 2, 1, 1, 0, __pyx_k__ellipkm1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ellipkm1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2544; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ellipkm1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2544; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_doc = __pyx_k_49;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[1]) = ((void *)((char *)__pyx_k__eval_chebyc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[3]) = ((void *)((char *)__pyx_k__eval_chebyc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[5]) = ((void *)((char *)__pyx_k__eval_chebyc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[7]) = ((void *)((char *)__pyx_k__eval_chebyc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[9]) = ((void *)((char *)__pyx_k__eval_chebyc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[11]) = ((void *)((char *)__pyx_k__eval_chebyc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_types, 6, 2, 1, 0, __pyx_k__eval_chebyc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2596; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_chebyc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2596; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_doc = __pyx_k_50;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[1]) = ((void *)((char *)__pyx_k__eval_chebys));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[3]) = ((void *)((char *)__pyx_k__eval_chebys));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[5]) = ((void *)((char *)__pyx_k__eval_chebys));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[7]) = ((void *)((char *)__pyx_k__eval_chebys));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[9]) = ((void *)((char *)__pyx_k__eval_chebys));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebys);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[11]) = ((void *)((char *)__pyx_k__eval_chebys));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_types, 6, 2, 1, 0, __pyx_k__eval_chebys, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebys_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2648; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_chebys, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2648; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_doc = __pyx_k_51;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[1]) = ((void *)((char *)__pyx_k__eval_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[3]) = ((void *)((char *)__pyx_k__eval_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[5]) = ((void *)((char *)__pyx_k__eval_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[7]) = ((void *)((char *)__pyx_k__eval_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[9]) = ((void *)((char *)__pyx_k__eval_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[11]) = ((void *)((char *)__pyx_k__eval_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_types, 6, 2, 1, 0, __pyx_k__eval_chebyt, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyt_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2703; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_chebyt, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2703; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_doc = __pyx_k_52;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[1]) = ((void *)((char *)__pyx_k__eval_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[3]) = ((void *)((char *)__pyx_k__eval_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[5]) = ((void *)((char *)__pyx_k__eval_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[7]) = ((void *)((char *)__pyx_k__eval_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[9]) = ((void *)((char *)__pyx_k__eval_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_chebyu);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[11]) = ((void *)((char *)__pyx_k__eval_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_types, 6, 2, 1, 0, __pyx_k__eval_chebyu, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_chebyu_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2755; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_chebyu, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2755; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_doc = __pyx_k_53;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ldd__As_lff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ldd__As_ldd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ffF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ddD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[4]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[12]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[13]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[14]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[15]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[16]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[17]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[18]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[19]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[20]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[21]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[22]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types[23]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[1]) = ((void *)((char *)__pyx_k__eval_gegenbauer));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[3]) = ((void *)((char *)__pyx_k__eval_gegenbauer));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[5]) = ((void *)((char *)__pyx_k__eval_gegenbauer));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[7]) = ((void *)((char *)__pyx_k__eval_gegenbauer));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[9]) = ((void *)((char *)__pyx_k__eval_gegenbauer));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_gegenbauer);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[11]) = ((void *)((char *)__pyx_k__eval_gegenbauer));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_types, 6, 3, 1, 0, __pyx_k__eval_gegenbauer, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_gegenbauer_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2813; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_gegenbauer, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2813; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_doc = __pyx_k_54;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ldd__As_lff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ldd__As_ldd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ffF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ddD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[4]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[12]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[13]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[14]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[15]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[16]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[17]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[18]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[19]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[20]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[21]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[22]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types[23]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[1]) = ((void *)((char *)__pyx_k__eval_genlaguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[3]) = ((void *)((char *)__pyx_k__eval_genlaguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[5]) = ((void *)((char *)__pyx_k__eval_genlaguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[7]) = ((void *)((char *)__pyx_k__eval_genlaguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[9]) = ((void *)((char *)__pyx_k__eval_genlaguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_genlaguerre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[11]) = ((void *)((char *)__pyx_k__eval_genlaguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_types, 6, 3, 1, 0, __pyx_k__eval_genlaguerre, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_genlaguerre_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2871; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_genlaguerre, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2871; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_doc = __pyx_k_55;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_hermite);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_ptr[1]) = ((void *)((char *)__pyx_k__eval_hermite));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_hermite);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_ptr[3]) = ((void *)((char *)__pyx_k__eval_hermite));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_types, 2, 2, 1, 0, __pyx_k__eval_hermite, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermite_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2895; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_hermite, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2895; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_doc = __pyx_k_56;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_hermitenorm);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_ptr[1]) = ((void *)((char *)__pyx_k__eval_hermitenorm));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_hermitenorm);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_ptr[3]) = ((void *)((char *)__pyx_k__eval_hermitenorm));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_types, 2, 2, 1, 0, __pyx_k__eval_hermitenorm, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_hermitenorm_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2919; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_hermitenorm, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2919; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_doc = __pyx_k_57;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_lddd__As_lfff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_lddd__As_lddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_fffF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_dddD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[5]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[12]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[13]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[14]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[15]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[16]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[17]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[18]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[19]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[20]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[21]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[22]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[23]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[24]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[25]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[26]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[27]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[28]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types[29]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[1]) = ((void *)((char *)__pyx_k__eval_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[3]) = ((void *)((char *)__pyx_k__eval_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[5]) = ((void *)((char *)__pyx_k__eval_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[7]) = ((void *)((char *)__pyx_k__eval_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[9]) = ((void *)((char *)__pyx_k__eval_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_jacobi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[11]) = ((void *)((char *)__pyx_k__eval_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_types, 6, 4, 1, 0, __pyx_k__eval_jacobi, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_jacobi_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2983; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_jacobi, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 2983; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_doc = __pyx_k_58;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[1]) = ((void *)((char *)__pyx_k__eval_laguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[3]) = ((void *)((char *)__pyx_k__eval_laguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[5]) = ((void *)((char *)__pyx_k__eval_laguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[7]) = ((void *)((char *)__pyx_k__eval_laguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[9]) = ((void *)((char *)__pyx_k__eval_laguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_laguerre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[11]) = ((void *)((char *)__pyx_k__eval_laguerre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_types, 6, 2, 1, 0, __pyx_k__eval_laguerre, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_laguerre_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3035; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_laguerre, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3035; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_doc = __pyx_k_59;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[1]) = ((void *)((char *)__pyx_k__eval_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[3]) = ((void *)((char *)__pyx_k__eval_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[5]) = ((void *)((char *)__pyx_k__eval_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[7]) = ((void *)((char *)__pyx_k__eval_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[9]) = ((void *)((char *)__pyx_k__eval_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_legendre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[11]) = ((void *)((char *)__pyx_k__eval_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_types, 6, 2, 1, 0, __pyx_k__eval_legendre, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_legendre_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3087; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_legendre, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3087; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_doc = __pyx_k_60;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[1]) = ((void *)((char *)__pyx_k__eval_sh_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[3]) = ((void *)((char *)__pyx_k__eval_sh_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[5]) = ((void *)((char *)__pyx_k__eval_sh_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[7]) = ((void *)((char *)__pyx_k__eval_sh_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[9]) = ((void *)((char *)__pyx_k__eval_sh_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyt);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[11]) = ((void *)((char *)__pyx_k__eval_sh_chebyt));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_types, 6, 2, 1, 0, __pyx_k__eval_sh_chebyt, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyt_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3139; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_sh_chebyt, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3139; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_doc = __pyx_k_61;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[1]) = ((void *)((char *)__pyx_k__eval_sh_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[3]) = ((void *)((char *)__pyx_k__eval_sh_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[5]) = ((void *)((char *)__pyx_k__eval_sh_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[7]) = ((void *)((char *)__pyx_k__eval_sh_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[9]) = ((void *)((char *)__pyx_k__eval_sh_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_chebyu);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[11]) = ((void *)((char *)__pyx_k__eval_sh_chebyu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_types, 6, 2, 1, 0, __pyx_k__eval_sh_chebyu, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_chebyu_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3191; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_sh_chebyu, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3191; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_doc = __pyx_k_62;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_lddd__As_lfff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_lddd__As_lddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_fffF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_dddD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[5]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[12]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[13]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[14]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[15]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[16]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[17]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[18]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[19]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[20]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[21]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[22]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[23]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[24]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[25]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[26]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[27]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[28]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types[29]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[1]) = ((void *)((char *)__pyx_k__eval_sh_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[3]) = ((void *)((char *)__pyx_k__eval_sh_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[5]) = ((void *)((char *)__pyx_k__eval_sh_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[7]) = ((void *)((char *)__pyx_k__eval_sh_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[9]) = ((void *)((char *)__pyx_k__eval_sh_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_jacobi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[11]) = ((void *)((char *)__pyx_k__eval_sh_jacobi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_types, 6, 4, 1, 0, __pyx_k__eval_sh_jacobi, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_jacobi_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3255; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_sh_jacobi, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3255; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_doc = __pyx_k_63;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ld__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_loops[4]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_loops[5]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[10]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[11]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[16]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types[17]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[0]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[1]) = ((void *)((char *)__pyx_k__eval_sh_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[2]) = ((void *)__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre_l);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[3]) = ((void *)((char *)__pyx_k__eval_sh_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[4]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[5]) = ((void *)((char *)__pyx_k__eval_sh_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[6]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[7]) = ((void *)((char *)__pyx_k__eval_sh_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[8]) = ((void *)__pyx_fuse_0__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[9]) = ((void *)((char *)__pyx_k__eval_sh_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[10]) = ((void *)__pyx_fuse_1__pyx_f_5scipy_7special_15orthogonal_eval_eval_sh_legendre);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[11]) = ((void *)((char *)__pyx_k__eval_sh_legendre));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[6]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_data[4]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[8]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_data[5]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_ptr[10]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_types, 6, 2, 1, 0, __pyx_k__eval_sh_legendre, __pyx_v_5scipy_7special_7_ufuncs_ufunc_eval_sh_legendre_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3307; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__eval_sh_legendre, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3307; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_doc = __pyx_k_64;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_F_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_D_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types[6]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types[7]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[0]) = ((void *)exp1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[1]) = ((void *)((char *)__pyx_k__exp1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[2]) = ((void *)exp1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[3]) = ((void *)((char *)__pyx_k__exp1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[4]) = ((void *)cexp1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[5]) = ((void *)((char *)__pyx_k__exp1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[6]) = ((void *)cexp1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[7]) = ((void *)((char *)__pyx_k__exp1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_types, 4, 1, 1, 0, __pyx_k__exp1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3340; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__exp1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3340; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_doc = __pyx_k_65;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_ptr[0]) = ((void *)exp10);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_ptr[1]) = ((void *)((char *)__pyx_k__exp10));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_ptr[2]) = ((void *)exp10);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_ptr[3]) = ((void *)((char *)__pyx_k__exp10));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_types, 2, 1, 1, 0, __pyx_k__exp10, __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp10_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3360; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__exp10, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3360; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_doc = __pyx_k_66;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_ptr[0]) = ((void *)exp2);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_ptr[1]) = ((void *)((char *)__pyx_k__exp2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_ptr[2]) = ((void *)exp2);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_ptr[3]) = ((void *)((char *)__pyx_k__exp2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_types, 2, 1, 1, 0, __pyx_k__exp2, __pyx_v_5scipy_7special_7_ufuncs_ufunc_exp2_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3380; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__exp2, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3380; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_doc = __pyx_k_67;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_F_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_D_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types[6]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types[7]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[0]) = ((void *)expi_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[1]) = ((void *)((char *)__pyx_k__expi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[2]) = ((void *)expi_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[3]) = ((void *)((char *)__pyx_k__expi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[4]) = ((void *)cexpi_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[5]) = ((void *)((char *)__pyx_k__expi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[6]) = ((void *)cexpi_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[7]) = ((void *)((char *)__pyx_k__expi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_types, 4, 1, 1, 0, __pyx_k__expi, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expi_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3414; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__expi, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3414; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_doc = __pyx_k_68;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_f_f__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_g_g__As_g_g);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_types[4]) = ((char)NPY_LONGDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_types[5]) = ((char)NPY_LONGDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[0]) = ((void *)expitf);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[1]) = ((void *)((char *)__pyx_k__expit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[2]) = ((void *)expit);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[3]) = ((void *)((char *)__pyx_k__expit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[4]) = ((void *)expitl);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[5]) = ((void *)((char *)__pyx_k__expit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_ptr[4]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_types, 3, 1, 1, 0, __pyx_k__expit, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expit_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3460; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__expit, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3460; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_doc = __pyx_k_69;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_ptr[0]) = ((void *)expm1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_ptr[1]) = ((void *)((char *)__pyx_k__expm1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_ptr[2]) = ((void *)expm1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_ptr[3]) = ((void *)((char *)__pyx_k__expm1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_types, 2, 1, 1, 0, __pyx_k__expm1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expm1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3480; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__expm1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3480; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_doc = __pyx_k_70;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[0]) = ((void *)expn);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[1]) = ((void *)((char *)__pyx_k__expn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[2]) = ((void *)expn);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[3]) = ((void *)((char *)__pyx_k__expn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_expn_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[5]) = ((void *)((char *)__pyx_k__expn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_expn_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[7]) = ((void *)((char *)__pyx_k__expn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_types, 4, 2, 1, 0, __pyx_k__expn, __pyx_v_5scipy_7special_7_ufuncs_ufunc_expn_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3517; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__expn, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3517; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_doc = __pyx_k_71;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_ptr[0]) = ((void *)fdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_ptr[1]) = ((void *)((char *)__pyx_k__fdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_ptr[2]) = ((void *)fdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_ptr[3]) = ((void *)((char *)__pyx_k__fdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_types, 2, 3, 1, 0, __pyx_k__fdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3545; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__fdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3545; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_doc = __pyx_k_72;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_ptr[0]) = ((void *)fdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_ptr[1]) = ((void *)((char *)__pyx_k__fdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_ptr[2]) = ((void *)fdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_ptr[3]) = ((void *)((char *)__pyx_k__fdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_types, 2, 3, 1, 0, __pyx_k__fdtrc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtrc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3569; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__fdtrc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3569; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_doc = __pyx_k_73;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_ptr[0]) = ((void *)fdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_ptr[1]) = ((void *)((char *)__pyx_k__fdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_ptr[2]) = ((void *)fdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_ptr[3]) = ((void *)((char *)__pyx_k__fdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_types, 2, 3, 1, 0, __pyx_k__fdtri, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtri_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3594; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__fdtri, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3594; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_doc = __pyx_k_74;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_ptr[0]) = ((void *)cdff4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_ptr[1]) = ((void *)((char *)__pyx_k__fdtridfd));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_ptr[2]) = ((void *)cdff4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_ptr[3]) = ((void *)((char *)__pyx_k__fdtridfd));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_types, 2, 3, 1, 0, __pyx_k__fdtridfd, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fdtridfd_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3619; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__fdtridfd, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3619; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_doc = __pyx_k_75;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_f_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_d_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DD_As_F_FF);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_D_DD_As_D_DD);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[6]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[7]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[8]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[9]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[10]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types[11]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[0]) = ((void *)fresnl);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[1]) = ((void *)((char *)__pyx_k__fresnel));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[2]) = ((void *)fresnl);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[3]) = ((void *)((char *)__pyx_k__fresnel));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[4]) = ((void *)cfresnl_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[5]) = ((void *)((char *)__pyx_k__fresnel));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[6]) = ((void *)cfresnl_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[7]) = ((void *)((char *)__pyx_k__fresnel));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_types, 4, 1, 2, 0, __pyx_k__fresnel, __pyx_v_5scipy_7special_7_ufuncs_ufunc_fresnel_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3657; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__fresnel, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3657; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_doc = __pyx_k_76;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_F_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_D_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types[6]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types[7]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[0]) = ((void *)Gamma);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[1]) = ((void *)((char *)__pyx_k__gamma));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[2]) = ((void *)Gamma);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[3]) = ((void *)((char *)__pyx_k__gamma));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[4]) = ((void *)cgamma_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[5]) = ((void *)((char *)__pyx_k__gamma));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[6]) = ((void *)cgamma_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[7]) = ((void *)((char *)__pyx_k__gamma));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_types, 4, 1, 1, 0, __pyx_k__gamma, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gamma_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3691; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gamma, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3691; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_doc = __pyx_k_77;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_ptr[0]) = ((void *)igam);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_ptr[1]) = ((void *)((char *)__pyx_k__gammainc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_ptr[2]) = ((void *)igam);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_ptr[3]) = ((void *)((char *)__pyx_k__gammainc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_types, 2, 2, 1, 0, __pyx_k__gammainc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3715; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gammainc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3715; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_doc = __pyx_k_78;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_ptr[0]) = ((void *)igamc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_ptr[1]) = ((void *)((char *)__pyx_k__gammaincc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_ptr[2]) = ((void *)igamc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_ptr[3]) = ((void *)((char *)__pyx_k__gammaincc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_types, 2, 2, 1, 0, __pyx_k__gammaincc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3739; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gammaincc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3739; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_doc = __pyx_k_79;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_ptr[0]) = ((void *)igami);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_ptr[1]) = ((void *)((char *)__pyx_k__gammainccinv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_ptr[2]) = ((void *)igami);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_ptr[3]) = ((void *)((char *)__pyx_k__gammainccinv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_types, 2, 2, 1, 0, __pyx_k__gammainccinv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammainccinv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3761; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gammainccinv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3761; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_doc = __pyx_k_80;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_ptr[0]) = ((void *)gammaincinv);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_ptr[1]) = ((void *)((char *)__pyx_k__gammaincinv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_ptr[2]) = ((void *)gammaincinv);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_ptr[3]) = ((void *)((char *)__pyx_k__gammaincinv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_types, 2, 2, 1, 0, __pyx_k__gammaincinv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaincinv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3783; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gammaincinv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3783; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_doc = __pyx_k_81;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_F_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_D_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types[6]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types[7]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[0]) = ((void *)lgam);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[1]) = ((void *)((char *)__pyx_k__gammaln));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[2]) = ((void *)lgam);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[3]) = ((void *)((char *)__pyx_k__gammaln));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[4]) = ((void *)clngamma_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[5]) = ((void *)((char *)__pyx_k__gammaln));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[6]) = ((void *)clngamma_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[7]) = ((void *)((char *)__pyx_k__gammaln));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_types, 4, 1, 1, 0, __pyx_k__gammaln, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammaln_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3820; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gammaln, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3820; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_doc = __pyx_k_82;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_ptr[0]) = ((void *)gammasgn);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_ptr[1]) = ((void *)((char *)__pyx_k__gammasgn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_ptr[2]) = ((void *)gammasgn);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_ptr[3]) = ((void *)((char *)__pyx_k__gammasgn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_types, 2, 1, 1, 0, __pyx_k__gammasgn, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gammasgn_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3844; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gammasgn, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3844; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_doc = __pyx_k_83;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_ptr[0]) = ((void *)gdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_ptr[1]) = ((void *)((char *)__pyx_k__gdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_ptr[2]) = ((void *)gdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_ptr[3]) = ((void *)((char *)__pyx_k__gdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_types, 2, 3, 1, 0, __pyx_k__gdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3870; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3870; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_doc = __pyx_k_84;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_ptr[0]) = ((void *)gdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_ptr[1]) = ((void *)((char *)__pyx_k__gdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_ptr[2]) = ((void *)gdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_ptr[3]) = ((void *)((char *)__pyx_k__gdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_types, 2, 3, 1, 0, __pyx_k__gdtrc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3895; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gdtrc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3895; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_ptr[0]) = ((void *)cdfgam4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_ptr[1]) = ((void *)((char *)__pyx_k__gdtria));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_ptr[2]) = ((void *)cdfgam4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_ptr[3]) = ((void *)((char *)__pyx_k__gdtria));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_types, 2, 3, 1, 0, __pyx_k__gdtria, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtria_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3919; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gdtria, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3919; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_ptr[0]) = ((void *)cdfgam3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_ptr[1]) = ((void *)((char *)__pyx_k__gdtrib));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_ptr[2]) = ((void *)cdfgam3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_ptr[3]) = ((void *)((char *)__pyx_k__gdtrib));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_types, 2, 3, 1, 0, __pyx_k__gdtrib, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrib_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3943; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gdtrib, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3943; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_ptr[0]) = ((void *)cdfgam2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_ptr[1]) = ((void *)((char *)__pyx_k__gdtrix));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_ptr[2]) = ((void *)cdfgam2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_ptr[3]) = ((void *)((char *)__pyx_k__gdtrix));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_types, 2, 3, 1, 0, __pyx_k__gdtrix, __pyx_v_5scipy_7special_7_ufuncs_ufunc_gdtrix_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3967; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__gdtrix, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3967; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_doc = __pyx_k_85;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_types[1]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_types[2]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_types[4]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_types[5]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_ptr[0]) = ((void *)cbesh_wrap1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_ptr[1]) = ((void *)((char *)__pyx_k__hankel1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_ptr[2]) = ((void *)cbesh_wrap1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_ptr[3]) = ((void *)((char *)__pyx_k__hankel1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_types, 2, 2, 1, 0, __pyx_k__hankel1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3989; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hankel1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 3989; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_doc = __pyx_k_86;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_types[1]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_types[2]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_types[4]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_types[5]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_ptr[0]) = ((void *)cbesh_wrap1_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_ptr[1]) = ((void *)((char *)__pyx_k__hankel1e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_ptr[2]) = ((void *)cbesh_wrap1_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_ptr[3]) = ((void *)((char *)__pyx_k__hankel1e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_types, 2, 2, 1, 0, __pyx_k__hankel1e, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel1e_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4013; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hankel1e, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4013; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_doc = __pyx_k_87;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_types[1]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_types[2]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_types[4]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_types[5]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_ptr[0]) = ((void *)cbesh_wrap2);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_ptr[1]) = ((void *)((char *)__pyx_k__hankel2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_ptr[2]) = ((void *)cbesh_wrap2);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_ptr[3]) = ((void *)((char *)__pyx_k__hankel2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_types, 2, 2, 1, 0, __pyx_k__hankel2, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4035; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hankel2, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4035; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_doc = __pyx_k_88;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_types[1]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_types[2]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_types[4]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_types[5]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_ptr[0]) = ((void *)cbesh_wrap2_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_ptr[1]) = ((void *)((char *)__pyx_k__hankel2e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_ptr[2]) = ((void *)cbesh_wrap2_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_ptr[3]) = ((void *)((char *)__pyx_k__hankel2e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_types, 2, 2, 1, 0, __pyx_k__hankel2e, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hankel2e_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4059; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hankel2e, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4059; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_doc = __pyx_k_89;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ffF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_ddD__As_ddD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[6]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[7]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[14]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types[15]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[0]) = ((void *)hyp1f1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[1]) = ((void *)((char *)__pyx_k__hyp1f1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[2]) = ((void *)chyp1f1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[3]) = ((void *)((char *)__pyx_k__hyp1f1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[4]) = ((void *)hyp1f1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[5]) = ((void *)((char *)__pyx_k__hyp1f1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[6]) = ((void *)chyp1f1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[7]) = ((void *)((char *)__pyx_k__hyp1f1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_types, 4, 3, 1, 0, __pyx_k__hyp1f1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4100; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hyp1f1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4100; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_doc = __pyx_k_90;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_ptr[0]) = ((void *)onef2);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_ptr[1]) = ((void *)((char *)__pyx_k__hyp1f2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_ptr[2]) = ((void *)onef2);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_ptr[3]) = ((void *)((char *)__pyx_k__hyp1f2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_types, 2, 4, 2, 0, __pyx_k__hyp1f2, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp1f2_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4128; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hyp1f2, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4128; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_doc = __pyx_k_91;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddi_d_As_fffl_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddi_d_As_dddl_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[15]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[16]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[17]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[18]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[19]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[20]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[21]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[22]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types[23]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[0]) = ((void *)hyp2f0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[1]) = ((void *)((char *)__pyx_k__hyp2f0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[2]) = ((void *)__pyx_f_5scipy_7special_7_legacy_hyp2f0_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[3]) = ((void *)((char *)__pyx_k__hyp2f0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[4]) = ((void *)hyp2f0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[5]) = ((void *)((char *)__pyx_k__hyp2f0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_hyp2f0_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[7]) = ((void *)((char *)__pyx_k__hyp2f0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_types, 4, 4, 2, 0, __pyx_k__hyp2f0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4177; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hyp2f0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4177; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_doc = __pyx_k_92;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_fffF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dddD__As_dddD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[8]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[9]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[16]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[17]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[18]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types[19]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[0]) = ((void *)hyp2f1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[1]) = ((void *)((char *)__pyx_k__hyp2f1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[2]) = ((void *)chyp2f1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[3]) = ((void *)((char *)__pyx_k__hyp2f1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[4]) = ((void *)hyp2f1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[5]) = ((void *)((char *)__pyx_k__hyp2f1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[6]) = ((void *)chyp2f1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[7]) = ((void *)((char *)__pyx_k__hyp2f1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_types, 4, 4, 1, 0, __pyx_k__hyp2f1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp2f1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4222; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hyp2f1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4222; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_doc = __pyx_k_93;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_ptr[0]) = ((void *)threef0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_ptr[1]) = ((void *)((char *)__pyx_k__hyp3f0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_ptr[2]) = ((void *)threef0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_ptr[3]) = ((void *)((char *)__pyx_k__hyp3f0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_types, 2, 4, 2, 0, __pyx_k__hyp3f0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyp3f0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4250; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hyp3f0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4250; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_doc = __pyx_k_94;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_ptr[0]) = ((void *)hypU_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_ptr[1]) = ((void *)((char *)__pyx_k__hyperu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_ptr[2]) = ((void *)hypU_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_ptr[3]) = ((void *)((char *)__pyx_k__hyperu));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_types, 2, 3, 1, 0, __pyx_k__hyperu, __pyx_v_5scipy_7special_7_ufuncs_ufunc_hyperu_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4275; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__hyperu, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4275; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_doc = __pyx_k_95;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_ptr[0]) = ((void *)i0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_ptr[1]) = ((void *)((char *)__pyx_k__i0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_ptr[2]) = ((void *)i0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_ptr[3]) = ((void *)((char *)__pyx_k__i0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_types, 2, 1, 1, 0, __pyx_k__i0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4295; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__i0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4295; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_doc = __pyx_k_96;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_ptr[0]) = ((void *)i0e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_ptr[1]) = ((void *)((char *)__pyx_k__i0e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_ptr[2]) = ((void *)i0e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_ptr[3]) = ((void *)((char *)__pyx_k__i0e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_types, 2, 1, 1, 0, __pyx_k__i0e, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i0e_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4316; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__i0e, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4316; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_doc = __pyx_k_97;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_ptr[0]) = ((void *)i1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_ptr[1]) = ((void *)((char *)__pyx_k__i1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_ptr[2]) = ((void *)i1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_ptr[3]) = ((void *)((char *)__pyx_k__i1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_types, 2, 1, 1, 0, __pyx_k__i1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4336; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__i1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4336; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_doc = __pyx_k_98;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_ptr[0]) = ((void *)i1e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_ptr[1]) = ((void *)((char *)__pyx_k__i1e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_ptr[2]) = ((void *)i1e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_ptr[3]) = ((void *)((char *)__pyx_k__i1e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_types, 2, 1, 1, 0, __pyx_k__i1e, __pyx_v_5scipy_7special_7_ufuncs_ufunc_i1e_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4357; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__i1e, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4357; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_doc = __pyx_k_99;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_f_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_d_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_ptr[0]) = ((void *)it2i0k0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_ptr[1]) = ((void *)((char *)__pyx_k__it2i0k0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_ptr[2]) = ((void *)it2i0k0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_ptr[3]) = ((void *)((char *)__pyx_k__it2i0k0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_types, 2, 1, 2, 0, __pyx_k__it2i0k0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2i0k0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4380; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__it2i0k0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4380; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_doc = __pyx_k_100;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_f_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_d_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_ptr[0]) = ((void *)it2j0y0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_ptr[1]) = ((void *)((char *)__pyx_k__it2j0y0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_ptr[2]) = ((void *)it2j0y0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_ptr[3]) = ((void *)((char *)__pyx_k__it2j0y0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_types, 2, 1, 2, 0, __pyx_k__it2j0y0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2j0y0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4403; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__it2j0y0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4403; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_doc = __pyx_k_101;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_ptr[0]) = ((void *)it2struve0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_ptr[1]) = ((void *)((char *)__pyx_k__it2struve0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_ptr[2]) = ((void *)it2struve0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_ptr[3]) = ((void *)((char *)__pyx_k__it2struve0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_types, 2, 1, 1, 0, __pyx_k__it2struve0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_it2struve0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4424; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__it2struve0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4424; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_doc = __pyx_k_102;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_f_ffff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dddd_As_d_dddd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_ptr[0]) = ((void *)itairy_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_ptr[1]) = ((void *)((char *)__pyx_k__itairy));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_ptr[2]) = ((void *)itairy_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_ptr[3]) = ((void *)((char *)__pyx_k__itairy));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_types, 2, 1, 4, 0, __pyx_k__itairy, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itairy_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4451; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__itairy, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4451; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_doc = __pyx_k_103;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_f_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_d_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_ptr[0]) = ((void *)it1i0k0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_ptr[1]) = ((void *)((char *)__pyx_k__iti0k0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_ptr[2]) = ((void *)it1i0k0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_ptr[3]) = ((void *)((char *)__pyx_k__iti0k0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_types, 2, 1, 2, 0, __pyx_k__iti0k0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_iti0k0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4474; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__iti0k0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4474; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_doc = __pyx_k_104;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_f_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_d_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_ptr[0]) = ((void *)it1j0y0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_ptr[1]) = ((void *)((char *)__pyx_k__itj0y0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_ptr[2]) = ((void *)it1j0y0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_ptr[3]) = ((void *)((char *)__pyx_k__itj0y0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_types, 2, 1, 2, 0, __pyx_k__itj0y0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itj0y0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4497; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__itj0y0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4497; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_doc = __pyx_k_105;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_ptr[0]) = ((void *)itmodstruve0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_ptr[1]) = ((void *)((char *)__pyx_k__itmodstruve0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_ptr[2]) = ((void *)itmodstruve0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_ptr[3]) = ((void *)((char *)__pyx_k__itmodstruve0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_types, 2, 1, 1, 0, __pyx_k__itmodstruve0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itmodstruve0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4518; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__itmodstruve0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4518; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_doc = __pyx_k_106;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_ptr[0]) = ((void *)itstruve0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_ptr[1]) = ((void *)((char *)__pyx_k__itstruve0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_ptr[2]) = ((void *)itstruve0_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_ptr[3]) = ((void *)((char *)__pyx_k__itstruve0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_types, 2, 1, 1, 0, __pyx_k__itstruve0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_itstruve0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4539; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__itstruve0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4539; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_doc = __pyx_k_107;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[10]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types[11]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[0]) = ((void *)iv);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[1]) = ((void *)((char *)__pyx_k__iv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[2]) = ((void *)cbesi_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[3]) = ((void *)((char *)__pyx_k__iv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[4]) = ((void *)iv);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[5]) = ((void *)((char *)__pyx_k__iv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[6]) = ((void *)cbesi_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[7]) = ((void *)((char *)__pyx_k__iv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_types, 4, 2, 1, 0, __pyx_k__iv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_iv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4576; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__iv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4576; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_doc = __pyx_k_108;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[10]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types[11]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[0]) = ((void *)cbesi_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[1]) = ((void *)((char *)__pyx_k__ive));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[2]) = ((void *)cbesi_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[3]) = ((void *)((char *)__pyx_k__ive));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[4]) = ((void *)cbesi_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[5]) = ((void *)((char *)__pyx_k__ive));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[6]) = ((void *)cbesi_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[7]) = ((void *)((char *)__pyx_k__ive));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_types, 4, 2, 1, 0, __pyx_k__ive, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ive_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4613; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ive, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4613; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_doc = __pyx_k_109;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_ptr[0]) = ((void *)j0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_ptr[1]) = ((void *)((char *)__pyx_k__j0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_ptr[2]) = ((void *)j0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_ptr[3]) = ((void *)((char *)__pyx_k__j0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_types, 2, 1, 1, 0, __pyx_k__j0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_j0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4633; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__j0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4633; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_doc = __pyx_k_110;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_ptr[0]) = ((void *)j1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_ptr[1]) = ((void *)((char *)__pyx_k__j1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_ptr[2]) = ((void *)j1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_ptr[3]) = ((void *)((char *)__pyx_k__j1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_types, 2, 1, 1, 0, __pyx_k__j1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_j1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4653; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__j1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4653; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_doc = __pyx_k_111;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[10]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types[11]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[0]) = ((void *)jv);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[1]) = ((void *)((char *)__pyx_k__jv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[2]) = ((void *)cbesj_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[3]) = ((void *)((char *)__pyx_k__jv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[4]) = ((void *)jv);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[5]) = ((void *)((char *)__pyx_k__jv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[6]) = ((void *)cbesj_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[7]) = ((void *)((char *)__pyx_k__jv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_types, 4, 2, 1, 0, __pyx_k__jv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_jv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4689; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__jv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4689; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_doc = __pyx_k_112;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[10]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types[11]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[0]) = ((void *)cbesj_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[1]) = ((void *)((char *)__pyx_k__jve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[2]) = ((void *)cbesj_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[3]) = ((void *)((char *)__pyx_k__jve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[4]) = ((void *)cbesj_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[5]) = ((void *)((char *)__pyx_k__jve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[6]) = ((void *)cbesj_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[7]) = ((void *)((char *)__pyx_k__jve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_types, 4, 2, 1, 0, __pyx_k__jve, __pyx_v_5scipy_7special_7_ufuncs_ufunc_jve_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4726; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__jve, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4726; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_doc = __pyx_k_113;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_ptr[0]) = ((void *)k0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_ptr[1]) = ((void *)((char *)__pyx_k__k0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_ptr[2]) = ((void *)k0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_ptr[3]) = ((void *)((char *)__pyx_k__k0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_types, 2, 1, 1, 0, __pyx_k__k0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4747; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__k0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4747; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_doc = __pyx_k_114;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_ptr[0]) = ((void *)k0e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_ptr[1]) = ((void *)((char *)__pyx_k__k0e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_ptr[2]) = ((void *)k0e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_ptr[3]) = ((void *)((char *)__pyx_k__k0e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_types, 2, 1, 1, 0, __pyx_k__k0e, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k0e_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4768; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__k0e, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4768; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_doc = __pyx_k_115;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_ptr[0]) = ((void *)k1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_ptr[1]) = ((void *)((char *)__pyx_k__k1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_ptr[2]) = ((void *)k1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_ptr[3]) = ((void *)((char *)__pyx_k__k1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_types, 2, 1, 1, 0, __pyx_k__k1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4789; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__k1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4789; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_doc = __pyx_k_116;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_ptr[0]) = ((void *)k1e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_ptr[1]) = ((void *)((char *)__pyx_k__k1e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_ptr[2]) = ((void *)k1e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_ptr[3]) = ((void *)((char *)__pyx_k__k1e));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_types, 2, 1, 1, 0, __pyx_k__k1e, __pyx_v_5scipy_7special_7_ufuncs_ufunc_k1e_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4810; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__k1e, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4810; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_doc = __pyx_k_117;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_ptr[0]) = ((void *)kei_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_ptr[1]) = ((void *)((char *)__pyx_k__kei));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_ptr[2]) = ((void *)kei_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_ptr[3]) = ((void *)((char *)__pyx_k__kei));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_types, 2, 1, 1, 0, __pyx_k__kei, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kei_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4830; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__kei, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4830; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_doc = __pyx_k_118;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_ptr[0]) = ((void *)keip_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_ptr[1]) = ((void *)((char *)__pyx_k__keip));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_ptr[2]) = ((void *)keip_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_ptr[3]) = ((void *)((char *)__pyx_k__keip));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_types, 2, 1, 1, 0, __pyx_k__keip, __pyx_v_5scipy_7special_7_ufuncs_ufunc_keip_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4850; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__keip, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4850; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_doc = __pyx_k_119;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DDDD_As_f_FFFF);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DDDD_As_d_DDDD);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[1]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[2]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[3]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[6]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[7]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[8]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types[9]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_ptr[0]) = ((void *)kelvin_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_ptr[1]) = ((void *)((char *)__pyx_k__kelvin));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_ptr[2]) = ((void *)kelvin_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_ptr[3]) = ((void *)((char *)__pyx_k__kelvin));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_types, 2, 1, 4, 0, __pyx_k__kelvin, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kelvin_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4880; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__kelvin, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4880; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_doc = __pyx_k_120;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_ptr[0]) = ((void *)ker_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_ptr[1]) = ((void *)((char *)__pyx_k__ker));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_ptr[2]) = ((void *)ker_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_ptr[3]) = ((void *)((char *)__pyx_k__ker));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_types, 2, 1, 1, 0, __pyx_k__ker, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ker_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4900; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ker, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4900; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_doc = __pyx_k_121;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_ptr[0]) = ((void *)kerp_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_ptr[1]) = ((void *)((char *)__pyx_k__kerp));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_ptr[2]) = ((void *)kerp_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_ptr[3]) = ((void *)((char *)__pyx_k__kerp));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_types, 2, 1, 1, 0, __pyx_k__kerp, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kerp_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4920; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__kerp, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4920; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_doc = __pyx_k_122;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[0]) = ((void *)kn);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[1]) = ((void *)((char *)__pyx_k__kn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[2]) = ((void *)kn);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[3]) = ((void *)((char *)__pyx_k__kn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_kn_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[5]) = ((void *)((char *)__pyx_k__kn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_kn_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[7]) = ((void *)((char *)__pyx_k__kn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_types, 4, 2, 1, 0, __pyx_k__kn, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kn_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4957; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__kn, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4957; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_doc = __pyx_k_123;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_ptr[0]) = ((void *)kolmogi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_ptr[1]) = ((void *)((char *)__pyx_k__kolmogi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_ptr[2]) = ((void *)kolmogi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_ptr[3]) = ((void *)((char *)__pyx_k__kolmogi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_types, 2, 1, 1, 0, __pyx_k__kolmogi, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogi_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4977; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__kolmogi, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 4977; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_doc = __pyx_k_124;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_ptr[0]) = ((void *)kolmogorov);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_ptr[1]) = ((void *)((char *)__pyx_k__kolmogorov));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_ptr[2]) = ((void *)kolmogorov);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_ptr[3]) = ((void *)((char *)__pyx_k__kolmogorov));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_types, 2, 1, 1, 0, __pyx_k__kolmogorov, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kolmogorov_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5001; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__kolmogorov, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5001; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_doc = __pyx_k_125;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[10]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types[11]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[0]) = ((void *)cbesk_wrap_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[1]) = ((void *)((char *)__pyx_k__kv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[2]) = ((void *)cbesk_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[3]) = ((void *)((char *)__pyx_k__kv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[4]) = ((void *)cbesk_wrap_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[5]) = ((void *)((char *)__pyx_k__kv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[6]) = ((void *)cbesk_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[7]) = ((void *)((char *)__pyx_k__kv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_types, 4, 2, 1, 0, __pyx_k__kv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5038; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__kv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5038; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_doc = __pyx_k_126;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[10]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types[11]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[0]) = ((void *)cbesk_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[1]) = ((void *)((char *)__pyx_k__kve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[2]) = ((void *)cbesk_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[3]) = ((void *)((char *)__pyx_k__kve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[4]) = ((void *)cbesk_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[5]) = ((void *)((char *)__pyx_k__kve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[6]) = ((void *)cbesk_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[7]) = ((void *)((char *)__pyx_k__kve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_types, 4, 2, 1, 0, __pyx_k__kve, __pyx_v_5scipy_7special_7_ufuncs_ufunc_kve_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5075; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__kve, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5075; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_doc = __pyx_k_127;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_ptr[0]) = ((void *)log1p);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_ptr[1]) = ((void *)((char *)__pyx_k__log1p));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_ptr[2]) = ((void *)log1p);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_ptr[3]) = ((void *)((char *)__pyx_k__log1p));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_types, 2, 1, 1, 0, __pyx_k__log1p, __pyx_v_5scipy_7special_7_ufuncs_ufunc_log1p_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5095; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__log1p, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5095; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_doc = __pyx_k_128;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_ptr[0]) = ((void *)log_ndtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_ptr[1]) = ((void *)((char *)__pyx_k__log_ndtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_ptr[2]) = ((void *)log_ndtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_ptr[3]) = ((void *)((char *)__pyx_k__log_ndtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_types, 2, 1, 1, 0, __pyx_k__log_ndtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_log_ndtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5117; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__log_ndtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5117; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_doc = __pyx_k_129;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_f_f__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_g_g__As_g_g);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_types[4]) = ((char)NPY_LONGDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_types[5]) = ((char)NPY_LONGDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[0]) = ((void *)logitf);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[1]) = ((void *)((char *)__pyx_k__logit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[2]) = ((void *)logit);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[3]) = ((void *)((char *)__pyx_k__logit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[4]) = ((void *)logitl);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[5]) = ((void *)((char *)__pyx_k__logit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_ptr[4]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_types, 3, 1, 1, 0, __pyx_k__logit, __pyx_v_5scipy_7special_7_ufuncs_ufunc_logit_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5164; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__logit, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5164; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_doc = __pyx_k_130;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_ptr[0]) = ((void *)pmv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_ptr[1]) = ((void *)((char *)__pyx_k__lpmv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_ptr[2]) = ((void *)pmv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_ptr[3]) = ((void *)((char *)__pyx_k__lpmv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_types, 2, 3, 1, 0, __pyx_k__lpmv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_lpmv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5189; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__lpmv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5189; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_doc = __pyx_k_131;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_ptr[0]) = ((void *)cem_cva_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_ptr[1]) = ((void *)((char *)__pyx_k__mathieu_a));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_ptr[2]) = ((void *)cem_cva_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_ptr[3]) = ((void *)((char *)__pyx_k__mathieu_a));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_types, 2, 2, 1, 0, __pyx_k__mathieu_a, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_a_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5212; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__mathieu_a, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5212; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_doc = __pyx_k_132;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_ptr[0]) = ((void *)sem_cva_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_ptr[1]) = ((void *)((char *)__pyx_k__mathieu_b));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_ptr[2]) = ((void *)sem_cva_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_ptr[3]) = ((void *)((char *)__pyx_k__mathieu_b));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_types, 2, 2, 1, 0, __pyx_k__mathieu_b, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_b_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5235; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__mathieu_b, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5235; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_doc = __pyx_k_133;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_fff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_ddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_ptr[0]) = ((void *)cem_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_ptr[1]) = ((void *)((char *)__pyx_k__mathieu_cem));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_ptr[2]) = ((void *)cem_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_ptr[3]) = ((void *)((char *)__pyx_k__mathieu_cem));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_types, 2, 3, 2, 0, __pyx_k__mathieu_cem, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_cem_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5263; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__mathieu_cem, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5263; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_doc = __pyx_k_134;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_fff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_ddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_ptr[0]) = ((void *)mcm1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_ptr[1]) = ((void *)((char *)__pyx_k__mathieu_modcem1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_ptr[2]) = ((void *)mcm1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_ptr[3]) = ((void *)((char *)__pyx_k__mathieu_modcem1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_types, 2, 3, 2, 0, __pyx_k__mathieu_modcem1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5291; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__mathieu_modcem1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5291; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_doc = __pyx_k_135;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_fff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_ddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_ptr[0]) = ((void *)mcm2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_ptr[1]) = ((void *)((char *)__pyx_k__mathieu_modcem2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_ptr[2]) = ((void *)mcm2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_ptr[3]) = ((void *)((char *)__pyx_k__mathieu_modcem2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_types, 2, 3, 2, 0, __pyx_k__mathieu_modcem2, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modcem2_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5319; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__mathieu_modcem2, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5319; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_doc = __pyx_k_136;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_fff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_ddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_ptr[0]) = ((void *)msm1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_ptr[1]) = ((void *)((char *)__pyx_k__mathieu_modsem1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_ptr[2]) = ((void *)msm1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_ptr[3]) = ((void *)((char *)__pyx_k__mathieu_modsem1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_types, 2, 3, 2, 0, __pyx_k__mathieu_modsem1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5347; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__mathieu_modsem1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5347; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_doc = __pyx_k_137;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_fff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_ddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_ptr[0]) = ((void *)msm2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_ptr[1]) = ((void *)((char *)__pyx_k__mathieu_modsem2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_ptr[2]) = ((void *)msm2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_ptr[3]) = ((void *)((char *)__pyx_k__mathieu_modsem2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_types, 2, 3, 2, 0, __pyx_k__mathieu_modsem2, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_modsem2_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5375; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__mathieu_modsem2, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5375; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_doc = __pyx_k_138;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_fff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddd_dd_As_ddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_ptr[0]) = ((void *)sem_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_ptr[1]) = ((void *)((char *)__pyx_k__mathieu_sem));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_ptr[2]) = ((void *)sem_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_ptr[3]) = ((void *)((char *)__pyx_k__mathieu_sem));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_types, 2, 3, 2, 0, __pyx_k__mathieu_sem, __pyx_v_5scipy_7special_7_ufuncs_ufunc_mathieu_sem_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5403; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__mathieu_sem, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5403; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_doc = __pyx_k_139;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DD_As_f_FF);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DD_As_d_DD);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_types[1]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_types[2]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_types[4]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_types[5]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_ptr[0]) = ((void *)modified_fresnel_minus_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_ptr[1]) = ((void *)((char *)__pyx_k__modfresnelm));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_ptr[2]) = ((void *)modified_fresnel_minus_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_ptr[3]) = ((void *)((char *)__pyx_k__modfresnelm));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_types, 2, 1, 2, 0, __pyx_k__modfresnelm, __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelm_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5426; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__modfresnelm, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5426; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_doc = __pyx_k_140;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DD_As_f_FF);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_DD_As_d_DD);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_types[1]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_types[2]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_types[4]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_types[5]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_ptr[0]) = ((void *)modified_fresnel_plus_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_ptr[1]) = ((void *)((char *)__pyx_k__modfresnelp));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_ptr[2]) = ((void *)modified_fresnel_plus_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_ptr[3]) = ((void *)((char *)__pyx_k__modfresnelp));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_types, 2, 1, 2, 0, __pyx_k__modfresnelp, __pyx_v_5scipy_7special_7_ufuncs_ufunc_modfresnelp_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5449; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__modfresnelp, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5449; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_doc = __pyx_k_141;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_ptr[0]) = ((void *)modstruve_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_ptr[1]) = ((void *)((char *)__pyx_k__modstruve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_ptr[2]) = ((void *)modstruve_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_ptr[3]) = ((void *)((char *)__pyx_k__modstruve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_types, 2, 2, 1, 0, __pyx_k__modstruve, __pyx_v_5scipy_7special_7_ufuncs_ufunc_modstruve_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5473; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__modstruve, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5473; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_doc = __pyx_k_142;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_llf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_lld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[1]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[4]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[5]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[0]) = ((void *)nbdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[1]) = ((void *)((char *)__pyx_k__nbdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[2]) = ((void *)nbdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[3]) = ((void *)((char *)__pyx_k__nbdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_nbdtr_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[5]) = ((void *)((char *)__pyx_k__nbdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_nbdtr_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[7]) = ((void *)((char *)__pyx_k__nbdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_types, 4, 3, 1, 0, __pyx_k__nbdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5516; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nbdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5516; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_doc = __pyx_k_143;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_llf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_lld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[1]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[4]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[5]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[0]) = ((void *)nbdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[1]) = ((void *)((char *)__pyx_k__nbdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[2]) = ((void *)nbdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[3]) = ((void *)((char *)__pyx_k__nbdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_nbdtrc_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[5]) = ((void *)((char *)__pyx_k__nbdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_nbdtrc_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[7]) = ((void *)((char *)__pyx_k__nbdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_types, 4, 3, 1, 0, __pyx_k__nbdtrc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5557; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nbdtrc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5557; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_doc = __pyx_k_144;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_llf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_iid__As_lld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[1]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[4]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[5]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[9]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[10]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[11]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[14]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types[15]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[0]) = ((void *)nbdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[1]) = ((void *)((char *)__pyx_k__nbdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[2]) = ((void *)nbdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[3]) = ((void *)((char *)__pyx_k__nbdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_nbdtri_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[5]) = ((void *)((char *)__pyx_k__nbdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_nbdtri_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[7]) = ((void *)((char *)__pyx_k__nbdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_types, 4, 3, 1, 0, __pyx_k__nbdtri, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtri_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5597; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nbdtri, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5597; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_doc = __pyx_k_145;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_ptr[0]) = ((void *)cdfnbn2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_ptr[1]) = ((void *)((char *)__pyx_k__nbdtrik));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_ptr[2]) = ((void *)cdfnbn2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_ptr[3]) = ((void *)((char *)__pyx_k__nbdtrik));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_types, 2, 3, 1, 0, __pyx_k__nbdtrik, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrik_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5621; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nbdtrik, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5621; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_doc = __pyx_k_146;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_ptr[0]) = ((void *)cdfnbn3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_ptr[1]) = ((void *)((char *)__pyx_k__nbdtrin));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_ptr[2]) = ((void *)cdfnbn3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_ptr[3]) = ((void *)((char *)__pyx_k__nbdtrin));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_types, 2, 3, 1, 0, __pyx_k__nbdtrin, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nbdtrin_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5645; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nbdtrin, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5645; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_ptr[0]) = ((void *)cdffnc1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_ptr[1]) = ((void *)((char *)__pyx_k__ncfdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_ptr[2]) = ((void *)cdffnc1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_ptr[3]) = ((void *)((char *)__pyx_k__ncfdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_types, 2, 4, 1, 0, __pyx_k__ncfdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5671; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ncfdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5671; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_ptr[0]) = ((void *)cdffnc2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_ptr[1]) = ((void *)((char *)__pyx_k__ncfdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_ptr[2]) = ((void *)cdffnc2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_ptr[3]) = ((void *)((char *)__pyx_k__ncfdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_types, 2, 4, 1, 0, __pyx_k__ncfdtri, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtri_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5697; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ncfdtri, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5697; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_ptr[0]) = ((void *)cdffnc4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_ptr[1]) = ((void *)((char *)__pyx_k__ncfdtridfd));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_ptr[2]) = ((void *)cdffnc4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_ptr[3]) = ((void *)((char *)__pyx_k__ncfdtridfd));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_types, 2, 4, 1, 0, __pyx_k__ncfdtridfd, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfd_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5723; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ncfdtridfd, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5723; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_ptr[0]) = ((void *)cdffnc3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_ptr[1]) = ((void *)((char *)__pyx_k__ncfdtridfn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_ptr[2]) = ((void *)cdffnc3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_ptr[3]) = ((void *)((char *)__pyx_k__ncfdtridfn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_types, 2, 4, 1, 0, __pyx_k__ncfdtridfn, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtridfn_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5749; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ncfdtridfn, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5749; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_ffff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd__As_dddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_ptr[0]) = ((void *)cdffnc5_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_ptr[1]) = ((void *)((char *)__pyx_k__ncfdtrinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_ptr[2]) = ((void *)cdffnc5_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_ptr[3]) = ((void *)((char *)__pyx_k__ncfdtrinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_types, 2, 4, 1, 0, __pyx_k__ncfdtrinc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ncfdtrinc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5775; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ncfdtrinc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5775; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_ptr[0]) = ((void *)cdftnc1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_ptr[1]) = ((void *)((char *)__pyx_k__nctdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_ptr[2]) = ((void *)cdftnc1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_ptr[3]) = ((void *)((char *)__pyx_k__nctdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_types, 2, 3, 1, 0, __pyx_k__nctdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5799; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nctdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5799; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_ptr[0]) = ((void *)cdftnc3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_ptr[1]) = ((void *)((char *)__pyx_k__nctdtridf));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_ptr[2]) = ((void *)cdftnc3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_ptr[3]) = ((void *)((char *)__pyx_k__nctdtridf));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_types, 2, 3, 1, 0, __pyx_k__nctdtridf, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtridf_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5823; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nctdtridf, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5823; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_ptr[0]) = ((void *)cdftnc4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_ptr[1]) = ((void *)((char *)__pyx_k__nctdtrinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_ptr[2]) = ((void *)cdftnc4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_ptr[3]) = ((void *)((char *)__pyx_k__nctdtrinc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_types, 2, 3, 1, 0, __pyx_k__nctdtrinc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrinc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5847; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nctdtrinc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5847; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_ptr[0]) = ((void *)cdftnc2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_ptr[1]) = ((void *)((char *)__pyx_k__nctdtrit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_ptr[2]) = ((void *)cdftnc2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_ptr[3]) = ((void *)((char *)__pyx_k__nctdtrit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_types, 2, 3, 1, 0, __pyx_k__nctdtrit, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nctdtrit_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5871; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nctdtrit, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5871; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_doc = __pyx_k_147;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_ptr[0]) = ((void *)ndtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_ptr[1]) = ((void *)((char *)__pyx_k__ndtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_ptr[2]) = ((void *)ndtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_ptr[3]) = ((void *)((char *)__pyx_k__ndtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_types, 2, 1, 1, 0, __pyx_k__ndtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5893; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ndtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5893; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_doc = __pyx_k_148;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_ptr[0]) = ((void *)ndtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_ptr[1]) = ((void *)((char *)__pyx_k__ndtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_ptr[2]) = ((void *)ndtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_ptr[3]) = ((void *)((char *)__pyx_k__ndtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_types, 2, 1, 1, 0, __pyx_k__ndtri, __pyx_v_5scipy_7special_7_ufuncs_ufunc_ndtri_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__ndtri, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_ptr[0]) = ((void *)cdfnor3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_ptr[1]) = ((void *)((char *)__pyx_k__nrdtrimn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_ptr[2]) = ((void *)cdfnor3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_ptr[3]) = ((void *)((char *)__pyx_k__nrdtrimn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_types, 2, 3, 1, 0, __pyx_k__nrdtrimn, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrimn_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5939; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nrdtrimn, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5939; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_ptr[0]) = ((void *)cdfnor4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_ptr[1]) = ((void *)((char *)__pyx_k__nrdtrisd));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_ptr[2]) = ((void *)cdfnor4_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_ptr[3]) = ((void *)((char *)__pyx_k__nrdtrisd));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_types, 2, 3, 1, 0, __pyx_k__nrdtrisd, __pyx_v_5scipy_7special_7_ufuncs_ufunc_nrdtrisd_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5963; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__nrdtrisd, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5963; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_doc = __pyx_k_149;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_ptr[0]) = ((void *)oblate_aswfa_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_ptr[1]) = ((void *)((char *)__pyx_k__obl_ang1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_ptr[2]) = ((void *)oblate_aswfa_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_ptr[3]) = ((void *)((char *)__pyx_k__obl_ang1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_types, 2, 4, 2, 0, __pyx_k__obl_ang1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5993; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__obl_ang1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 5993; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_doc = __pyx_k_150;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_fffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_ddddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_ptr[0]) = ((void *)oblate_aswfa_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_ptr[1]) = ((void *)((char *)__pyx_k__obl_ang1_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_ptr[2]) = ((void *)oblate_aswfa_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_ptr[3]) = ((void *)((char *)__pyx_k__obl_ang1_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_types, 2, 5, 2, 0, __pyx_k__obl_ang1_cv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_ang1_cv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6026; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__obl_ang1_cv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6026; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_doc = __pyx_k_151;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_ptr[0]) = ((void *)oblate_segv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_ptr[1]) = ((void *)((char *)__pyx_k__obl_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_ptr[2]) = ((void *)oblate_segv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_ptr[3]) = ((void *)((char *)__pyx_k__obl_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_types, 2, 3, 1, 0, __pyx_k__obl_cv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_cv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6051; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__obl_cv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6051; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_doc = __pyx_k_152;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_ptr[0]) = ((void *)oblate_radial1_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_ptr[1]) = ((void *)((char *)__pyx_k__obl_rad1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_ptr[2]) = ((void *)oblate_radial1_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_ptr[3]) = ((void *)((char *)__pyx_k__obl_rad1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_types, 2, 4, 2, 0, __pyx_k__obl_rad1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6081; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__obl_rad1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6081; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_doc = __pyx_k_153;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_fffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_ddddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_ptr[0]) = ((void *)oblate_radial1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_ptr[1]) = ((void *)((char *)__pyx_k__obl_rad1_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_ptr[2]) = ((void *)oblate_radial1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_ptr[3]) = ((void *)((char *)__pyx_k__obl_rad1_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_types, 2, 5, 2, 0, __pyx_k__obl_rad1_cv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad1_cv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6114; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__obl_rad1_cv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6114; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_doc = __pyx_k_154;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_ptr[0]) = ((void *)oblate_radial2_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_ptr[1]) = ((void *)((char *)__pyx_k__obl_rad2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_ptr[2]) = ((void *)oblate_radial2_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_ptr[3]) = ((void *)((char *)__pyx_k__obl_rad2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_types, 2, 4, 2, 0, __pyx_k__obl_rad2, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6144; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__obl_rad2, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6144; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_doc = __pyx_k_155;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_fffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_ddddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_ptr[0]) = ((void *)oblate_radial2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_ptr[1]) = ((void *)((char *)__pyx_k__obl_rad2_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_ptr[2]) = ((void *)oblate_radial2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_ptr[3]) = ((void *)((char *)__pyx_k__obl_rad2_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_types, 2, 5, 2, 0, __pyx_k__obl_rad2_cv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_obl_rad2_cv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6177; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__obl_rad2_cv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6177; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_doc = __pyx_k_156;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_ff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_dd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_ptr[0]) = ((void *)pbdv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_ptr[1]) = ((void *)((char *)__pyx_k__pbdv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_ptr[2]) = ((void *)pbdv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_ptr[3]) = ((void *)((char *)__pyx_k__pbdv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_types, 2, 2, 2, 0, __pyx_k__pbdv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbdv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6202; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pbdv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6202; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_doc = __pyx_k_157;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_ff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_dd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_ptr[0]) = ((void *)pbvv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_ptr[1]) = ((void *)((char *)__pyx_k__pbvv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_ptr[2]) = ((void *)pbvv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_ptr[3]) = ((void *)((char *)__pyx_k__pbvv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_types, 2, 2, 2, 0, __pyx_k__pbvv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbvv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6227; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pbvv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6227; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_doc = __pyx_k_158;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_ff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_dd_dd_As_dd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_ptr[0]) = ((void *)pbwa_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_ptr[1]) = ((void *)((char *)__pyx_k__pbwa));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_ptr[2]) = ((void *)pbwa_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_ptr[3]) = ((void *)((char *)__pyx_k__pbwa));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_types, 2, 2, 2, 0, __pyx_k__pbwa, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pbwa_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6253; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pbwa, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6253; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_doc = __pyx_k_159;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[0]) = ((void *)pdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[1]) = ((void *)((char *)__pyx_k__pdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[2]) = ((void *)pdtr);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[3]) = ((void *)((char *)__pyx_k__pdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_pdtr_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[5]) = ((void *)((char *)__pyx_k__pdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_pdtr_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[7]) = ((void *)((char *)__pyx_k__pdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_types, 4, 2, 1, 0, __pyx_k__pdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6291; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6291; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_doc = __pyx_k_160;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[0]) = ((void *)pdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[1]) = ((void *)((char *)__pyx_k__pdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[2]) = ((void *)pdtrc);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[3]) = ((void *)((char *)__pyx_k__pdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_pdtrc_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[5]) = ((void *)((char *)__pyx_k__pdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_pdtrc_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[7]) = ((void *)((char *)__pyx_k__pdtrc));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_types, 4, 2, 1, 0, __pyx_k__pdtrc, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrc_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pdtrc, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_doc = __pyx_k_161;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[0]) = ((void *)pdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[1]) = ((void *)((char *)__pyx_k__pdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[2]) = ((void *)pdtri);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[3]) = ((void *)((char *)__pyx_k__pdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_pdtri_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[5]) = ((void *)((char *)__pyx_k__pdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_pdtri_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[7]) = ((void *)((char *)__pyx_k__pdtri));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_types, 4, 2, 1, 0, __pyx_k__pdtri, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtri_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6368; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pdtri, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6368; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_doc = __pyx_k_162;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_ptr[0]) = ((void *)cdfpoi2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_ptr[1]) = ((void *)((char *)__pyx_k__pdtrik));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_ptr[2]) = ((void *)cdfpoi2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_ptr[3]) = ((void *)((char *)__pyx_k__pdtrik));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_types, 2, 2, 1, 0, __pyx_k__pdtrik, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pdtrik_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6390; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pdtrik, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6390; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_doc = __pyx_k_163;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_ptr[0]) = ((void *)prolate_aswfa_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_ptr[1]) = ((void *)((char *)__pyx_k__pro_ang1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_ptr[2]) = ((void *)prolate_aswfa_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_ptr[3]) = ((void *)((char *)__pyx_k__pro_ang1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_types, 2, 4, 2, 0, __pyx_k__pro_ang1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6420; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pro_ang1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6420; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_doc = __pyx_k_164;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_fffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_ddddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_ptr[0]) = ((void *)prolate_aswfa_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_ptr[1]) = ((void *)((char *)__pyx_k__pro_ang1_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_ptr[2]) = ((void *)prolate_aswfa_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_ptr[3]) = ((void *)((char *)__pyx_k__pro_ang1_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_types, 2, 5, 2, 0, __pyx_k__pro_ang1_cv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_ang1_cv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6453; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pro_ang1_cv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6453; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_doc = __pyx_k_165;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_ptr[0]) = ((void *)prolate_segv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_ptr[1]) = ((void *)((char *)__pyx_k__pro_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_ptr[2]) = ((void *)prolate_segv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_ptr[3]) = ((void *)((char *)__pyx_k__pro_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_types, 2, 3, 1, 0, __pyx_k__pro_cv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_cv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6478; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pro_cv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6478; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_doc = __pyx_k_166;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_ptr[0]) = ((void *)prolate_radial1_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_ptr[1]) = ((void *)((char *)__pyx_k__pro_rad1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_ptr[2]) = ((void *)prolate_radial1_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_ptr[3]) = ((void *)((char *)__pyx_k__pro_rad1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_types, 2, 4, 2, 0, __pyx_k__pro_rad1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6508; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pro_rad1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6508; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_doc = __pyx_k_167;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_fffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_ddddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_ptr[0]) = ((void *)prolate_radial1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_ptr[1]) = ((void *)((char *)__pyx_k__pro_rad1_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_ptr[2]) = ((void *)prolate_radial1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_ptr[3]) = ((void *)((char *)__pyx_k__pro_rad1_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_types, 2, 5, 2, 0, __pyx_k__pro_rad1_cv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad1_cv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6541; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pro_rad1_cv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6541; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_doc = __pyx_k_168;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_ffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dddd_d_As_dddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_ptr[0]) = ((void *)prolate_radial2_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_ptr[1]) = ((void *)((char *)__pyx_k__pro_rad2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_ptr[2]) = ((void *)prolate_radial2_nocv_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_ptr[3]) = ((void *)((char *)__pyx_k__pro_rad2));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_types, 2, 4, 2, 0, __pyx_k__pro_rad2, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6571; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pro_rad2, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6571; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_doc = __pyx_k_169;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_fffff_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_ddddd_dd_As_ddddd_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[4]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[5]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[12]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types[13]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_ptr[0]) = ((void *)prolate_radial2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_ptr[1]) = ((void *)((char *)__pyx_k__pro_rad2_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_ptr[2]) = ((void *)prolate_radial2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_ptr[3]) = ((void *)((char *)__pyx_k__pro_rad2_cv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_types, 2, 5, 2, 0, __pyx_k__pro_rad2_cv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_pro_rad2_cv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6604; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__pro_rad2_cv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6604; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_doc = __pyx_k_170;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_F_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_D_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types[6]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types[7]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[0]) = ((void *)psi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[1]) = ((void *)((char *)__pyx_k__psi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[2]) = ((void *)psi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[3]) = ((void *)((char *)__pyx_k__psi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[4]) = ((void *)cpsi_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[5]) = ((void *)((char *)__pyx_k__psi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[6]) = ((void *)cpsi_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[7]) = ((void *)((char *)__pyx_k__psi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_types, 4, 1, 1, 0, __pyx_k__psi, __pyx_v_5scipy_7special_7_ufuncs_ufunc_psi_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6637; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__psi, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6637; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_doc = __pyx_k_171;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_fff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_ddd__As_ddd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_ptr[0]) = ((void *)radian);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_ptr[1]) = ((void *)((char *)__pyx_k__radian));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_ptr[2]) = ((void *)radian);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_ptr[3]) = ((void *)((char *)__pyx_k__radian));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_types, 2, 3, 1, 0, __pyx_k__radian, __pyx_v_5scipy_7special_7_ufuncs_ufunc_radian_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6662; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__radian, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6662; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_doc = __pyx_k_172;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_F_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_D__As_D_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types[6]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types[7]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[0]) = ((void *)rgamma);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[1]) = ((void *)((char *)__pyx_k__rgamma));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[2]) = ((void *)rgamma);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[3]) = ((void *)((char *)__pyx_k__rgamma));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[4]) = ((void *)crgamma_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[5]) = ((void *)((char *)__pyx_k__rgamma));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[6]) = ((void *)crgamma_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[7]) = ((void *)((char *)__pyx_k__rgamma));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_types, 4, 1, 1, 0, __pyx_k__rgamma, __pyx_v_5scipy_7special_7_ufuncs_ufunc_rgamma_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6694; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__rgamma, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6694; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_round_doc = __pyx_k_173;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_ptr[0]) = ((void *)round);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_ptr[1]) = ((void *)((char *)__pyx_k__round));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_ptr[2]) = ((void *)round);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_ptr[3]) = ((void *)((char *)__pyx_k__round));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_round_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_round_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_round_types, 2, 1, 1, 0, __pyx_k__round, __pyx_v_5scipy_7special_7_ufuncs_ufunc_round_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6716; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__round, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6716; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_doc = __pyx_k_174;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_f_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_d_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_ptr[0]) = ((void *)shichi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_ptr[1]) = ((void *)((char *)__pyx_k__shichi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_ptr[2]) = ((void *)shichi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_ptr[3]) = ((void *)((char *)__pyx_k__shichi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_types, 2, 1, 2, 0, __pyx_k__shichi, __pyx_v_5scipy_7special_7_ufuncs_ufunc_shichi_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6740; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__shichi, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6740; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_doc = __pyx_k_175;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_f_ff);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_i_d_dd_As_d_dd);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_ptr[0]) = ((void *)sici);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_ptr[1]) = ((void *)((char *)__pyx_k__sici));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_ptr[2]) = ((void *)sici);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_ptr[3]) = ((void *)((char *)__pyx_k__sici));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_types, 2, 1, 2, 0, __pyx_k__sici, __pyx_v_5scipy_7special_7_ufuncs_ufunc_sici_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6764; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__sici, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6764; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_doc = __pyx_k_176;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_ptr[0]) = ((void *)sindg);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_ptr[1]) = ((void *)((char *)__pyx_k__sindg));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_ptr[2]) = ((void *)sindg);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_ptr[3]) = ((void *)((char *)__pyx_k__sindg));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_types, 2, 1, 1, 0, __pyx_k__sindg, __pyx_v_5scipy_7special_7_ufuncs_ufunc_sindg_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6784; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__sindg, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6784; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_doc = __pyx_k_177;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[0]) = ((void *)smirnov);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[1]) = ((void *)((char *)__pyx_k__smirnov));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[2]) = ((void *)smirnov);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[3]) = ((void *)((char *)__pyx_k__smirnov));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_smirnov_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[5]) = ((void *)((char *)__pyx_k__smirnov));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_smirnov_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[7]) = ((void *)((char *)__pyx_k__smirnov));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_types, 4, 2, 1, 0, __pyx_k__smirnov, __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnov_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6824; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__smirnov, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6824; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_doc = __pyx_k_178;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[0]) = ((void *)smirnovi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[1]) = ((void *)((char *)__pyx_k__smirnovi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[2]) = ((void *)smirnovi);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[3]) = ((void *)((char *)__pyx_k__smirnovi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_smirnovi_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[5]) = ((void *)((char *)__pyx_k__smirnovi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_smirnovi_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[7]) = ((void *)((char *)__pyx_k__smirnovi));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_types, 4, 2, 1, 0, __pyx_k__smirnovi, __pyx_v_5scipy_7special_7_ufuncs_ufunc_smirnovi_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6860; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__smirnovi, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6860; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_doc = __pyx_k_179;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_ptr[0]) = ((void *)spence);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_ptr[1]) = ((void *)((char *)__pyx_k__spence));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_ptr[2]) = ((void *)spence);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_ptr[3]) = ((void *)((char *)__pyx_k__spence));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_types, 2, 1, 1, 0, __pyx_k__spence, __pyx_v_5scipy_7special_7_ufuncs_ufunc_spence_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6881; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__spence, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6881; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_doc = __pyx_k_180;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_ptr[0]) = ((void *)cdft1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_ptr[1]) = ((void *)((char *)__pyx_k__stdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_ptr[2]) = ((void *)cdft1_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_ptr[3]) = ((void *)((char *)__pyx_k__stdtr));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_types, 2, 2, 1, 0, __pyx_k__stdtr, __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtr_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6906; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__stdtr, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6906; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_doc = __pyx_k_181;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_ptr[0]) = ((void *)cdft3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_ptr[1]) = ((void *)((char *)__pyx_k__stdtridf));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_ptr[2]) = ((void *)cdft3_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_ptr[3]) = ((void *)((char *)__pyx_k__stdtridf));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_types, 2, 2, 1, 0, __pyx_k__stdtridf, __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtridf_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6928; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__stdtridf, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6928; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_doc = __pyx_k_182;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_ptr[0]) = ((void *)cdft2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_ptr[1]) = ((void *)((char *)__pyx_k__stdtrit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_ptr[2]) = ((void *)cdft2_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_ptr[3]) = ((void *)((char *)__pyx_k__stdtrit));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_types, 2, 2, 1, 0, __pyx_k__stdtrit, __pyx_v_5scipy_7special_7_ufuncs_ufunc_stdtrit_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6950; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__stdtrit, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6950; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_doc = __pyx_k_183;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_ptr[0]) = ((void *)struve_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_ptr[1]) = ((void *)((char *)__pyx_k__struve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_ptr[2]) = ((void *)struve_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_ptr[3]) = ((void *)((char *)__pyx_k__struve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_types, 2, 2, 1, 0, __pyx_k__struve, __pyx_v_5scipy_7special_7_ufuncs_ufunc_struve_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6973; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__struve, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6973; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_doc = __pyx_k_184;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_ptr[0]) = ((void *)tandg);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_ptr[1]) = ((void *)((char *)__pyx_k__tandg));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_ptr[2]) = ((void *)tandg);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_ptr[3]) = ((void *)((char *)__pyx_k__tandg));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_types, 2, 1, 1, 0, __pyx_k__tandg, __pyx_v_5scipy_7special_7_ufuncs_ufunc_tandg_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6993; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__tandg, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 6993; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_doc = __pyx_k_24;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_ptr[0]) = ((void *)tukeylambdacdf);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_ptr[1]) = ((void *)((char *)__pyx_k__tklmbda));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_ptr[2]) = ((void *)tukeylambdacdf);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_ptr[3]) = ((void *)((char *)__pyx_k__tklmbda));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_types, 2, 2, 1, 0, __pyx_k__tklmbda, __pyx_v_5scipy_7special_7_ufuncs_ufunc_tklmbda_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7015; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__tklmbda, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7015; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_doc = __pyx_k_185;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_ptr[0]) = ((void *)y0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_ptr[1]) = ((void *)((char *)__pyx_k__y0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_ptr[2]) = ((void *)y0);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_ptr[3]) = ((void *)((char *)__pyx_k__y0));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_types, 2, 1, 1, 0, __pyx_k__y0, __pyx_v_5scipy_7special_7_ufuncs_ufunc_y0_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7035; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__y0, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7035; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_doc = __pyx_k_186;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_ptr[0]) = ((void *)y1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_ptr[1]) = ((void *)((char *)__pyx_k__y1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_ptr[2]) = ((void *)y1);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_ptr[3]) = ((void *)((char *)__pyx_k__y1));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_types, 2, 1, 1, 0, __pyx_k__y1, __pyx_v_5scipy_7special_7_ufuncs_ufunc_y1_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7055; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__y1, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7055; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_doc = __pyx_k_187;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_lf_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_id__As_ld_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[0]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[3]) = ((char)NPY_LONG);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[6]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[7]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[8]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[10]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types[11]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[0]) = ((void *)yn);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[1]) = ((void *)((char *)__pyx_k__yn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[2]) = ((void *)yn);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[3]) = ((void *)((char *)__pyx_k__yn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[4]) = ((void *)__pyx_f_5scipy_7special_7_legacy_yn_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[5]) = ((void *)((char *)__pyx_k__yn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[6]) = ((void *)__pyx_f_5scipy_7special_7_legacy_yn_unsafe);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[7]) = ((void *)((char *)__pyx_k__yn));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_types, 4, 2, 1, 0, __pyx_k__yn, __pyx_v_5scipy_7special_7_ufuncs_ufunc_yn_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7092; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__yn, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7092; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_doc = __pyx_k_188;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[10]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types[11]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[0]) = ((void *)yv);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[1]) = ((void *)((char *)__pyx_k__yv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[2]) = ((void *)cbesy_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[3]) = ((void *)((char *)__pyx_k__yv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[4]) = ((void *)yv);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[5]) = ((void *)((char *)__pyx_k__yv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[6]) = ((void *)cbesy_wrap);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[7]) = ((void *)((char *)__pyx_k__yv));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_types, 4, 2, 1, 0, __pyx_k__yv, __pyx_v_5scipy_7special_7_ufuncs_ufunc_yv_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7129; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__yv, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7129; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_doc = __pyx_k_189;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_fF_F);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_loops[2]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_loops[3]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_D_dD__As_dD_D);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[3]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[4]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[5]) = ((char)NPY_CFLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[6]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[7]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[8]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[9]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[10]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types[11]) = ((char)NPY_CDOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[0]) = ((void *)cbesy_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[1]) = ((void *)((char *)__pyx_k__yve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[2]) = ((void *)cbesy_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[3]) = ((void *)((char *)__pyx_k__yve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[4]) = ((void *)cbesy_wrap_e_real);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[5]) = ((void *)((char *)__pyx_k__yve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[6]) = ((void *)cbesy_wrap_e);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[7]) = ((void *)((char *)__pyx_k__yve));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[2]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_data[2]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[4]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_data[3]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_ptr[6]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_types, 4, 2, 1, 0, __pyx_k__yve, __pyx_v_5scipy_7special_7_ufuncs_ufunc_yve_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7166; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__yve, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7166; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_doc = __pyx_k_190;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_ff_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_dd__As_dd_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_types[2]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_types[4]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_types[5]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_ptr[0]) = ((void *)zeta);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_ptr[1]) = ((void *)((char *)__pyx_k__zeta));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_ptr[2]) = ((void *)zeta);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_ptr[3]) = ((void *)((char *)__pyx_k__zeta));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_types, 2, 2, 1, 0, __pyx_k__zeta, __pyx_v_5scipy_7special_7_ufuncs_ufunc_zeta_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7189; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__zeta, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7189; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_doc = __pyx_k_191;

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_loops[0]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_f_f);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_loops[1]) = ((PyUFuncGenericFunction)__pyx_f_5scipy_7special_7_ufuncs_loop_d_d__As_d_d);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_types[0]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_types[1]) = ((char)NPY_FLOAT);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_types[2]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_types[3]) = ((char)NPY_DOUBLE);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_ptr[0]) = ((void *)zetac);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_ptr[1]) = ((void *)((char *)__pyx_k__zetac));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_ptr[2]) = ((void *)zetac);

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_ptr[3]) = ((void *)((char *)__pyx_k__zetac));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_data[0]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_ptr[0]));

  
  (__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_data[1]) = (&(__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_ptr[2]));

  
  __pyx_t_1 = PyUFunc_FromFuncAndData(__pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_loops, __pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_data, __pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_types, 2, 1, 1, 0, __pyx_k__zetac, __pyx_v_5scipy_7special_7_ufuncs_ufunc_zetac_doc, 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7209; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__zetac, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7209; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_t_1 = PyCFunction_NewEx(&__pyx_mdef_5scipy_7special_7_ufuncs_1_errprint, NULL, __pyx_n_s_195); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7228; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s___errprint, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7228; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__jv); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7257; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__jn, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 7257; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_t_1 = PyDict_New(); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_1));
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s____test__, ((PyObject *)__pyx_t_1)) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(((PyObject *)__pyx_t_1)); __pyx_t_1 = 0;

  
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  if (__pyx_m) {
    __Pyx_AddTraceback("init scipy.special._ufuncs", __pyx_clineno, __pyx_lineno, __pyx_filename);
    Py_DECREF(__pyx_m); __pyx_m = 0;
  } else if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_ImportError, "init scipy.special._ufuncs");
  }
  __pyx_L0:;
  __Pyx_RefNannyFinishContext();
  #if PY_MAJOR_VERSION < 3
  return;
  #else
  return __pyx_m;
  #endif
}


#if CYTHON_REFNANNY
static __Pyx_RefNannyAPIStruct *__Pyx_RefNannyImportAPI(const char *modname) {
    PyObject *m = NULL, *p = NULL;
    void *r = NULL;
    m = PyImport_ImportModule((char *)modname);
    if (!m) goto end;
    p = PyObject_GetAttrString(m, (char *)"RefNannyAPI");
    if (!p) goto end;
    r = PyLong_AsVoidPtr(p);
end:
    Py_XDECREF(p);
    Py_XDECREF(m);
    return (__Pyx_RefNannyAPIStruct *)r;
}
#endif 

static PyObject *__Pyx_GetName(PyObject *dict, PyObject *name) {
    PyObject *result;
    result = PyObject_GetAttr(dict, name);
    if (!result) {
        if (dict != __pyx_b) {
            PyErr_Clear();
            result = PyObject_GetAttr(__pyx_b, name);
        }
        if (!result) {
            PyErr_SetObject(PyExc_NameError, name);
        }
    }
    return result;
}

static void __Pyx_RaiseDoubleKeywordsError(
    const char* func_name,
    PyObject* kw_name)
{
    PyErr_Format(PyExc_TypeError,
        #if PY_MAJOR_VERSION >= 3
        "%s() got multiple values for keyword argument '%U'", func_name, kw_name);
        #else
        "%s() got multiple values for keyword argument '%s'", func_name,
        PyString_AsString(kw_name));
        #endif
}

static int __Pyx_ParseOptionalKeywords(
    PyObject *kwds,
    PyObject **argnames[],
    PyObject *kwds2,
    PyObject *values[],
    Py_ssize_t num_pos_args,
    const char* function_name)
{
    PyObject *key = 0, *value = 0;
    Py_ssize_t pos = 0;
    PyObject*** name;
    PyObject*** first_kw_arg = argnames + num_pos_args;
    while (PyDict_Next(kwds, &pos, &key, &value)) {
        name = first_kw_arg;
        while (*name && (**name != key)) name++;
        if (*name) {
            values[name-argnames] = value;
            continue;
        }
        name = first_kw_arg;
        #if PY_MAJOR_VERSION < 3
        if (likely(PyString_CheckExact(key)) || likely(PyString_Check(key))) {
            while (*name) {
                if ((CYTHON_COMPILING_IN_PYPY || PyString_GET_SIZE(**name) == PyString_GET_SIZE(key))
                        && _PyString_Eq(**name, key)) {
                    values[name-argnames] = value;
                    break;
                }
                name++;
            }
            if (*name) continue;
            else {
                PyObject*** argname = argnames;
                while (argname != first_kw_arg) {
                    if ((**argname == key) || (
                            (CYTHON_COMPILING_IN_PYPY || PyString_GET_SIZE(**argname) == PyString_GET_SIZE(key))
                             && _PyString_Eq(**argname, key))) {
                        goto arg_passed_twice;
                    }
                    argname++;
                }
            }
        } else
        #endif
        if (likely(PyUnicode_Check(key))) {
            while (*name) {
                int cmp = (**name == key) ? 0 :
                #if !CYTHON_COMPILING_IN_PYPY && PY_MAJOR_VERSION >= 3
                    (PyUnicode_GET_SIZE(**name) != PyUnicode_GET_SIZE(key)) ? 1 :
                #endif
                    PyUnicode_Compare(**name, key);
                if (cmp < 0 && unlikely(PyErr_Occurred())) goto bad;
                if (cmp == 0) {
                    values[name-argnames] = value;
                    break;
                }
                name++;
            }
            if (*name) continue;
            else {
                PyObject*** argname = argnames;
                while (argname != first_kw_arg) {
                    int cmp = (**argname == key) ? 0 :
                    #if !CYTHON_COMPILING_IN_PYPY && PY_MAJOR_VERSION >= 3
                        (PyUnicode_GET_SIZE(**argname) != PyUnicode_GET_SIZE(key)) ? 1 :
                    #endif
                        PyUnicode_Compare(**argname, key);
                    if (cmp < 0 && unlikely(PyErr_Occurred())) goto bad;
                    if (cmp == 0) goto arg_passed_twice;
                    argname++;
                }
            }
        } else
            goto invalid_keyword_type;
        if (kwds2) {
            if (unlikely(PyDict_SetItem(kwds2, key, value))) goto bad;
        } else {
            goto invalid_keyword;
        }
    }
    return 0;
arg_passed_twice:
    __Pyx_RaiseDoubleKeywordsError(function_name, key);
    goto bad;
invalid_keyword_type:
    PyErr_Format(PyExc_TypeError,
        "%s() keywords must be strings", function_name);
    goto bad;
invalid_keyword:
    PyErr_Format(PyExc_TypeError,
    #if PY_MAJOR_VERSION < 3
        "%s() got an unexpected keyword argument '%s'",
        function_name, PyString_AsString(key));
    #else
        "%s() got an unexpected keyword argument '%U'",
        function_name, key);
    #endif
bad:
    return -1;
}

static void __Pyx_RaiseArgtupleInvalid(
    const char* func_name,
    int exact,
    Py_ssize_t num_min,
    Py_ssize_t num_max,
    Py_ssize_t num_found)
{
    Py_ssize_t num_expected;
    const char *more_or_less;
    if (num_found < num_min) {
        num_expected = num_min;
        more_or_less = "at least";
    } else {
        num_expected = num_max;
        more_or_less = "at most";
    }
    if (exact) {
        more_or_less = "exactly";
    }
    PyErr_Format(PyExc_TypeError,
                 "%s() takes %s %" CYTHON_FORMAT_SSIZE_T "d positional argument%s (%" CYTHON_FORMAT_SSIZE_T "d given)",
                 func_name, more_or_less, num_expected,
                 (num_expected == 1) ? "" : "s", num_found);
}

static CYTHON_INLINE long __Pyx_mod_long(long a, long b) {
    long r = a % b;
    r += ((r != 0) & ((r ^ b) < 0)) * b;
    return r;
}

static CYTHON_INLINE long __Pyx_div_long(long a, long b) {
    long q = a / b;
    long r = a - q*b;
    q -= ((r != 0) & ((r ^ b) < 0));
    return q;
}

static CYTHON_INLINE void __Pyx_ErrRestore(PyObject *type, PyObject *value, PyObject *tb) {
#if CYTHON_COMPILING_IN_CPYTHON
    PyObject *tmp_type, *tmp_value, *tmp_tb;
    PyThreadState *tstate = PyThreadState_GET();
    tmp_type = tstate->curexc_type;
    tmp_value = tstate->curexc_value;
    tmp_tb = tstate->curexc_traceback;
    tstate->curexc_type = type;
    tstate->curexc_value = value;
    tstate->curexc_traceback = tb;
    Py_XDECREF(tmp_type);
    Py_XDECREF(tmp_value);
    Py_XDECREF(tmp_tb);
#else
    PyErr_Restore(type, value, tb);
#endif
}
static CYTHON_INLINE void __Pyx_ErrFetch(PyObject **type, PyObject **value, PyObject **tb) {
#if CYTHON_COMPILING_IN_CPYTHON
    PyThreadState *tstate = PyThreadState_GET();
    *type = tstate->curexc_type;
    *value = tstate->curexc_value;
    *tb = tstate->curexc_traceback;
    tstate->curexc_type = 0;
    tstate->curexc_value = 0;
    tstate->curexc_traceback = 0;
#else
    PyErr_Fetch(type, value, tb);
#endif
}

#if PY_MAJOR_VERSION < 3
static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb,
                        CYTHON_UNUSED PyObject *cause) {
    Py_XINCREF(type);
    if (!value || value == Py_None)
        value = NULL;
    else
        Py_INCREF(value);
    if (!tb || tb == Py_None)
        tb = NULL;
    else {
        Py_INCREF(tb);
        if (!PyTraceBack_Check(tb)) {
            PyErr_SetString(PyExc_TypeError,
                "raise: arg 3 must be a traceback or None");
            goto raise_error;
        }
    }
    #if PY_VERSION_HEX < 0x02050000
    if (PyClass_Check(type)) {
    #else
    if (PyType_Check(type)) {
    #endif
#if CYTHON_COMPILING_IN_PYPY
        if (!value) {
            Py_INCREF(Py_None);
            value = Py_None;
        }
#endif
        PyErr_NormalizeException(&type, &value, &tb);
    } else {
        if (value) {
            PyErr_SetString(PyExc_TypeError,
                "instance exception may not have a separate value");
            goto raise_error;
        }
        value = type;
        #if PY_VERSION_HEX < 0x02050000
            if (PyInstance_Check(type)) {
                type = (PyObject*) ((PyInstanceObject*)type)->in_class;
                Py_INCREF(type);
            }
            else {
                type = 0;
                PyErr_SetString(PyExc_TypeError,
                    "raise: exception must be an old-style class or instance");
                goto raise_error;
            }
        #else
            type = (PyObject*) Py_TYPE(type);
            Py_INCREF(type);
            if (!PyType_IsSubtype((PyTypeObject *)type, (PyTypeObject *)PyExc_BaseException)) {
                PyErr_SetString(PyExc_TypeError,
                    "raise: exception class must be a subclass of BaseException");
                goto raise_error;
            }
        #endif
    }
    __Pyx_ErrRestore(type, value, tb);
    return;
raise_error:
    Py_XDECREF(value);
    Py_XDECREF(type);
    Py_XDECREF(tb);
    return;
}
#else 
static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb, PyObject *cause) {
    PyObject* owned_instance = NULL;
    if (tb == Py_None) {
        tb = 0;
    } else if (tb && !PyTraceBack_Check(tb)) {
        PyErr_SetString(PyExc_TypeError,
            "raise: arg 3 must be a traceback or None");
        goto bad;
    }
    if (value == Py_None)
        value = 0;
    if (PyExceptionInstance_Check(type)) {
        if (value) {
            PyErr_SetString(PyExc_TypeError,
                "instance exception may not have a separate value");
            goto bad;
        }
        value = type;
        type = (PyObject*) Py_TYPE(value);
    } else if (PyExceptionClass_Check(type)) {
        PyObject *args;
        if (!value)
            args = PyTuple_New(0);
        else if (PyTuple_Check(value)) {
            Py_INCREF(value);
            args = value;
        }
        else
            args = PyTuple_Pack(1, value);
        if (!args)
            goto bad;
        owned_instance = PyEval_CallObject(type, args);
        Py_DECREF(args);
        if (!owned_instance)
            goto bad;
        value = owned_instance;
        if (!PyExceptionInstance_Check(value)) {
            PyErr_Format(PyExc_TypeError,
                         "calling %R should have returned an instance of "
                         "BaseException, not %R",
                         type, Py_TYPE(value));
            goto bad;
        }
    } else {
        PyErr_SetString(PyExc_TypeError,
            "raise: exception class must be a subclass of BaseException");
        goto bad;
    }
    if (cause && cause != Py_None) {
        PyObject *fixed_cause;
        if (PyExceptionClass_Check(cause)) {
            fixed_cause = PyObject_CallObject(cause, NULL);
            if (fixed_cause == NULL)
                goto bad;
        }
        else if (PyExceptionInstance_Check(cause)) {
            fixed_cause = cause;
            Py_INCREF(fixed_cause);
        }
        else {
            PyErr_SetString(PyExc_TypeError,
                            "exception causes must derive from "
                            "BaseException");
            goto bad;
        }
        PyException_SetCause(value, fixed_cause);
    }
    PyErr_SetObject(type, value);
    if (tb) {
        PyThreadState *tstate = PyThreadState_GET();
        PyObject* tmp_tb = tstate->curexc_traceback;
        if (tb != tmp_tb) {
            Py_INCREF(tb);
            tstate->curexc_traceback = tb;
            Py_XDECREF(tmp_tb);
        }
    }
bad:
    Py_XDECREF(owned_instance);
    return;
}
#endif

static CYTHON_INLINE void __Pyx_RaiseTooManyValuesError(Py_ssize_t expected) {
    PyErr_Format(PyExc_ValueError,
                 "too many values to unpack (expected %" CYTHON_FORMAT_SSIZE_T "d)", expected);
}

static CYTHON_INLINE void __Pyx_RaiseNeedMoreValuesError(Py_ssize_t index) {
    PyErr_Format(PyExc_ValueError,
                 "need more than %" CYTHON_FORMAT_SSIZE_T "d value%s to unpack",
                 index, (index == 1) ? "" : "s");
}

static CYTHON_INLINE void __Pyx_RaiseNoneNotIterableError(void) {
    PyErr_SetString(PyExc_TypeError, "'NoneType' object is not iterable");
}

static CYTHON_INLINE int __Pyx_IterFinish(void) {
#if CYTHON_COMPILING_IN_CPYTHON
    PyThreadState *tstate = PyThreadState_GET();
    PyObject* exc_type = tstate->curexc_type;
    if (unlikely(exc_type)) {
        if (likely(exc_type == PyExc_StopIteration) || PyErr_GivenExceptionMatches(exc_type, PyExc_StopIteration)) {
            PyObject *exc_value, *exc_tb;
            exc_value = tstate->curexc_value;
            exc_tb = tstate->curexc_traceback;
            tstate->curexc_type = 0;
            tstate->curexc_value = 0;
            tstate->curexc_traceback = 0;
            Py_DECREF(exc_type);
            Py_XDECREF(exc_value);
            Py_XDECREF(exc_tb);
            return 0;
        } else {
            return -1;
        }
    }
    return 0;
#else
    if (unlikely(PyErr_Occurred())) {
        if (likely(PyErr_ExceptionMatches(PyExc_StopIteration))) {
            PyErr_Clear();
            return 0;
        } else {
            return -1;
        }
    }
    return 0;
#endif
}

static int __Pyx_IternextUnpackEndCheck(PyObject *retval, Py_ssize_t expected) {
    if (unlikely(retval)) {
        Py_DECREF(retval);
        __Pyx_RaiseTooManyValuesError(expected);
        return -1;
    } else {
        return __Pyx_IterFinish();
    }
    return 0;
}

static CYTHON_INLINE int __Pyx_TypeTest(PyObject *obj, PyTypeObject *type) {
    if (unlikely(!type)) {
        PyErr_Format(PyExc_SystemError, "Missing type object");
        return 0;
    }
    if (likely(PyObject_TypeCheck(obj, type)))
        return 1;
    PyErr_Format(PyExc_TypeError, "Cannot convert %.200s to %.200s",
                 Py_TYPE(obj)->tp_name, type->tp_name);
    return 0;
}

#if CYTHON_CCOMPLEX
  #ifdef __cplusplus
    static CYTHON_INLINE __pyx_t_double_complex __pyx_t_double_complex_from_parts(double x, double y) {
      return ::std::complex< double >(x, y);
    }
  #else
    static CYTHON_INLINE __pyx_t_double_complex __pyx_t_double_complex_from_parts(double x, double y) {
      return x + y*(__pyx_t_double_complex)_Complex_I;
    }
  #endif
#else
    static CYTHON_INLINE __pyx_t_double_complex __pyx_t_double_complex_from_parts(double x, double y) {
      __pyx_t_double_complex z;
      z.real = x;
      z.imag = y;
      return z;
    }
#endif

#if CYTHON_CCOMPLEX
#else
    static CYTHON_INLINE int __Pyx_c_eq(__pyx_t_double_complex a, __pyx_t_double_complex b) {
       return (a.real == b.real) && (a.imag == b.imag);
    }
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_sum(__pyx_t_double_complex a, __pyx_t_double_complex b) {
        __pyx_t_double_complex z;
        z.real = a.real + b.real;
        z.imag = a.imag + b.imag;
        return z;
    }
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_diff(__pyx_t_double_complex a, __pyx_t_double_complex b) {
        __pyx_t_double_complex z;
        z.real = a.real - b.real;
        z.imag = a.imag - b.imag;
        return z;
    }
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_prod(__pyx_t_double_complex a, __pyx_t_double_complex b) {
        __pyx_t_double_complex z;
        z.real = a.real * b.real - a.imag * b.imag;
        z.imag = a.real * b.imag + a.imag * b.real;
        return z;
    }
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_quot(__pyx_t_double_complex a, __pyx_t_double_complex b) {
        __pyx_t_double_complex z;
        double denom = b.real * b.real + b.imag * b.imag;
        z.real = (a.real * b.real + a.imag * b.imag) / denom;
        z.imag = (a.imag * b.real - a.real * b.imag) / denom;
        return z;
    }
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_neg(__pyx_t_double_complex a) {
        __pyx_t_double_complex z;
        z.real = -a.real;
        z.imag = -a.imag;
        return z;
    }
    static CYTHON_INLINE int __Pyx_c_is_zero(__pyx_t_double_complex a) {
       return (a.real == 0) && (a.imag == 0);
    }
    static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_conj(__pyx_t_double_complex a) {
        __pyx_t_double_complex z;
        z.real =  a.real;
        z.imag = -a.imag;
        return z;
    }
    #if 1
        static CYTHON_INLINE double __Pyx_c_abs(__pyx_t_double_complex z) {
          #if !defined(HAVE_HYPOT) || defined(_MSC_VER)
            return sqrt(z.real*z.real + z.imag*z.imag);
          #else
            return hypot(z.real, z.imag);
          #endif
        }
        static CYTHON_INLINE __pyx_t_double_complex __Pyx_c_pow(__pyx_t_double_complex a, __pyx_t_double_complex b) {
            __pyx_t_double_complex z;
            double r, lnr, theta, z_r, z_theta;
            if (b.imag == 0 && b.real == (int)b.real) {
                if (b.real < 0) {
                    double denom = a.real * a.real + a.imag * a.imag;
                    a.real = a.real / denom;
                    a.imag = -a.imag / denom;
                    b.real = -b.real;
                }
                switch ((int)b.real) {
                    case 0:
                        z.real = 1;
                        z.imag = 0;
                        return z;
                    case 1:
                        return a;
                    case 2:
                        z = __Pyx_c_prod(a, a);
                        return __Pyx_c_prod(a, a);
                    case 3:
                        z = __Pyx_c_prod(a, a);
                        return __Pyx_c_prod(z, a);
                    case 4:
                        z = __Pyx_c_prod(a, a);
                        return __Pyx_c_prod(z, z);
                }
            }
            if (a.imag == 0) {
                if (a.real == 0) {
                    return a;
                }
                r = a.real;
                theta = 0;
            } else {
                r = __Pyx_c_abs(a);
                theta = atan2(a.imag, a.real);
            }
            lnr = log(r);
            z_r = exp(lnr * b.real - theta * b.imag);
            z_theta = theta * b.real + lnr * b.imag;
            z.real = z_r * cos(z_theta);
            z.imag = z_r * sin(z_theta);
            return z;
        }
    #endif
#endif

static CYTHON_INLINE PyObject *__Pyx_PyInt_to_py_Py_intptr_t(Py_intptr_t val) {
    const Py_intptr_t neg_one = (Py_intptr_t)-1, const_zero = (Py_intptr_t)0;
    const int is_unsigned = const_zero < neg_one;
    if ((sizeof(Py_intptr_t) == sizeof(char))  ||
        (sizeof(Py_intptr_t) == sizeof(short))) {
        return PyInt_FromLong((long)val);
    } else if ((sizeof(Py_intptr_t) == sizeof(int)) ||
               (sizeof(Py_intptr_t) == sizeof(long))) {
        if (is_unsigned)
            return PyLong_FromUnsignedLong((unsigned long)val);
        else
            return PyInt_FromLong((long)val);
    } else if (sizeof(Py_intptr_t) == sizeof(PY_LONG_LONG)) {
        if (is_unsigned)
            return PyLong_FromUnsignedLongLong((unsigned PY_LONG_LONG)val);
        else
            return PyLong_FromLongLong((PY_LONG_LONG)val);
    } else {
        int one = 1; int little = (int)*(unsigned char *)&one;
        unsigned char *bytes = (unsigned char *)&val;
        return _PyLong_FromByteArray(bytes, sizeof(Py_intptr_t),
                                     little, !is_unsigned);
    }
}

static CYTHON_INLINE Py_intptr_t __Pyx_PyInt_from_py_Py_intptr_t(PyObject* x) {
    const Py_intptr_t neg_one = (Py_intptr_t)-1, const_zero = (Py_intptr_t)0;
    const int is_unsigned = const_zero < neg_one;
    if (sizeof(Py_intptr_t) == sizeof(char)) {
        if (is_unsigned)
            return (Py_intptr_t)__Pyx_PyInt_AsUnsignedChar(x);
        else
            return (Py_intptr_t)__Pyx_PyInt_AsSignedChar(x);
    } else if (sizeof(Py_intptr_t) == sizeof(short)) {
        if (is_unsigned)
            return (Py_intptr_t)__Pyx_PyInt_AsUnsignedShort(x);
        else
            return (Py_intptr_t)__Pyx_PyInt_AsSignedShort(x);
    } else if (sizeof(Py_intptr_t) == sizeof(int)) {
        if (is_unsigned)
            return (Py_intptr_t)__Pyx_PyInt_AsUnsignedInt(x);
        else
            return (Py_intptr_t)__Pyx_PyInt_AsSignedInt(x);
    } else if (sizeof(Py_intptr_t) == sizeof(long)) {
        if (is_unsigned)
            return (Py_intptr_t)__Pyx_PyInt_AsUnsignedLong(x);
        else
            return (Py_intptr_t)__Pyx_PyInt_AsSignedLong(x);
    } else if (sizeof(Py_intptr_t) == sizeof(PY_LONG_LONG)) {
        if (is_unsigned)
            return (Py_intptr_t)__Pyx_PyInt_AsUnsignedLongLong(x);
        else
            return (Py_intptr_t)__Pyx_PyInt_AsSignedLongLong(x);
    }  else {
        #if CYTHON_COMPILING_IN_PYPY && !defined(_PyLong_AsByteArray)
        PyErr_SetString(PyExc_RuntimeError,
                        "_PyLong_AsByteArray() not available in PyPy, cannot convert large numbers");
        #else
        Py_intptr_t val;
        PyObject *v = __Pyx_PyNumber_Int(x);
        #if PY_VERSION_HEX < 0x03000000
        if (likely(v) && !PyLong_Check(v)) {
            PyObject *tmp = v;
            v = PyNumber_Long(tmp);
            Py_DECREF(tmp);
        }
        #endif
        if (likely(v)) {
            int one = 1; int is_little = (int)*(unsigned char *)&one;
            unsigned char *bytes = (unsigned char *)&val;
            int ret = _PyLong_AsByteArray((PyLongObject *)v,
                                          bytes, sizeof(val),
                                          is_little, !is_unsigned);
            Py_DECREF(v);
            if (likely(!ret))
                return val;
        }
        #endif
        return (Py_intptr_t)-1;
    }
}

#if CYTHON_CCOMPLEX
  #ifdef __cplusplus
    static CYTHON_INLINE __pyx_t_float_complex __pyx_t_float_complex_from_parts(float x, float y) {
      return ::std::complex< float >(x, y);
    }
  #else
    static CYTHON_INLINE __pyx_t_float_complex __pyx_t_float_complex_from_parts(float x, float y) {
      return x + y*(__pyx_t_float_complex)_Complex_I;
    }
  #endif
#else
    static CYTHON_INLINE __pyx_t_float_complex __pyx_t_float_complex_from_parts(float x, float y) {
      __pyx_t_float_complex z;
      z.real = x;
      z.imag = y;
      return z;
    }
#endif

#if CYTHON_CCOMPLEX
#else
    static CYTHON_INLINE int __Pyx_c_eqf(__pyx_t_float_complex a, __pyx_t_float_complex b) {
       return (a.real == b.real) && (a.imag == b.imag);
    }
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_sumf(__pyx_t_float_complex a, __pyx_t_float_complex b) {
        __pyx_t_float_complex z;
        z.real = a.real + b.real;
        z.imag = a.imag + b.imag;
        return z;
    }
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_difff(__pyx_t_float_complex a, __pyx_t_float_complex b) {
        __pyx_t_float_complex z;
        z.real = a.real - b.real;
        z.imag = a.imag - b.imag;
        return z;
    }
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_prodf(__pyx_t_float_complex a, __pyx_t_float_complex b) {
        __pyx_t_float_complex z;
        z.real = a.real * b.real - a.imag * b.imag;
        z.imag = a.real * b.imag + a.imag * b.real;
        return z;
    }
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_quotf(__pyx_t_float_complex a, __pyx_t_float_complex b) {
        __pyx_t_float_complex z;
        float denom = b.real * b.real + b.imag * b.imag;
        z.real = (a.real * b.real + a.imag * b.imag) / denom;
        z.imag = (a.imag * b.real - a.real * b.imag) / denom;
        return z;
    }
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_negf(__pyx_t_float_complex a) {
        __pyx_t_float_complex z;
        z.real = -a.real;
        z.imag = -a.imag;
        return z;
    }
    static CYTHON_INLINE int __Pyx_c_is_zerof(__pyx_t_float_complex a) {
       return (a.real == 0) && (a.imag == 0);
    }
    static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_conjf(__pyx_t_float_complex a) {
        __pyx_t_float_complex z;
        z.real =  a.real;
        z.imag = -a.imag;
        return z;
    }
    #if 1
        static CYTHON_INLINE float __Pyx_c_absf(__pyx_t_float_complex z) {
          #if !defined(HAVE_HYPOT) || defined(_MSC_VER)
            return sqrtf(z.real*z.real + z.imag*z.imag);
          #else
            return hypotf(z.real, z.imag);
          #endif
        }
        static CYTHON_INLINE __pyx_t_float_complex __Pyx_c_powf(__pyx_t_float_complex a, __pyx_t_float_complex b) {
            __pyx_t_float_complex z;
            float r, lnr, theta, z_r, z_theta;
            if (b.imag == 0 && b.real == (int)b.real) {
                if (b.real < 0) {
                    float denom = a.real * a.real + a.imag * a.imag;
                    a.real = a.real / denom;
                    a.imag = -a.imag / denom;
                    b.real = -b.real;
                }
                switch ((int)b.real) {
                    case 0:
                        z.real = 1;
                        z.imag = 0;
                        return z;
                    case 1:
                        return a;
                    case 2:
                        z = __Pyx_c_prodf(a, a);
                        return __Pyx_c_prodf(a, a);
                    case 3:
                        z = __Pyx_c_prodf(a, a);
                        return __Pyx_c_prodf(z, a);
                    case 4:
                        z = __Pyx_c_prodf(a, a);
                        return __Pyx_c_prodf(z, z);
                }
            }
            if (a.imag == 0) {
                if (a.real == 0) {
                    return a;
                }
                r = a.real;
                theta = 0;
            } else {
                r = __Pyx_c_absf(a);
                theta = atan2f(a.imag, a.real);
            }
            lnr = logf(r);
            z_r = expf(lnr * b.real - theta * b.imag);
            z_theta = theta * b.real + lnr * b.imag;
            z.real = z_r * cosf(z_theta);
            z.imag = z_r * sinf(z_theta);
            return z;
        }
    #endif
#endif

static CYTHON_INLINE long __Pyx_pow_long(long b, long e) {
    long t = b;
    switch (e) {
        case 3:
            t *= b;
        case 2:
            t *= b;
        case 1:
            return t;
        case 0:
            return 1;
    }
    if (unlikely(e<0)) return 0;
    t = 1;
    while (likely(e)) {
        t *= (b * (e&1)) | ((~e)&1);    
        b *= b;
        e >>= 1;
    }
    return t;
}

static CYTHON_INLINE unsigned char __Pyx_PyInt_AsUnsignedChar(PyObject* x) {
    const unsigned char neg_one = (unsigned char)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(unsigned char) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(unsigned char)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to unsigned char" :
                    "value too large to convert to unsigned char");
            }
            return (unsigned char)-1;
        }
        return (unsigned char)val;
    }
    return (unsigned char)__Pyx_PyInt_AsUnsignedLong(x);
}

static CYTHON_INLINE unsigned short __Pyx_PyInt_AsUnsignedShort(PyObject* x) {
    const unsigned short neg_one = (unsigned short)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(unsigned short) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(unsigned short)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to unsigned short" :
                    "value too large to convert to unsigned short");
            }
            return (unsigned short)-1;
        }
        return (unsigned short)val;
    }
    return (unsigned short)__Pyx_PyInt_AsUnsignedLong(x);
}

static CYTHON_INLINE unsigned int __Pyx_PyInt_AsUnsignedInt(PyObject* x) {
    const unsigned int neg_one = (unsigned int)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(unsigned int) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(unsigned int)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to unsigned int" :
                    "value too large to convert to unsigned int");
            }
            return (unsigned int)-1;
        }
        return (unsigned int)val;
    }
    return (unsigned int)__Pyx_PyInt_AsUnsignedLong(x);
}

static CYTHON_INLINE char __Pyx_PyInt_AsChar(PyObject* x) {
    const char neg_one = (char)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(char) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(char)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to char" :
                    "value too large to convert to char");
            }
            return (char)-1;
        }
        return (char)val;
    }
    return (char)__Pyx_PyInt_AsLong(x);
}

static CYTHON_INLINE short __Pyx_PyInt_AsShort(PyObject* x) {
    const short neg_one = (short)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(short) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(short)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to short" :
                    "value too large to convert to short");
            }
            return (short)-1;
        }
        return (short)val;
    }
    return (short)__Pyx_PyInt_AsLong(x);
}

static CYTHON_INLINE int __Pyx_PyInt_AsInt(PyObject* x) {
    const int neg_one = (int)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(int) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(int)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to int" :
                    "value too large to convert to int");
            }
            return (int)-1;
        }
        return (int)val;
    }
    return (int)__Pyx_PyInt_AsLong(x);
}

static CYTHON_INLINE signed char __Pyx_PyInt_AsSignedChar(PyObject* x) {
    const signed char neg_one = (signed char)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(signed char) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(signed char)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to signed char" :
                    "value too large to convert to signed char");
            }
            return (signed char)-1;
        }
        return (signed char)val;
    }
    return (signed char)__Pyx_PyInt_AsSignedLong(x);
}

static CYTHON_INLINE signed short __Pyx_PyInt_AsSignedShort(PyObject* x) {
    const signed short neg_one = (signed short)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(signed short) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(signed short)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to signed short" :
                    "value too large to convert to signed short");
            }
            return (signed short)-1;
        }
        return (signed short)val;
    }
    return (signed short)__Pyx_PyInt_AsSignedLong(x);
}

static CYTHON_INLINE signed int __Pyx_PyInt_AsSignedInt(PyObject* x) {
    const signed int neg_one = (signed int)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(signed int) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(signed int)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to signed int" :
                    "value too large to convert to signed int");
            }
            return (signed int)-1;
        }
        return (signed int)val;
    }
    return (signed int)__Pyx_PyInt_AsSignedLong(x);
}

static CYTHON_INLINE int __Pyx_PyInt_AsLongDouble(PyObject* x) {
    const int neg_one = (int)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
    if (sizeof(int) < sizeof(long)) {
        long val = __Pyx_PyInt_AsLong(x);
        if (unlikely(val != (long)(int)val)) {
            if (!unlikely(val == -1 && PyErr_Occurred())) {
                PyErr_SetString(PyExc_OverflowError,
                    (is_unsigned && unlikely(val < 0)) ?
                    "can't convert negative value to int" :
                    "value too large to convert to int");
            }
            return (int)-1;
        }
        return (int)val;
    }
    return (int)__Pyx_PyInt_AsLong(x);
}

static CYTHON_INLINE unsigned long __Pyx_PyInt_AsUnsignedLong(PyObject* x) {
    const unsigned long neg_one = (unsigned long)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
#if PY_VERSION_HEX < 0x03000000
    if (likely(PyInt_Check(x))) {
        long val = PyInt_AS_LONG(x);
        if (is_unsigned && unlikely(val < 0)) {
            PyErr_SetString(PyExc_OverflowError,
                            "can't convert negative value to unsigned long");
            return (unsigned long)-1;
        }
        return (unsigned long)val;
    } else
#endif
    if (likely(PyLong_Check(x))) {
        if (is_unsigned) {
            if (unlikely(Py_SIZE(x) < 0)) {
                PyErr_SetString(PyExc_OverflowError,
                                "can't convert negative value to unsigned long");
                return (unsigned long)-1;
            }
            return (unsigned long)PyLong_AsUnsignedLong(x);
        } else {
            return (unsigned long)PyLong_AsLong(x);
        }
    } else {
        unsigned long val;
        PyObject *tmp = __Pyx_PyNumber_Int(x);
        if (!tmp) return (unsigned long)-1;
        val = __Pyx_PyInt_AsUnsignedLong(tmp);
        Py_DECREF(tmp);
        return val;
    }
}

static CYTHON_INLINE unsigned PY_LONG_LONG __Pyx_PyInt_AsUnsignedLongLong(PyObject* x) {
    const unsigned PY_LONG_LONG neg_one = (unsigned PY_LONG_LONG)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
#if PY_VERSION_HEX < 0x03000000
    if (likely(PyInt_Check(x))) {
        long val = PyInt_AS_LONG(x);
        if (is_unsigned && unlikely(val < 0)) {
            PyErr_SetString(PyExc_OverflowError,
                            "can't convert negative value to unsigned PY_LONG_LONG");
            return (unsigned PY_LONG_LONG)-1;
        }
        return (unsigned PY_LONG_LONG)val;
    } else
#endif
    if (likely(PyLong_Check(x))) {
        if (is_unsigned) {
            if (unlikely(Py_SIZE(x) < 0)) {
                PyErr_SetString(PyExc_OverflowError,
                                "can't convert negative value to unsigned PY_LONG_LONG");
                return (unsigned PY_LONG_LONG)-1;
            }
            return (unsigned PY_LONG_LONG)PyLong_AsUnsignedLongLong(x);
        } else {
            return (unsigned PY_LONG_LONG)PyLong_AsLongLong(x);
        }
    } else {
        unsigned PY_LONG_LONG val;
        PyObject *tmp = __Pyx_PyNumber_Int(x);
        if (!tmp) return (unsigned PY_LONG_LONG)-1;
        val = __Pyx_PyInt_AsUnsignedLongLong(tmp);
        Py_DECREF(tmp);
        return val;
    }
}

static CYTHON_INLINE long __Pyx_PyInt_AsLong(PyObject* x) {
    const long neg_one = (long)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
#if PY_VERSION_HEX < 0x03000000
    if (likely(PyInt_Check(x))) {
        long val = PyInt_AS_LONG(x);
        if (is_unsigned && unlikely(val < 0)) {
            PyErr_SetString(PyExc_OverflowError,
                            "can't convert negative value to long");
            return (long)-1;
        }
        return (long)val;
    } else
#endif
    if (likely(PyLong_Check(x))) {
        if (is_unsigned) {
            if (unlikely(Py_SIZE(x) < 0)) {
                PyErr_SetString(PyExc_OverflowError,
                                "can't convert negative value to long");
                return (long)-1;
            }
            return (long)PyLong_AsUnsignedLong(x);
        } else {
            return (long)PyLong_AsLong(x);
        }
    } else {
        long val;
        PyObject *tmp = __Pyx_PyNumber_Int(x);
        if (!tmp) return (long)-1;
        val = __Pyx_PyInt_AsLong(tmp);
        Py_DECREF(tmp);
        return val;
    }
}

static CYTHON_INLINE PY_LONG_LONG __Pyx_PyInt_AsLongLong(PyObject* x) {
    const PY_LONG_LONG neg_one = (PY_LONG_LONG)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
#if PY_VERSION_HEX < 0x03000000
    if (likely(PyInt_Check(x))) {
        long val = PyInt_AS_LONG(x);
        if (is_unsigned && unlikely(val < 0)) {
            PyErr_SetString(PyExc_OverflowError,
                            "can't convert negative value to PY_LONG_LONG");
            return (PY_LONG_LONG)-1;
        }
        return (PY_LONG_LONG)val;
    } else
#endif
    if (likely(PyLong_Check(x))) {
        if (is_unsigned) {
            if (unlikely(Py_SIZE(x) < 0)) {
                PyErr_SetString(PyExc_OverflowError,
                                "can't convert negative value to PY_LONG_LONG");
                return (PY_LONG_LONG)-1;
            }
            return (PY_LONG_LONG)PyLong_AsUnsignedLongLong(x);
        } else {
            return (PY_LONG_LONG)PyLong_AsLongLong(x);
        }
    } else {
        PY_LONG_LONG val;
        PyObject *tmp = __Pyx_PyNumber_Int(x);
        if (!tmp) return (PY_LONG_LONG)-1;
        val = __Pyx_PyInt_AsLongLong(tmp);
        Py_DECREF(tmp);
        return val;
    }
}

static CYTHON_INLINE signed long __Pyx_PyInt_AsSignedLong(PyObject* x) {
    const signed long neg_one = (signed long)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
#if PY_VERSION_HEX < 0x03000000
    if (likely(PyInt_Check(x))) {
        long val = PyInt_AS_LONG(x);
        if (is_unsigned && unlikely(val < 0)) {
            PyErr_SetString(PyExc_OverflowError,
                            "can't convert negative value to signed long");
            return (signed long)-1;
        }
        return (signed long)val;
    } else
#endif
    if (likely(PyLong_Check(x))) {
        if (is_unsigned) {
            if (unlikely(Py_SIZE(x) < 0)) {
                PyErr_SetString(PyExc_OverflowError,
                                "can't convert negative value to signed long");
                return (signed long)-1;
            }
            return (signed long)PyLong_AsUnsignedLong(x);
        } else {
            return (signed long)PyLong_AsLong(x);
        }
    } else {
        signed long val;
        PyObject *tmp = __Pyx_PyNumber_Int(x);
        if (!tmp) return (signed long)-1;
        val = __Pyx_PyInt_AsSignedLong(tmp);
        Py_DECREF(tmp);
        return val;
    }
}

static CYTHON_INLINE signed PY_LONG_LONG __Pyx_PyInt_AsSignedLongLong(PyObject* x) {
    const signed PY_LONG_LONG neg_one = (signed PY_LONG_LONG)-1, const_zero = 0;
    const int is_unsigned = neg_one > const_zero;
#if PY_VERSION_HEX < 0x03000000
    if (likely(PyInt_Check(x))) {
        long val = PyInt_AS_LONG(x);
        if (is_unsigned && unlikely(val < 0)) {
            PyErr_SetString(PyExc_OverflowError,
                            "can't convert negative value to signed PY_LONG_LONG");
            return (signed PY_LONG_LONG)-1;
        }
        return (signed PY_LONG_LONG)val;
    } else
#endif
    if (likely(PyLong_Check(x))) {
        if (is_unsigned) {
            if (unlikely(Py_SIZE(x) < 0)) {
                PyErr_SetString(PyExc_OverflowError,
                                "can't convert negative value to signed PY_LONG_LONG");
                return (signed PY_LONG_LONG)-1;
            }
            return (signed PY_LONG_LONG)PyLong_AsUnsignedLongLong(x);
        } else {
            return (signed PY_LONG_LONG)PyLong_AsLongLong(x);
        }
    } else {
        signed PY_LONG_LONG val;
        PyObject *tmp = __Pyx_PyNumber_Int(x);
        if (!tmp) return (signed PY_LONG_LONG)-1;
        val = __Pyx_PyInt_AsSignedLongLong(tmp);
        Py_DECREF(tmp);
        return val;
    }
}

static void __Pyx_WriteUnraisable(const char *name, CYTHON_UNUSED int clineno,
                                  CYTHON_UNUSED int lineno, CYTHON_UNUSED const char *filename) {
    PyObject *old_exc, *old_val, *old_tb;
    PyObject *ctx;
    __Pyx_ErrFetch(&old_exc, &old_val, &old_tb);
    #if PY_MAJOR_VERSION < 3
    ctx = PyString_FromString(name);
    #else
    ctx = PyUnicode_FromString(name);
    #endif
    __Pyx_ErrRestore(old_exc, old_val, old_tb);
    if (!ctx) {
        PyErr_WriteUnraisable(Py_None);
    } else {
        PyErr_WriteUnraisable(ctx);
        Py_DECREF(ctx);
    }
}

static int __Pyx_check_binary_version(void) {
    char ctversion[4], rtversion[4];
    PyOS_snprintf(ctversion, 4, "%d.%d", PY_MAJOR_VERSION, PY_MINOR_VERSION);
    PyOS_snprintf(rtversion, 4, "%s", Py_GetVersion());
    if (ctversion[0] != rtversion[0] || ctversion[2] != rtversion[2]) {
        char message[200];
        PyOS_snprintf(message, sizeof(message),
                      "compiletime version %s of module '%.100s' "
                      "does not match runtime version %s",
                      ctversion, __Pyx_MODULE_NAME, rtversion);
        #if PY_VERSION_HEX < 0x02050000
        return PyErr_Warn(NULL, message);
        #else
        return PyErr_WarnEx(NULL, message, 1);
        #endif
    }
    return 0;
}

#ifndef __PYX_HAVE_RT_ImportModule
#define __PYX_HAVE_RT_ImportModule
static PyObject *__Pyx_ImportModule(const char *name) {
    PyObject *py_name = 0;
    PyObject *py_module = 0;
    py_name = __Pyx_PyIdentifier_FromString(name);
    if (!py_name)
        goto bad;
    py_module = PyImport_Import(py_name);
    Py_DECREF(py_name);
    return py_module;
bad:
    Py_XDECREF(py_name);
    return 0;
}
#endif

#ifndef __PYX_HAVE_RT_ImportType
#define __PYX_HAVE_RT_ImportType
static PyTypeObject *__Pyx_ImportType(const char *module_name, const char *class_name,
    size_t size, int strict)
{
    PyObject *py_module = 0;
    PyObject *result = 0;
    PyObject *py_name = 0;
    char warning[200];
    py_module = __Pyx_ImportModule(module_name);
    if (!py_module)
        goto bad;
    py_name = __Pyx_PyIdentifier_FromString(class_name);
    if (!py_name)
        goto bad;
    result = PyObject_GetAttr(py_module, py_name);
    Py_DECREF(py_name);
    py_name = 0;
    Py_DECREF(py_module);
    py_module = 0;
    if (!result)
        goto bad;
    if (!PyType_Check(result)) {
        PyErr_Format(PyExc_TypeError,
            "%s.%s is not a type object",
            module_name, class_name);
        goto bad;
    }
    if (!strict && (size_t)((PyTypeObject *)result)->tp_basicsize > size) {
        PyOS_snprintf(warning, sizeof(warning),
            "%s.%s size changed, may indicate binary incompatibility",
            module_name, class_name);
        #if PY_VERSION_HEX < 0x02050000
        if (PyErr_Warn(NULL, warning) < 0) goto bad;
        #else
        if (PyErr_WarnEx(NULL, warning, 0) < 0) goto bad;
        #endif
    }
    else if ((size_t)((PyTypeObject *)result)->tp_basicsize != size) {
        PyErr_Format(PyExc_ValueError,
            "%s.%s has the wrong size, try recompiling",
            module_name, class_name);
        goto bad;
    }
    return (PyTypeObject *)result;
bad:
    Py_XDECREF(py_module);
    Py_XDECREF(result);
    return NULL;
}
#endif

static int __pyx_bisect_code_objects(__Pyx_CodeObjectCacheEntry* entries, int count, int code_line) {
    int start = 0, mid = 0, end = count - 1;
    if (end >= 0 && code_line > entries[end].code_line) {
        return count;
    }
    while (start < end) {
        mid = (start + end) / 2;
        if (code_line < entries[mid].code_line) {
            end = mid;
        } else if (code_line > entries[mid].code_line) {
             start = mid + 1;
        } else {
            return mid;
        }
    }
    if (code_line <= entries[mid].code_line) {
        return mid;
    } else {
        return mid + 1;
    }
}
static PyCodeObject *__pyx_find_code_object(int code_line) {
    PyCodeObject* code_object;
    int pos;
    if (unlikely(!code_line) || unlikely(!__pyx_code_cache.entries)) {
        return NULL;
    }
    pos = __pyx_bisect_code_objects(__pyx_code_cache.entries, __pyx_code_cache.count, code_line);
    if (unlikely(pos >= __pyx_code_cache.count) || unlikely(__pyx_code_cache.entries[pos].code_line != code_line)) {
        return NULL;
    }
    code_object = __pyx_code_cache.entries[pos].code_object;
    Py_INCREF(code_object);
    return code_object;
}
static void __pyx_insert_code_object(int code_line, PyCodeObject* code_object) {
    int pos, i;
    __Pyx_CodeObjectCacheEntry* entries = __pyx_code_cache.entries;
    if (unlikely(!code_line)) {
        return;
    }
    if (unlikely(!entries)) {
        entries = (__Pyx_CodeObjectCacheEntry*)PyMem_Malloc(64*sizeof(__Pyx_CodeObjectCacheEntry));
        if (likely(entries)) {
            __pyx_code_cache.entries = entries;
            __pyx_code_cache.max_count = 64;
            __pyx_code_cache.count = 1;
            entries[0].code_line = code_line;
            entries[0].code_object = code_object;
            Py_INCREF(code_object);
        }
        return;
    }
    pos = __pyx_bisect_code_objects(__pyx_code_cache.entries, __pyx_code_cache.count, code_line);
    if ((pos < __pyx_code_cache.count) && unlikely(__pyx_code_cache.entries[pos].code_line == code_line)) {
        PyCodeObject* tmp = entries[pos].code_object;
        entries[pos].code_object = code_object;
        Py_DECREF(tmp);
        return;
    }
    if (__pyx_code_cache.count == __pyx_code_cache.max_count) {
        int new_max = __pyx_code_cache.max_count + 64;
        entries = (__Pyx_CodeObjectCacheEntry*)PyMem_Realloc(
            __pyx_code_cache.entries, new_max*sizeof(__Pyx_CodeObjectCacheEntry));
        if (unlikely(!entries)) {
            return;
        }
        __pyx_code_cache.entries = entries;
        __pyx_code_cache.max_count = new_max;
    }
    for (i=__pyx_code_cache.count; i>pos; i--) {
        entries[i] = entries[i-1];
    }
    entries[pos].code_line = code_line;
    entries[pos].code_object = code_object;
    __pyx_code_cache.count++;
    Py_INCREF(code_object);
}

#include "compile.h"
#include "frameobject.h"
#include "traceback.h"
static PyCodeObject* __Pyx_CreateCodeObjectForTraceback(
            const char *funcname, int c_line,
            int py_line, const char *filename) {
    PyCodeObject *py_code = 0;
    PyObject *py_srcfile = 0;
    PyObject *py_funcname = 0;
    #if PY_MAJOR_VERSION < 3
    py_srcfile = PyString_FromString(filename);
    #else
    py_srcfile = PyUnicode_FromString(filename);
    #endif
    if (!py_srcfile) goto bad;
    if (c_line) {
        #if PY_MAJOR_VERSION < 3
        py_funcname = PyString_FromFormat( "%s (%s:%d)", funcname, __pyx_cfilenm, c_line);
        #else
        py_funcname = PyUnicode_FromFormat( "%s (%s:%d)", funcname, __pyx_cfilenm, c_line);
        #endif
    }
    else {
        #if PY_MAJOR_VERSION < 3
        py_funcname = PyString_FromString(funcname);
        #else
        py_funcname = PyUnicode_FromString(funcname);
        #endif
    }
    if (!py_funcname) goto bad;
    py_code = __Pyx_PyCode_New(
        0,            
        0,            
        0,            
        0,            
        0,            
        __pyx_empty_bytes, 
        __pyx_empty_tuple, 
        __pyx_empty_tuple, 
        __pyx_empty_tuple, 
        __pyx_empty_tuple, 
        __pyx_empty_tuple, 
        py_srcfile,   
        py_funcname,  
        py_line,      
        __pyx_empty_bytes  
    );
    Py_DECREF(py_srcfile);
    Py_DECREF(py_funcname);
    return py_code;
bad:
    Py_XDECREF(py_srcfile);
    Py_XDECREF(py_funcname);
    return NULL;
}
static void __Pyx_AddTraceback(const char *funcname, int c_line,
                               int py_line, const char *filename) {
    PyCodeObject *py_code = 0;
    PyObject *py_globals = 0;
    PyFrameObject *py_frame = 0;
    py_code = __pyx_find_code_object(c_line ? c_line : py_line);
    if (!py_code) {
        py_code = __Pyx_CreateCodeObjectForTraceback(
            funcname, c_line, py_line, filename);
        if (!py_code) goto bad;
        __pyx_insert_code_object(c_line ? c_line : py_line, py_code);
    }
    py_globals = PyModule_GetDict(__pyx_m);
    if (!py_globals) goto bad;
    py_frame = PyFrame_New(
        PyThreadState_GET(), 
        py_code,             
        py_globals,          
        0                    
    );
    if (!py_frame) goto bad;
    py_frame->f_lineno = py_line;
    PyTraceBack_Here(py_frame);
bad:
    Py_XDECREF(py_code);
    Py_XDECREF(py_frame);
}

static int __Pyx_InitStrings(__Pyx_StringTabEntry *t) {
    while (t->p) {
        #if PY_MAJOR_VERSION < 3
        if (t->is_unicode) {
            *t->p = PyUnicode_DecodeUTF8(t->s, t->n - 1, NULL);
        } else if (t->intern) {
            *t->p = PyString_InternFromString(t->s);
        } else {
            *t->p = PyString_FromStringAndSize(t->s, t->n - 1);
        }
        #else  
        if (t->is_unicode | t->is_str) {
            if (t->intern) {
                *t->p = PyUnicode_InternFromString(t->s);
            } else if (t->encoding) {
                *t->p = PyUnicode_Decode(t->s, t->n - 1, t->encoding, NULL);
            } else {
                *t->p = PyUnicode_FromStringAndSize(t->s, t->n - 1);
            }
        } else {
            *t->p = PyBytes_FromStringAndSize(t->s, t->n - 1);
        }
        #endif
        if (!*t->p)
            return -1;
        ++t;
    }
    return 0;
}




static CYTHON_INLINE int __Pyx_PyObject_IsTrue(PyObject* x) {
   int is_true = x == Py_True;
   if (is_true | (x == Py_False) | (x == Py_None)) return is_true;
   else return PyObject_IsTrue(x);
}

static CYTHON_INLINE PyObject* __Pyx_PyNumber_Int(PyObject* x) {
  PyNumberMethods *m;
  const char *name = NULL;
  PyObject *res = NULL;
#if PY_VERSION_HEX < 0x03000000
  if (PyInt_Check(x) || PyLong_Check(x))
#else
  if (PyLong_Check(x))
#endif
    return Py_INCREF(x), x;
  m = Py_TYPE(x)->tp_as_number;
#if PY_VERSION_HEX < 0x03000000
  if (m && m->nb_int) {
    name = "int";
    res = PyNumber_Int(x);
  }
  else if (m && m->nb_long) {
    name = "long";
    res = PyNumber_Long(x);
  }
#else
  if (m && m->nb_int) {
    name = "int";
    res = PyNumber_Long(x);
  }
#endif
  if (res) {
#if PY_VERSION_HEX < 0x03000000
    if (!PyInt_Check(res) && !PyLong_Check(res)) {
#else
    if (!PyLong_Check(res)) {
#endif
      PyErr_Format(PyExc_TypeError,
                   "__%s__ returned non-%s (type %.200s)",
                   name, name, Py_TYPE(res)->tp_name);
      Py_DECREF(res);
      return NULL;
    }
  }
  else if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_TypeError,
                    "an integer is required");
  }
  return res;
}

static CYTHON_INLINE Py_ssize_t __Pyx_PyIndex_AsSsize_t(PyObject* b) {
  Py_ssize_t ival;
  PyObject* x = PyNumber_Index(b);
  if (!x) return -1;
  ival = PyInt_AsSsize_t(x);
  Py_DECREF(x);
  return ival;
}

static CYTHON_INLINE PyObject * __Pyx_PyInt_FromSize_t(size_t ival) {
#if PY_VERSION_HEX < 0x02050000
   if (ival <= LONG_MAX)
       return PyInt_FromLong((long)ival);
   else {
       unsigned char *bytes = (unsigned char *) &ival;
       int one = 1; int little = (int)*(unsigned char*)&one;
       return _PyLong_FromByteArray(bytes, sizeof(size_t), little, 0);
   }
#else
   return PyInt_FromSize_t(ival);
#endif
}

static CYTHON_INLINE size_t __Pyx_PyInt_AsSize_t(PyObject* x) {
   unsigned PY_LONG_LONG val = __Pyx_PyInt_AsUnsignedLongLong(x);
   if (unlikely(val == (unsigned PY_LONG_LONG)-1 && PyErr_Occurred())) {
       return (size_t)-1;
   } else if (unlikely(val != (unsigned PY_LONG_LONG)(size_t)val)) {
       PyErr_SetString(PyExc_OverflowError,
                       "value too large to convert to size_t");
       return (size_t)-1;
   }
   return (size_t)val;
}


#endif 
