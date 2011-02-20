

#define PY_SSIZE_T_CLEAN
#include "Python.h"
#ifndef Py_PYTHON_H
    #error Python headers needed to compile C extensions, please install development version of Python.
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

#if PY_VERSION_HEX < 0x02040000
  #define METH_COEXIST 0
  #define PyDict_CheckExact(op) (Py_TYPE(op) == &PyDict_Type)
  #define PyDict_Contains(d,o)   PySequence_Contains(d,o)
#endif

#if PY_VERSION_HEX < 0x02050000
  typedef int Py_ssize_t;
  #define PY_SSIZE_T_MAX INT_MAX
  #define PY_SSIZE_T_MIN INT_MIN
  #define PY_FORMAT_SIZE_T ""
  #define PyInt_FromSsize_t(z) PyInt_FromLong(z)
  #define PyInt_AsSsize_t(o)   PyInt_AsLong(o)
  #define PyNumber_Index(o)    PyNumber_Int(o)
  #define PyIndex_Check(o)     PyNumber_Check(o)
  #define PyErr_WarnEx(category, message, stacklevel) PyErr_Warn(category, message)
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

#endif

#if PY_MAJOR_VERSION < 3
  #define __Pyx_BUILTIN_MODULE_NAME "__builtin__"
#else
  #define __Pyx_BUILTIN_MODULE_NAME "builtins"
#endif

#if PY_MAJOR_VERSION >= 3
  #define Py_TPFLAGS_CHECKTYPES 0
  #define Py_TPFLAGS_HAVE_INDEX 0
#endif

#if (PY_VERSION_HEX < 0x02060000) || (PY_MAJOR_VERSION >= 3)
  #define Py_TPFLAGS_HAVE_NEWBUFFER 0
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
  #define PySet_Check(obj)             PyObject_TypeCheck(obj, &PySet_Type)
  #define PyFrozenSet_Check(obj)       PyObject_TypeCheck(obj, &PyFrozenSet_Type)
#endif

#ifndef PySet_CheckExact
#  define PySet_CheckExact(obj)          (Py_TYPE(obj) == &PySet_Type)
#endif

#if PY_MAJOR_VERSION >= 3
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
  #define PyBoolObject PyLongObject
#endif


#if PY_MAJOR_VERSION >= 3
  #define __Pyx_PyNumber_Divide(x,y)         PyNumber_TrueDivide(x,y)
  #define __Pyx_PyNumber_InPlaceDivide(x,y)  PyNumber_InPlaceTrueDivide(x,y)
#else
  #define __Pyx_PyNumber_Divide(x,y)         PyNumber_Divide(x,y)
  #define __Pyx_PyNumber_InPlaceDivide(x,y)  PyNumber_InPlaceDivide(x,y)
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

#ifdef __cplusplus
#define __PYX_EXTERN_C extern "C"
#else
#define __PYX_EXTERN_C extern
#endif

#if defined(WIN32) || defined(MS_WINDOWS)
#define _USE_MATH_DEFINES
#endif
#include <math.h>
#define __PYX_HAVE_API__scipy__spatial__qhull
#include "numpy/ndarrayobject.h"
#include "stdio.h"
#include "stdlib.h"
#include "numpy/arrayobject.h"
#include "numpy/ufuncobject.h"
#include "math.h"
#include "qhull/src/qset.h"
#include "qhull/src/qhull.h"
#include "qhull_blas.h"


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
# elif defined(__ICC) || defined(__INTEL_COMPILER)
#   define CYTHON_UNUSED __attribute__ ((__unused__)) 
# else
#   define CYTHON_UNUSED 
# endif
#endif

typedef struct {PyObject **p; char *s; const long n; const char* encoding; const char is_unicode; const char is_str; const char intern; } __Pyx_StringTabEntry; 




#define __Pyx_PyBytes_FromUString(s) PyBytes_FromString((char*)s)
#define __Pyx_PyBytes_AsUString(s)   ((unsigned char*) PyBytes_AsString(s))

#define __Pyx_PyBool_FromLong(b) ((b) ? (Py_INCREF(Py_True), Py_True) : (Py_INCREF(Py_False), Py_False))
static CYTHON_INLINE int __Pyx_PyObject_IsTrue(PyObject*);
static CYTHON_INLINE PyObject* __Pyx_PyNumber_Int(PyObject* x);

static CYTHON_INLINE Py_ssize_t __Pyx_PyIndex_AsSsize_t(PyObject*);
static CYTHON_INLINE PyObject * __Pyx_PyInt_FromSize_t(size_t);
static CYTHON_INLINE size_t __Pyx_PyInt_AsSize_t(PyObject*);

#define __pyx_PyFloat_AsDouble(x) (PyFloat_CheckExact(x) ? PyFloat_AS_DOUBLE(x) : PyFloat_AsDouble(x))


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
  "qhull.pyx",
  "numpy.pxd",
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

typedef npy_intp __pyx_t_5numpy_intp_t;

typedef npy_uintp __pyx_t_5numpy_uintp_t;

typedef npy_ulong __pyx_t_5numpy_uint_t;

typedef npy_ulonglong __pyx_t_5numpy_ulong_t;

typedef npy_double __pyx_t_5numpy_float_t;

typedef npy_double __pyx_t_5numpy_double_t;

typedef npy_longdouble __pyx_t_5numpy_longdouble_t;

#if CYTHON_CCOMPLEX
  #ifdef __cplusplus
    typedef ::std::complex< float > __pyx_t_float_complex;
  #else
    typedef float _Complex __pyx_t_float_complex;
  #endif
#else
    typedef struct { float real, imag; } __pyx_t_float_complex;
#endif

#if CYTHON_CCOMPLEX
  #ifdef __cplusplus
    typedef ::std::complex< double > __pyx_t_double_complex;
  #else
    typedef double _Complex __pyx_t_double_complex;
  #endif
#else
    typedef struct { double real, imag; } __pyx_t_double_complex;
#endif



typedef npy_cfloat __pyx_t_5numpy_cfloat_t;

typedef npy_cdouble __pyx_t_5numpy_cdouble_t;

typedef npy_clongdouble __pyx_t_5numpy_clongdouble_t;

typedef npy_cdouble __pyx_t_5numpy_complex_t;



typedef struct {
  int ndim;
  int npoints;
  int nsimplex;
  double *points;
  int *vertices;
  int *neighbors;
  double *equations;
  double *transform;
  int *vertex_to_simplex;
  double paraboloid_scale;
  double paraboloid_shift;
  double *max_bound;
  double *min_bound;
} __pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t;



typedef struct {
  __pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *info;
  int index;
  int vertex;
  int vertex2;
  int triangle;
  int start_triangle;
  int start_index;
  int restart;
} __pyx_t_5scipy_7spatial_5qhull_RidgeIter2D_t;



struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D {
  PyObject_HEAD
  __pyx_t_5scipy_7spatial_5qhull_RidgeIter2D_t it;
  PyObject *delaunay;
  __pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t info;
};

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
  static __Pyx_RefNannyAPIStruct * __Pyx_RefNannyImportAPI(const char *modname) {
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
  #define __Pyx_RefNannySetupContext(name)           void *__pyx_refnanny = __Pyx_RefNanny->SetupContext((name), __LINE__, __FILE__)
  #define __Pyx_RefNannyFinishContext()           __Pyx_RefNanny->FinishContext(&__pyx_refnanny)
  #define __Pyx_INCREF(r) __Pyx_RefNanny->INCREF(__pyx_refnanny, (PyObject *)(r), __LINE__)
  #define __Pyx_DECREF(r) __Pyx_RefNanny->DECREF(__pyx_refnanny, (PyObject *)(r), __LINE__)
  #define __Pyx_GOTREF(r) __Pyx_RefNanny->GOTREF(__pyx_refnanny, (PyObject *)(r), __LINE__)
  #define __Pyx_GIVEREF(r) __Pyx_RefNanny->GIVEREF(__pyx_refnanny, (PyObject *)(r), __LINE__)
  #define __Pyx_XDECREF(r) do { if((r) != NULL) {__Pyx_DECREF(r);} } while(0)
#else
  #define __Pyx_RefNannySetupContext(name)
  #define __Pyx_RefNannyFinishContext()
  #define __Pyx_INCREF(r) Py_INCREF(r)
  #define __Pyx_DECREF(r) Py_DECREF(r)
  #define __Pyx_GOTREF(r)
  #define __Pyx_GIVEREF(r)
  #define __Pyx_XDECREF(r) Py_XDECREF(r)
#endif 
#define __Pyx_XGIVEREF(r) do { if((r) != NULL) {__Pyx_GIVEREF(r);} } while(0)
#define __Pyx_XGOTREF(r) do { if((r) != NULL) {__Pyx_GOTREF(r);} } while(0)

static PyObject *__Pyx_GetName(PyObject *dict, PyObject *name); 

static int __Pyx_ArgTypeTest(PyObject *obj, PyTypeObject *type, int none_allowed,
    const char *name, int exact); 


struct __Pyx_StructField_;

typedef struct {
  const char* name; 
  struct __Pyx_StructField_* fields;
  size_t size;     
  char typegroup; 
} __Pyx_TypeInfo;

typedef struct __Pyx_StructField_ {
  __Pyx_TypeInfo* type;
  const char* name;
  size_t offset;
} __Pyx_StructField;

typedef struct {
  __Pyx_StructField* field;
  size_t parent_offset;
} __Pyx_BufFmt_StackElem;


static CYTHON_INLINE int  __Pyx_GetBufferAndValidate(Py_buffer* buf, PyObject* obj, __Pyx_TypeInfo* dtype, int flags, int nd, int cast, __Pyx_BufFmt_StackElem* stack);
static CYTHON_INLINE void __Pyx_SafeReleaseBuffer(Py_buffer* info);

static CYTHON_INLINE int __Pyx_TypeTest(PyObject *obj, PyTypeObject *type); 

static void __Pyx_RaiseBufferFallbackError(void); 

static CYTHON_INLINE void __Pyx_RaiseNeedMoreValuesError(Py_ssize_t index);

static CYTHON_INLINE void __Pyx_RaiseTooManyValuesError(Py_ssize_t expected);

static PyObject *__Pyx_UnpackItem(PyObject *, Py_ssize_t index); 
static int __Pyx_EndUnpack(PyObject *, Py_ssize_t expected); 

static CYTHON_INLINE void __Pyx_ErrRestore(PyObject *type, PyObject *value, PyObject *tb); 
static CYTHON_INLINE void __Pyx_ErrFetch(PyObject **type, PyObject **value, PyObject **tb); 

static void __Pyx_RaiseDoubleKeywordsError(
    const char* func_name, PyObject* kw_name); 

static void __Pyx_RaiseArgtupleInvalid(const char* func_name, int exact,
    Py_ssize_t num_min, Py_ssize_t num_max, Py_ssize_t num_found); 

static int __Pyx_ParseOptionalKeywords(PyObject *kwds, PyObject **argnames[],     PyObject *kwds2, PyObject *values[], Py_ssize_t num_pos_args,     const char* function_name); 
static void __Pyx_RaiseBufferIndexError(int axis); 
#define __Pyx_BufPtrStrided1d(type, buf, i0, s0) (type)((char*)buf + i0 * s0)
#define __Pyx_BufPtrStrided2d(type, buf, i0, s0, i1, s1) (type)((char*)buf + i0 * s0 + i1 * s1)
#define __Pyx_BufPtrStrided3d(type, buf, i0, s0, i1, s1, i2, s2) (type)((char*)buf + i0 * s0 + i1 * s1 + i2 * s2)


static CYTHON_INLINE PyObject *__Pyx_GetItemInt_Generic(PyObject *o, PyObject* j) {
    PyObject *r;
    if (!j) return NULL;
    r = PyObject_GetItem(o, j);
    Py_DECREF(j);
    return r;
}


#define __Pyx_GetItemInt_List(o, i, size, to_py_func) (((size) <= sizeof(Py_ssize_t)) ? \
                                                    __Pyx_GetItemInt_List_Fast(o, i) : \
                                                    __Pyx_GetItemInt_Generic(o, to_py_func(i)))

static CYTHON_INLINE PyObject *__Pyx_GetItemInt_List_Fast(PyObject *o, Py_ssize_t i) {
    if (likely(o != Py_None)) {
        if (likely((0 <= i) & (i < PyList_GET_SIZE(o)))) {
            PyObject *r = PyList_GET_ITEM(o, i);
            Py_INCREF(r);
            return r;
        }
        else if ((-PyList_GET_SIZE(o) <= i) & (i < 0)) {
            PyObject *r = PyList_GET_ITEM(o, PyList_GET_SIZE(o) + i);
            Py_INCREF(r);
            return r;
        }
    }
    return __Pyx_GetItemInt_Generic(o, PyInt_FromSsize_t(i));
}

#define __Pyx_GetItemInt_Tuple(o, i, size, to_py_func) (((size) <= sizeof(Py_ssize_t)) ? \
                                                    __Pyx_GetItemInt_Tuple_Fast(o, i) : \
                                                    __Pyx_GetItemInt_Generic(o, to_py_func(i)))

static CYTHON_INLINE PyObject *__Pyx_GetItemInt_Tuple_Fast(PyObject *o, Py_ssize_t i) {
    if (likely(o != Py_None)) {
        if (likely((0 <= i) & (i < PyTuple_GET_SIZE(o)))) {
            PyObject *r = PyTuple_GET_ITEM(o, i);
            Py_INCREF(r);
            return r;
        }
        else if ((-PyTuple_GET_SIZE(o) <= i) & (i < 0)) {
            PyObject *r = PyTuple_GET_ITEM(o, PyTuple_GET_SIZE(o) + i);
            Py_INCREF(r);
            return r;
        }
    }
    return __Pyx_GetItemInt_Generic(o, PyInt_FromSsize_t(i));
}


#define __Pyx_GetItemInt(o, i, size, to_py_func) (((size) <= sizeof(Py_ssize_t)) ? \
                                                    __Pyx_GetItemInt_Fast(o, i) : \
                                                    __Pyx_GetItemInt_Generic(o, to_py_func(i)))

static CYTHON_INLINE PyObject *__Pyx_GetItemInt_Fast(PyObject *o, Py_ssize_t i) {
    PyObject *r;
    if (PyList_CheckExact(o) && ((0 <= i) & (i < PyList_GET_SIZE(o)))) {
        r = PyList_GET_ITEM(o, i);
        Py_INCREF(r);
    }
    else if (PyTuple_CheckExact(o) && ((0 <= i) & (i < PyTuple_GET_SIZE(o)))) {
        r = PyTuple_GET_ITEM(o, i);
        Py_INCREF(r);
    }
    else if (Py_TYPE(o)->tp_as_sequence && Py_TYPE(o)->tp_as_sequence->sq_item && (likely(i >= 0))) {
        r = PySequence_GetItem(o, i);
    }
    else {
        r = __Pyx_GetItemInt_Generic(o, PyInt_FromSsize_t(i));
    }
    return r;
}

static CYTHON_INLINE void __Pyx_RaiseNoneNotIterableError(void);

static void __Pyx_UnpackTupleError(PyObject *, Py_ssize_t index); 
#if PY_MAJOR_VERSION < 3
static int __Pyx_GetBuffer(PyObject *obj, Py_buffer *view, int flags);
static void __Pyx_ReleaseBuffer(Py_buffer *view);
#else
#define __Pyx_GetBuffer PyObject_GetBuffer
#define __Pyx_ReleaseBuffer PyBuffer_Release
#endif

Py_ssize_t __Pyx_zeros[] = {0, 0, 0};
Py_ssize_t __Pyx_minusones[] = {-1, -1, -1};

static PyObject *__Pyx_Import(PyObject *name, PyObject *from_list); 

static PyObject *__Pyx_CreateClass(PyObject *bases, PyObject *dict, PyObject *name, const char *modname); 

static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb); 

static CYTHON_INLINE PyObject *__Pyx_PyInt_to_py_flagT(flagT);

static CYTHON_INLINE PyObject *__Pyx_PyInt_to_py_Py_intptr_t(Py_intptr_t);

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
    
  #else
    #define __Pyx_c_is_zerof(z) ((z)==0)
    #define __Pyx_c_conjf(z)    (conjf(z))
    
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
    
  #else
    #define __Pyx_c_is_zero(z) ((z)==0)
    #define __Pyx_c_conj(z)    (conj(z))
    
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
    
#endif

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

static void __Pyx_WriteUnraisable(const char *name); 

static int __Pyx_ExportFunction(const char *name, void (*f)(void), const char *sig); 

static PyTypeObject *__Pyx_ImportType(const char *module_name, const char *class_name, long size, int strict);  

static PyObject *__Pyx_ImportModule(const char *name); 

static void __Pyx_AddTraceback(const char *funcname); 

static int __Pyx_InitStrings(__Pyx_StringTabEntry *t); 














static PyTypeObject *__pyx_ptype_5numpy_dtype = 0;
static PyTypeObject *__pyx_ptype_5numpy_flatiter = 0;
static PyTypeObject *__pyx_ptype_5numpy_broadcast = 0;
static PyTypeObject *__pyx_ptype_5numpy_ndarray = 0;
static PyTypeObject *__pyx_ptype_5numpy_ufunc = 0;
static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew1(PyObject *); 
static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew2(PyObject *, PyObject *); 
static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew3(PyObject *, PyObject *, PyObject *); 
static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew4(PyObject *, PyObject *, PyObject *, PyObject *); 
static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew5(PyObject *, PyObject *, PyObject *, PyObject *, PyObject *); 
static CYTHON_INLINE char *__pyx_f_5numpy__util_dtypestring(PyArray_Descr *, char *, char *, int *); 
static CYTHON_INLINE void __pyx_f_5numpy_set_array_base(PyArrayObject *, PyObject *); 
static CYTHON_INLINE PyObject *__pyx_f_5numpy_get_array_base(PyArrayObject *); 




static PyTypeObject *__pyx_ptype_5scipy_7spatial_5qhull_RidgeIter2D = 0;
static void __pyx_f_5scipy_7spatial_5qhull__get_delaunay_info(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, PyObject *, int, int); 
static int __pyx_f_5scipy_7spatial_5qhull__barycentric_inside(int, double *, double *, double *, double); 
static void __pyx_f_5scipy_7spatial_5qhull__barycentric_coordinate_single(int, double *, double *, double *, int); 
static void __pyx_f_5scipy_7spatial_5qhull__barycentric_coordinates(int, double *, double *, double *); 
static void __pyx_f_5scipy_7spatial_5qhull__lift_point(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double *); 
static double __pyx_f_5scipy_7spatial_5qhull__distplane(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, int, double *); 
static int __pyx_f_5scipy_7spatial_5qhull__is_point_fully_outside(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double); 
static int __pyx_f_5scipy_7spatial_5qhull__find_simplex_bruteforce(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double *, double); 
static int __pyx_f_5scipy_7spatial_5qhull__find_simplex_directed(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double *, int *, double); 
static int __pyx_f_5scipy_7spatial_5qhull__find_simplex(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double *, int *, double); 
static void __pyx_f_5scipy_7spatial_5qhull__RidgeIter2D_init(__pyx_t_5scipy_7spatial_5qhull_RidgeIter2D_t *, __pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, int); 
static void __pyx_f_5scipy_7spatial_5qhull__RidgeIter2D_next(__pyx_t_5scipy_7spatial_5qhull_RidgeIter2D_t *); 
static __Pyx_TypeInfo __Pyx_TypeInfo_nn___pyx_t_5numpy_double_t = { "double_t", NULL, sizeof(__pyx_t_5numpy_double_t), 'R' };
static __Pyx_TypeInfo __Pyx_TypeInfo_nn_npy_int = { "npy_int", NULL, sizeof(npy_int), 'I' };
#define __Pyx_MODULE_NAME "scipy.spatial.qhull"
int __pyx_module_is_main_scipy__spatial__qhull = 0;


static PyObject *__pyx_builtin_object;
static PyObject *__pyx_builtin_property;
static PyObject *__pyx_builtin_ValueError;
static PyObject *__pyx_builtin_RuntimeError;
static PyObject *__pyx_builtin_xrange;
static PyObject *__pyx_builtin_StopIteration;
static PyObject *__pyx_builtin_range;
static char __pyx_k_1[] = "qhull d Qz Qbb Qt";
static char __pyx_k_2[] = "No points to triangulate";
static char __pyx_k_3[] = "Need at least 2-D data to triangulate";
static char __pyx_k_4[] = "Qhull error";
static char __pyx_k_5[] = "_qhull_get_facet_array";
static char __pyx_k_6[] = "qhull: did not free %d bytes (%d pieces)";
static char __pyx_k_7[] = "non-simplical facet encountered";
static char __pyx_k_8[] = "RidgeIter2D supports only 2-D";
static char __pyx_k_9[] = "_get_barycentric_transforms";
static char __pyx_k_11[] = "wrong dimensionality in xi";
static char __pyx_k_12[] = "xi has different dimensionality than triangulation";
static char __pyx_k_13[] = "ndarray is not C contiguous";
static char __pyx_k_14[] = "ndarray is not Fortran contiguous";
static char __pyx_k_15[] = "Non-native byte order not supported";
static char __pyx_k_16[] = "unknown dtype code in numpy.pxd (%d)";
static char __pyx_k_17[] = "Format string allocated too short, see comment in numpy.pxd";
static char __pyx_k_18[] = "Format string allocated too short.";
static char __pyx_k_19[] = "\nWrappers for Qhull triangulation, plus some additional N-D geometry utilities\n\n.. versionadded:: 0.9\n\n";
static char __pyx_k_20[] = "\n    Delaunay(points)\n\n    Delaunay tesselation in N dimensions\n\n    .. versionadded:: 0.9\n\n    Parameters\n    ----------\n    points : ndarray of floats, shape (npoints, ndim)\n        Coordinates of points to triangulate\n\n    Attributes\n    ----------\n    points : ndarray of double, shape (npoints, ndim)\n        Points in the triangulation\n    vertices : ndarray of ints, shape (nsimplex, ndim+1)\n        Indices of vertices forming simplices in the triangulation\n    neighbors : ndarray of ints, shape (nsimplex, ndim+1)\n        Indices of neighbor simplices for each simplex.\n        The kth neighbor is opposite to the kth vertex.\n        For simplices at the boundary, -1 denotes no neighbor.\n    equations : ndarray of double, shape (nsimplex, ndim+2)\n        [normal, offset] forming the hyperplane equation of the facet\n        on the paraboloid. (See [Qhull]_ documentation for more.)\n    paraboloid_scale, paraboloid_shift : float\n        Scale and shift for the extra paraboloid dimension.\n        (See [Qhull]_ documentation for more.)\n    transform : ndarray of double, shape (nsimplex, ndim+1, ndim)\n        Affine transform from ``x`` to the barycentric coordinates ``c``.\n        This is defined by::\n\n            T c = x - r\n\n        At vertex ``j``, ``c_j = 1`` and the other coordinates zero.\n\n        For simplex ``i``, ``transform[i,:ndim,:ndim]`` contains\n        inverse of the matrix ``T``, and ``transform[i,ndim,:]``\n        contains the vector ``r``.\n    vertex_to_simplex : ndarray of int, shape (npoints,)\n        Lookup array, from a vertex, to some simplex which it is a part of.\n    convex_hull : ndarray of int, shape (nfaces, ndim)\n        Vertices of facets forming the convex hull of the point set.\n        The array contains the indices of the points belonging to\n        the (N-1)-dimensional facets that form the convex hull\n        of the triangulation.\n\n    Notes\n    -----\n    The tesselation is computed usi""ng the Qhull libary [Qhull]_.\n\n    References\n    ----------\n\n    .. [Qhull] http://www.qhull.org/\n\n    ";
static char __pyx_k_21[] = "_construct_delaunay (line 134)";
static char __pyx_k_22[] = "_qhull_get_facet_array (line 197)";
static char __pyx_k_23[] = "_get_barycentric_transforms (line 278)";
static char __pyx_k_24[] = "Delaunay.transform (line 934)";
static char __pyx_k_25[] = "Delaunay.vertex_to_simplex (line 957)";
static char __pyx_k_26[] = "Delaunay.convex_hull (line 987)";
static char __pyx_k_27[] = "Delaunay.find_simplex (line 1033)";
static char __pyx_k_28[] = "Delaunay.plane_distance (line 1103)";
static char __pyx_k_29[] = "Delaunay.lift_points (line 1136)";
static char __pyx_k_30[] = "tsearch (line 1151)";
static char __pyx_k__B[] = "B";
static char __pyx_k__H[] = "H";
static char __pyx_k__I[] = "I";
static char __pyx_k__L[] = "L";
static char __pyx_k__O[] = "O";
static char __pyx_k__Q[] = "Q";
static char __pyx_k__b[] = "b";
static char __pyx_k__d[] = "d";
static char __pyx_k__e[] = "e";
static char __pyx_k__f[] = "f";
static char __pyx_k__g[] = "g";
static char __pyx_k__h[] = "h";
static char __pyx_k__i[] = "i";
static char __pyx_k__l[] = "l";
static char __pyx_k__p[] = "p";
static char __pyx_k__q[] = "q";
static char __pyx_k__x[] = "x";
static char __pyx_k__Zd[] = "Zd";
static char __pyx_k__Zf[] = "Zf";
static char __pyx_k__Zg[] = "Zg";
static char __pyx_k__id[] = "id";
static char __pyx_k__it[] = "it";
static char __pyx_k__np[] = "np";
static char __pyx_k__xi[] = "xi";
static char __pyx_k__buf[] = "buf";
static char __pyx_k__eps[] = "eps";
static char __pyx_k__max[] = "max";
static char __pyx_k__min[] = "min";
static char __pyx_k__nan[] = "nan";
static char __pyx_k__obj[] = "obj";
static char __pyx_k__sum[] = "sum";
static char __pyx_k__tri[] = "tri";
static char __pyx_k__Lock[] = "Lock";
static char __pyx_k__axis[] = "axis";
static char __pyx_k__base[] = "base";
static char __pyx_k__data[] = "data";
static char __pyx_k__fill[] = "fill";
static char __pyx_k__info[] = "info";
static char __pyx_k__intc[] = "intc";
static char __pyx_k__ndim[] = "ndim";
static char __pyx_k__next[] = "next";
static char __pyx_k__prod[] = "prod";
static char __pyx_k__self[] = "self";
static char __pyx_k__descr[] = "descr";
static char __pyx_k__dtype[] = "dtype";
static char __pyx_k__empty[] = "empty";
static char __pyx_k__finfo[] = "finfo";
static char __pyx_k__index[] = "index";
static char __pyx_k__names[] = "names";
static char __pyx_k__numpy[] = "numpy";
static char __pyx_k__point[] = "point";
static char __pyx_k__range[] = "range";
static char __pyx_k__shape[] = "shape";
static char __pyx_k__zeros[] = "zeros";
static char __pyx_k__astype[] = "astype";
static char __pyx_k__double[] = "double";
static char __pyx_k__fields[] = "fields";
static char __pyx_k__format[] = "format";
static char __pyx_k__normal[] = "normal";
static char __pyx_k__object[] = "object";
static char __pyx_k__offset[] = "offset";
static char __pyx_k__points[] = "points";
static char __pyx_k__resize[] = "resize";
static char __pyx_k__vertex[] = "vertex";
static char __pyx_k__xrange[] = "xrange";
static char __pyx_k____all__[] = "__all__";
static char __pyx_k__acquire[] = "acquire";
static char __pyx_k__ivertex[] = "ivertex";
static char __pyx_k__npoints[] = "npoints";
static char __pyx_k__release[] = "release";
static char __pyx_k__reshape[] = "reshape";
static char __pyx_k__restart[] = "restart";
static char __pyx_k__strides[] = "strides";
static char __pyx_k__tsearch[] = "tsearch";
static char __pyx_k__vertex2[] = "vertex2";
static char __pyx_k__Delaunay[] = "Delaunay";
static char __pyx_k____init__[] = "__init__";
static char __pyx_k____main__[] = "__main__";
static char __pyx_k____test__[] = "__test__";
static char __pyx_k__delaunay[] = "delaunay";
static char __pyx_k__facet_id[] = "facet_id";
static char __pyx_k__itemsize[] = "itemsize";
static char __pyx_k__last_low[] = "last_low";
static char __pyx_k__nsimplex[] = "nsimplex";
static char __pyx_k__property[] = "property";
static char __pyx_k__readonly[] = "readonly";
static char __pyx_k__triangle[] = "triangle";
static char __pyx_k__type_num[] = "type_num";
static char __pyx_k__vertices[] = "vertices";
static char __pyx_k__NOerrexit[] = "NOerrexit";
static char __pyx_k__SCALElast[] = "SCALElast";
static char __pyx_k__byteorder[] = "byteorder";
static char __pyx_k__equations[] = "equations";
static char __pyx_k__last_high[] = "last_high";
static char __pyx_k__max_bound[] = "max_bound";
static char __pyx_k__min_bound[] = "min_bound";
static char __pyx_k__neighbors[] = "neighbors";
static char __pyx_k__numpoints[] = "numpoints";
static char __pyx_k__threading[] = "threading";
static char __pyx_k__transform[] = "transform";
static char __pyx_k__ValueError[] = "ValueError";
static char __pyx_k___transform[] = "_transform";
static char __pyx_k__asanyarray[] = "asanyarray";
static char __pyx_k__bruteforce[] = "bruteforce";
static char __pyx_k__facet_list[] = "facet_list";
static char __pyx_k__simplicial[] = "simplicial";
static char __pyx_k__suboffsets[] = "suboffsets";
static char __pyx_k___qhull_lock[] = "_qhull_lock";
static char __pyx_k__convex_hull[] = "convex_hull";
static char __pyx_k__lift_points[] = "lift_points";
static char __pyx_k__start_index[] = "start_index";
static char __pyx_k__RuntimeError[] = "RuntimeError";
static char __pyx_k__find_simplex[] = "find_simplex";
static char __pyx_k__last_newhigh[] = "last_newhigh";
static char __pyx_k__StopIteration[] = "StopIteration";
static char __pyx_k__upperdelaunay[] = "upperdelaunay";
static char __pyx_k__plane_distance[] = "plane_distance";
static char __pyx_k__start_triangle[] = "start_triangle";
static char __pyx_k__paraboloid_scale[] = "paraboloid_scale";
static char __pyx_k__paraboloid_shift[] = "paraboloid_shift";
static char __pyx_k__ascontiguousarray[] = "ascontiguousarray";
static char __pyx_k__vertex_to_simplex[] = "vertex_to_simplex";
static char __pyx_k___vertex_to_simplex[] = "_vertex_to_simplex";
static char __pyx_k___construct_delaunay[] = "_construct_delaunay";
static PyObject *__pyx_kp_s_11;
static PyObject *__pyx_kp_s_12;
static PyObject *__pyx_kp_u_13;
static PyObject *__pyx_kp_u_14;
static PyObject *__pyx_kp_u_15;
static PyObject *__pyx_kp_u_16;
static PyObject *__pyx_kp_u_17;
static PyObject *__pyx_kp_u_18;
static PyObject *__pyx_kp_s_2;
static PyObject *__pyx_kp_s_20;
static PyObject *__pyx_kp_u_21;
static PyObject *__pyx_kp_u_22;
static PyObject *__pyx_kp_u_23;
static PyObject *__pyx_kp_u_24;
static PyObject *__pyx_kp_u_25;
static PyObject *__pyx_kp_u_26;
static PyObject *__pyx_kp_u_27;
static PyObject *__pyx_kp_u_28;
static PyObject *__pyx_kp_u_29;
static PyObject *__pyx_kp_s_3;
static PyObject *__pyx_kp_u_30;
static PyObject *__pyx_kp_s_4;
static PyObject *__pyx_n_s_5;
static PyObject *__pyx_kp_s_6;
static PyObject *__pyx_kp_s_7;
static PyObject *__pyx_kp_s_8;
static PyObject *__pyx_n_s_9;
static PyObject *__pyx_n_s__Delaunay;
static PyObject *__pyx_n_s__Lock;
static PyObject *__pyx_n_s__NOerrexit;
static PyObject *__pyx_n_s__RuntimeError;
static PyObject *__pyx_n_s__SCALElast;
static PyObject *__pyx_n_s__StopIteration;
static PyObject *__pyx_n_s__ValueError;
static PyObject *__pyx_n_s____all__;
static PyObject *__pyx_n_s____init__;
static PyObject *__pyx_n_s____main__;
static PyObject *__pyx_n_s____test__;
static PyObject *__pyx_n_s___construct_delaunay;
static PyObject *__pyx_n_s___qhull_lock;
static PyObject *__pyx_n_s___transform;
static PyObject *__pyx_n_s___vertex_to_simplex;
static PyObject *__pyx_n_s__acquire;
static PyObject *__pyx_n_s__asanyarray;
static PyObject *__pyx_n_s__ascontiguousarray;
static PyObject *__pyx_n_s__astype;
static PyObject *__pyx_n_s__axis;
static PyObject *__pyx_n_s__base;
static PyObject *__pyx_n_s__bruteforce;
static PyObject *__pyx_n_s__buf;
static PyObject *__pyx_n_s__byteorder;
static PyObject *__pyx_n_s__convex_hull;
static PyObject *__pyx_n_s__data;
static PyObject *__pyx_n_s__delaunay;
static PyObject *__pyx_n_s__descr;
static PyObject *__pyx_n_s__double;
static PyObject *__pyx_n_s__dtype;
static PyObject *__pyx_n_s__e;
static PyObject *__pyx_n_s__empty;
static PyObject *__pyx_n_s__eps;
static PyObject *__pyx_n_s__equations;
static PyObject *__pyx_n_s__facet_id;
static PyObject *__pyx_n_s__facet_list;
static PyObject *__pyx_n_s__fields;
static PyObject *__pyx_n_s__fill;
static PyObject *__pyx_n_s__find_simplex;
static PyObject *__pyx_n_s__finfo;
static PyObject *__pyx_n_s__format;
static PyObject *__pyx_n_s__id;
static PyObject *__pyx_n_s__index;
static PyObject *__pyx_n_s__info;
static PyObject *__pyx_n_s__intc;
static PyObject *__pyx_n_s__it;
static PyObject *__pyx_n_s__itemsize;
static PyObject *__pyx_n_s__ivertex;
static PyObject *__pyx_n_s__last_high;
static PyObject *__pyx_n_s__last_low;
static PyObject *__pyx_n_s__last_newhigh;
static PyObject *__pyx_n_s__lift_points;
static PyObject *__pyx_n_s__max;
static PyObject *__pyx_n_s__max_bound;
static PyObject *__pyx_n_s__min;
static PyObject *__pyx_n_s__min_bound;
static PyObject *__pyx_n_s__names;
static PyObject *__pyx_n_s__nan;
static PyObject *__pyx_n_s__ndim;
static PyObject *__pyx_n_s__neighbors;
static PyObject *__pyx_n_s__next;
static PyObject *__pyx_n_s__normal;
static PyObject *__pyx_n_s__np;
static PyObject *__pyx_n_s__npoints;
static PyObject *__pyx_n_s__nsimplex;
static PyObject *__pyx_n_s__numpoints;
static PyObject *__pyx_n_s__numpy;
static PyObject *__pyx_n_s__obj;
static PyObject *__pyx_n_s__object;
static PyObject *__pyx_n_s__offset;
static PyObject *__pyx_n_s__p;
static PyObject *__pyx_n_s__paraboloid_scale;
static PyObject *__pyx_n_s__paraboloid_shift;
static PyObject *__pyx_n_s__plane_distance;
static PyObject *__pyx_n_s__point;
static PyObject *__pyx_n_s__points;
static PyObject *__pyx_n_s__prod;
static PyObject *__pyx_n_s__property;
static PyObject *__pyx_n_s__range;
static PyObject *__pyx_n_s__readonly;
static PyObject *__pyx_n_s__release;
static PyObject *__pyx_n_s__reshape;
static PyObject *__pyx_n_s__resize;
static PyObject *__pyx_n_s__restart;
static PyObject *__pyx_n_s__self;
static PyObject *__pyx_n_s__shape;
static PyObject *__pyx_n_s__simplicial;
static PyObject *__pyx_n_s__start_index;
static PyObject *__pyx_n_s__start_triangle;
static PyObject *__pyx_n_s__strides;
static PyObject *__pyx_n_s__suboffsets;
static PyObject *__pyx_n_s__sum;
static PyObject *__pyx_n_s__threading;
static PyObject *__pyx_n_s__transform;
static PyObject *__pyx_n_s__tri;
static PyObject *__pyx_n_s__triangle;
static PyObject *__pyx_n_s__tsearch;
static PyObject *__pyx_n_s__type_num;
static PyObject *__pyx_n_s__upperdelaunay;
static PyObject *__pyx_n_s__vertex;
static PyObject *__pyx_n_s__vertex2;
static PyObject *__pyx_n_s__vertex_to_simplex;
static PyObject *__pyx_n_s__vertices;
static PyObject *__pyx_n_s__x;
static PyObject *__pyx_n_s__xi;
static PyObject *__pyx_n_s__xrange;
static PyObject *__pyx_n_s__zeros;
static PyObject *__pyx_int_0;
static PyObject *__pyx_int_1;
static PyObject *__pyx_int_2;
static PyObject *__pyx_int_neg_1;
static PyObject *__pyx_int_10;
static PyObject *__pyx_int_15;
static PyObject *__pyx_k_10;



static PyObject *__pyx_pf_5scipy_7spatial_5qhull__construct_delaunay(PyObject *__pyx_self, PyObject *__pyx_v_points); 
static char __pyx_doc_5scipy_7spatial_5qhull__construct_delaunay[] = "\n    Perform Delaunay triangulation of the given set of points.\n\n    ";
static PyObject *__pyx_pf_5scipy_7spatial_5qhull__construct_delaunay(PyObject *__pyx_self, PyObject *__pyx_v_points) {
  char *__pyx_v_options;
  int __pyx_v_curlong;
  int __pyx_v_totlong;
  int __pyx_v_dim;
  int __pyx_v_numpoints;
  int __pyx_v_exitcode;
  double __pyx_v_paraboloid_scale;
  double __pyx_v_paraboloid_shift;
  PyObject *__pyx_v_vertices;
  PyObject *__pyx_v_neighbors;
  PyObject *__pyx_v_equations;
  Py_buffer __pyx_bstruct_points;
  Py_ssize_t __pyx_bstride_0_points = 0;
  Py_ssize_t __pyx_bstride_1_points = 0;
  Py_ssize_t __pyx_bshape_0_points = 0;
  Py_ssize_t __pyx_bshape_1_points = 0;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  PyArrayObject *__pyx_t_4 = NULL;
  int __pyx_t_5;
  PyObject *__pyx_t_6 = NULL;
  PyObject *__pyx_t_7 = NULL;
  PyObject *__pyx_t_8 = NULL;
  int __pyx_t_9;
  realT __pyx_t_10;
  PyObject *__pyx_t_11 = NULL;
  PyObject *__pyx_t_12 = NULL;
  int __pyx_t_13;
  int __pyx_t_14;
  __Pyx_RefNannySetupContext("_construct_delaunay");
  __pyx_self = __pyx_self;
  __Pyx_INCREF((PyObject *)__pyx_v_points);
  __pyx_v_vertices = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_neighbors = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_equations = Py_None; __Pyx_INCREF(Py_None);
  __pyx_bstruct_points.buf = NULL;
  if (unlikely(!__Pyx_ArgTypeTest(((PyObject *)__pyx_v_points), __pyx_ptype_5numpy_ndarray, 1, "points", 0))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 134; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_points, (PyObject*)__pyx_v_points, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 134; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_bstride_0_points = __pyx_bstruct_points.strides[0]; __pyx_bstride_1_points = __pyx_bstruct_points.strides[1];
  __pyx_bshape_0_points = __pyx_bstruct_points.shape[0]; __pyx_bshape_1_points = __pyx_bstruct_points.shape[1];

  
  __pyx_v_options = __pyx_k_1;

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 153; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__ascontiguousarray); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 153; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 153; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_INCREF(__pyx_v_points);
  PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_v_points);
  __Pyx_GIVEREF(__pyx_v_points);
  __pyx_t_3 = PyObject_Call(__pyx_t_2, __pyx_t_1, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 153; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (!(likely(((__pyx_t_3) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_3, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 153; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_4 = __pyx_t_3;
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_points);
    __pyx_t_5 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_points, (PyObject*)__pyx_t_4, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_5 < 0)) {
      PyErr_Fetch(&__pyx_t_6, &__pyx_t_7, &__pyx_t_8);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_points, (PyObject*)__pyx_v_points, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_6); Py_XDECREF(__pyx_t_7); Py_XDECREF(__pyx_t_8);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_6, __pyx_t_7, __pyx_t_8);
      }
    }
    __pyx_bstride_0_points = __pyx_bstruct_points.strides[0]; __pyx_bstride_1_points = __pyx_bstruct_points.strides[1];
    __pyx_bshape_0_points = __pyx_bstruct_points.shape[0]; __pyx_bshape_1_points = __pyx_bstruct_points.shape[1];
    if (unlikely(__pyx_t_5 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 153; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_4 = 0;
  __Pyx_DECREF(__pyx_v_points);
  __pyx_v_points = __pyx_t_3;
  __pyx_t_3 = 0;

  
  __pyx_v_numpoints = (((PyArrayObject *)__pyx_v_points)->dimensions[0]);

  
  __pyx_v_dim = (((PyArrayObject *)__pyx_v_points)->dimensions[1]);

  
  __pyx_t_9 = (__pyx_v_numpoints <= 0);
  if (__pyx_t_9) {

    
    __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 158; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __Pyx_INCREF(((PyObject *)__pyx_kp_s_2));
    PyTuple_SET_ITEM(__pyx_t_3, 0, ((PyObject *)__pyx_kp_s_2));
    __Pyx_GIVEREF(((PyObject *)__pyx_kp_s_2));
    __pyx_t_1 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_3, NULL); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 158; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    __Pyx_Raise(__pyx_t_1, 0, 0);
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    {__pyx_filename = __pyx_f[0]; __pyx_lineno = 158; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L5;
  }
  __pyx_L5:;

  
  __pyx_t_9 = (__pyx_v_dim < 2);
  if (__pyx_t_9) {

    
    __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 161; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __Pyx_INCREF(((PyObject *)__pyx_kp_s_3));
    PyTuple_SET_ITEM(__pyx_t_1, 0, ((PyObject *)__pyx_kp_s_3));
    __Pyx_GIVEREF(((PyObject *)__pyx_kp_s_3));
    __pyx_t_3 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_1, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 161; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    __Pyx_Raise(__pyx_t_3, 0, 0);
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    {__pyx_filename = __pyx_f[0]; __pyx_lineno = 161; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L6;
  }
  __pyx_L6:;

  
  __pyx_t_3 = __Pyx_GetName(__pyx_m, __pyx_n_s___qhull_lock); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 163; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_1 = PyObject_GetAttr(__pyx_t_3, __pyx_n_s__acquire); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 163; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyObject_Call(__pyx_t_1, ((PyObject *)__pyx_empty_tuple), NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 163; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;

  
   {

    
    qh_qh.NOerrexit = 1;

    
    __pyx_v_exitcode = qh_new_qhull(__pyx_v_dim, __pyx_v_numpoints, ((realT *)((PyArrayObject *)__pyx_v_points)->data), 0, __pyx_v_options, NULL, stderr);

    
     {

      
      __pyx_t_9 = (__pyx_v_exitcode != 0);
      if (__pyx_t_9) {

        
        __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 170; __pyx_clineno = __LINE__; goto __pyx_L11;}
        __Pyx_GOTREF(__pyx_t_3);
        __Pyx_INCREF(((PyObject *)__pyx_kp_s_4));
        PyTuple_SET_ITEM(__pyx_t_3, 0, ((PyObject *)__pyx_kp_s_4));
        __Pyx_GIVEREF(((PyObject *)__pyx_kp_s_4));
        __pyx_t_1 = PyObject_Call(__pyx_builtin_RuntimeError, __pyx_t_3, NULL); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 170; __pyx_clineno = __LINE__; goto __pyx_L11;}
        __Pyx_GOTREF(__pyx_t_1);
        __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
        __Pyx_Raise(__pyx_t_1, 0, 0);
        __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 170; __pyx_clineno = __LINE__; goto __pyx_L11;}
        goto __pyx_L13;
      }
      __pyx_L13:;

      
      qh_triangulate();

      
      if (qh_qh.SCALElast) {

        
        __pyx_t_10 = (qh_qh.last_high - qh_qh.last_low);
        if (unlikely(__pyx_t_10 == 0)) {
          PyErr_Format(PyExc_ZeroDivisionError, "float division");
          {__pyx_filename = __pyx_f[0]; __pyx_lineno = 175; __pyx_clineno = __LINE__; goto __pyx_L11;}
        }
        __pyx_v_paraboloid_scale = (qh_qh.last_newhigh / __pyx_t_10);

        
        __pyx_v_paraboloid_shift = ((-qh_qh.last_low) * __pyx_v_paraboloid_scale);
        goto __pyx_L14;
      }
       {

        
        __pyx_v_paraboloid_scale = 1.0;

        
        __pyx_v_paraboloid_shift = 0.0;
      }
      __pyx_L14:;

      
      __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s_5); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 183; __pyx_clineno = __LINE__; goto __pyx_L11;}
      __Pyx_GOTREF(__pyx_t_1);
      __pyx_t_3 = PyInt_FromLong(__pyx_v_dim); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 183; __pyx_clineno = __LINE__; goto __pyx_L11;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_2 = PyInt_FromLong(__pyx_v_numpoints); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 183; __pyx_clineno = __LINE__; goto __pyx_L11;}
      __Pyx_GOTREF(__pyx_t_2);
      __pyx_t_11 = PyTuple_New(2); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 183; __pyx_clineno = __LINE__; goto __pyx_L11;}
      __Pyx_GOTREF(__pyx_t_11);
      PyTuple_SET_ITEM(__pyx_t_11, 0, __pyx_t_3);
      __Pyx_GIVEREF(__pyx_t_3);
      PyTuple_SET_ITEM(__pyx_t_11, 1, __pyx_t_2);
      __Pyx_GIVEREF(__pyx_t_2);
      __pyx_t_3 = 0;
      __pyx_t_2 = 0;
      __pyx_t_2 = PyObject_Call(__pyx_t_1, __pyx_t_11, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 183; __pyx_clineno = __LINE__; goto __pyx_L11;}
      __Pyx_GOTREF(__pyx_t_2);
      __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
      __Pyx_DECREF(__pyx_t_11); __pyx_t_11 = 0;
      if (PyTuple_CheckExact(__pyx_t_2) && likely(PyTuple_GET_SIZE(__pyx_t_2) == 3)) {
        PyObject* tuple = __pyx_t_2;
        __pyx_t_11 = PyTuple_GET_ITEM(tuple, 0); __Pyx_INCREF(__pyx_t_11);
        __pyx_t_1 = PyTuple_GET_ITEM(tuple, 1); __Pyx_INCREF(__pyx_t_1);
        __pyx_t_3 = PyTuple_GET_ITEM(tuple, 2); __Pyx_INCREF(__pyx_t_3);

        
        __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
        __Pyx_DECREF(__pyx_v_vertices);
        __pyx_v_vertices = __pyx_t_11;
        __pyx_t_11 = 0;
        __Pyx_DECREF(__pyx_v_neighbors);
        __pyx_v_neighbors = __pyx_t_1;
        __pyx_t_1 = 0;
        __Pyx_DECREF(__pyx_v_equations);
        __pyx_v_equations = __pyx_t_3;
        __pyx_t_3 = 0;
      } else {
        __pyx_t_12 = PyObject_GetIter(__pyx_t_2); if (unlikely(!__pyx_t_12)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 182; __pyx_clineno = __LINE__; goto __pyx_L11;}
        __Pyx_GOTREF(__pyx_t_12);
        __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
        __pyx_t_11 = __Pyx_UnpackItem(__pyx_t_12, 0); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 182; __pyx_clineno = __LINE__; goto __pyx_L11;}
        __Pyx_GOTREF(__pyx_t_11);
        __pyx_t_1 = __Pyx_UnpackItem(__pyx_t_12, 1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 182; __pyx_clineno = __LINE__; goto __pyx_L11;}
        __Pyx_GOTREF(__pyx_t_1);
        __pyx_t_3 = __Pyx_UnpackItem(__pyx_t_12, 2); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 182; __pyx_clineno = __LINE__; goto __pyx_L11;}
        __Pyx_GOTREF(__pyx_t_3);
        if (__Pyx_EndUnpack(__pyx_t_12, 3) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 182; __pyx_clineno = __LINE__; goto __pyx_L11;}
        __Pyx_DECREF(__pyx_t_12); __pyx_t_12 = 0;
        __Pyx_DECREF(__pyx_v_vertices);
        __pyx_v_vertices = __pyx_t_11;
        __pyx_t_11 = 0;
        __Pyx_DECREF(__pyx_v_neighbors);
        __pyx_v_neighbors = __pyx_t_1;
        __pyx_t_1 = 0;
        __Pyx_DECREF(__pyx_v_equations);
        __pyx_v_equations = __pyx_t_3;
        __pyx_t_3 = 0;
      }

      
      __Pyx_XDECREF(__pyx_r);

      
      __pyx_t_2 = PyFloat_FromDouble(__pyx_v_paraboloid_scale); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 186; __pyx_clineno = __LINE__; goto __pyx_L11;}
      __Pyx_GOTREF(__pyx_t_2);
      __pyx_t_3 = PyFloat_FromDouble(__pyx_v_paraboloid_shift); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 186; __pyx_clineno = __LINE__; goto __pyx_L11;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_1 = PyTuple_New(5); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 185; __pyx_clineno = __LINE__; goto __pyx_L11;}
      __Pyx_GOTREF(__pyx_t_1);
      __Pyx_INCREF(__pyx_v_vertices);
      PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_v_vertices);
      __Pyx_GIVEREF(__pyx_v_vertices);
      __Pyx_INCREF(__pyx_v_neighbors);
      PyTuple_SET_ITEM(__pyx_t_1, 1, __pyx_v_neighbors);
      __Pyx_GIVEREF(__pyx_v_neighbors);
      __Pyx_INCREF(__pyx_v_equations);
      PyTuple_SET_ITEM(__pyx_t_1, 2, __pyx_v_equations);
      __Pyx_GIVEREF(__pyx_v_equations);
      PyTuple_SET_ITEM(__pyx_t_1, 3, __pyx_t_2);
      __Pyx_GIVEREF(__pyx_t_2);
      PyTuple_SET_ITEM(__pyx_t_1, 4, __pyx_t_3);
      __Pyx_GIVEREF(__pyx_t_3);
      __pyx_t_2 = 0;
      __pyx_t_3 = 0;
      __pyx_r = __pyx_t_1;
      __pyx_t_1 = 0;
      goto __pyx_L10;
    }

    
     {
      int __pyx_why;
      PyObject *__pyx_exc_type, *__pyx_exc_value, *__pyx_exc_tb;
      int __pyx_exc_lineno;
      __pyx_exc_type = 0; __pyx_exc_value = 0; __pyx_exc_tb = 0; __pyx_exc_lineno = 0;
      __pyx_why = 0; goto __pyx_L12;
      __pyx_L10: __pyx_exc_type = 0; __pyx_exc_value = 0; __pyx_exc_tb = 0; __pyx_exc_lineno = 0;
      __pyx_why = 3; goto __pyx_L12;
      __pyx_L11: {
        __pyx_why = 4;
        __Pyx_XDECREF(__pyx_t_12); __pyx_t_12 = 0;
        __Pyx_XDECREF(__pyx_t_11); __pyx_t_11 = 0;
        __Pyx_XDECREF(__pyx_t_2); __pyx_t_2 = 0;
        __Pyx_XDECREF(__pyx_t_3); __pyx_t_3 = 0;
        __Pyx_XDECREF(__pyx_t_1); __pyx_t_1 = 0;
        __Pyx_ErrFetch(&__pyx_exc_type, &__pyx_exc_value, &__pyx_exc_tb);
        __pyx_exc_lineno = __pyx_lineno;
        goto __pyx_L12;
      }
      __pyx_L12:;
      qh_freeqhull(0);

      
      qh_memfreeshort((&__pyx_v_curlong), (&__pyx_v_totlong));

      
      __pyx_t_9 = (__pyx_v_curlong != 0);
      if (!__pyx_t_9) {
        __pyx_t_13 = (__pyx_v_totlong != 0);
        __pyx_t_14 = __pyx_t_13;
      } else {
        __pyx_t_14 = __pyx_t_9;
      }
      if (__pyx_t_14) {

        
        __pyx_t_1 = PyInt_FromLong(__pyx_v_totlong); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 192; __pyx_clineno = __LINE__; goto __pyx_L15_error;}
        __Pyx_GOTREF(__pyx_t_1);
        __pyx_t_3 = PyInt_FromLong(__pyx_v_curlong); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 192; __pyx_clineno = __LINE__; goto __pyx_L15_error;}
        __Pyx_GOTREF(__pyx_t_3);
        __pyx_t_2 = PyTuple_New(2); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 192; __pyx_clineno = __LINE__; goto __pyx_L15_error;}
        __Pyx_GOTREF(__pyx_t_2);
        PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_t_1);
        __Pyx_GIVEREF(__pyx_t_1);
        PyTuple_SET_ITEM(__pyx_t_2, 1, __pyx_t_3);
        __Pyx_GIVEREF(__pyx_t_3);
        __pyx_t_1 = 0;
        __pyx_t_3 = 0;
        __pyx_t_3 = PyNumber_Remainder(((PyObject *)__pyx_kp_s_6), __pyx_t_2); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 191; __pyx_clineno = __LINE__; goto __pyx_L15_error;}
        __Pyx_GOTREF(((PyObject *)__pyx_t_3));
        __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
        __pyx_t_2 = PyTuple_New(1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 191; __pyx_clineno = __LINE__; goto __pyx_L15_error;}
        __Pyx_GOTREF(__pyx_t_2);
        PyTuple_SET_ITEM(__pyx_t_2, 0, ((PyObject *)__pyx_t_3));
        __Pyx_GIVEREF(((PyObject *)__pyx_t_3));
        __pyx_t_3 = 0;
        __pyx_t_3 = PyObject_Call(__pyx_builtin_RuntimeError, __pyx_t_2, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 191; __pyx_clineno = __LINE__; goto __pyx_L15_error;}
        __Pyx_GOTREF(__pyx_t_3);
        __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
        __Pyx_Raise(__pyx_t_3, 0, 0);
        __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 191; __pyx_clineno = __LINE__; goto __pyx_L15_error;}
        goto __pyx_L16;
      }
      __pyx_L16:;
      goto __pyx_L17;
      __pyx_L15_error:;
      if (__pyx_why == 4) {
        Py_XDECREF(__pyx_exc_type);
        Py_XDECREF(__pyx_exc_value);
        Py_XDECREF(__pyx_exc_tb);
      }
      goto __pyx_L8;
      __pyx_L17:;
      switch (__pyx_why) {
        case 3: goto __pyx_L7;
        case 4: {
          __Pyx_ErrRestore(__pyx_exc_type, __pyx_exc_value, __pyx_exc_tb);
          __pyx_lineno = __pyx_exc_lineno;
          __pyx_exc_type = 0;
          __pyx_exc_value = 0;
          __pyx_exc_tb = 0;
          goto __pyx_L8;
        }
      }
    }
  }

  
   {
    int __pyx_why;
    PyObject *__pyx_exc_type, *__pyx_exc_value, *__pyx_exc_tb;
    int __pyx_exc_lineno;
    __pyx_exc_type = 0; __pyx_exc_value = 0; __pyx_exc_tb = 0; __pyx_exc_lineno = 0;
    __pyx_why = 0; goto __pyx_L9;
    __pyx_L7: __pyx_exc_type = 0; __pyx_exc_value = 0; __pyx_exc_tb = 0; __pyx_exc_lineno = 0;
    __pyx_why = 3; goto __pyx_L9;
    __pyx_L8: {
      __pyx_why = 4;
      __Pyx_XDECREF(__pyx_t_12); __pyx_t_12 = 0;
      __Pyx_XDECREF(__pyx_t_11); __pyx_t_11 = 0;
      __Pyx_XDECREF(__pyx_t_1); __pyx_t_1 = 0;
      __Pyx_XDECREF(__pyx_t_2); __pyx_t_2 = 0;
      __Pyx_XDECREF(__pyx_t_3); __pyx_t_3 = 0;
      __Pyx_ErrFetch(&__pyx_exc_type, &__pyx_exc_value, &__pyx_exc_tb);
      __pyx_exc_lineno = __pyx_lineno;
      goto __pyx_L9;
    }
    __pyx_L9:;
    __pyx_t_3 = __Pyx_GetName(__pyx_m, __pyx_n_s___qhull_lock); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 194; __pyx_clineno = __LINE__; goto __pyx_L18_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __pyx_t_2 = PyObject_GetAttr(__pyx_t_3, __pyx_n_s__release); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 194; __pyx_clineno = __LINE__; goto __pyx_L18_error;}
    __Pyx_GOTREF(__pyx_t_2);
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    __pyx_t_3 = PyObject_Call(__pyx_t_2, ((PyObject *)__pyx_empty_tuple), NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 194; __pyx_clineno = __LINE__; goto __pyx_L18_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    goto __pyx_L19;
    __pyx_L18_error:;
    if (__pyx_why == 4) {
      Py_XDECREF(__pyx_exc_type);
      Py_XDECREF(__pyx_exc_value);
      Py_XDECREF(__pyx_exc_tb);
    }
    goto __pyx_L1_error;
    __pyx_L19:;
    switch (__pyx_why) {
      case 3: goto __pyx_L0;
      case 4: {
        __Pyx_ErrRestore(__pyx_exc_type, __pyx_exc_value, __pyx_exc_tb);
        __pyx_lineno = __pyx_exc_lineno;
        __pyx_exc_type = 0;
        __pyx_exc_value = 0;
        __pyx_exc_tb = 0;
        goto __pyx_L1_error;
      }
    }
  }

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_11);
  __Pyx_XDECREF(__pyx_t_12);
  { PyObject *__pyx_type, *__pyx_value, *__pyx_tb;
    __Pyx_ErrFetch(&__pyx_type, &__pyx_value, &__pyx_tb);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_points);
  __Pyx_ErrRestore(__pyx_type, __pyx_value, __pyx_tb);}
  __Pyx_AddTraceback("scipy.spatial.qhull._construct_delaunay");
  __pyx_r = NULL;
  goto __pyx_L2;
  __pyx_L0:;
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_points);
  __pyx_L2:;
  __Pyx_DECREF(__pyx_v_vertices);
  __Pyx_DECREF(__pyx_v_neighbors);
  __Pyx_DECREF(__pyx_v_equations);
  __Pyx_DECREF((PyObject *)__pyx_v_points);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull__qhull_get_facet_array(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); 
static char __pyx_doc_5scipy_7spatial_5qhull__qhull_get_facet_array[] = "\n    Return array of simplical facets currently in Qhull.\n\n    Returns\n    -------\n    vertices : array of int, shape (nfacets, ndim+1)\n        Indices of coordinates of vertices forming the simplical facets\n    neighbors : array of int, shape (nfacets, ndim)\n        Indices of neighboring facets.  The kth neighbor is opposite\n        the kth vertex, and the first neighbor is the horizon facet\n        for the first vertex.\n\n        Facets extending to infinity are denoted with index -1.\n\n    ";
static PyObject *__pyx_pf_5scipy_7spatial_5qhull__qhull_get_facet_array(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  int __pyx_v_ndim;
  int __pyx_v_numpoints;
  facetT *__pyx_v_facet;
  facetT *__pyx_v_neighbor;
  vertexT *__pyx_v_vertex;
  int __pyx_v_i;
  int __pyx_v_j;
  int __pyx_v_point;
  PyArrayObject *__pyx_v_vertices;
  PyArrayObject *__pyx_v_neighbors;
  PyArrayObject *__pyx_v_equations;
  PyArrayObject *__pyx_v_id_map;
  Py_buffer __pyx_bstruct_neighbors;
  Py_ssize_t __pyx_bstride_0_neighbors = 0;
  Py_ssize_t __pyx_bstride_1_neighbors = 0;
  Py_ssize_t __pyx_bshape_0_neighbors = 0;
  Py_ssize_t __pyx_bshape_1_neighbors = 0;
  Py_buffer __pyx_bstruct_id_map;
  Py_ssize_t __pyx_bstride_0_id_map = 0;
  Py_ssize_t __pyx_bshape_0_id_map = 0;
  Py_buffer __pyx_bstruct_vertices;
  Py_ssize_t __pyx_bstride_0_vertices = 0;
  Py_ssize_t __pyx_bstride_1_vertices = 0;
  Py_ssize_t __pyx_bshape_0_vertices = 0;
  Py_ssize_t __pyx_bshape_1_vertices = 0;
  Py_buffer __pyx_bstruct_equations;
  Py_ssize_t __pyx_bstride_0_equations = 0;
  Py_ssize_t __pyx_bstride_1_equations = 0;
  Py_ssize_t __pyx_bshape_0_equations = 0;
  Py_ssize_t __pyx_bshape_1_equations = 0;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  PyArrayObject *__pyx_t_6 = NULL;
  int __pyx_t_7;
  PyObject *__pyx_t_8 = NULL;
  PyObject *__pyx_t_9 = NULL;
  PyObject *__pyx_t_10 = NULL;
  int __pyx_t_11;
  int __pyx_t_12;
  unsigned int __pyx_t_13;
  PyArrayObject *__pyx_t_14 = NULL;
  PyArrayObject *__pyx_t_15 = NULL;
  PyArrayObject *__pyx_t_16 = NULL;
  long __pyx_t_17;
  int __pyx_t_18;
  int __pyx_t_19;
  int __pyx_t_20;
  unsigned int __pyx_t_21;
  int __pyx_t_22;
  int __pyx_t_23;
  int __pyx_t_24;
  int __pyx_t_25;
  static PyObject **__pyx_pyargnames[] = {&__pyx_n_s__ndim,&__pyx_n_s__numpoints,0};
  __Pyx_RefNannySetupContext("_qhull_get_facet_array");
  __pyx_self = __pyx_self;
  if (unlikely(__pyx_kwds)) {
    Py_ssize_t kw_args = PyDict_Size(__pyx_kwds);
    PyObject* values[2] = {0,0};
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  2: values[1] = PyTuple_GET_ITEM(__pyx_args, 1);
      case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
      case  0: break;
      default: goto __pyx_L5_argtuple_error;
    }
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  0:
      values[0] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__ndim);
      if (likely(values[0])) kw_args--;
      else goto __pyx_L5_argtuple_error;
      case  1:
      values[1] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__numpoints);
      if (likely(values[1])) kw_args--;
      else {
        __Pyx_RaiseArgtupleInvalid("_qhull_get_facet_array", 1, 2, 2, 1); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 197; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
      }
    }
    if (unlikely(kw_args > 0)) {
      if (unlikely(__Pyx_ParseOptionalKeywords(__pyx_kwds, __pyx_pyargnames, 0, values, PyTuple_GET_SIZE(__pyx_args), "_qhull_get_facet_array") < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 197; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    }
    __pyx_v_ndim = __Pyx_PyInt_AsInt(values[0]); if (unlikely((__pyx_v_ndim == (int)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 197; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    __pyx_v_numpoints = __Pyx_PyInt_AsInt(values[1]); if (unlikely((__pyx_v_numpoints == (int)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 197; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  } else if (PyTuple_GET_SIZE(__pyx_args) != 2) {
    goto __pyx_L5_argtuple_error;
  } else {
    __pyx_v_ndim = __Pyx_PyInt_AsInt(PyTuple_GET_ITEM(__pyx_args, 0)); if (unlikely((__pyx_v_ndim == (int)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 197; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    __pyx_v_numpoints = __Pyx_PyInt_AsInt(PyTuple_GET_ITEM(__pyx_args, 1)); if (unlikely((__pyx_v_numpoints == (int)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 197; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L5_argtuple_error:;
  __Pyx_RaiseArgtupleInvalid("_qhull_get_facet_array", 1, 2, 2, PyTuple_GET_SIZE(__pyx_args)); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 197; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  __pyx_L3_error:;
  __Pyx_AddTraceback("scipy.spatial.qhull._qhull_get_facet_array");
  __Pyx_RefNannyFinishContext();
  return NULL;
  __pyx_L4_argument_unpacking_done:;
  __pyx_v_vertices = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_neighbors = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_equations = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_id_map = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_bstruct_vertices.buf = NULL;
  __pyx_bstruct_neighbors.buf = NULL;
  __pyx_bstruct_equations.buf = NULL;
  __pyx_bstruct_id_map.buf = NULL;

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__empty); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyLong_FromUnsignedLong(qh_qh.facet_id); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_1);
  __Pyx_GIVEREF(__pyx_t_1);
  __pyx_t_1 = 0;
  __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_t_3);
  __Pyx_GIVEREF(__pyx_t_3);
  __pyx_t_3 = 0;
  __pyx_t_3 = PyDict_New(); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_3));
  __pyx_t_4 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_5 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__intc); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyDict_SetItem(__pyx_t_3, ((PyObject *)__pyx_n_s__dtype), __pyx_t_5) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __pyx_t_5 = PyEval_CallObjectWithKeywords(__pyx_t_2, __pyx_t_1, ((PyObject *)__pyx_t_3)); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_3)); __pyx_t_3 = 0;
  if (!(likely(((__pyx_t_5) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_5, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_6 = ((PyArrayObject *)__pyx_t_5);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_id_map);
    __pyx_t_7 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_id_map, (PyObject*)__pyx_t_6, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 1, 0, __pyx_stack);
    if (unlikely(__pyx_t_7 < 0)) {
      PyErr_Fetch(&__pyx_t_8, &__pyx_t_9, &__pyx_t_10);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_id_map, (PyObject*)__pyx_v_id_map, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 1, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_8); Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_10);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_8, __pyx_t_9, __pyx_t_10);
      }
    }
    __pyx_bstride_0_id_map = __pyx_bstruct_id_map.strides[0];
    __pyx_bshape_0_id_map = __pyx_bstruct_id_map.shape[0];
    if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 223; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_6 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_id_map));
  __pyx_v_id_map = ((PyArrayObject *)__pyx_t_5);
  __pyx_t_5 = 0;

  
  __pyx_t_5 = PyObject_GetAttr(((PyObject *)__pyx_v_id_map), __pyx_n_s__fill); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 224; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 224; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_INCREF(__pyx_int_neg_1);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_int_neg_1);
  __Pyx_GIVEREF(__pyx_int_neg_1);
  __pyx_t_1 = PyObject_Call(__pyx_t_5, __pyx_t_3, NULL); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 224; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_v_facet = qh_qh.facet_list;

  
  __pyx_v_j = 0;

  
  while (1) {
    if ((__pyx_v_facet != 0)) {
      __pyx_t_11 = (__pyx_v_facet->next != 0);
    } else {
      __pyx_t_11 = (__pyx_v_facet != 0);
    }
    if (!__pyx_t_11) break;

    
    if (__pyx_v_facet->simplicial) {
      __pyx_t_11 = (!__pyx_v_facet->upperdelaunay);
      __pyx_t_12 = __pyx_t_11;
    } else {
      __pyx_t_12 = __pyx_v_facet->simplicial;
    }
    if (__pyx_t_12) {

      
      __pyx_t_13 = __pyx_v_facet->id;
      __pyx_t_7 = -1;
      if (unlikely(__pyx_t_13 >= (size_t)__pyx_bshape_0_id_map)) __pyx_t_7 = 0;
      if (unlikely(__pyx_t_7 != -1)) {
        __Pyx_RaiseBufferIndexError(__pyx_t_7);
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 231; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      *__Pyx_BufPtrStrided1d(npy_int *, __pyx_bstruct_id_map.buf, __pyx_t_13, __pyx_bstride_0_id_map) = __pyx_v_j;

      
      __pyx_v_j += 1;
      goto __pyx_L8;
    }
    __pyx_L8:;

    
    __pyx_v_facet = __pyx_v_facet->next;
  }

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__zeros); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyInt_FromLong(__pyx_v_j); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_5 = PyInt_FromLong((__pyx_v_ndim + 1)); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_2 = PyTuple_New(2); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_t_1);
  __Pyx_GIVEREF(__pyx_t_1);
  PyTuple_SET_ITEM(__pyx_t_2, 1, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_1 = 0;
  __pyx_t_5 = 0;
  __pyx_t_5 = PyTuple_New(1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  PyTuple_SET_ITEM(__pyx_t_5, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyDict_New(); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_2));
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_4 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__intc); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (PyDict_SetItem(__pyx_t_2, ((PyObject *)__pyx_n_s__dtype), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyEval_CallObjectWithKeywords(__pyx_t_3, __pyx_t_5, ((PyObject *)__pyx_t_2)); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_2)); __pyx_t_2 = 0;
  if (!(likely(((__pyx_t_4) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_4, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_14 = ((PyArrayObject *)__pyx_t_4);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
    __pyx_t_7 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_vertices, (PyObject*)__pyx_t_14, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_7 < 0)) {
      PyErr_Fetch(&__pyx_t_10, &__pyx_t_9, &__pyx_t_8);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_vertices, (PyObject*)__pyx_v_vertices, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_10); Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_8);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_10, __pyx_t_9, __pyx_t_8);
      }
    }
    __pyx_bstride_0_vertices = __pyx_bstruct_vertices.strides[0]; __pyx_bstride_1_vertices = __pyx_bstruct_vertices.strides[1];
    __pyx_bshape_0_vertices = __pyx_bstruct_vertices.shape[0]; __pyx_bshape_1_vertices = __pyx_bstruct_vertices.shape[1];
    if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 236; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_14 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_vertices));
  __pyx_v_vertices = ((PyArrayObject *)__pyx_t_4);
  __pyx_t_4 = 0;

  
  __pyx_t_4 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__zeros); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyInt_FromLong(__pyx_v_j); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_5 = PyInt_FromLong((__pyx_v_ndim + 1)); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_3 = PyTuple_New(2); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_4);
  __Pyx_GIVEREF(__pyx_t_4);
  PyTuple_SET_ITEM(__pyx_t_3, 1, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_4 = 0;
  __pyx_t_5 = 0;
  __pyx_t_5 = PyTuple_New(1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  PyTuple_SET_ITEM(__pyx_t_5, 0, __pyx_t_3);
  __Pyx_GIVEREF(__pyx_t_3);
  __pyx_t_3 = 0;
  __pyx_t_3 = PyDict_New(); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_3));
  __pyx_t_4 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_1 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__intc); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyDict_SetItem(__pyx_t_3, ((PyObject *)__pyx_n_s__dtype), __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyEval_CallObjectWithKeywords(__pyx_t_2, __pyx_t_5, ((PyObject *)__pyx_t_3)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_3)); __pyx_t_3 = 0;
  if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_15 = ((PyArrayObject *)__pyx_t_1);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_neighbors);
    __pyx_t_7 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_neighbors, (PyObject*)__pyx_t_15, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_7 < 0)) {
      PyErr_Fetch(&__pyx_t_8, &__pyx_t_9, &__pyx_t_10);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_neighbors, (PyObject*)__pyx_v_neighbors, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_8); Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_10);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_8, __pyx_t_9, __pyx_t_10);
      }
    }
    __pyx_bstride_0_neighbors = __pyx_bstruct_neighbors.strides[0]; __pyx_bstride_1_neighbors = __pyx_bstruct_neighbors.strides[1];
    __pyx_bshape_0_neighbors = __pyx_bstruct_neighbors.shape[0]; __pyx_bshape_1_neighbors = __pyx_bstruct_neighbors.shape[1];
    if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 237; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_15 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_neighbors));
  __pyx_v_neighbors = ((PyArrayObject *)__pyx_t_1);
  __pyx_t_1 = 0;

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__zeros); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyInt_FromLong(__pyx_v_j); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_5 = PyInt_FromLong((__pyx_v_ndim + 2)); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_2 = PyTuple_New(2); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_t_1);
  __Pyx_GIVEREF(__pyx_t_1);
  PyTuple_SET_ITEM(__pyx_t_2, 1, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_1 = 0;
  __pyx_t_5 = 0;
  __pyx_t_5 = PyTuple_New(1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  PyTuple_SET_ITEM(__pyx_t_5, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyDict_New(); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_2));
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_4 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__double); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (PyDict_SetItem(__pyx_t_2, ((PyObject *)__pyx_n_s__dtype), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyEval_CallObjectWithKeywords(__pyx_t_3, __pyx_t_5, ((PyObject *)__pyx_t_2)); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_2)); __pyx_t_2 = 0;
  if (!(likely(((__pyx_t_4) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_4, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_16 = ((PyArrayObject *)__pyx_t_4);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_equations);
    __pyx_t_7 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_equations, (PyObject*)__pyx_t_16, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_7 < 0)) {
      PyErr_Fetch(&__pyx_t_10, &__pyx_t_9, &__pyx_t_8);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_equations, (PyObject*)__pyx_v_equations, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_10); Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_8);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_10, __pyx_t_9, __pyx_t_8);
      }
    }
    __pyx_bstride_0_equations = __pyx_bstruct_equations.strides[0]; __pyx_bstride_1_equations = __pyx_bstruct_equations.strides[1];
    __pyx_bshape_0_equations = __pyx_bstruct_equations.shape[0]; __pyx_bshape_1_equations = __pyx_bstruct_equations.shape[1];
    if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 238; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_16 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_equations));
  __pyx_v_equations = ((PyArrayObject *)__pyx_t_4);
  __pyx_t_4 = 0;

  
  __pyx_v_facet = qh_qh.facet_list;

  
  __pyx_v_j = 0;

  
  while (1) {
    if ((__pyx_v_facet != 0)) {
      __pyx_t_12 = (__pyx_v_facet->next != 0);
    } else {
      __pyx_t_12 = (__pyx_v_facet != 0);
    }
    if (!__pyx_t_12) break;

    
    __pyx_t_12 = (!__pyx_v_facet->simplicial);
    if (__pyx_t_12) {

      
      __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 245; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_4);
      __Pyx_INCREF(((PyObject *)__pyx_kp_s_7));
      PyTuple_SET_ITEM(__pyx_t_4, 0, ((PyObject *)__pyx_kp_s_7));
      __Pyx_GIVEREF(((PyObject *)__pyx_kp_s_7));
      __pyx_t_2 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_4, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 245; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_2);
      __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
      __Pyx_Raise(__pyx_t_2, 0, 0);
      __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
      {__pyx_filename = __pyx_f[0]; __pyx_lineno = 245; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      goto __pyx_L11;
    }
    __pyx_L11:;

    
    if (__pyx_v_facet->upperdelaunay) {

      
      __pyx_v_facet = __pyx_v_facet->next;

      
      goto __pyx_L9_continue;
      goto __pyx_L12;
    }
    __pyx_L12:;

    
    __pyx_t_17 = (__pyx_v_ndim + 1);
    for (__pyx_t_7 = 0; __pyx_t_7 < __pyx_t_17; __pyx_t_7+=1) {
      __pyx_v_i = __pyx_t_7;

      
      __pyx_v_vertex = ((vertexT *)(__pyx_v_facet->vertices->e[__pyx_v_i]).p);

      
      __pyx_v_point = qh_pointid(__pyx_v_vertex->point);

      
      __pyx_t_18 = __pyx_v_j;
      __pyx_t_19 = __pyx_v_i;
      __pyx_t_20 = -1;
      if (__pyx_t_18 < 0) {
        __pyx_t_18 += __pyx_bshape_0_vertices;
        if (unlikely(__pyx_t_18 < 0)) __pyx_t_20 = 0;
      } else if (unlikely(__pyx_t_18 >= __pyx_bshape_0_vertices)) __pyx_t_20 = 0;
      if (__pyx_t_19 < 0) {
        __pyx_t_19 += __pyx_bshape_1_vertices;
        if (unlikely(__pyx_t_19 < 0)) __pyx_t_20 = 1;
      } else if (unlikely(__pyx_t_19 >= __pyx_bshape_1_vertices)) __pyx_t_20 = 1;
      if (unlikely(__pyx_t_20 != -1)) {
        __Pyx_RaiseBufferIndexError(__pyx_t_20);
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 255; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      *__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_18, __pyx_bstride_0_vertices, __pyx_t_19, __pyx_bstride_1_vertices) = __pyx_v_point;
    }

    
    __pyx_t_17 = (__pyx_v_ndim + 1);
    for (__pyx_t_7 = 0; __pyx_t_7 < __pyx_t_17; __pyx_t_7+=1) {
      __pyx_v_i = __pyx_t_7;

      
      __pyx_v_neighbor = ((facetT *)(__pyx_v_facet->neighbors->e[__pyx_v_i]).p);

      
      __pyx_t_21 = __pyx_v_neighbor->id;
      __pyx_t_20 = -1;
      if (unlikely(__pyx_t_21 >= (size_t)__pyx_bshape_0_id_map)) __pyx_t_20 = 0;
      if (unlikely(__pyx_t_20 != -1)) {
        __Pyx_RaiseBufferIndexError(__pyx_t_20);
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 260; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      __pyx_t_20 = __pyx_v_j;
      __pyx_t_22 = __pyx_v_i;
      __pyx_t_23 = -1;
      if (__pyx_t_20 < 0) {
        __pyx_t_20 += __pyx_bshape_0_neighbors;
        if (unlikely(__pyx_t_20 < 0)) __pyx_t_23 = 0;
      } else if (unlikely(__pyx_t_20 >= __pyx_bshape_0_neighbors)) __pyx_t_23 = 0;
      if (__pyx_t_22 < 0) {
        __pyx_t_22 += __pyx_bshape_1_neighbors;
        if (unlikely(__pyx_t_22 < 0)) __pyx_t_23 = 1;
      } else if (unlikely(__pyx_t_22 >= __pyx_bshape_1_neighbors)) __pyx_t_23 = 1;
      if (unlikely(__pyx_t_23 != -1)) {
        __Pyx_RaiseBufferIndexError(__pyx_t_23);
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 260; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      *__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_neighbors.buf, __pyx_t_20, __pyx_bstride_0_neighbors, __pyx_t_22, __pyx_bstride_1_neighbors) = (*__Pyx_BufPtrStrided1d(npy_int *, __pyx_bstruct_id_map.buf, __pyx_t_21, __pyx_bstride_0_id_map));
    }

    
    __pyx_t_17 = (__pyx_v_ndim + 1);
    for (__pyx_t_7 = 0; __pyx_t_7 < __pyx_t_17; __pyx_t_7+=1) {
      __pyx_v_i = __pyx_t_7;

      
      __pyx_t_23 = __pyx_v_j;
      __pyx_t_24 = __pyx_v_i;
      __pyx_t_25 = -1;
      if (__pyx_t_23 < 0) {
        __pyx_t_23 += __pyx_bshape_0_equations;
        if (unlikely(__pyx_t_23 < 0)) __pyx_t_25 = 0;
      } else if (unlikely(__pyx_t_23 >= __pyx_bshape_0_equations)) __pyx_t_25 = 0;
      if (__pyx_t_24 < 0) {
        __pyx_t_24 += __pyx_bshape_1_equations;
        if (unlikely(__pyx_t_24 < 0)) __pyx_t_25 = 1;
      } else if (unlikely(__pyx_t_24 >= __pyx_bshape_1_equations)) __pyx_t_25 = 1;
      if (unlikely(__pyx_t_25 != -1)) {
        __Pyx_RaiseBufferIndexError(__pyx_t_25);
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 264; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      *__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_equations.buf, __pyx_t_23, __pyx_bstride_0_equations, __pyx_t_24, __pyx_bstride_1_equations) = (__pyx_v_facet->normal[__pyx_v_i]);
    }

    
    __pyx_t_7 = __pyx_v_j;
    __pyx_t_17 = (__pyx_v_ndim + 1);
    __pyx_t_25 = -1;
    if (__pyx_t_7 < 0) {
      __pyx_t_7 += __pyx_bshape_0_equations;
      if (unlikely(__pyx_t_7 < 0)) __pyx_t_25 = 0;
    } else if (unlikely(__pyx_t_7 >= __pyx_bshape_0_equations)) __pyx_t_25 = 0;
    if (__pyx_t_17 < 0) {
      __pyx_t_17 += __pyx_bshape_1_equations;
      if (unlikely(__pyx_t_17 < 0)) __pyx_t_25 = 1;
    } else if (unlikely(__pyx_t_17 >= __pyx_bshape_1_equations)) __pyx_t_25 = 1;
    if (unlikely(__pyx_t_25 != -1)) {
      __Pyx_RaiseBufferIndexError(__pyx_t_25);
      {__pyx_filename = __pyx_f[0]; __pyx_lineno = 265; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    }
    *__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_equations.buf, __pyx_t_7, __pyx_bstride_0_equations, __pyx_t_17, __pyx_bstride_1_equations) = __pyx_v_facet->offset;

    
    __pyx_v_j += 1;

    
    __pyx_v_facet = __pyx_v_facet->next;
    __pyx_L9_continue:;
  }

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_2 = PyTuple_New(3); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 270; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_INCREF(((PyObject *)__pyx_v_vertices));
  PyTuple_SET_ITEM(__pyx_t_2, 0, ((PyObject *)__pyx_v_vertices));
  __Pyx_GIVEREF(((PyObject *)__pyx_v_vertices));
  __Pyx_INCREF(((PyObject *)__pyx_v_neighbors));
  PyTuple_SET_ITEM(__pyx_t_2, 1, ((PyObject *)__pyx_v_neighbors));
  __Pyx_GIVEREF(((PyObject *)__pyx_v_neighbors));
  __Pyx_INCREF(((PyObject *)__pyx_v_equations));
  PyTuple_SET_ITEM(__pyx_t_2, 2, ((PyObject *)__pyx_v_equations));
  __Pyx_GIVEREF(((PyObject *)__pyx_v_equations));
  __pyx_r = __pyx_t_2;
  __pyx_t_2 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  { PyObject *__pyx_type, *__pyx_value, *__pyx_tb;
    __Pyx_ErrFetch(&__pyx_type, &__pyx_value, &__pyx_tb);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_neighbors);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_id_map);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_equations);
  __Pyx_ErrRestore(__pyx_type, __pyx_value, __pyx_tb);}
  __Pyx_AddTraceback("scipy.spatial.qhull._qhull_get_facet_array");
  __pyx_r = NULL;
  goto __pyx_L2;
  __pyx_L0:;
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_neighbors);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_id_map);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_equations);
  __pyx_L2:;
  __Pyx_DECREF((PyObject *)__pyx_v_vertices);
  __Pyx_DECREF((PyObject *)__pyx_v_neighbors);
  __Pyx_DECREF((PyObject *)__pyx_v_equations);
  __Pyx_DECREF((PyObject *)__pyx_v_id_map);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull__get_barycentric_transforms(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); 
static char __pyx_doc_5scipy_7spatial_5qhull__get_barycentric_transforms[] = "\n    Compute barycentric affine coordinate transformations for given\n    simplices.\n\n    Returns\n    -------\n    Tinvs : array, shape (nsimplex, ndim+1, ndim)\n        Barycentric transforms for each simplex.\n\n        Tinvs[i,:ndim,:ndim] contains inverse of the matrix ``T``,\n        and Tinvs[i,ndim,:] contains the vector ``r_n`` (see below).\n\n    Notes\n    -----\n    Barycentric transform from ``x`` to ``c`` is defined by::\n\n        T c = x - r_n\n\n    where the ``r_1, ..., r_n`` are the vertices of the simplex.\n    The matrix ``T`` is defined by the condition::\n\n        T e_j = r_j - r_n\n\n    where ``e_j`` is the unit axis vector, e.g, ``e_2 = [0,1,0,0,...]``\n    This implies that ``T_ij = (r_j - r_n)_i``.\n\n    For the barycentric transforms, we need to compute the inverse\n    matrix ``T^-1`` and store the vectors ``r_n`` for each vertex.\n    These are stacked into the `Tinvs` returned.\n\n    ";
static PyObject *__pyx_pf_5scipy_7spatial_5qhull__get_barycentric_transforms(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyArrayObject *__pyx_v_points = 0;
  PyArrayObject *__pyx_v_vertices = 0;
  PyArrayObject *__pyx_v_T;
  PyArrayObject *__pyx_v_Tinvs;
  int __pyx_v_ivertex;
  int __pyx_v_i;
  int __pyx_v_j;
  int __pyx_v_n;
  int __pyx_v_nrhs;
  int __pyx_v_lda;
  int __pyx_v_ldb;
  int __pyx_v_info;
  int __pyx_v_ipiv[(NPY_MAXDIMS + 1)];
  int __pyx_v_ndim;
  int __pyx_v_nvertex;
  double __pyx_v_nan;
  double __pyx_v_x1;
  double __pyx_v_x2;
  double __pyx_v_x3;
  double __pyx_v_y1;
  double __pyx_v_y2;
  double __pyx_v_y3;
  double __pyx_v_det;
  Py_buffer __pyx_bstruct_Tinvs;
  Py_ssize_t __pyx_bstride_0_Tinvs = 0;
  Py_ssize_t __pyx_bstride_1_Tinvs = 0;
  Py_ssize_t __pyx_bstride_2_Tinvs = 0;
  Py_ssize_t __pyx_bshape_0_Tinvs = 0;
  Py_ssize_t __pyx_bshape_1_Tinvs = 0;
  Py_ssize_t __pyx_bshape_2_Tinvs = 0;
  Py_buffer __pyx_bstruct_T;
  Py_ssize_t __pyx_bstride_0_T = 0;
  Py_ssize_t __pyx_bstride_1_T = 0;
  Py_ssize_t __pyx_bshape_0_T = 0;
  Py_ssize_t __pyx_bshape_1_T = 0;
  Py_buffer __pyx_bstruct_vertices;
  Py_ssize_t __pyx_bstride_0_vertices = 0;
  Py_ssize_t __pyx_bstride_1_vertices = 0;
  Py_ssize_t __pyx_bshape_0_vertices = 0;
  Py_ssize_t __pyx_bshape_1_vertices = 0;
  Py_buffer __pyx_bstruct_points;
  Py_ssize_t __pyx_bstride_0_points = 0;
  Py_ssize_t __pyx_bstride_1_points = 0;
  Py_ssize_t __pyx_bshape_0_points = 0;
  Py_ssize_t __pyx_bshape_1_points = 0;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  double __pyx_t_3;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  PyObject *__pyx_t_6 = NULL;
  PyArrayObject *__pyx_t_7 = NULL;
  int __pyx_t_8;
  PyObject *__pyx_t_9 = NULL;
  PyObject *__pyx_t_10 = NULL;
  PyObject *__pyx_t_11 = NULL;
  PyArrayObject *__pyx_t_12 = NULL;
  int __pyx_t_13;
  int __pyx_t_14;
  int __pyx_t_15;
  long __pyx_t_16;
  npy_int __pyx_t_17;
  long __pyx_t_18;
  int __pyx_t_19;
  long __pyx_t_20;
  npy_int __pyx_t_21;
  long __pyx_t_22;
  int __pyx_t_23;
  long __pyx_t_24;
  npy_int __pyx_t_25;
  long __pyx_t_26;
  int __pyx_t_27;
  long __pyx_t_28;
  npy_int __pyx_t_29;
  long __pyx_t_30;
  int __pyx_t_31;
  long __pyx_t_32;
  npy_int __pyx_t_33;
  long __pyx_t_34;
  int __pyx_t_35;
  long __pyx_t_36;
  npy_int __pyx_t_37;
  long __pyx_t_38;
  int __pyx_t_39;
  long __pyx_t_40;
  long __pyx_t_41;
  int __pyx_t_42;
  long __pyx_t_43;
  long __pyx_t_44;
  int __pyx_t_45;
  long __pyx_t_46;
  long __pyx_t_47;
  int __pyx_t_48;
  long __pyx_t_49;
  long __pyx_t_50;
  int __pyx_t_51;
  long __pyx_t_52;
  long __pyx_t_53;
  int __pyx_t_54;
  long __pyx_t_55;
  long __pyx_t_56;
  int __pyx_t_57;
  int __pyx_t_58;
  int __pyx_t_59;
  int __pyx_t_60;
  npy_int __pyx_t_61;
  int __pyx_t_62;
  int __pyx_t_63;
  int __pyx_t_64;
  int __pyx_t_65;
  int __pyx_t_66;
  int __pyx_t_67;
  int __pyx_t_68;
  int __pyx_t_69;
  npy_int __pyx_t_70;
  int __pyx_t_71;
  int __pyx_t_72;
  int __pyx_t_73;
  int __pyx_t_74;
  int __pyx_t_75;
  int __pyx_t_76;
  int __pyx_t_77;
  long __pyx_t_78;
  int __pyx_t_79;
  int __pyx_t_80;
  int __pyx_t_81;
  int __pyx_t_82;
  static PyObject **__pyx_pyargnames[] = {&__pyx_n_s__points,&__pyx_n_s__vertices,0};
  __Pyx_RefNannySetupContext("_get_barycentric_transforms");
  __pyx_self = __pyx_self;
  if (unlikely(__pyx_kwds)) {
    Py_ssize_t kw_args = PyDict_Size(__pyx_kwds);
    PyObject* values[2] = {0,0};
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  2: values[1] = PyTuple_GET_ITEM(__pyx_args, 1);
      case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
      case  0: break;
      default: goto __pyx_L5_argtuple_error;
    }
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  0:
      values[0] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__points);
      if (likely(values[0])) kw_args--;
      else goto __pyx_L5_argtuple_error;
      case  1:
      values[1] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__vertices);
      if (likely(values[1])) kw_args--;
      else {
        __Pyx_RaiseArgtupleInvalid("_get_barycentric_transforms", 1, 2, 2, 1); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 278; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
      }
    }
    if (unlikely(kw_args > 0)) {
      if (unlikely(__Pyx_ParseOptionalKeywords(__pyx_kwds, __pyx_pyargnames, 0, values, PyTuple_GET_SIZE(__pyx_args), "_get_barycentric_transforms") < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 278; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    }
    __pyx_v_points = ((PyArrayObject *)values[0]);
    __pyx_v_vertices = ((PyArrayObject *)values[1]);
  } else if (PyTuple_GET_SIZE(__pyx_args) != 2) {
    goto __pyx_L5_argtuple_error;
  } else {
    __pyx_v_points = ((PyArrayObject *)PyTuple_GET_ITEM(__pyx_args, 0));
    __pyx_v_vertices = ((PyArrayObject *)PyTuple_GET_ITEM(__pyx_args, 1));
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L5_argtuple_error:;
  __Pyx_RaiseArgtupleInvalid("_get_barycentric_transforms", 1, 2, 2, PyTuple_GET_SIZE(__pyx_args)); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 278; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  __pyx_L3_error:;
  __Pyx_AddTraceback("scipy.spatial.qhull._get_barycentric_transforms");
  __Pyx_RefNannyFinishContext();
  return NULL;
  __pyx_L4_argument_unpacking_done:;
  __pyx_v_T = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_Tinvs = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_bstruct_T.buf = NULL;
  __pyx_bstruct_Tinvs.buf = NULL;
  __pyx_bstruct_points.buf = NULL;
  __pyx_bstruct_vertices.buf = NULL;
  if (unlikely(!__Pyx_ArgTypeTest(((PyObject *)__pyx_v_points), __pyx_ptype_5numpy_ndarray, 1, "points", 0))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 278; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (unlikely(!__Pyx_ArgTypeTest(((PyObject *)__pyx_v_vertices), __pyx_ptype_5numpy_ndarray, 1, "vertices", 0))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 279; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_points, (PyObject*)__pyx_v_points, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 278; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_bstride_0_points = __pyx_bstruct_points.strides[0]; __pyx_bstride_1_points = __pyx_bstruct_points.strides[1];
  __pyx_bshape_0_points = __pyx_bstruct_points.shape[0]; __pyx_bshape_1_points = __pyx_bstruct_points.shape[1];
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_vertices, (PyObject*)__pyx_v_vertices, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 278; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_bstride_0_vertices = __pyx_bstruct_vertices.strides[0]; __pyx_bstride_1_vertices = __pyx_bstruct_vertices.strides[1];
  __pyx_bshape_0_vertices = __pyx_bstruct_vertices.shape[0]; __pyx_bshape_1_vertices = __pyx_bstruct_vertices.shape[1];

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 324; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__nan); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 324; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_3 = __pyx_PyFloat_AsDouble(__pyx_t_2); if (unlikely((__pyx_t_3 == (double)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 324; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_v_nan = __pyx_t_3;

  
  __pyx_v_ndim = (__pyx_v_points->dimensions[1]);

  
  __pyx_v_nvertex = (__pyx_v_vertices->dimensions[0]);

  
  __pyx_t_2 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_1 = PyObject_GetAttr(__pyx_t_2, __pyx_n_s__zeros); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_t_2 = PyInt_FromLong(__pyx_v_ndim); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_4 = PyInt_FromLong(__pyx_v_ndim); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_5 = PyTuple_New(2); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  PyTuple_SET_ITEM(__pyx_t_5, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_5, 1, __pyx_t_4);
  __Pyx_GIVEREF(__pyx_t_4);
  __pyx_t_2 = 0;
  __pyx_t_4 = 0;
  __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  PyTuple_SET_ITEM(__pyx_t_4, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_5 = 0;
  __pyx_t_5 = PyDict_New(); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_5));
  __pyx_t_2 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_6 = PyObject_GetAttr(__pyx_t_2, __pyx_n_s__double); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  if (PyDict_SetItem(__pyx_t_5, ((PyObject *)__pyx_n_s__dtype), __pyx_t_6) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;
  __pyx_t_6 = PyEval_CallObjectWithKeywords(__pyx_t_1, __pyx_t_4, ((PyObject *)__pyx_t_5)); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_5)); __pyx_t_5 = 0;
  if (!(likely(((__pyx_t_6) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_6, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_7 = ((PyArrayObject *)__pyx_t_6);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_T);
    __pyx_t_8 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_T, (PyObject*)__pyx_t_7, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_8 < 0)) {
      PyErr_Fetch(&__pyx_t_9, &__pyx_t_10, &__pyx_t_11);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_T, (PyObject*)__pyx_v_T, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_10); Py_XDECREF(__pyx_t_11);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_9, __pyx_t_10, __pyx_t_11);
      }
    }
    __pyx_bstride_0_T = __pyx_bstruct_T.strides[0]; __pyx_bstride_1_T = __pyx_bstruct_T.strides[1];
    __pyx_bshape_0_T = __pyx_bstruct_T.shape[0]; __pyx_bshape_1_T = __pyx_bstruct_T.shape[1];
    if (unlikely(__pyx_t_8 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 328; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_7 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_T));
  __pyx_v_T = ((PyArrayObject *)__pyx_t_6);
  __pyx_t_6 = 0;

  
  __pyx_t_6 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __pyx_t_5 = PyObject_GetAttr(__pyx_t_6, __pyx_n_s__zeros); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;
  __pyx_t_6 = PyInt_FromLong(__pyx_v_nvertex); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __pyx_t_4 = PyInt_FromLong((__pyx_v_ndim + 1)); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_1 = PyInt_FromLong(__pyx_v_ndim); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyTuple_New(3); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_t_6);
  __Pyx_GIVEREF(__pyx_t_6);
  PyTuple_SET_ITEM(__pyx_t_2, 1, __pyx_t_4);
  __Pyx_GIVEREF(__pyx_t_4);
  PyTuple_SET_ITEM(__pyx_t_2, 2, __pyx_t_1);
  __Pyx_GIVEREF(__pyx_t_1);
  __pyx_t_6 = 0;
  __pyx_t_4 = 0;
  __pyx_t_1 = 0;
  __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyDict_New(); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_2));
  __pyx_t_4 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_6 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__double); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyDict_SetItem(__pyx_t_2, ((PyObject *)__pyx_n_s__dtype), __pyx_t_6) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;
  __pyx_t_6 = PyEval_CallObjectWithKeywords(__pyx_t_5, __pyx_t_1, ((PyObject *)__pyx_t_2)); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_2)); __pyx_t_2 = 0;
  if (!(likely(((__pyx_t_6) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_6, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_12 = ((PyArrayObject *)__pyx_t_6);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_Tinvs);
    __pyx_t_8 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_Tinvs, (PyObject*)__pyx_t_12, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 3, 0, __pyx_stack);
    if (unlikely(__pyx_t_8 < 0)) {
      PyErr_Fetch(&__pyx_t_11, &__pyx_t_10, &__pyx_t_9);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_Tinvs, (PyObject*)__pyx_v_Tinvs, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 3, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_11); Py_XDECREF(__pyx_t_10); Py_XDECREF(__pyx_t_9);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_11, __pyx_t_10, __pyx_t_9);
      }
    }
    __pyx_bstride_0_Tinvs = __pyx_bstruct_Tinvs.strides[0]; __pyx_bstride_1_Tinvs = __pyx_bstruct_Tinvs.strides[1]; __pyx_bstride_2_Tinvs = __pyx_bstruct_Tinvs.strides[2];
    __pyx_bshape_0_Tinvs = __pyx_bstruct_Tinvs.shape[0]; __pyx_bshape_1_Tinvs = __pyx_bstruct_Tinvs.shape[1]; __pyx_bshape_2_Tinvs = __pyx_bstruct_Tinvs.shape[2];
    if (unlikely(__pyx_t_8 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 329; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_12 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_Tinvs));
  __pyx_v_Tinvs = ((PyArrayObject *)__pyx_t_6);
  __pyx_t_6 = 0;

  
  __pyx_t_8 = __pyx_v_nvertex;
  for (__pyx_t_13 = 0; __pyx_t_13 < __pyx_t_8; __pyx_t_13+=1) {
    __pyx_v_ivertex = __pyx_t_13;

    
    __pyx_t_14 = (__pyx_v_ndim == 2);
    if (__pyx_t_14) {

      
      __pyx_t_15 = __pyx_v_ivertex;
      __pyx_t_16 = 0;
      if (__pyx_t_15 < 0) __pyx_t_15 += __pyx_bshape_0_vertices;
      if (__pyx_t_16 < 0) __pyx_t_16 += __pyx_bshape_1_vertices;
      __pyx_t_17 = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_15, __pyx_bstride_0_vertices, __pyx_t_16, __pyx_bstride_1_vertices));
      __pyx_t_18 = 0;
      if (__pyx_t_17 < 0) __pyx_t_17 += __pyx_bshape_0_points;
      if (__pyx_t_18 < 0) __pyx_t_18 += __pyx_bshape_1_points;
      __pyx_v_x1 = (*__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_points.buf, __pyx_t_17, __pyx_bstride_0_points, __pyx_t_18, __pyx_bstride_1_points));

      
      __pyx_t_19 = __pyx_v_ivertex;
      __pyx_t_20 = 1;
      if (__pyx_t_19 < 0) __pyx_t_19 += __pyx_bshape_0_vertices;
      if (__pyx_t_20 < 0) __pyx_t_20 += __pyx_bshape_1_vertices;
      __pyx_t_21 = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_19, __pyx_bstride_0_vertices, __pyx_t_20, __pyx_bstride_1_vertices));
      __pyx_t_22 = 0;
      if (__pyx_t_21 < 0) __pyx_t_21 += __pyx_bshape_0_points;
      if (__pyx_t_22 < 0) __pyx_t_22 += __pyx_bshape_1_points;
      __pyx_v_x2 = (*__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_points.buf, __pyx_t_21, __pyx_bstride_0_points, __pyx_t_22, __pyx_bstride_1_points));

      
      __pyx_t_23 = __pyx_v_ivertex;
      __pyx_t_24 = 2;
      if (__pyx_t_23 < 0) __pyx_t_23 += __pyx_bshape_0_vertices;
      if (__pyx_t_24 < 0) __pyx_t_24 += __pyx_bshape_1_vertices;
      __pyx_t_25 = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_23, __pyx_bstride_0_vertices, __pyx_t_24, __pyx_bstride_1_vertices));
      __pyx_t_26 = 0;
      if (__pyx_t_25 < 0) __pyx_t_25 += __pyx_bshape_0_points;
      if (__pyx_t_26 < 0) __pyx_t_26 += __pyx_bshape_1_points;
      __pyx_v_x3 = (*__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_points.buf, __pyx_t_25, __pyx_bstride_0_points, __pyx_t_26, __pyx_bstride_1_points));

      
      __pyx_t_27 = __pyx_v_ivertex;
      __pyx_t_28 = 0;
      if (__pyx_t_27 < 0) __pyx_t_27 += __pyx_bshape_0_vertices;
      if (__pyx_t_28 < 0) __pyx_t_28 += __pyx_bshape_1_vertices;
      __pyx_t_29 = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_27, __pyx_bstride_0_vertices, __pyx_t_28, __pyx_bstride_1_vertices));
      __pyx_t_30 = 1;
      if (__pyx_t_29 < 0) __pyx_t_29 += __pyx_bshape_0_points;
      if (__pyx_t_30 < 0) __pyx_t_30 += __pyx_bshape_1_points;
      __pyx_v_y1 = (*__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_points.buf, __pyx_t_29, __pyx_bstride_0_points, __pyx_t_30, __pyx_bstride_1_points));

      
      __pyx_t_31 = __pyx_v_ivertex;
      __pyx_t_32 = 1;
      if (__pyx_t_31 < 0) __pyx_t_31 += __pyx_bshape_0_vertices;
      if (__pyx_t_32 < 0) __pyx_t_32 += __pyx_bshape_1_vertices;
      __pyx_t_33 = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_31, __pyx_bstride_0_vertices, __pyx_t_32, __pyx_bstride_1_vertices));
      __pyx_t_34 = 1;
      if (__pyx_t_33 < 0) __pyx_t_33 += __pyx_bshape_0_points;
      if (__pyx_t_34 < 0) __pyx_t_34 += __pyx_bshape_1_points;
      __pyx_v_y2 = (*__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_points.buf, __pyx_t_33, __pyx_bstride_0_points, __pyx_t_34, __pyx_bstride_1_points));

      
      __pyx_t_35 = __pyx_v_ivertex;
      __pyx_t_36 = 2;
      if (__pyx_t_35 < 0) __pyx_t_35 += __pyx_bshape_0_vertices;
      if (__pyx_t_36 < 0) __pyx_t_36 += __pyx_bshape_1_vertices;
      __pyx_t_37 = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_35, __pyx_bstride_0_vertices, __pyx_t_36, __pyx_bstride_1_vertices));
      __pyx_t_38 = 1;
      if (__pyx_t_37 < 0) __pyx_t_37 += __pyx_bshape_0_points;
      if (__pyx_t_38 < 0) __pyx_t_38 += __pyx_bshape_1_points;
      __pyx_v_y3 = (*__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_points.buf, __pyx_t_37, __pyx_bstride_0_points, __pyx_t_38, __pyx_bstride_1_points));

      
      __pyx_v_x1 -= __pyx_v_x3;

      
      __pyx_v_x2 -= __pyx_v_x3;

      
      __pyx_v_y1 -= __pyx_v_y3;

      
      __pyx_v_y2 -= __pyx_v_y3;

      
      __pyx_v_det = ((__pyx_v_x1 * __pyx_v_y2) - (__pyx_v_x2 * __pyx_v_y1));

      
      __pyx_t_14 = (__pyx_v_det == 0.0);
      if (__pyx_t_14) {

        
        __pyx_v_info = 1;
        goto __pyx_L9;
      }
       {

        
        __pyx_v_info = 0;

        
        if (unlikely(__pyx_v_det == 0)) {
          PyErr_Format(PyExc_ZeroDivisionError, "float division");
          {__pyx_filename = __pyx_f[0]; __pyx_lineno = 359; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        }
        __pyx_t_39 = __pyx_v_ivertex;
        __pyx_t_40 = 0;
        __pyx_t_41 = 0;
        if (__pyx_t_39 < 0) __pyx_t_39 += __pyx_bshape_0_Tinvs;
        if (__pyx_t_40 < 0) __pyx_t_40 += __pyx_bshape_1_Tinvs;
        if (__pyx_t_41 < 0) __pyx_t_41 += __pyx_bshape_2_Tinvs;
        *__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_39, __pyx_bstride_0_Tinvs, __pyx_t_40, __pyx_bstride_1_Tinvs, __pyx_t_41, __pyx_bstride_2_Tinvs) = (__pyx_v_y2 / __pyx_v_det);

        
        __pyx_t_3 = (-__pyx_v_x2);
        if (unlikely(__pyx_v_det == 0)) {
          PyErr_Format(PyExc_ZeroDivisionError, "float division");
          {__pyx_filename = __pyx_f[0]; __pyx_lineno = 360; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        }
        __pyx_t_42 = __pyx_v_ivertex;
        __pyx_t_43 = 0;
        __pyx_t_44 = 1;
        if (__pyx_t_42 < 0) __pyx_t_42 += __pyx_bshape_0_Tinvs;
        if (__pyx_t_43 < 0) __pyx_t_43 += __pyx_bshape_1_Tinvs;
        if (__pyx_t_44 < 0) __pyx_t_44 += __pyx_bshape_2_Tinvs;
        *__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_42, __pyx_bstride_0_Tinvs, __pyx_t_43, __pyx_bstride_1_Tinvs, __pyx_t_44, __pyx_bstride_2_Tinvs) = (__pyx_t_3 / __pyx_v_det);

        
        __pyx_t_3 = (-__pyx_v_y1);
        if (unlikely(__pyx_v_det == 0)) {
          PyErr_Format(PyExc_ZeroDivisionError, "float division");
          {__pyx_filename = __pyx_f[0]; __pyx_lineno = 361; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        }
        __pyx_t_45 = __pyx_v_ivertex;
        __pyx_t_46 = 1;
        __pyx_t_47 = 0;
        if (__pyx_t_45 < 0) __pyx_t_45 += __pyx_bshape_0_Tinvs;
        if (__pyx_t_46 < 0) __pyx_t_46 += __pyx_bshape_1_Tinvs;
        if (__pyx_t_47 < 0) __pyx_t_47 += __pyx_bshape_2_Tinvs;
        *__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_45, __pyx_bstride_0_Tinvs, __pyx_t_46, __pyx_bstride_1_Tinvs, __pyx_t_47, __pyx_bstride_2_Tinvs) = (__pyx_t_3 / __pyx_v_det);

        
        if (unlikely(__pyx_v_det == 0)) {
          PyErr_Format(PyExc_ZeroDivisionError, "float division");
          {__pyx_filename = __pyx_f[0]; __pyx_lineno = 362; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        }
        __pyx_t_48 = __pyx_v_ivertex;
        __pyx_t_49 = 1;
        __pyx_t_50 = 1;
        if (__pyx_t_48 < 0) __pyx_t_48 += __pyx_bshape_0_Tinvs;
        if (__pyx_t_49 < 0) __pyx_t_49 += __pyx_bshape_1_Tinvs;
        if (__pyx_t_50 < 0) __pyx_t_50 += __pyx_bshape_2_Tinvs;
        *__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_48, __pyx_bstride_0_Tinvs, __pyx_t_49, __pyx_bstride_1_Tinvs, __pyx_t_50, __pyx_bstride_2_Tinvs) = (__pyx_v_x1 / __pyx_v_det);

        
        __pyx_t_51 = __pyx_v_ivertex;
        __pyx_t_52 = 2;
        __pyx_t_53 = 0;
        if (__pyx_t_51 < 0) __pyx_t_51 += __pyx_bshape_0_Tinvs;
        if (__pyx_t_52 < 0) __pyx_t_52 += __pyx_bshape_1_Tinvs;
        if (__pyx_t_53 < 0) __pyx_t_53 += __pyx_bshape_2_Tinvs;
        *__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_51, __pyx_bstride_0_Tinvs, __pyx_t_52, __pyx_bstride_1_Tinvs, __pyx_t_53, __pyx_bstride_2_Tinvs) = __pyx_v_x3;

        
        __pyx_t_54 = __pyx_v_ivertex;
        __pyx_t_55 = 2;
        __pyx_t_56 = 1;
        if (__pyx_t_54 < 0) __pyx_t_54 += __pyx_bshape_0_Tinvs;
        if (__pyx_t_55 < 0) __pyx_t_55 += __pyx_bshape_1_Tinvs;
        if (__pyx_t_56 < 0) __pyx_t_56 += __pyx_bshape_2_Tinvs;
        *__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_54, __pyx_bstride_0_Tinvs, __pyx_t_55, __pyx_bstride_1_Tinvs, __pyx_t_56, __pyx_bstride_2_Tinvs) = __pyx_v_y3;
      }
      __pyx_L9:;
      goto __pyx_L8;
    }
     {

      
      __pyx_t_57 = __pyx_v_ndim;
      for (__pyx_t_58 = 0; __pyx_t_58 < __pyx_t_57; __pyx_t_58+=1) {
        __pyx_v_i = __pyx_t_58;

        
        __pyx_t_59 = __pyx_v_ivertex;
        __pyx_t_60 = __pyx_v_ndim;
        if (__pyx_t_59 < 0) __pyx_t_59 += __pyx_bshape_0_vertices;
        if (__pyx_t_60 < 0) __pyx_t_60 += __pyx_bshape_1_vertices;
        __pyx_t_61 = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_59, __pyx_bstride_0_vertices, __pyx_t_60, __pyx_bstride_1_vertices));
        __pyx_t_62 = __pyx_v_i;
        if (__pyx_t_61 < 0) __pyx_t_61 += __pyx_bshape_0_points;
        if (__pyx_t_62 < 0) __pyx_t_62 += __pyx_bshape_1_points;
        __pyx_t_63 = __pyx_v_ivertex;
        __pyx_t_64 = __pyx_v_ndim;
        __pyx_t_65 = __pyx_v_i;
        if (__pyx_t_63 < 0) __pyx_t_63 += __pyx_bshape_0_Tinvs;
        if (__pyx_t_64 < 0) __pyx_t_64 += __pyx_bshape_1_Tinvs;
        if (__pyx_t_65 < 0) __pyx_t_65 += __pyx_bshape_2_Tinvs;
        *__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_63, __pyx_bstride_0_Tinvs, __pyx_t_64, __pyx_bstride_1_Tinvs, __pyx_t_65, __pyx_bstride_2_Tinvs) = (*__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_points.buf, __pyx_t_61, __pyx_bstride_0_points, __pyx_t_62, __pyx_bstride_1_points));

        
        __pyx_t_66 = __pyx_v_ndim;
        for (__pyx_t_67 = 0; __pyx_t_67 < __pyx_t_66; __pyx_t_67+=1) {
          __pyx_v_j = __pyx_t_67;

          
          __pyx_t_68 = __pyx_v_ivertex;
          __pyx_t_69 = __pyx_v_j;
          if (__pyx_t_68 < 0) __pyx_t_68 += __pyx_bshape_0_vertices;
          if (__pyx_t_69 < 0) __pyx_t_69 += __pyx_bshape_1_vertices;
          __pyx_t_70 = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_68, __pyx_bstride_0_vertices, __pyx_t_69, __pyx_bstride_1_vertices));
          __pyx_t_71 = __pyx_v_i;
          if (__pyx_t_70 < 0) __pyx_t_70 += __pyx_bshape_0_points;
          if (__pyx_t_71 < 0) __pyx_t_71 += __pyx_bshape_1_points;

          
          __pyx_t_72 = __pyx_v_ivertex;
          __pyx_t_73 = __pyx_v_ndim;
          __pyx_t_74 = __pyx_v_i;
          if (__pyx_t_72 < 0) __pyx_t_72 += __pyx_bshape_0_Tinvs;
          if (__pyx_t_73 < 0) __pyx_t_73 += __pyx_bshape_1_Tinvs;
          if (__pyx_t_74 < 0) __pyx_t_74 += __pyx_bshape_2_Tinvs;

          
          __pyx_t_75 = __pyx_v_i;
          __pyx_t_76 = __pyx_v_j;
          if (__pyx_t_75 < 0) __pyx_t_75 += __pyx_bshape_0_T;
          if (__pyx_t_76 < 0) __pyx_t_76 += __pyx_bshape_1_T;
          *__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_T.buf, __pyx_t_75, __pyx_bstride_0_T, __pyx_t_76, __pyx_bstride_1_T) = ((*__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_points.buf, __pyx_t_70, __pyx_bstride_0_points, __pyx_t_71, __pyx_bstride_1_points)) - (*__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_72, __pyx_bstride_0_Tinvs, __pyx_t_73, __pyx_bstride_1_Tinvs, __pyx_t_74, __pyx_bstride_2_Tinvs)));
        }

        
        __pyx_t_66 = __pyx_v_ivertex;
        __pyx_t_67 = __pyx_v_i;
        __pyx_t_77 = __pyx_v_i;
        if (__pyx_t_66 < 0) __pyx_t_66 += __pyx_bshape_0_Tinvs;
        if (__pyx_t_67 < 0) __pyx_t_67 += __pyx_bshape_1_Tinvs;
        if (__pyx_t_77 < 0) __pyx_t_77 += __pyx_bshape_2_Tinvs;
        *__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_66, __pyx_bstride_0_Tinvs, __pyx_t_67, __pyx_bstride_1_Tinvs, __pyx_t_77, __pyx_bstride_2_Tinvs) = 1.0;
      }

      
      __pyx_v_n = __pyx_v_ndim;

      
      __pyx_v_nrhs = __pyx_v_ndim;

      
      __pyx_v_lda = __pyx_v_ndim;

      
      __pyx_v_ldb = __pyx_v_ndim;

      
      qh_dgesv((&__pyx_v_n), (&__pyx_v_nrhs), ((double *)__pyx_v_T->data), (&__pyx_v_lda), __pyx_v_ipiv, (((double *)__pyx_v_Tinvs->data) + ((__pyx_v_ndim * (__pyx_v_ndim + 1)) * __pyx_v_ivertex)), (&__pyx_v_ldb), (&__pyx_v_info));
    }
    __pyx_L8:;

    
    __pyx_t_14 = (__pyx_v_info != 0);
    if (__pyx_t_14) {

      
      __pyx_t_78 = (__pyx_v_ndim + 1);
      for (__pyx_t_57 = 0; __pyx_t_57 < __pyx_t_78; __pyx_t_57+=1) {
        __pyx_v_i = __pyx_t_57;

        
        __pyx_t_58 = __pyx_v_ndim;
        for (__pyx_t_79 = 0; __pyx_t_79 < __pyx_t_58; __pyx_t_79+=1) {
          __pyx_v_j = __pyx_t_79;

          
          __pyx_t_80 = __pyx_v_ivertex;
          __pyx_t_81 = __pyx_v_i;
          __pyx_t_82 = __pyx_v_j;
          if (__pyx_t_80 < 0) __pyx_t_80 += __pyx_bshape_0_Tinvs;
          if (__pyx_t_81 < 0) __pyx_t_81 += __pyx_bshape_1_Tinvs;
          if (__pyx_t_82 < 0) __pyx_t_82 += __pyx_bshape_2_Tinvs;
          *__Pyx_BufPtrStrided3d(__pyx_t_5numpy_double_t *, __pyx_bstruct_Tinvs.buf, __pyx_t_80, __pyx_bstride_0_Tinvs, __pyx_t_81, __pyx_bstride_1_Tinvs, __pyx_t_82, __pyx_bstride_2_Tinvs) = __pyx_v_nan;
        }
      }
      goto __pyx_L14;
    }
    __pyx_L14:;
  }

  
  __Pyx_XDECREF(__pyx_r);
  __Pyx_INCREF(((PyObject *)__pyx_v_Tinvs));
  __pyx_r = ((PyObject *)__pyx_v_Tinvs);
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_XDECREF(__pyx_t_6);
  { PyObject *__pyx_type, *__pyx_value, *__pyx_tb;
    __Pyx_ErrFetch(&__pyx_type, &__pyx_value, &__pyx_tb);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_Tinvs);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_T);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_points);
  __Pyx_ErrRestore(__pyx_type, __pyx_value, __pyx_tb);}
  __Pyx_AddTraceback("scipy.spatial.qhull._get_barycentric_transforms");
  __pyx_r = NULL;
  goto __pyx_L2;
  __pyx_L0:;
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_Tinvs);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_T);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_points);
  __pyx_L2:;
  __Pyx_DECREF((PyObject *)__pyx_v_T);
  __Pyx_DECREF((PyObject *)__pyx_v_Tinvs);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static  int __pyx_f_5scipy_7spatial_5qhull__barycentric_inside(int __pyx_v_ndim, double *__pyx_v_transform, double *__pyx_v_x, double *__pyx_v_c, double __pyx_v_eps) {
  int __pyx_v_i;
  int __pyx_v_j;
  int __pyx_r;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  int __pyx_t_4;
  double __pyx_t_5;
  int __pyx_t_6;
  int __pyx_t_7;

  
  (__pyx_v_c[__pyx_v_ndim]) = 1.0;

  
  __pyx_t_1 = __pyx_v_ndim;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    (__pyx_v_c[__pyx_v_i]) = 0.0;

    
    __pyx_t_3 = __pyx_v_ndim;
    for (__pyx_t_4 = 0; __pyx_t_4 < __pyx_t_3; __pyx_t_4+=1) {
      __pyx_v_j = __pyx_t_4;

      
      (__pyx_v_c[__pyx_v_i]) += ((__pyx_v_transform[((__pyx_v_ndim * __pyx_v_i) + __pyx_v_j)]) * ((__pyx_v_x[__pyx_v_j]) - (__pyx_v_transform[((__pyx_v_ndim * __pyx_v_ndim) + __pyx_v_j)])));
    }

    
    (__pyx_v_c[__pyx_v_ndim]) -= (__pyx_v_c[__pyx_v_i]);

    
    __pyx_t_5 = (__pyx_v_c[__pyx_v_i]);
    __pyx_t_6 = ((-__pyx_v_eps) <= __pyx_t_5);
    if (__pyx_t_6) {
      __pyx_t_6 = (__pyx_t_5 <= (1.0 + __pyx_v_eps));
    }
    __pyx_t_7 = (!__pyx_t_6);
    if (__pyx_t_7) {

      
      __pyx_r = 0;
      goto __pyx_L0;
      goto __pyx_L7;
    }
    __pyx_L7:;
  }

  
  __pyx_t_5 = (__pyx_v_c[__pyx_v_ndim]);
  __pyx_t_7 = ((-__pyx_v_eps) <= __pyx_t_5);
  if (__pyx_t_7) {
    __pyx_t_7 = (__pyx_t_5 <= (1.0 + __pyx_v_eps));
  }
  __pyx_t_6 = (!__pyx_t_7);
  if (__pyx_t_6) {

    
    __pyx_r = 0;
    goto __pyx_L0;
    goto __pyx_L8;
  }
  __pyx_L8:;

  
  __pyx_r = 1;
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static  void __pyx_f_5scipy_7spatial_5qhull__barycentric_coordinate_single(int __pyx_v_ndim, double *__pyx_v_transform, double *__pyx_v_x, double *__pyx_v_c, int __pyx_v_i) {
  int __pyx_v_j;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;

  
  __pyx_t_1 = (__pyx_v_i == __pyx_v_ndim);
  if (__pyx_t_1) {

    
    (__pyx_v_c[__pyx_v_ndim]) = 1.0;

    
    __pyx_t_2 = __pyx_v_ndim;
    for (__pyx_t_3 = 0; __pyx_t_3 < __pyx_t_2; __pyx_t_3+=1) {
      __pyx_v_j = __pyx_t_3;

      
      (__pyx_v_c[__pyx_v_ndim]) -= (__pyx_v_c[__pyx_v_j]);
    }
    goto __pyx_L3;
  }
   {

    
    (__pyx_v_c[__pyx_v_i]) = 0.0;

    
    __pyx_t_2 = __pyx_v_ndim;
    for (__pyx_t_3 = 0; __pyx_t_3 < __pyx_t_2; __pyx_t_3+=1) {
      __pyx_v_j = __pyx_t_3;

      
      (__pyx_v_c[__pyx_v_i]) += ((__pyx_v_transform[((__pyx_v_ndim * __pyx_v_i) + __pyx_v_j)]) * ((__pyx_v_x[__pyx_v_j]) - (__pyx_v_transform[((__pyx_v_ndim * __pyx_v_ndim) + __pyx_v_j)])));
    }
  }
  __pyx_L3:;

}



static  void __pyx_f_5scipy_7spatial_5qhull__barycentric_coordinates(int __pyx_v_ndim, double *__pyx_v_transform, double *__pyx_v_x, double *__pyx_v_c) {
  int __pyx_v_i;
  int __pyx_v_j;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  int __pyx_t_4;

  
  (__pyx_v_c[__pyx_v_ndim]) = 1.0;

  
  __pyx_t_1 = __pyx_v_ndim;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    (__pyx_v_c[__pyx_v_i]) = 0.0;

    
    __pyx_t_3 = __pyx_v_ndim;
    for (__pyx_t_4 = 0; __pyx_t_4 < __pyx_t_3; __pyx_t_4+=1) {
      __pyx_v_j = __pyx_t_4;

      
      (__pyx_v_c[__pyx_v_i]) += ((__pyx_v_transform[((__pyx_v_ndim * __pyx_v_i) + __pyx_v_j)]) * ((__pyx_v_x[__pyx_v_j]) - (__pyx_v_transform[((__pyx_v_ndim * __pyx_v_ndim) + __pyx_v_j)])));
    }

    
    (__pyx_v_c[__pyx_v_ndim]) -= (__pyx_v_c[__pyx_v_i]);
  }

}



static  void __pyx_f_5scipy_7spatial_5qhull__lift_point(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *__pyx_v_d, double *__pyx_v_x, double *__pyx_v_z) {
  int __pyx_v_i;
  int __pyx_t_1;
  int __pyx_t_2;

  
  (__pyx_v_z[__pyx_v_d->ndim]) = 0.0;

  
  __pyx_t_1 = __pyx_v_d->ndim;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    (__pyx_v_z[__pyx_v_i]) = (__pyx_v_x[__pyx_v_i]);

    
    (__pyx_v_z[__pyx_v_d->ndim]) += pow((__pyx_v_x[__pyx_v_i]), 2.0);
  }

  
  (__pyx_v_z[__pyx_v_d->ndim]) *= __pyx_v_d->paraboloid_scale;

  
  (__pyx_v_z[__pyx_v_d->ndim]) += __pyx_v_d->paraboloid_shift;

}



static  double __pyx_f_5scipy_7spatial_5qhull__distplane(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *__pyx_v_d, int __pyx_v_isimplex, double *__pyx_v_point) {
  double __pyx_v_dist;
  int __pyx_v_k;
  double __pyx_r;
  long __pyx_t_1;
  int __pyx_t_2;

  
  __pyx_v_dist = (__pyx_v_d->equations[(((__pyx_v_isimplex * (__pyx_v_d->ndim + 2)) + __pyx_v_d->ndim) + 1)]);

  
  __pyx_t_1 = (__pyx_v_d->ndim + 1);
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_k = __pyx_t_2;

    
    __pyx_v_dist += ((__pyx_v_d->equations[((__pyx_v_isimplex * (__pyx_v_d->ndim + 2)) + __pyx_v_k)]) * (__pyx_v_point[__pyx_v_k]));
  }

  
  __pyx_r = __pyx_v_dist;
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static  void __pyx_f_5scipy_7spatial_5qhull__RidgeIter2D_init(__pyx_t_5scipy_7spatial_5qhull_RidgeIter2D_t *__pyx_v_it, __pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *__pyx_v_d, int __pyx_v_vertex) {
  int __pyx_v_k;
  int __pyx_v_ivertex;
  int __pyx_v_start;
  int __pyx_t_1;
  int __pyx_t_2;

  
  __pyx_v_start = 0;

  
  __pyx_v_it->info = __pyx_v_d;

  
  __pyx_v_it->vertex = __pyx_v_vertex;

  
  __pyx_v_it->triangle = (__pyx_v_d->vertex_to_simplex[__pyx_v_vertex]);

  
  __pyx_v_it->start_triangle = __pyx_v_it->triangle;

  
  __pyx_v_it->restart = 0;

  
  __pyx_t_1 = (__pyx_v_it->triangle != -1);
  if (__pyx_t_1) {

    
    for (__pyx_t_2 = 0; __pyx_t_2 < 3; __pyx_t_2+=1) {
      __pyx_v_k = __pyx_t_2;

      
      __pyx_v_ivertex = (__pyx_v_it->info->vertices[((__pyx_v_it->triangle * 3) + __pyx_v_k)]);

      
      __pyx_t_1 = (__pyx_v_ivertex != __pyx_v_vertex);
      if (__pyx_t_1) {

        
        __pyx_v_it->vertex2 = __pyx_v_ivertex;

        
        __pyx_v_it->index = __pyx_v_k;

        
        __pyx_v_it->start_index = __pyx_v_k;

        
        goto __pyx_L5_break;
        goto __pyx_L6;
      }
      __pyx_L6:;
    }
    __pyx_L5_break:;
    goto __pyx_L3;
  }
   {

    
    __pyx_v_it->start_index = -1;

    
    __pyx_v_it->index = -1;
  }
  __pyx_L3:;

}



static  void __pyx_f_5scipy_7spatial_5qhull__RidgeIter2D_next(__pyx_t_5scipy_7spatial_5qhull_RidgeIter2D_t *__pyx_v_it) {
  int __pyx_v_itri;
  int __pyx_v_k;
  int __pyx_v_ivertex;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  int __pyx_t_4;

  
  if (__pyx_v_it->restart) {

    
    __pyx_t_1 = (__pyx_v_it->start_index == -1);
    if (__pyx_t_1) {

      
      __pyx_v_it->index = -1;

      
      goto __pyx_L0;
      goto __pyx_L4;
    }
    __pyx_L4:;

    
    __pyx_v_it->triangle = __pyx_v_it->start_triangle;

    
    for (__pyx_t_2 = 0; __pyx_t_2 < 3; __pyx_t_2+=1) {
      __pyx_v_k = __pyx_t_2;

      
      __pyx_v_ivertex = (__pyx_v_it->info->vertices[((__pyx_v_it->triangle * 3) + __pyx_v_k)]);

      
      __pyx_t_1 = (__pyx_v_ivertex != __pyx_v_it->vertex);
      if (__pyx_t_1) {
        __pyx_t_3 = (__pyx_v_k != __pyx_v_it->start_index);
        __pyx_t_4 = __pyx_t_3;
      } else {
        __pyx_t_4 = __pyx_t_1;
      }
      if (__pyx_t_4) {

        
        __pyx_v_it->index = __pyx_v_k;

        
        __pyx_v_it->vertex2 = __pyx_v_ivertex;

        
        goto __pyx_L6_break;
        goto __pyx_L7;
      }
      __pyx_L7:;
    }
    __pyx_L6_break:;

    
    __pyx_v_it->start_index = -1;

    
    __pyx_v_it->restart = 0;

    
    __pyx_t_4 = ((__pyx_v_it->info->neighbors[((__pyx_v_it->triangle * 3) + __pyx_v_it->index)]) == -1);
    if (__pyx_t_4) {

      
      __pyx_v_it->index = -1;

      
      goto __pyx_L0;
      goto __pyx_L8;
    }
     {

      
      __pyx_f_5scipy_7spatial_5qhull__RidgeIter2D_next(__pyx_v_it);

      
      __pyx_t_4 = (__pyx_v_it->index == -1);
      if (__pyx_t_4) {

        
        goto __pyx_L0;
        goto __pyx_L9;
      }
      __pyx_L9:;
    }
    __pyx_L8:;
    goto __pyx_L3;
  }
  __pyx_L3:;

  
  __pyx_v_itri = (__pyx_v_it->info->neighbors[((__pyx_v_it->triangle * 3) + __pyx_v_it->index)]);

  
  __pyx_t_4 = (__pyx_v_itri == -1);
  if (__pyx_t_4) {

    
    for (__pyx_t_2 = 0; __pyx_t_2 < 3; __pyx_t_2+=1) {
      __pyx_v_k = __pyx_t_2;

      
      __pyx_v_ivertex = (__pyx_v_it->info->vertices[((__pyx_v_it->triangle * 3) + __pyx_v_k)]);

      
      __pyx_t_4 = (__pyx_v_ivertex != __pyx_v_it->vertex);
      if (__pyx_t_4) {
        __pyx_t_1 = (__pyx_v_k != __pyx_v_it->index);
        __pyx_t_3 = __pyx_t_1;
      } else {
        __pyx_t_3 = __pyx_t_4;
      }
      if (__pyx_t_3) {

        
        __pyx_v_it->index = __pyx_v_k;

        
        __pyx_v_it->vertex2 = __pyx_v_ivertex;

        
        goto __pyx_L12_break;
        goto __pyx_L13;
      }
      __pyx_L13:;
    }
    __pyx_L12_break:;

    
    __pyx_v_it->restart = 1;

    
    goto __pyx_L0;
    goto __pyx_L10;
  }
  __pyx_L10:;

  
  for (__pyx_t_2 = 0; __pyx_t_2 < 3; __pyx_t_2+=1) {
    __pyx_v_k = __pyx_t_2;

    
    __pyx_v_ivertex = (__pyx_v_it->info->vertices[((__pyx_v_itri * 3) + __pyx_v_k)]);

    
    __pyx_t_3 = ((__pyx_v_it->info->neighbors[((__pyx_v_itri * 3) + __pyx_v_k)]) != __pyx_v_it->triangle);
    if (__pyx_t_3) {

      
      __pyx_t_4 = (__pyx_v_ivertex != __pyx_v_it->vertex);
      __pyx_t_1 = __pyx_t_4;
    } else {
      __pyx_t_1 = __pyx_t_3;
    }
    if (__pyx_t_1) {

      
      __pyx_v_it->index = __pyx_v_k;

      
      __pyx_v_it->vertex2 = __pyx_v_ivertex;

      
      goto __pyx_L15_break;
      goto __pyx_L16;
    }
    __pyx_L16:;
  }
  __pyx_L15_break:;

  
  __pyx_v_it->triangle = __pyx_v_itri;

  
  __pyx_t_1 = (__pyx_v_it->triangle == __pyx_v_it->start_triangle);
  if (__pyx_t_1) {

    
    __pyx_v_it->index = -1;

    
    goto __pyx_L0;
    goto __pyx_L17;
  }
  __pyx_L17:;

  __pyx_L0:;
}



static int __pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___init__(PyObject *__pyx_v_self, PyObject *__pyx_args, PyObject *__pyx_kwds); 
static int __pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___init__(PyObject *__pyx_v_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyObject *__pyx_v_delaunay = 0;
  PyObject *__pyx_v_ivertex = 0;
  int __pyx_r;
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  int __pyx_t_3;
  int __pyx_t_4;
  static PyObject **__pyx_pyargnames[] = {&__pyx_n_s__delaunay,&__pyx_n_s__ivertex,0};
  __Pyx_RefNannySetupContext("__init__");
  if (unlikely(__pyx_kwds)) {
    Py_ssize_t kw_args = PyDict_Size(__pyx_kwds);
    PyObject* values[2] = {0,0};
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  2: values[1] = PyTuple_GET_ITEM(__pyx_args, 1);
      case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
      case  0: break;
      default: goto __pyx_L5_argtuple_error;
    }
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  0:
      values[0] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__delaunay);
      if (likely(values[0])) kw_args--;
      else goto __pyx_L5_argtuple_error;
      case  1:
      values[1] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__ivertex);
      if (likely(values[1])) kw_args--;
      else {
        __Pyx_RaiseArgtupleInvalid("__init__", 1, 2, 2, 1); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 598; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
      }
    }
    if (unlikely(kw_args > 0)) {
      if (unlikely(__Pyx_ParseOptionalKeywords(__pyx_kwds, __pyx_pyargnames, 0, values, PyTuple_GET_SIZE(__pyx_args), "__init__") < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 598; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    }
    __pyx_v_delaunay = values[0];
    __pyx_v_ivertex = values[1];
  } else if (PyTuple_GET_SIZE(__pyx_args) != 2) {
    goto __pyx_L5_argtuple_error;
  } else {
    __pyx_v_delaunay = PyTuple_GET_ITEM(__pyx_args, 0);
    __pyx_v_ivertex = PyTuple_GET_ITEM(__pyx_args, 1);
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L5_argtuple_error:;
  __Pyx_RaiseArgtupleInvalid("__init__", 1, 2, 2, PyTuple_GET_SIZE(__pyx_args)); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 598; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  __pyx_L3_error:;
  __Pyx_AddTraceback("scipy.spatial.qhull.RidgeIter2D.__init__");
  __Pyx_RefNannyFinishContext();
  return -1;
  __pyx_L4_argument_unpacking_done:;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_delaunay, __pyx_n_s__ndim); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 599; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyObject_RichCompare(__pyx_t_1, __pyx_int_2, Py_NE); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 599; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_3 = __Pyx_PyObject_IsTrue(__pyx_t_2); if (unlikely(__pyx_t_3 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 599; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  if (__pyx_t_3) {

    
    __pyx_t_2 = PyTuple_New(1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 600; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_2);
    __Pyx_INCREF(((PyObject *)__pyx_kp_s_8));
    PyTuple_SET_ITEM(__pyx_t_2, 0, ((PyObject *)__pyx_kp_s_8));
    __Pyx_GIVEREF(((PyObject *)__pyx_kp_s_8));
    __pyx_t_1 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_2, NULL); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 600; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
    __Pyx_Raise(__pyx_t_1, 0, 0);
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    {__pyx_filename = __pyx_f[0]; __pyx_lineno = 600; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L6;
  }
  __pyx_L6:;

  
  __Pyx_INCREF(__pyx_v_delaunay);
  __Pyx_GIVEREF(__pyx_v_delaunay);
  __Pyx_GOTREF(((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->delaunay);
  __Pyx_DECREF(((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->delaunay);
  ((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->delaunay = __pyx_v_delaunay;

  
  __pyx_f_5scipy_7spatial_5qhull__get_delaunay_info((&((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->info), __pyx_v_delaunay, 0, 1);

  
  __pyx_t_4 = __Pyx_PyInt_AsInt(__pyx_v_ivertex); if (unlikely((__pyx_t_4 == (int)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 603; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_f_5scipy_7spatial_5qhull__RidgeIter2D_init((&((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->it), (&((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->info), __pyx_t_4);

  __pyx_r = 0;
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_AddTraceback("scipy.spatial.qhull.RidgeIter2D.__init__");
  __pyx_r = -1;
  __pyx_L0:;
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___iter__(PyObject *__pyx_v_self); 
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___iter__(PyObject *__pyx_v_self) {
  PyObject *__pyx_r = NULL;
  __Pyx_RefNannySetupContext("__iter__");

  
  __Pyx_XDECREF(__pyx_r);
  __Pyx_INCREF(__pyx_v_self);
  __pyx_r = __pyx_v_self;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___next__(PyObject *__pyx_v_self); 
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___next__(PyObject *__pyx_v_self) {
  PyObject *__pyx_v_ret;
  PyObject *__pyx_r = NULL;
  int __pyx_t_1;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  PyObject *__pyx_t_6 = NULL;
  __Pyx_RefNannySetupContext("__next__");
  __pyx_v_ret = ((PyObject *)Py_None); __Pyx_INCREF(Py_None);

  
  __pyx_t_1 = (((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->it.index == -1);
  if (__pyx_t_1) {

    
    __pyx_t_2 = PyObject_Call(__pyx_builtin_StopIteration, ((PyObject *)__pyx_empty_tuple), NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 610; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_2);
    __Pyx_Raise(__pyx_t_2, 0, 0);
    __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
    {__pyx_filename = __pyx_f[0]; __pyx_lineno = 610; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L5;
  }
  __pyx_L5:;

  
  __pyx_t_2 = PyInt_FromLong(((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->it.vertex); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 611; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_3 = PyInt_FromLong(((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->it.vertex2); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 611; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_4 = PyInt_FromLong(((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->it.index); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 611; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_5 = PyInt_FromLong(((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->it.triangle); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 611; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_6 = PyTuple_New(4); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 611; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  PyTuple_SET_ITEM(__pyx_t_6, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_6, 1, __pyx_t_3);
  __Pyx_GIVEREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_6, 2, __pyx_t_4);
  __Pyx_GIVEREF(__pyx_t_4);
  PyTuple_SET_ITEM(__pyx_t_6, 3, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_2 = 0;
  __pyx_t_3 = 0;
  __pyx_t_4 = 0;
  __pyx_t_5 = 0;
  if (!(likely(PyTuple_CheckExact(__pyx_t_6))||(PyErr_Format(PyExc_TypeError, "Expected tuple, got %.200s", Py_TYPE(__pyx_t_6)->tp_name), 0))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 611; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(((PyObject *)__pyx_v_ret));
  __pyx_v_ret = ((PyObject *)__pyx_t_6);
  __pyx_t_6 = 0;

  
  __pyx_f_5scipy_7spatial_5qhull__RidgeIter2D_next((&((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)__pyx_v_self)->it));

  
  __Pyx_XDECREF(__pyx_r);
  __Pyx_INCREF(((PyObject *)__pyx_v_ret));
  __pyx_r = ((PyObject *)__pyx_v_ret);
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_XDECREF(__pyx_t_6);
  __Pyx_AddTraceback("scipy.spatial.qhull.RidgeIter2D.__next__");
  __pyx_r = NULL;
  __pyx_L0:;
  __Pyx_DECREF(__pyx_v_ret);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static  int __pyx_f_5scipy_7spatial_5qhull__is_point_fully_outside(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *__pyx_v_d, double *__pyx_v_x, double __pyx_v_eps) {
  int __pyx_v_i;
  int __pyx_r;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  int __pyx_t_4;
  int __pyx_t_5;

  
  __pyx_t_1 = __pyx_v_d->ndim;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_i = __pyx_t_2;

    
    __pyx_t_3 = ((__pyx_v_x[__pyx_v_i]) < ((__pyx_v_d->min_bound[__pyx_v_i]) - __pyx_v_eps));
    if (!__pyx_t_3) {
      __pyx_t_4 = ((__pyx_v_x[__pyx_v_i]) > ((__pyx_v_d->max_bound[__pyx_v_i]) + __pyx_v_eps));
      __pyx_t_5 = __pyx_t_4;
    } else {
      __pyx_t_5 = __pyx_t_3;
    }
    if (__pyx_t_5) {

      
      __pyx_r = 1;
      goto __pyx_L0;
      goto __pyx_L5;
    }
    __pyx_L5:;
  }

  
  __pyx_r = 0;
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static  int __pyx_f_5scipy_7spatial_5qhull__find_simplex_bruteforce(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *__pyx_v_d, double *__pyx_v_c, double *__pyx_v_x, double __pyx_v_eps) {
  int __pyx_v_inside;
  int __pyx_v_isimplex;
  int __pyx_r;
  int __pyx_t_1;
  int __pyx_t_2;

  
  __pyx_t_1 = __pyx_f_5scipy_7spatial_5qhull__is_point_fully_outside(__pyx_v_d, __pyx_v_x, __pyx_v_eps);
  if (__pyx_t_1) {

    
    __pyx_r = -1;
    goto __pyx_L0;
    goto __pyx_L3;
  }
  __pyx_L3:;

  
  __pyx_t_1 = __pyx_v_d->nsimplex;
  for (__pyx_t_2 = 0; __pyx_t_2 < __pyx_t_1; __pyx_t_2+=1) {
    __pyx_v_isimplex = __pyx_t_2;

    
    __pyx_v_inside = __pyx_f_5scipy_7spatial_5qhull__barycentric_inside(__pyx_v_d->ndim, (__pyx_v_d->transform + ((__pyx_v_isimplex * __pyx_v_d->ndim) * (__pyx_v_d->ndim + 1))), __pyx_v_x, __pyx_v_c, __pyx_v_eps);

    
    if (__pyx_v_inside) {

      
      __pyx_r = __pyx_v_isimplex;
      goto __pyx_L0;
      goto __pyx_L6;
    }
    __pyx_L6:;
  }

  
  __pyx_r = -1;
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static  int __pyx_f_5scipy_7spatial_5qhull__find_simplex_directed(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *__pyx_v_d, double *__pyx_v_c, double *__pyx_v_x, int *__pyx_v_start, double __pyx_v_eps) {
  int __pyx_v_k;
  int __pyx_v_m;
  int __pyx_v_ndim;
  int __pyx_v_inside;
  int __pyx_v_isimplex;
  double *__pyx_v_transform;
  int __pyx_r;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  long __pyx_t_4;
  int __pyx_t_5;

  
  __pyx_v_ndim = __pyx_v_d->ndim;

  
  __pyx_v_isimplex = (__pyx_v_start[0]);

  
  __pyx_t_1 = (__pyx_v_isimplex < 0);
  if (!__pyx_t_1) {
    __pyx_t_2 = (__pyx_v_isimplex >= __pyx_v_d->nsimplex);
    __pyx_t_3 = __pyx_t_2;
  } else {
    __pyx_t_3 = __pyx_t_1;
  }
  if (__pyx_t_3) {

    
    __pyx_v_isimplex = 0;
    goto __pyx_L3;
  }
  __pyx_L3:;

  
  while (1) {
    __pyx_t_3 = (__pyx_v_isimplex != -1);
    if (!__pyx_t_3) break;

    
    __pyx_v_transform = (__pyx_v_d->transform + ((__pyx_v_isimplex * __pyx_v_ndim) * (__pyx_v_ndim + 1)));

    
    __pyx_v_inside = 1;

    
    __pyx_t_4 = (__pyx_v_ndim + 1);
    for (__pyx_t_5 = 0; __pyx_t_5 < __pyx_t_4; __pyx_t_5+=1) {
      __pyx_v_k = __pyx_t_5;

      
      __pyx_f_5scipy_7spatial_5qhull__barycentric_coordinate_single(__pyx_v_ndim, __pyx_v_transform, __pyx_v_x, __pyx_v_c, __pyx_v_k);

      
      __pyx_t_3 = ((__pyx_v_c[__pyx_v_k]) < (-__pyx_v_eps));
      if (__pyx_t_3) {

        
        __pyx_v_m = (__pyx_v_d->neighbors[(((__pyx_v_ndim + 1) * __pyx_v_isimplex) + __pyx_v_k)]);

        
        __pyx_t_3 = (__pyx_v_m == -1);
        if (__pyx_t_3) {

          
          (__pyx_v_start[0]) = __pyx_v_isimplex;

          
          __pyx_r = -1;
          goto __pyx_L0;
          goto __pyx_L9;
        }
        __pyx_L9:;

        
        __pyx_v_isimplex = __pyx_v_m;

        
        __pyx_v_inside = -1;

        
        goto __pyx_L7_break;
        goto __pyx_L8;
      }

      
      __pyx_t_3 = ((__pyx_v_c[__pyx_v_k]) <= (1.0 + __pyx_v_eps));
      if (__pyx_t_3) {
        goto __pyx_L8;
      }
       {

        
        __pyx_v_inside = 0;
      }
      __pyx_L8:;
    }
    __pyx_L7_break:;

    
    switch (__pyx_v_inside) {
      case -1:

      
      goto __pyx_L4_continue;
      break;

      
      case 1:

      
      goto __pyx_L5_break;
      break;
      default:

      
      __pyx_v_isimplex = __pyx_f_5scipy_7spatial_5qhull__find_simplex_bruteforce(__pyx_v_d, __pyx_v_c, __pyx_v_x, __pyx_v_eps);

      
      goto __pyx_L5_break;
      break;
    }
    __pyx_L4_continue:;
  }
  __pyx_L5_break:;

  
  (__pyx_v_start[0]) = __pyx_v_isimplex;

  
  __pyx_r = __pyx_v_isimplex;
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static  int __pyx_f_5scipy_7spatial_5qhull__find_simplex(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *__pyx_v_d, double *__pyx_v_c, double *__pyx_v_x, int *__pyx_v_start, double __pyx_v_eps) {
  int __pyx_v_isimplex;
  int __pyx_v_k;
  int __pyx_v_ineigh;
  int __pyx_v_ndim;
  double __pyx_v_z[(NPY_MAXDIMS + 1)];
  double __pyx_v_best_dist;
  double __pyx_v_dist;
  int __pyx_v_changed;
  int __pyx_r;
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  int __pyx_t_4;
  long __pyx_t_5;

  
  __pyx_t_1 = __pyx_f_5scipy_7spatial_5qhull__is_point_fully_outside(__pyx_v_d, __pyx_v_x, __pyx_v_eps);
  if (__pyx_t_1) {

    
    __pyx_r = -1;
    goto __pyx_L0;
    goto __pyx_L3;
  }
  __pyx_L3:;

  
  __pyx_t_2 = (__pyx_v_d->nsimplex <= 0);
  if (__pyx_t_2) {

    
    __pyx_r = -1;
    goto __pyx_L0;
    goto __pyx_L4;
  }
  __pyx_L4:;

  
  __pyx_v_ndim = __pyx_v_d->ndim;

  
  __pyx_v_isimplex = (__pyx_v_start[0]);

  
  __pyx_t_2 = (__pyx_v_isimplex < 0);
  if (!__pyx_t_2) {
    __pyx_t_3 = (__pyx_v_isimplex >= __pyx_v_d->nsimplex);
    __pyx_t_4 = __pyx_t_3;
  } else {
    __pyx_t_4 = __pyx_t_2;
  }
  if (__pyx_t_4) {

    
    __pyx_v_isimplex = 0;
    goto __pyx_L5;
  }
  __pyx_L5:;

  
  __pyx_f_5scipy_7spatial_5qhull__lift_point(__pyx_v_d, __pyx_v_x, __pyx_v_z);

  
  __pyx_v_best_dist = __pyx_f_5scipy_7spatial_5qhull__distplane(__pyx_v_d, __pyx_v_isimplex, __pyx_v_z);

  
  __pyx_v_changed = 1;

  
  while (1) {
    if (!__pyx_v_changed) break;

    
    __pyx_t_4 = (__pyx_v_best_dist > 0.0);
    if (__pyx_t_4) {

      
      goto __pyx_L7_break;
      goto __pyx_L8;
    }
    __pyx_L8:;

    
    __pyx_v_changed = 0;

    
    __pyx_t_5 = (__pyx_v_ndim + 1);
    for (__pyx_t_1 = 0; __pyx_t_1 < __pyx_t_5; __pyx_t_1+=1) {
      __pyx_v_k = __pyx_t_1;

      
      __pyx_v_ineigh = (__pyx_v_d->neighbors[(((__pyx_v_ndim + 1) * __pyx_v_isimplex) + __pyx_v_k)]);

      
      __pyx_t_4 = (__pyx_v_ineigh == -1);
      if (__pyx_t_4) {

        
        goto __pyx_L9_continue;
        goto __pyx_L11;
      }
      __pyx_L11:;

      
      __pyx_v_dist = __pyx_f_5scipy_7spatial_5qhull__distplane(__pyx_v_d, __pyx_v_ineigh, __pyx_v_z);

      
      __pyx_t_4 = (__pyx_v_dist > (__pyx_v_best_dist + (__pyx_v_eps * (1.0 + fabs(__pyx_v_best_dist)))));
      if (__pyx_t_4) {

        
        __pyx_v_isimplex = __pyx_v_ineigh;

        
        __pyx_v_best_dist = __pyx_v_dist;

        
        __pyx_v_changed = 1;
        goto __pyx_L12;
      }
      __pyx_L12:;
      __pyx_L9_continue:;
    }
  }
  __pyx_L7_break:;

  
  (__pyx_v_start[0]) = __pyx_v_isimplex;

  
  __pyx_r = __pyx_f_5scipy_7spatial_5qhull__find_simplex_directed(__pyx_v_d, __pyx_v_c, __pyx_v_x, __pyx_v_start, __pyx_v_eps);
  goto __pyx_L0;

  __pyx_r = 0;
  __pyx_L0:;
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay___init__(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); 
static PyMethodDef __pyx_mdef_5scipy_7spatial_5qhull_8Delaunay___init__ = {__Pyx_NAMESTR("__init__"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull_8Delaunay___init__, METH_VARARGS|METH_KEYWORDS, __Pyx_DOCSTR(0)};
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay___init__(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyObject *__pyx_v_self = 0;
  PyObject *__pyx_v_points = 0;
  PyObject *__pyx_v_vertices;
  PyObject *__pyx_v_neighbors;
  PyObject *__pyx_v_equations;
  PyObject *__pyx_v_paraboloid_scale;
  PyObject *__pyx_v_paraboloid_shift;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  PyObject *__pyx_t_6 = NULL;
  PyObject *__pyx_t_7 = NULL;
  static PyObject **__pyx_pyargnames[] = {&__pyx_n_s__self,&__pyx_n_s__points,0};
  __Pyx_RefNannySetupContext("__init__");
  __pyx_self = __pyx_self;
  if (unlikely(__pyx_kwds)) {
    Py_ssize_t kw_args = PyDict_Size(__pyx_kwds);
    PyObject* values[2] = {0,0};
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  2: values[1] = PyTuple_GET_ITEM(__pyx_args, 1);
      case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
      case  0: break;
      default: goto __pyx_L5_argtuple_error;
    }
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  0:
      values[0] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__self);
      if (likely(values[0])) kw_args--;
      else goto __pyx_L5_argtuple_error;
      case  1:
      values[1] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__points);
      if (likely(values[1])) kw_args--;
      else {
        __Pyx_RaiseArgtupleInvalid("__init__", 1, 2, 2, 1); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 914; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
      }
    }
    if (unlikely(kw_args > 0)) {
      if (unlikely(__Pyx_ParseOptionalKeywords(__pyx_kwds, __pyx_pyargnames, 0, values, PyTuple_GET_SIZE(__pyx_args), "__init__") < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 914; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    }
    __pyx_v_self = values[0];
    __pyx_v_points = values[1];
  } else if (PyTuple_GET_SIZE(__pyx_args) != 2) {
    goto __pyx_L5_argtuple_error;
  } else {
    __pyx_v_self = PyTuple_GET_ITEM(__pyx_args, 0);
    __pyx_v_points = PyTuple_GET_ITEM(__pyx_args, 1);
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L5_argtuple_error:;
  __Pyx_RaiseArgtupleInvalid("__init__", 1, 2, 2, PyTuple_GET_SIZE(__pyx_args)); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 914; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  __pyx_L3_error:;
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.__init__");
  __Pyx_RefNannyFinishContext();
  return NULL;
  __pyx_L4_argument_unpacking_done:;
  __Pyx_INCREF(__pyx_v_points);
  __pyx_v_vertices = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_neighbors = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_equations = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_paraboloid_scale = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_paraboloid_shift = Py_None; __Pyx_INCREF(Py_None);

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__ascontiguousarray); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_INCREF(__pyx_v_points);
  PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_v_points);
  __Pyx_GIVEREF(__pyx_v_points);
  __pyx_t_3 = PyObject_Call(__pyx_t_2, __pyx_t_1, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyObject_GetAttr(__pyx_t_3, __pyx_n_s__astype); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_3, __pyx_n_s__double); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_Call(__pyx_t_1, __pyx_t_3, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 915; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(__pyx_v_points);
  __pyx_v_points = __pyx_t_2;
  __pyx_t_2 = 0;

  
  __pyx_t_2 = __Pyx_GetName(__pyx_m, __pyx_n_s___construct_delaunay); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 917; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 917; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_INCREF(__pyx_v_points);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_v_points);
  __Pyx_GIVEREF(__pyx_v_points);
  __pyx_t_1 = PyObject_Call(__pyx_t_2, __pyx_t_3, NULL); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 917; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyTuple_CheckExact(__pyx_t_1) && likely(PyTuple_GET_SIZE(__pyx_t_1) == 5)) {
    PyObject* tuple = __pyx_t_1;
    __pyx_t_3 = PyTuple_GET_ITEM(tuple, 0); __Pyx_INCREF(__pyx_t_3);
    __pyx_t_2 = PyTuple_GET_ITEM(tuple, 1); __Pyx_INCREF(__pyx_t_2);
    __pyx_t_4 = PyTuple_GET_ITEM(tuple, 2); __Pyx_INCREF(__pyx_t_4);
    __pyx_t_5 = PyTuple_GET_ITEM(tuple, 3); __Pyx_INCREF(__pyx_t_5);
    __pyx_t_6 = PyTuple_GET_ITEM(tuple, 4); __Pyx_INCREF(__pyx_t_6);

    
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    __Pyx_DECREF(__pyx_v_vertices);
    __pyx_v_vertices = __pyx_t_3;
    __pyx_t_3 = 0;
    __Pyx_DECREF(__pyx_v_neighbors);
    __pyx_v_neighbors = __pyx_t_2;
    __pyx_t_2 = 0;
    __Pyx_DECREF(__pyx_v_equations);
    __pyx_v_equations = __pyx_t_4;
    __pyx_t_4 = 0;
    __Pyx_DECREF(__pyx_v_paraboloid_scale);
    __pyx_v_paraboloid_scale = __pyx_t_5;
    __pyx_t_5 = 0;
    __Pyx_DECREF(__pyx_v_paraboloid_shift);
    __pyx_v_paraboloid_shift = __pyx_t_6;
    __pyx_t_6 = 0;
  } else {
    __pyx_t_7 = PyObject_GetIter(__pyx_t_1); if (unlikely(!__pyx_t_7)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 916; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_7);
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    __pyx_t_3 = __Pyx_UnpackItem(__pyx_t_7, 0); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 916; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __pyx_t_2 = __Pyx_UnpackItem(__pyx_t_7, 1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 916; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_2);
    __pyx_t_4 = __Pyx_UnpackItem(__pyx_t_7, 2); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 916; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __pyx_t_5 = __Pyx_UnpackItem(__pyx_t_7, 3); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 916; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_5);
    __pyx_t_6 = __Pyx_UnpackItem(__pyx_t_7, 4); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 916; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_6);
    if (__Pyx_EndUnpack(__pyx_t_7, 5) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 916; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(__pyx_t_7); __pyx_t_7 = 0;
    __Pyx_DECREF(__pyx_v_vertices);
    __pyx_v_vertices = __pyx_t_3;
    __pyx_t_3 = 0;
    __Pyx_DECREF(__pyx_v_neighbors);
    __pyx_v_neighbors = __pyx_t_2;
    __pyx_t_2 = 0;
    __Pyx_DECREF(__pyx_v_equations);
    __pyx_v_equations = __pyx_t_4;
    __pyx_t_4 = 0;
    __Pyx_DECREF(__pyx_v_paraboloid_scale);
    __pyx_v_paraboloid_scale = __pyx_t_5;
    __pyx_t_5 = 0;
    __Pyx_DECREF(__pyx_v_paraboloid_shift);
    __pyx_v_paraboloid_shift = __pyx_t_6;
    __pyx_t_6 = 0;
  }

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_points, __pyx_n_s__shape); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 919; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_6 = __Pyx_GetItemInt(__pyx_t_1, 1, sizeof(long), PyInt_FromLong); if (!__pyx_t_6) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 919; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__ndim, __pyx_t_6) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 919; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;

  
  __pyx_t_6 = PyObject_GetAttr(__pyx_v_points, __pyx_n_s__shape); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 920; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __pyx_t_1 = __Pyx_GetItemInt(__pyx_t_6, 0, sizeof(long), PyInt_FromLong); if (!__pyx_t_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 920; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__npoints, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 920; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_vertices, __pyx_n_s__shape); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 921; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_6 = __Pyx_GetItemInt(__pyx_t_1, 0, sizeof(long), PyInt_FromLong); if (!__pyx_t_6) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 921; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__nsimplex, __pyx_t_6) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 921; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;

  
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__points, __pyx_v_points) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 922; __pyx_clineno = __LINE__; goto __pyx_L1_error;}

  
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__vertices, __pyx_v_vertices) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 923; __pyx_clineno = __LINE__; goto __pyx_L1_error;}

  
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__neighbors, __pyx_v_neighbors) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 924; __pyx_clineno = __LINE__; goto __pyx_L1_error;}

  
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__equations, __pyx_v_equations) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 925; __pyx_clineno = __LINE__; goto __pyx_L1_error;}

  
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__paraboloid_scale, __pyx_v_paraboloid_scale) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 926; __pyx_clineno = __LINE__; goto __pyx_L1_error;}

  
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__paraboloid_shift, __pyx_v_paraboloid_shift) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 927; __pyx_clineno = __LINE__; goto __pyx_L1_error;}

  
  __pyx_t_6 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__points); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 928; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __pyx_t_1 = PyObject_GetAttr(__pyx_t_6, __pyx_n_s__min); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 928; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;
  __pyx_t_6 = PyDict_New(); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 928; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_6));
  if (PyDict_SetItem(__pyx_t_6, ((PyObject *)__pyx_n_s__axis), __pyx_int_0) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 928; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_5 = PyEval_CallObjectWithKeywords(__pyx_t_1, ((PyObject *)__pyx_empty_tuple), ((PyObject *)__pyx_t_6)); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 928; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_6)); __pyx_t_6 = 0;
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__min_bound, __pyx_t_5) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 928; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;

  
  __pyx_t_5 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__points); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 929; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_6 = PyObject_GetAttr(__pyx_t_5, __pyx_n_s__max); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 929; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_6);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __pyx_t_5 = PyDict_New(); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 929; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_5));
  if (PyDict_SetItem(__pyx_t_5, ((PyObject *)__pyx_n_s__axis), __pyx_int_0) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 929; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_1 = PyEval_CallObjectWithKeywords(__pyx_t_6, ((PyObject *)__pyx_empty_tuple), ((PyObject *)__pyx_t_5)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 929; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_5)); __pyx_t_5 = 0;
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s__max_bound, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 929; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s___transform, Py_None) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 930; __pyx_clineno = __LINE__; goto __pyx_L1_error;}

  
  if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s___vertex_to_simplex, Py_None) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 931; __pyx_clineno = __LINE__; goto __pyx_L1_error;}

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_XDECREF(__pyx_t_6);
  __Pyx_XDECREF(__pyx_t_7);
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.__init__");
  __pyx_r = NULL;
  __pyx_L0:;
  __Pyx_DECREF(__pyx_v_vertices);
  __Pyx_DECREF(__pyx_v_neighbors);
  __Pyx_DECREF(__pyx_v_equations);
  __Pyx_DECREF(__pyx_v_paraboloid_scale);
  __Pyx_DECREF(__pyx_v_paraboloid_shift);
  __Pyx_DECREF(__pyx_v_points);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_transform(PyObject *__pyx_self, PyObject *__pyx_v_self); 
static char __pyx_doc_5scipy_7spatial_5qhull_8Delaunay_transform[] = "\n        Affine transform from ``x`` to the barycentric coordinates ``c``.\n\n        :type: ndarray of double, shape (nsimplex, ndim+1, ndim)\n\n        This is defined by::\n\n            T c = x - r\n\n        At vertex ``j``, ``c_j = 1`` and the other coordinates zero.\n\n        For simplex ``i``, ``transform[i,:ndim,:ndim]`` contains\n        inverse of the matrix ``T``, and ``transform[i,ndim,:]``\n        contains the vector ``r``.\n\n        ";
static PyMethodDef __pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_transform = {__Pyx_NAMESTR("transform"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_transform, METH_O, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull_8Delaunay_transform)};
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_transform(PyObject *__pyx_self, PyObject *__pyx_v_self) {
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  int __pyx_t_2;
  PyObject *__pyx_t_3 = NULL;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  __Pyx_RefNannySetupContext("transform");
  __pyx_self = __pyx_self;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s___transform); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 951; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = (__pyx_t_1 == Py_None);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (__pyx_t_2) {

    
    __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s_9); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 952; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __pyx_t_3 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__points); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 952; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);

    
    __pyx_t_4 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__vertices); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 953; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __pyx_t_5 = PyTuple_New(2); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 952; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_5);
    PyTuple_SET_ITEM(__pyx_t_5, 0, __pyx_t_3);
    __Pyx_GIVEREF(__pyx_t_3);
    PyTuple_SET_ITEM(__pyx_t_5, 1, __pyx_t_4);
    __Pyx_GIVEREF(__pyx_t_4);
    __pyx_t_3 = 0;
    __pyx_t_4 = 0;
    __pyx_t_4 = PyObject_Call(__pyx_t_1, __pyx_t_5, NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 952; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;

    
    if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s___transform, __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 952; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
    goto __pyx_L5;
  }
  __pyx_L5:;

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_4 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s___transform); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 954; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_r = __pyx_t_4;
  __pyx_t_4 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.transform");
  __pyx_r = NULL;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_vertex_to_simplex(PyObject *__pyx_self, PyObject *__pyx_v_self); 
static char __pyx_doc_5scipy_7spatial_5qhull_8Delaunay_vertex_to_simplex[] = "\n        Lookup array, from a vertex, to some simplex which it is a part of.\n\n        :type: ndarray of int, shape (npoints,)\n        ";
static PyMethodDef __pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_vertex_to_simplex = {__Pyx_NAMESTR("vertex_to_simplex"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_vertex_to_simplex, METH_O, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull_8Delaunay_vertex_to_simplex)};
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_vertex_to_simplex(PyObject *__pyx_self, PyObject *__pyx_v_self) {
  int __pyx_v_isimplex;
  int __pyx_v_k;
  int __pyx_v_ivertex;
  int __pyx_v_nsimplex;
  int __pyx_v_ndim;
  PyArrayObject *__pyx_v_vertices;
  PyArrayObject *__pyx_v_arr;
  Py_buffer __pyx_bstruct_arr;
  Py_ssize_t __pyx_bstride_0_arr = 0;
  Py_ssize_t __pyx_bshape_0_arr = 0;
  Py_buffer __pyx_bstruct_vertices;
  Py_ssize_t __pyx_bstride_0_vertices = 0;
  Py_ssize_t __pyx_bstride_1_vertices = 0;
  Py_ssize_t __pyx_bshape_0_vertices = 0;
  Py_ssize_t __pyx_bshape_1_vertices = 0;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  int __pyx_t_2;
  PyObject *__pyx_t_3 = NULL;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  PyObject *__pyx_t_6 = NULL;
  PyArrayObject *__pyx_t_7 = NULL;
  int __pyx_t_8;
  PyObject *__pyx_t_9 = NULL;
  PyObject *__pyx_t_10 = NULL;
  PyObject *__pyx_t_11 = NULL;
  PyArrayObject *__pyx_t_12 = NULL;
  int __pyx_t_13;
  long __pyx_t_14;
  int __pyx_t_15;
  int __pyx_t_16;
  int __pyx_t_17;
  int __pyx_t_18;
  int __pyx_t_19;
  int __pyx_t_20;
  __Pyx_RefNannySetupContext("vertex_to_simplex");
  __pyx_self = __pyx_self;
  __pyx_v_vertices = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_arr = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_bstruct_vertices.buf = NULL;
  __pyx_bstruct_arr.buf = NULL;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s___vertex_to_simplex); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 967; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = (__pyx_t_1 == Py_None);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (__pyx_t_2) {

    
    __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __pyx_t_3 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__empty); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__npoints); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    PyTuple_SET_ITEM(__pyx_t_4, 0, __pyx_t_1);
    __Pyx_GIVEREF(__pyx_t_1);
    __pyx_t_1 = 0;
    __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_t_4);
    __Pyx_GIVEREF(__pyx_t_4);
    __pyx_t_4 = 0;
    __pyx_t_4 = PyDict_New(); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(((PyObject *)__pyx_t_4));
    __pyx_t_5 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_5);
    __pyx_t_6 = PyObject_GetAttr(__pyx_t_5, __pyx_n_s__intc); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_6);
    __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
    if (PyDict_SetItem(__pyx_t_4, ((PyObject *)__pyx_n_s__dtype), __pyx_t_6) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;
    __pyx_t_6 = PyEval_CallObjectWithKeywords(__pyx_t_3, __pyx_t_1, ((PyObject *)__pyx_t_4)); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_6);
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    __Pyx_DECREF(((PyObject *)__pyx_t_4)); __pyx_t_4 = 0;
    if (PyObject_SetAttr(__pyx_v_self, __pyx_n_s___vertex_to_simplex, __pyx_t_6) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 968; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;

    
    __pyx_t_6 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s___vertex_to_simplex); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 969; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_6);
    __pyx_t_4 = PyObject_GetAttr(__pyx_t_6, __pyx_n_s__fill); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 969; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;
    __pyx_t_6 = PyTuple_New(1); if (unlikely(!__pyx_t_6)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 969; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_6);
    __Pyx_INCREF(__pyx_int_neg_1);
    PyTuple_SET_ITEM(__pyx_t_6, 0, __pyx_int_neg_1);
    __Pyx_GIVEREF(__pyx_int_neg_1);
    __pyx_t_1 = PyObject_Call(__pyx_t_4, __pyx_t_6, NULL); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 969; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
    __Pyx_DECREF(__pyx_t_6); __pyx_t_6 = 0;
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

    
    __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s___vertex_to_simplex); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 971; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 971; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __pyx_t_7 = ((PyArrayObject *)__pyx_t_1);
    {
      __Pyx_BufFmt_StackElem __pyx_stack[1];
      __Pyx_SafeReleaseBuffer(&__pyx_bstruct_arr);
      __pyx_t_8 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_t_7, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 1, 0, __pyx_stack);
      if (unlikely(__pyx_t_8 < 0)) {
        PyErr_Fetch(&__pyx_t_9, &__pyx_t_10, &__pyx_t_11);
        if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_v_arr, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 1, 0, __pyx_stack) == -1)) {
          Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_10); Py_XDECREF(__pyx_t_11);
          __Pyx_RaiseBufferFallbackError();
        } else {
          PyErr_Restore(__pyx_t_9, __pyx_t_10, __pyx_t_11);
        }
      }
      __pyx_bstride_0_arr = __pyx_bstruct_arr.strides[0];
      __pyx_bshape_0_arr = __pyx_bstruct_arr.shape[0];
      if (unlikely(__pyx_t_8 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 971; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    }
    __pyx_t_7 = 0;
    __Pyx_DECREF(((PyObject *)__pyx_v_arr));
    __pyx_v_arr = ((PyArrayObject *)__pyx_t_1);
    __pyx_t_1 = 0;

    
    __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__vertices); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 972; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 972; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __pyx_t_12 = ((PyArrayObject *)__pyx_t_1);
    {
      __Pyx_BufFmt_StackElem __pyx_stack[1];
      __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
      __pyx_t_8 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_vertices, (PyObject*)__pyx_t_12, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack);
      if (unlikely(__pyx_t_8 < 0)) {
        PyErr_Fetch(&__pyx_t_11, &__pyx_t_10, &__pyx_t_9);
        if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_vertices, (PyObject*)__pyx_v_vertices, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
          Py_XDECREF(__pyx_t_11); Py_XDECREF(__pyx_t_10); Py_XDECREF(__pyx_t_9);
          __Pyx_RaiseBufferFallbackError();
        } else {
          PyErr_Restore(__pyx_t_11, __pyx_t_10, __pyx_t_9);
        }
      }
      __pyx_bstride_0_vertices = __pyx_bstruct_vertices.strides[0]; __pyx_bstride_1_vertices = __pyx_bstruct_vertices.strides[1];
      __pyx_bshape_0_vertices = __pyx_bstruct_vertices.shape[0]; __pyx_bshape_1_vertices = __pyx_bstruct_vertices.shape[1];
      if (unlikely(__pyx_t_8 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 972; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    }
    __pyx_t_12 = 0;
    __Pyx_DECREF(((PyObject *)__pyx_v_vertices));
    __pyx_v_vertices = ((PyArrayObject *)__pyx_t_1);
    __pyx_t_1 = 0;

    
    __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__nsimplex); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 974; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __pyx_t_8 = __Pyx_PyInt_AsInt(__pyx_t_1); if (unlikely((__pyx_t_8 == (int)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 974; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    __pyx_v_nsimplex = __pyx_t_8;

    
    __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__ndim); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 975; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __pyx_t_8 = __Pyx_PyInt_AsInt(__pyx_t_1); if (unlikely((__pyx_t_8 == (int)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 975; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    __pyx_v_ndim = __pyx_t_8;

    
    __pyx_t_8 = __pyx_v_nsimplex;
    for (__pyx_t_13 = 0; __pyx_t_13 < __pyx_t_8; __pyx_t_13+=1) {
      __pyx_v_isimplex = __pyx_t_13;

      
      __pyx_t_14 = (__pyx_v_ndim + 1);
      for (__pyx_t_15 = 0; __pyx_t_15 < __pyx_t_14; __pyx_t_15+=1) {
        __pyx_v_k = __pyx_t_15;

        
        __pyx_t_16 = __pyx_v_isimplex;
        __pyx_t_17 = __pyx_v_k;
        __pyx_t_18 = -1;
        if (__pyx_t_16 < 0) {
          __pyx_t_16 += __pyx_bshape_0_vertices;
          if (unlikely(__pyx_t_16 < 0)) __pyx_t_18 = 0;
        } else if (unlikely(__pyx_t_16 >= __pyx_bshape_0_vertices)) __pyx_t_18 = 0;
        if (__pyx_t_17 < 0) {
          __pyx_t_17 += __pyx_bshape_1_vertices;
          if (unlikely(__pyx_t_17 < 0)) __pyx_t_18 = 1;
        } else if (unlikely(__pyx_t_17 >= __pyx_bshape_1_vertices)) __pyx_t_18 = 1;
        if (unlikely(__pyx_t_18 != -1)) {
          __Pyx_RaiseBufferIndexError(__pyx_t_18);
          {__pyx_filename = __pyx_f[0]; __pyx_lineno = 979; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        }
        __pyx_v_ivertex = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_16, __pyx_bstride_0_vertices, __pyx_t_17, __pyx_bstride_1_vertices));

        
        __pyx_t_18 = __pyx_v_ivertex;
        __pyx_t_19 = -1;
        if (__pyx_t_18 < 0) {
          __pyx_t_18 += __pyx_bshape_0_arr;
          if (unlikely(__pyx_t_18 < 0)) __pyx_t_19 = 0;
        } else if (unlikely(__pyx_t_18 >= __pyx_bshape_0_arr)) __pyx_t_19 = 0;
        if (unlikely(__pyx_t_19 != -1)) {
          __Pyx_RaiseBufferIndexError(__pyx_t_19);
          {__pyx_filename = __pyx_f[0]; __pyx_lineno = 980; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        }
        __pyx_t_2 = ((*__Pyx_BufPtrStrided1d(npy_int *, __pyx_bstruct_arr.buf, __pyx_t_18, __pyx_bstride_0_arr)) == -1);
        if (__pyx_t_2) {

          
          __pyx_t_19 = __pyx_v_ivertex;
          __pyx_t_20 = -1;
          if (__pyx_t_19 < 0) {
            __pyx_t_19 += __pyx_bshape_0_arr;
            if (unlikely(__pyx_t_19 < 0)) __pyx_t_20 = 0;
          } else if (unlikely(__pyx_t_19 >= __pyx_bshape_0_arr)) __pyx_t_20 = 0;
          if (unlikely(__pyx_t_20 != -1)) {
            __Pyx_RaiseBufferIndexError(__pyx_t_20);
            {__pyx_filename = __pyx_f[0]; __pyx_lineno = 981; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
          }
          *__Pyx_BufPtrStrided1d(npy_int *, __pyx_bstruct_arr.buf, __pyx_t_19, __pyx_bstride_0_arr) = __pyx_v_isimplex;
          goto __pyx_L10;
        }
        __pyx_L10:;
      }
    }
    goto __pyx_L5;
  }
  __pyx_L5:;

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s___vertex_to_simplex); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 983; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_XDECREF(__pyx_t_6);
  { PyObject *__pyx_type, *__pyx_value, *__pyx_tb;
    __Pyx_ErrFetch(&__pyx_type, &__pyx_value, &__pyx_tb);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_arr);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
  __Pyx_ErrRestore(__pyx_type, __pyx_value, __pyx_tb);}
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.vertex_to_simplex");
  __pyx_r = NULL;
  goto __pyx_L2;
  __pyx_L0:;
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_arr);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
  __pyx_L2:;
  __Pyx_DECREF((PyObject *)__pyx_v_vertices);
  __Pyx_DECREF((PyObject *)__pyx_v_arr);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_convex_hull(PyObject *__pyx_self, PyObject *__pyx_v_self); 
static char __pyx_doc_5scipy_7spatial_5qhull_8Delaunay_convex_hull[] = "\n        Vertices of facets forming the convex hull of the point set.\n\n        :type: ndarray of int, shape (nfaces, ndim)\n\n        The array contains the indices of the points\n        belonging to the (N-1)-dimensional facets that form the convex\n        hull of the triangulation.\n\n        ";
static PyMethodDef __pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_convex_hull = {__Pyx_NAMESTR("convex_hull"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_convex_hull, METH_O, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull_8Delaunay_convex_hull)};
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_convex_hull(PyObject *__pyx_self, PyObject *__pyx_v_self) {
  int __pyx_v_isimplex;
  int __pyx_v_k;
  int __pyx_v_j;
  int __pyx_v_ndim;
  int __pyx_v_nsimplex;
  int __pyx_v_m;
  int __pyx_v_msize;
  PyArrayObject *__pyx_v_arr;
  PyArrayObject *__pyx_v_neighbors;
  PyArrayObject *__pyx_v_vertices;
  PyObject *__pyx_v_out;
  Py_buffer __pyx_bstruct_neighbors;
  Py_ssize_t __pyx_bstride_0_neighbors = 0;
  Py_ssize_t __pyx_bstride_1_neighbors = 0;
  Py_ssize_t __pyx_bshape_0_neighbors = 0;
  Py_ssize_t __pyx_bshape_1_neighbors = 0;
  Py_buffer __pyx_bstruct_arr;
  Py_ssize_t __pyx_bstride_0_arr = 0;
  Py_ssize_t __pyx_bstride_1_arr = 0;
  Py_ssize_t __pyx_bshape_0_arr = 0;
  Py_ssize_t __pyx_bshape_1_arr = 0;
  Py_buffer __pyx_bstruct_vertices;
  Py_ssize_t __pyx_bstride_0_vertices = 0;
  Py_ssize_t __pyx_bstride_1_vertices = 0;
  Py_ssize_t __pyx_bshape_0_vertices = 0;
  Py_ssize_t __pyx_bshape_1_vertices = 0;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  PyArrayObject *__pyx_t_2 = NULL;
  int __pyx_t_3;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  PyObject *__pyx_t_6 = NULL;
  PyArrayObject *__pyx_t_7 = NULL;
  PyObject *__pyx_t_8 = NULL;
  PyObject *__pyx_t_9 = NULL;
  PyObject *__pyx_t_10 = NULL;
  PyObject *__pyx_t_11 = NULL;
  PyArrayObject *__pyx_t_12 = NULL;
  int __pyx_t_13;
  long __pyx_t_14;
  int __pyx_t_15;
  int __pyx_t_16;
  int __pyx_t_17;
  int __pyx_t_18;
  long __pyx_t_19;
  int __pyx_t_20;
  int __pyx_t_21;
  int __pyx_t_22;
  int __pyx_t_23;
  int __pyx_t_24;
  int __pyx_t_25;
  int __pyx_t_26;
  int __pyx_t_27;
  long __pyx_t_28;
  __Pyx_RefNannySetupContext("convex_hull");
  __pyx_self = __pyx_self;
  __pyx_v_arr = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_neighbors = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_vertices = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_out = Py_None; __Pyx_INCREF(Py_None);
  __pyx_bstruct_arr.buf = NULL;
  __pyx_bstruct_neighbors.buf = NULL;
  __pyx_bstruct_vertices.buf = NULL;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__neighbors); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1003; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1003; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_2 = ((PyArrayObject *)__pyx_t_1);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_neighbors);
    __pyx_t_3 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_neighbors, (PyObject*)__pyx_t_2, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_3 < 0)) {
      PyErr_Fetch(&__pyx_t_4, &__pyx_t_5, &__pyx_t_6);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_neighbors, (PyObject*)__pyx_v_neighbors, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_4); Py_XDECREF(__pyx_t_5); Py_XDECREF(__pyx_t_6);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_4, __pyx_t_5, __pyx_t_6);
      }
    }
    __pyx_bstride_0_neighbors = __pyx_bstruct_neighbors.strides[0]; __pyx_bstride_1_neighbors = __pyx_bstruct_neighbors.strides[1];
    __pyx_bshape_0_neighbors = __pyx_bstruct_neighbors.shape[0]; __pyx_bshape_1_neighbors = __pyx_bstruct_neighbors.shape[1];
    if (unlikely(__pyx_t_3 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1003; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_2 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_neighbors));
  __pyx_v_neighbors = ((PyArrayObject *)__pyx_t_1);
  __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__vertices); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1004; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1004; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_7 = ((PyArrayObject *)__pyx_t_1);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
    __pyx_t_3 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_vertices, (PyObject*)__pyx_t_7, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_3 < 0)) {
      PyErr_Fetch(&__pyx_t_6, &__pyx_t_5, &__pyx_t_4);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_vertices, (PyObject*)__pyx_v_vertices, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_6); Py_XDECREF(__pyx_t_5); Py_XDECREF(__pyx_t_4);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_6, __pyx_t_5, __pyx_t_4);
      }
    }
    __pyx_bstride_0_vertices = __pyx_bstruct_vertices.strides[0]; __pyx_bstride_1_vertices = __pyx_bstruct_vertices.strides[1];
    __pyx_bshape_0_vertices = __pyx_bstruct_vertices.shape[0]; __pyx_bshape_1_vertices = __pyx_bstruct_vertices.shape[1];
    if (unlikely(__pyx_t_3 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1004; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_7 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_vertices));
  __pyx_v_vertices = ((PyArrayObject *)__pyx_t_1);
  __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__ndim); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1005; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_3 = __Pyx_PyInt_AsInt(__pyx_t_1); if (unlikely((__pyx_t_3 == (int)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1005; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_v_ndim = __pyx_t_3;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__nsimplex); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1006; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_3 = __Pyx_PyInt_AsInt(__pyx_t_1); if (unlikely((__pyx_t_3 == (int)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1006; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_v_nsimplex = __pyx_t_3;

  
  __pyx_v_msize = 10;

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_8 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__empty); if (unlikely(!__pyx_t_8)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_8);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyInt_FromLong(__pyx_v_msize); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_9 = PyInt_FromLong(__pyx_v_ndim); if (unlikely(!__pyx_t_9)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_9);
  __pyx_t_10 = PyTuple_New(2); if (unlikely(!__pyx_t_10)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_10);
  PyTuple_SET_ITEM(__pyx_t_10, 0, __pyx_t_1);
  __Pyx_GIVEREF(__pyx_t_1);
  PyTuple_SET_ITEM(__pyx_t_10, 1, __pyx_t_9);
  __Pyx_GIVEREF(__pyx_t_9);
  __pyx_t_1 = 0;
  __pyx_t_9 = 0;
  __pyx_t_9 = PyTuple_New(1); if (unlikely(!__pyx_t_9)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_9);
  PyTuple_SET_ITEM(__pyx_t_9, 0, __pyx_t_10);
  __Pyx_GIVEREF(__pyx_t_10);
  __pyx_t_10 = 0;
  __pyx_t_10 = PyDict_New(); if (unlikely(!__pyx_t_10)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_10));
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_11 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__intc); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_11);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (PyDict_SetItem(__pyx_t_10, ((PyObject *)__pyx_n_s__dtype), __pyx_t_11) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_11); __pyx_t_11 = 0;
  __pyx_t_11 = PyEval_CallObjectWithKeywords(__pyx_t_8, __pyx_t_9, ((PyObject *)__pyx_t_10)); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1009; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_11);
  __Pyx_DECREF(__pyx_t_8); __pyx_t_8 = 0;
  __Pyx_DECREF(__pyx_t_9); __pyx_t_9 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_10)); __pyx_t_10 = 0;
  __Pyx_DECREF(__pyx_v_out);
  __pyx_v_out = __pyx_t_11;
  __pyx_t_11 = 0;

  
  if (!(likely(((__pyx_v_out) == Py_None) || likely(__Pyx_TypeTest(__pyx_v_out, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1010; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_12 = ((PyArrayObject *)__pyx_v_out);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_arr);
    __pyx_t_3 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_t_12, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_3 < 0)) {
      PyErr_Fetch(&__pyx_t_4, &__pyx_t_5, &__pyx_t_6);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_v_arr, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_4); Py_XDECREF(__pyx_t_5); Py_XDECREF(__pyx_t_6);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_4, __pyx_t_5, __pyx_t_6);
      }
    }
    __pyx_bstride_0_arr = __pyx_bstruct_arr.strides[0]; __pyx_bstride_1_arr = __pyx_bstruct_arr.strides[1];
    __pyx_bshape_0_arr = __pyx_bstruct_arr.shape[0]; __pyx_bshape_1_arr = __pyx_bstruct_arr.shape[1];
    if (unlikely(__pyx_t_3 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1010; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_12 = 0;
  __Pyx_INCREF(__pyx_v_out);
  __Pyx_DECREF(((PyObject *)__pyx_v_arr));
  __pyx_v_arr = ((PyArrayObject *)__pyx_v_out);

  
  __pyx_v_m = 0;

  
  __pyx_t_3 = __pyx_v_nsimplex;
  for (__pyx_t_13 = 0; __pyx_t_13 < __pyx_t_3; __pyx_t_13+=1) {
    __pyx_v_isimplex = __pyx_t_13;

    
    __pyx_t_14 = (__pyx_v_ndim + 1);
    for (__pyx_t_15 = 0; __pyx_t_15 < __pyx_t_14; __pyx_t_15+=1) {
      __pyx_v_k = __pyx_t_15;

      
      __pyx_t_16 = __pyx_v_isimplex;
      __pyx_t_17 = __pyx_v_k;
      if (__pyx_t_16 < 0) __pyx_t_16 += __pyx_bshape_0_neighbors;
      if (__pyx_t_17 < 0) __pyx_t_17 += __pyx_bshape_1_neighbors;
      __pyx_t_18 = ((*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_neighbors.buf, __pyx_t_16, __pyx_bstride_0_neighbors, __pyx_t_17, __pyx_bstride_1_neighbors)) == -1);
      if (__pyx_t_18) {

        
        __pyx_t_19 = (__pyx_v_ndim + 1);
        for (__pyx_t_20 = 0; __pyx_t_20 < __pyx_t_19; __pyx_t_20+=1) {
          __pyx_v_j = __pyx_t_20;

          
          __pyx_t_18 = (__pyx_v_j < __pyx_v_k);
          if (__pyx_t_18) {

            
            __pyx_t_21 = __pyx_v_isimplex;
            __pyx_t_22 = __pyx_v_j;
            if (__pyx_t_21 < 0) __pyx_t_21 += __pyx_bshape_0_vertices;
            if (__pyx_t_22 < 0) __pyx_t_22 += __pyx_bshape_1_vertices;
            __pyx_t_23 = __pyx_v_m;
            __pyx_t_24 = __pyx_v_j;
            if (__pyx_t_23 < 0) __pyx_t_23 += __pyx_bshape_0_arr;
            if (__pyx_t_24 < 0) __pyx_t_24 += __pyx_bshape_1_arr;
            *__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_arr.buf, __pyx_t_23, __pyx_bstride_0_arr, __pyx_t_24, __pyx_bstride_1_arr) = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_21, __pyx_bstride_0_vertices, __pyx_t_22, __pyx_bstride_1_vertices));
            goto __pyx_L12;
          }

          
          __pyx_t_18 = (__pyx_v_j > __pyx_v_k);
          if (__pyx_t_18) {

            
            __pyx_t_25 = __pyx_v_isimplex;
            __pyx_t_26 = __pyx_v_j;
            if (__pyx_t_25 < 0) __pyx_t_25 += __pyx_bshape_0_vertices;
            if (__pyx_t_26 < 0) __pyx_t_26 += __pyx_bshape_1_vertices;
            __pyx_t_27 = __pyx_v_m;
            __pyx_t_28 = (__pyx_v_j - 1);
            if (__pyx_t_27 < 0) __pyx_t_27 += __pyx_bshape_0_arr;
            if (__pyx_t_28 < 0) __pyx_t_28 += __pyx_bshape_1_arr;
            *__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_arr.buf, __pyx_t_27, __pyx_bstride_0_arr, __pyx_t_28, __pyx_bstride_1_arr) = (*__Pyx_BufPtrStrided2d(npy_int *, __pyx_bstruct_vertices.buf, __pyx_t_25, __pyx_bstride_0_vertices, __pyx_t_26, __pyx_bstride_1_vertices));
            goto __pyx_L12;
          }
          __pyx_L12:;
        }

        
        __pyx_v_m += 1;

        
        __pyx_t_18 = (__pyx_v_m >= __pyx_v_msize);
        if (__pyx_t_18) {

          
          __pyx_t_12 = ((PyArrayObject *)Py_None);
          {
            __Pyx_BufFmt_StackElem __pyx_stack[1];
            __Pyx_SafeReleaseBuffer(&__pyx_bstruct_arr);
            __pyx_t_20 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_t_12, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack);
            if (unlikely(__pyx_t_20 < 0)) {
              PyErr_Fetch(&__pyx_t_6, &__pyx_t_5, &__pyx_t_4);
              if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_v_arr, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack) == -1)) {
                Py_XDECREF(__pyx_t_6); Py_XDECREF(__pyx_t_5); Py_XDECREF(__pyx_t_4);
                __Pyx_RaiseBufferFallbackError();
              } else {
                PyErr_Restore(__pyx_t_6, __pyx_t_5, __pyx_t_4);
              }
            }
            __pyx_bstride_0_arr = __pyx_bstruct_arr.strides[0]; __pyx_bstride_1_arr = __pyx_bstruct_arr.strides[1];
            __pyx_bshape_0_arr = __pyx_bstruct_arr.shape[0]; __pyx_bshape_1_arr = __pyx_bstruct_arr.shape[1];
            if (unlikely(__pyx_t_20 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1024; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
          }
          __pyx_t_12 = 0;
          __Pyx_INCREF(Py_None);
          __Pyx_DECREF(((PyObject *)__pyx_v_arr));
          __pyx_v_arr = ((PyArrayObject *)Py_None);

          
          __pyx_v_msize = ((2 * __pyx_v_msize) + 1);

          
          __pyx_t_11 = PyObject_GetAttr(__pyx_v_out, __pyx_n_s__resize); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1026; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
          __Pyx_GOTREF(__pyx_t_11);
          __pyx_t_10 = PyInt_FromLong(__pyx_v_msize); if (unlikely(!__pyx_t_10)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1026; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
          __Pyx_GOTREF(__pyx_t_10);
          __pyx_t_9 = PyInt_FromLong(__pyx_v_ndim); if (unlikely(!__pyx_t_9)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1026; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
          __Pyx_GOTREF(__pyx_t_9);
          __pyx_t_8 = PyTuple_New(2); if (unlikely(!__pyx_t_8)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1026; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
          __Pyx_GOTREF(__pyx_t_8);
          PyTuple_SET_ITEM(__pyx_t_8, 0, __pyx_t_10);
          __Pyx_GIVEREF(__pyx_t_10);
          PyTuple_SET_ITEM(__pyx_t_8, 1, __pyx_t_9);
          __Pyx_GIVEREF(__pyx_t_9);
          __pyx_t_10 = 0;
          __pyx_t_9 = 0;
          __pyx_t_9 = PyObject_Call(__pyx_t_11, __pyx_t_8, NULL); if (unlikely(!__pyx_t_9)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1026; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
          __Pyx_GOTREF(__pyx_t_9);
          __Pyx_DECREF(__pyx_t_11); __pyx_t_11 = 0;
          __Pyx_DECREF(__pyx_t_8); __pyx_t_8 = 0;
          __Pyx_DECREF(__pyx_t_9); __pyx_t_9 = 0;

          
          if (!(likely(((__pyx_v_out) == Py_None) || likely(__Pyx_TypeTest(__pyx_v_out, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1027; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
          __pyx_t_12 = ((PyArrayObject *)__pyx_v_out);
          {
            __Pyx_BufFmt_StackElem __pyx_stack[1];
            __Pyx_SafeReleaseBuffer(&__pyx_bstruct_arr);
            __pyx_t_20 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_t_12, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack);
            if (unlikely(__pyx_t_20 < 0)) {
              PyErr_Fetch(&__pyx_t_4, &__pyx_t_5, &__pyx_t_6);
              if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_v_arr, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack) == -1)) {
                Py_XDECREF(__pyx_t_4); Py_XDECREF(__pyx_t_5); Py_XDECREF(__pyx_t_6);
                __Pyx_RaiseBufferFallbackError();
              } else {
                PyErr_Restore(__pyx_t_4, __pyx_t_5, __pyx_t_6);
              }
            }
            __pyx_bstride_0_arr = __pyx_bstruct_arr.strides[0]; __pyx_bstride_1_arr = __pyx_bstruct_arr.strides[1];
            __pyx_bshape_0_arr = __pyx_bstruct_arr.shape[0]; __pyx_bshape_1_arr = __pyx_bstruct_arr.shape[1];
            if (unlikely(__pyx_t_20 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1027; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
          }
          __pyx_t_12 = 0;
          __Pyx_INCREF(__pyx_v_out);
          __Pyx_DECREF(((PyObject *)__pyx_v_arr));
          __pyx_v_arr = ((PyArrayObject *)__pyx_v_out);
          goto __pyx_L13;
        }
        __pyx_L13:;
        goto __pyx_L9;
      }
      __pyx_L9:;
    }
  }

  
  __pyx_t_12 = ((PyArrayObject *)Py_None);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_arr);
    __pyx_t_3 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_t_12, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_3 < 0)) {
      PyErr_Fetch(&__pyx_t_6, &__pyx_t_5, &__pyx_t_4);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_arr, (PyObject*)__pyx_v_arr, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_6); Py_XDECREF(__pyx_t_5); Py_XDECREF(__pyx_t_4);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_6, __pyx_t_5, __pyx_t_4);
      }
    }
    __pyx_bstride_0_arr = __pyx_bstruct_arr.strides[0]; __pyx_bstride_1_arr = __pyx_bstruct_arr.strides[1];
    __pyx_bshape_0_arr = __pyx_bstruct_arr.shape[0]; __pyx_bshape_1_arr = __pyx_bstruct_arr.shape[1];
    if (unlikely(__pyx_t_3 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1029; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_12 = 0;
  __Pyx_INCREF(Py_None);
  __Pyx_DECREF(((PyObject *)__pyx_v_arr));
  __pyx_v_arr = ((PyArrayObject *)Py_None);

  
  __pyx_t_9 = PyObject_GetAttr(__pyx_v_out, __pyx_n_s__resize); if (unlikely(!__pyx_t_9)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1030; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_9);
  __pyx_t_8 = PyInt_FromLong(__pyx_v_m); if (unlikely(!__pyx_t_8)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1030; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_8);
  __pyx_t_11 = PyInt_FromLong(__pyx_v_ndim); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1030; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_11);
  __pyx_t_10 = PyTuple_New(2); if (unlikely(!__pyx_t_10)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1030; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_10);
  PyTuple_SET_ITEM(__pyx_t_10, 0, __pyx_t_8);
  __Pyx_GIVEREF(__pyx_t_8);
  PyTuple_SET_ITEM(__pyx_t_10, 1, __pyx_t_11);
  __Pyx_GIVEREF(__pyx_t_11);
  __pyx_t_8 = 0;
  __pyx_t_11 = 0;
  __pyx_t_11 = PyObject_Call(__pyx_t_9, __pyx_t_10, NULL); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1030; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_11);
  __Pyx_DECREF(__pyx_t_9); __pyx_t_9 = 0;
  __Pyx_DECREF(__pyx_t_10); __pyx_t_10 = 0;
  __Pyx_DECREF(__pyx_t_11); __pyx_t_11 = 0;

  
  __Pyx_XDECREF(__pyx_r);
  __Pyx_INCREF(__pyx_v_out);
  __pyx_r = __pyx_v_out;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_8);
  __Pyx_XDECREF(__pyx_t_9);
  __Pyx_XDECREF(__pyx_t_10);
  __Pyx_XDECREF(__pyx_t_11);
  { PyObject *__pyx_type, *__pyx_value, *__pyx_tb;
    __Pyx_ErrFetch(&__pyx_type, &__pyx_value, &__pyx_tb);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_neighbors);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_arr);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
  __Pyx_ErrRestore(__pyx_type, __pyx_value, __pyx_tb);}
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.convex_hull");
  __pyx_r = NULL;
  goto __pyx_L2;
  __pyx_L0:;
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_neighbors);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_arr);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
  __pyx_L2:;
  __Pyx_DECREF((PyObject *)__pyx_v_arr);
  __Pyx_DECREF((PyObject *)__pyx_v_neighbors);
  __Pyx_DECREF((PyObject *)__pyx_v_vertices);
  __Pyx_DECREF(__pyx_v_out);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_find_simplex(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); 
static char __pyx_doc_5scipy_7spatial_5qhull_8Delaunay_find_simplex[] = "\n        find_simplex(xi, bruteforce=False)\n\n        Find the simplices containing the given points.\n\n        Parameters\n        ----------\n        tri : DelaunayInfo\n            Delaunay triangulation\n        xi : ndarray of double, shape (..., ndim)\n            Points to locate\n        bruteforce : bool, optional\n            Whether to only perform a brute-force search\n\n        Returns\n        -------\n        i : ndarray of int, same shape as `xi`\n            Indices of simplices containing each point.\n            Points outside the triangulation get the value -1.\n\n        Notes\n        -----\n        This uses an algorithm adapted from Qhull's qh_findbestfacet,\n        which makes use of the connection between a convex hull and a\n        Delaunay triangulation. After finding the simplex closest to\n        the point in N+1 dimensions, the algorithm falls back to\n        directed search in N dimensions.\n\n        ";
static PyMethodDef __pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_find_simplex = {__Pyx_NAMESTR("find_simplex"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_find_simplex, METH_VARARGS|METH_KEYWORDS, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull_8Delaunay_find_simplex)};
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_find_simplex(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyObject *__pyx_v_self = 0;
  PyObject *__pyx_v_xi = 0;
  PyObject *__pyx_v_bruteforce = 0;
  __pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t __pyx_v_info;
  int __pyx_v_isimplex;
  double __pyx_v_c[NPY_MAXDIMS];
  double __pyx_v_eps;
  int __pyx_v_start;
  int __pyx_v_k;
  PyArrayObject *__pyx_v_x;
  PyArrayObject *__pyx_v_out_;
  PyObject *__pyx_v_xi_shape;
  PyObject *__pyx_v_out;
  Py_buffer __pyx_bstruct_out_;
  Py_ssize_t __pyx_bstride_0_out_ = 0;
  Py_ssize_t __pyx_bshape_0_out_ = 0;
  Py_buffer __pyx_bstruct_x;
  Py_ssize_t __pyx_bstride_0_x = 0;
  Py_ssize_t __pyx_bstride_1_x = 0;
  Py_ssize_t __pyx_bshape_0_x = 0;
  Py_ssize_t __pyx_bshape_1_x = 0;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  int __pyx_t_4;
  PyObject *__pyx_t_5 = NULL;
  PyArrayObject *__pyx_t_6 = NULL;
  int __pyx_t_7;
  PyObject *__pyx_t_8 = NULL;
  PyObject *__pyx_t_9 = NULL;
  PyObject *__pyx_t_10 = NULL;
  double __pyx_t_11;
  PyObject *__pyx_t_12 = NULL;
  PyArrayObject *__pyx_t_13 = NULL;
  npy_intp __pyx_t_14;
  int __pyx_t_15;
  int __pyx_t_16;
  int __pyx_t_17;
  static PyObject **__pyx_pyargnames[] = {&__pyx_n_s__self,&__pyx_n_s__xi,&__pyx_n_s__bruteforce,0};
  __Pyx_RefNannySetupContext("find_simplex");
  __pyx_self = __pyx_self;
  if (unlikely(__pyx_kwds)) {
    Py_ssize_t kw_args = PyDict_Size(__pyx_kwds);
    PyObject* values[3] = {0,0,0};
    values[2] = __pyx_k_10;
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  3: values[2] = PyTuple_GET_ITEM(__pyx_args, 2);
      case  2: values[1] = PyTuple_GET_ITEM(__pyx_args, 1);
      case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
      case  0: break;
      default: goto __pyx_L5_argtuple_error;
    }
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  0:
      values[0] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__self);
      if (likely(values[0])) kw_args--;
      else goto __pyx_L5_argtuple_error;
      case  1:
      values[1] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__xi);
      if (likely(values[1])) kw_args--;
      else {
        __Pyx_RaiseArgtupleInvalid("find_simplex", 0, 2, 3, 1); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1033; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
      }
      case  2:
      if (kw_args > 0) {
        PyObject* value = PyDict_GetItem(__pyx_kwds, __pyx_n_s__bruteforce);
        if (value) { values[2] = value; kw_args--; }
      }
    }
    if (unlikely(kw_args > 0)) {
      if (unlikely(__Pyx_ParseOptionalKeywords(__pyx_kwds, __pyx_pyargnames, 0, values, PyTuple_GET_SIZE(__pyx_args), "find_simplex") < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1033; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    }
    __pyx_v_self = values[0];
    __pyx_v_xi = values[1];
    __pyx_v_bruteforce = values[2];
  } else {
    __pyx_v_bruteforce = __pyx_k_10;
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  3:
      __pyx_v_bruteforce = PyTuple_GET_ITEM(__pyx_args, 2);
      case  2:
      __pyx_v_xi = PyTuple_GET_ITEM(__pyx_args, 1);
      __pyx_v_self = PyTuple_GET_ITEM(__pyx_args, 0);
      break;
      default: goto __pyx_L5_argtuple_error;
    }
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L5_argtuple_error:;
  __Pyx_RaiseArgtupleInvalid("find_simplex", 0, 2, 3, PyTuple_GET_SIZE(__pyx_args)); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1033; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  __pyx_L3_error:;
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.find_simplex");
  __Pyx_RefNannyFinishContext();
  return NULL;
  __pyx_L4_argument_unpacking_done:;
  __Pyx_INCREF(__pyx_v_xi);
  __pyx_v_x = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_out_ = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_xi_shape = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_out = Py_None; __Pyx_INCREF(Py_None);
  __pyx_bstruct_x.buf = NULL;
  __pyx_bstruct_out_.buf = NULL;

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1072; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__asanyarray); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1072; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1072; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_INCREF(__pyx_v_xi);
  PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_v_xi);
  __Pyx_GIVEREF(__pyx_v_xi);
  __pyx_t_3 = PyObject_Call(__pyx_t_2, __pyx_t_1, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1072; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_v_xi);
  __pyx_v_xi = __pyx_t_3;
  __pyx_t_3 = 0;

  
  __pyx_t_3 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__shape); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1074; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_1 = __Pyx_GetItemInt(__pyx_t_3, -1, sizeof(long), PyInt_FromLong); if (!__pyx_t_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1074; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__ndim); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1074; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_2 = PyObject_RichCompare(__pyx_t_1, __pyx_t_3, Py_NE); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1074; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_4 = __Pyx_PyObject_IsTrue(__pyx_t_2); if (unlikely(__pyx_t_4 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1074; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  if (__pyx_t_4) {

    
    __pyx_t_2 = PyTuple_New(1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1075; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_2);
    __Pyx_INCREF(((PyObject *)__pyx_kp_s_11));
    PyTuple_SET_ITEM(__pyx_t_2, 0, ((PyObject *)__pyx_kp_s_11));
    __Pyx_GIVEREF(((PyObject *)__pyx_kp_s_11));
    __pyx_t_3 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_2, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1075; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
    __Pyx_Raise(__pyx_t_3, 0, 0);
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1075; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L6;
  }
  __pyx_L6:;

  
  __pyx_t_3 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__shape); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1077; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_v_xi_shape);
  __pyx_v_xi_shape = __pyx_t_3;
  __pyx_t_3 = 0;

  
  __pyx_t_3 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__reshape); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_2 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_1 = PyObject_GetAttr(__pyx_t_2, __pyx_n_s__prod); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__shape); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_5 = PySequence_GetSlice(__pyx_t_2, 0, -1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_t_2 = PyTuple_New(1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_5 = 0;
  __pyx_t_5 = PyObject_Call(__pyx_t_1, __pyx_t_2, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__shape); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_1 = __Pyx_GetItemInt(__pyx_t_2, -1, sizeof(long), PyInt_FromLong); if (!__pyx_t_1) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_t_2 = PyTuple_New(2); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  PyTuple_SET_ITEM(__pyx_t_2, 1, __pyx_t_1);
  __Pyx_GIVEREF(__pyx_t_1);
  __pyx_t_5 = 0;
  __pyx_t_1 = 0;
  __pyx_t_1 = PyObject_Call(__pyx_t_3, __pyx_t_2, NULL); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1078; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_v_xi);
  __pyx_v_xi = __pyx_t_1;
  __pyx_t_1 = 0;

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__ascontiguousarray); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__astype); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_3 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_5 = PyObject_GetAttr(__pyx_t_3, __pyx_n_s__double); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_5 = 0;
  __pyx_t_5 = PyObject_Call(__pyx_t_1, __pyx_t_3, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_5 = 0;
  __pyx_t_5 = PyObject_Call(__pyx_t_2, __pyx_t_3, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (!(likely(((__pyx_t_5) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_5, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_6 = ((PyArrayObject *)__pyx_t_5);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_x);
    __pyx_t_7 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_x, (PyObject*)__pyx_t_6, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_7 < 0)) {
      PyErr_Fetch(&__pyx_t_8, &__pyx_t_9, &__pyx_t_10);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_x, (PyObject*)__pyx_v_x, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_8); Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_10);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_8, __pyx_t_9, __pyx_t_10);
      }
    }
    __pyx_bstride_0_x = __pyx_bstruct_x.strides[0]; __pyx_bstride_1_x = __pyx_bstruct_x.strides[1];
    __pyx_bshape_0_x = __pyx_bstruct_x.shape[0]; __pyx_bshape_1_x = __pyx_bstruct_x.shape[1];
    if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1079; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_6 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_x));
  __pyx_v_x = ((PyArrayObject *)__pyx_t_5);
  __pyx_t_5 = 0;

  
  __pyx_v_start = 0;

  
  __pyx_t_5 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1083; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_5, __pyx_n_s__finfo); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1083; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __pyx_t_5 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1083; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_5, __pyx_n_s__double); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1083; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __pyx_t_5 = PyTuple_New(1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1083; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  PyTuple_SET_ITEM(__pyx_t_5, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_Call(__pyx_t_3, __pyx_t_5, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1083; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __pyx_t_5 = PyObject_GetAttr(__pyx_t_2, __pyx_n_s__eps); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1083; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_t_2 = PyNumber_Multiply(__pyx_t_5, __pyx_int_10); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1083; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __pyx_t_11 = __pyx_PyFloat_AsDouble(__pyx_t_2); if (unlikely((__pyx_t_11 == (double)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1083; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_v_eps = __pyx_t_11;

  
  __pyx_t_2 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_5 = PyObject_GetAttr(__pyx_t_2, __pyx_n_s__zeros); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__shape); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_3 = __Pyx_GetItemInt(__pyx_t_2, 0, sizeof(long), PyInt_FromLong); if (!__pyx_t_3) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_t_2 = PyTuple_New(1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_t_3);
  __Pyx_GIVEREF(__pyx_t_3);
  __pyx_t_3 = 0;
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyDict_New(); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_2));
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_12 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__intc); if (unlikely(!__pyx_t_12)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_12);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (PyDict_SetItem(__pyx_t_2, ((PyObject *)__pyx_n_s__dtype), __pyx_t_12) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_12); __pyx_t_12 = 0;
  __pyx_t_12 = PyEval_CallObjectWithKeywords(__pyx_t_5, __pyx_t_3, ((PyObject *)__pyx_t_2)); if (unlikely(!__pyx_t_12)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1084; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_12);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_2)); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_v_out);
  __pyx_v_out = __pyx_t_12;
  __pyx_t_12 = 0;

  
  if (!(likely(((__pyx_v_out) == Py_None) || likely(__Pyx_TypeTest(__pyx_v_out, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1085; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_13 = ((PyArrayObject *)__pyx_v_out);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_out_);
    __pyx_t_7 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_out_, (PyObject*)__pyx_t_13, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 1, 0, __pyx_stack);
    if (unlikely(__pyx_t_7 < 0)) {
      PyErr_Fetch(&__pyx_t_10, &__pyx_t_9, &__pyx_t_8);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_out_, (PyObject*)__pyx_v_out_, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 1, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_10); Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_8);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_10, __pyx_t_9, __pyx_t_8);
      }
    }
    __pyx_bstride_0_out_ = __pyx_bstruct_out_.strides[0];
    __pyx_bshape_0_out_ = __pyx_bstruct_out_.shape[0];
    if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1085; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_13 = 0;
  __Pyx_INCREF(__pyx_v_out);
  __Pyx_DECREF(((PyObject *)__pyx_v_out_));
  __pyx_v_out_ = ((PyArrayObject *)__pyx_v_out);

  
  __pyx_f_5scipy_7spatial_5qhull__get_delaunay_info((&__pyx_v_info), __pyx_v_self, 1, 0);

  
  __pyx_t_4 = __Pyx_PyObject_IsTrue(__pyx_v_bruteforce); if (unlikely(__pyx_t_4 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1088; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__pyx_t_4) {

    
    __pyx_t_14 = (__pyx_v_x->dimensions[0]);
    for (__pyx_t_7 = 0; __pyx_t_7 < __pyx_t_14; __pyx_t_7+=1) {
      __pyx_v_k = __pyx_t_7;

      
      __pyx_v_isimplex = __pyx_f_5scipy_7spatial_5qhull__find_simplex_bruteforce((&__pyx_v_info), __pyx_v_c, (((double *)__pyx_v_x->data) + (__pyx_v_info.ndim * __pyx_v_k)), __pyx_v_eps);

      
      __pyx_t_15 = __pyx_v_k;
      __pyx_t_16 = -1;
      if (__pyx_t_15 < 0) {
        __pyx_t_15 += __pyx_bshape_0_out_;
        if (unlikely(__pyx_t_15 < 0)) __pyx_t_16 = 0;
      } else if (unlikely(__pyx_t_15 >= __pyx_bshape_0_out_)) __pyx_t_16 = 0;
      if (unlikely(__pyx_t_16 != -1)) {
        __Pyx_RaiseBufferIndexError(__pyx_t_16);
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1094; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      *__Pyx_BufPtrStrided1d(npy_int *, __pyx_bstruct_out_.buf, __pyx_t_15, __pyx_bstride_0_out_) = __pyx_v_isimplex;
    }
    goto __pyx_L7;
  }
   {

    
    __pyx_t_14 = (__pyx_v_x->dimensions[0]);
    for (__pyx_t_7 = 0; __pyx_t_7 < __pyx_t_14; __pyx_t_7+=1) {
      __pyx_v_k = __pyx_t_7;

      
      __pyx_v_isimplex = __pyx_f_5scipy_7spatial_5qhull__find_simplex((&__pyx_v_info), __pyx_v_c, (((double *)__pyx_v_x->data) + (__pyx_v_info.ndim * __pyx_v_k)), (&__pyx_v_start), __pyx_v_eps);

      
      __pyx_t_16 = __pyx_v_k;
      __pyx_t_17 = -1;
      if (__pyx_t_16 < 0) {
        __pyx_t_16 += __pyx_bshape_0_out_;
        if (unlikely(__pyx_t_16 < 0)) __pyx_t_17 = 0;
      } else if (unlikely(__pyx_t_16 >= __pyx_bshape_0_out_)) __pyx_t_17 = 0;
      if (unlikely(__pyx_t_17 != -1)) {
        __Pyx_RaiseBufferIndexError(__pyx_t_17);
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1099; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      *__Pyx_BufPtrStrided1d(npy_int *, __pyx_bstruct_out_.buf, __pyx_t_16, __pyx_bstride_0_out_) = __pyx_v_isimplex;
    }
  }
  __pyx_L7:;

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_12 = PyObject_GetAttr(__pyx_v_out, __pyx_n_s__reshape); if (unlikely(!__pyx_t_12)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1101; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_12);
  __pyx_t_2 = PySequence_GetSlice(__pyx_v_xi_shape, 0, -1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1101; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1101; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_Call(__pyx_t_12, __pyx_t_3, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1101; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_12); __pyx_t_12 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_r = __pyx_t_2;
  __pyx_t_2 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_XDECREF(__pyx_t_12);
  { PyObject *__pyx_type, *__pyx_value, *__pyx_tb;
    __Pyx_ErrFetch(&__pyx_type, &__pyx_value, &__pyx_tb);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_out_);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_x);
  __Pyx_ErrRestore(__pyx_type, __pyx_value, __pyx_tb);}
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.find_simplex");
  __pyx_r = NULL;
  goto __pyx_L2;
  __pyx_L0:;
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_out_);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_x);
  __pyx_L2:;
  __Pyx_DECREF((PyObject *)__pyx_v_x);
  __Pyx_DECREF((PyObject *)__pyx_v_out_);
  __Pyx_DECREF(__pyx_v_xi_shape);
  __Pyx_DECREF(__pyx_v_out);
  __Pyx_DECREF(__pyx_v_xi);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_plane_distance(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); 
static char __pyx_doc_5scipy_7spatial_5qhull_8Delaunay_plane_distance[] = "\n        plane_distance(xi)\n\n        Compute hyperplane distances to the point `xi` from all simplices.\n\n        ";
static PyMethodDef __pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_plane_distance = {__Pyx_NAMESTR("plane_distance"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_plane_distance, METH_VARARGS|METH_KEYWORDS, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull_8Delaunay_plane_distance)};
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_plane_distance(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyObject *__pyx_v_self = 0;
  PyObject *__pyx_v_xi = 0;
  PyArrayObject *__pyx_v_x;
  PyArrayObject *__pyx_v_out_;
  __pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t __pyx_v_info;
  double __pyx_v_z[(NPY_MAXDIMS + 1)];
  int __pyx_v_i;
  int __pyx_v_j;
  PyObject *__pyx_v_xi_shape;
  PyObject *__pyx_v_out;
  Py_buffer __pyx_bstruct_out_;
  Py_ssize_t __pyx_bstride_0_out_ = 0;
  Py_ssize_t __pyx_bstride_1_out_ = 0;
  Py_ssize_t __pyx_bshape_0_out_ = 0;
  Py_ssize_t __pyx_bshape_1_out_ = 0;
  Py_buffer __pyx_bstruct_x;
  Py_ssize_t __pyx_bstride_0_x = 0;
  Py_ssize_t __pyx_bstride_1_x = 0;
  Py_ssize_t __pyx_bshape_0_x = 0;
  Py_ssize_t __pyx_bshape_1_x = 0;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  int __pyx_t_4;
  PyObject *__pyx_t_5 = NULL;
  PyArrayObject *__pyx_t_6 = NULL;
  int __pyx_t_7;
  PyObject *__pyx_t_8 = NULL;
  PyObject *__pyx_t_9 = NULL;
  PyObject *__pyx_t_10 = NULL;
  PyObject *__pyx_t_11 = NULL;
  PyArrayObject *__pyx_t_12 = NULL;
  npy_intp __pyx_t_13;
  int __pyx_t_14;
  int __pyx_t_15;
  int __pyx_t_16;
  int __pyx_t_17;
  int __pyx_t_18;
  static PyObject **__pyx_pyargnames[] = {&__pyx_n_s__self,&__pyx_n_s__xi,0};
  __Pyx_RefNannySetupContext("plane_distance");
  __pyx_self = __pyx_self;
  if (unlikely(__pyx_kwds)) {
    Py_ssize_t kw_args = PyDict_Size(__pyx_kwds);
    PyObject* values[2] = {0,0};
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  2: values[1] = PyTuple_GET_ITEM(__pyx_args, 1);
      case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
      case  0: break;
      default: goto __pyx_L5_argtuple_error;
    }
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  0:
      values[0] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__self);
      if (likely(values[0])) kw_args--;
      else goto __pyx_L5_argtuple_error;
      case  1:
      values[1] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__xi);
      if (likely(values[1])) kw_args--;
      else {
        __Pyx_RaiseArgtupleInvalid("plane_distance", 1, 2, 2, 1); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1103; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
      }
    }
    if (unlikely(kw_args > 0)) {
      if (unlikely(__Pyx_ParseOptionalKeywords(__pyx_kwds, __pyx_pyargnames, 0, values, PyTuple_GET_SIZE(__pyx_args), "plane_distance") < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1103; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    }
    __pyx_v_self = values[0];
    __pyx_v_xi = values[1];
  } else if (PyTuple_GET_SIZE(__pyx_args) != 2) {
    goto __pyx_L5_argtuple_error;
  } else {
    __pyx_v_self = PyTuple_GET_ITEM(__pyx_args, 0);
    __pyx_v_xi = PyTuple_GET_ITEM(__pyx_args, 1);
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L5_argtuple_error:;
  __Pyx_RaiseArgtupleInvalid("plane_distance", 1, 2, 2, PyTuple_GET_SIZE(__pyx_args)); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1103; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  __pyx_L3_error:;
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.plane_distance");
  __Pyx_RefNannyFinishContext();
  return NULL;
  __pyx_L4_argument_unpacking_done:;
  __Pyx_INCREF(__pyx_v_xi);
  __pyx_v_x = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_out_ = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_xi_shape = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_out = Py_None; __Pyx_INCREF(Py_None);
  __pyx_bstruct_x.buf = NULL;
  __pyx_bstruct_out_.buf = NULL;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__shape); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1116; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = __Pyx_GetItemInt(__pyx_t_1, -1, sizeof(long), PyInt_FromLong); if (!__pyx_t_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1116; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__ndim); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1116; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_3 = PyObject_RichCompare(__pyx_t_2, __pyx_t_1, Py_NE); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1116; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_4 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_4 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1116; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (__pyx_t_4) {

    
    __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1117; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __Pyx_INCREF(((PyObject *)__pyx_kp_s_12));
    PyTuple_SET_ITEM(__pyx_t_3, 0, ((PyObject *)__pyx_kp_s_12));
    __Pyx_GIVEREF(((PyObject *)__pyx_kp_s_12));
    __pyx_t_1 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_3, NULL); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1117; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    __Pyx_Raise(__pyx_t_1, 0, 0);
    __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
    {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1117; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L6;
  }
  __pyx_L6:;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__shape); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1120; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_v_xi_shape);
  __pyx_v_xi_shape = __pyx_t_1;
  __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__reshape); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_3 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_3, __pyx_n_s__prod); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__shape); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_5 = PySequence_GetSlice(__pyx_t_3, 0, -1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_5 = 0;
  __pyx_t_5 = PyObject_Call(__pyx_t_2, __pyx_t_3, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__shape); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_2 = __Pyx_GetItemInt(__pyx_t_3, -1, sizeof(long), PyInt_FromLong); if (!__pyx_t_2) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyTuple_New(2); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  PyTuple_SET_ITEM(__pyx_t_3, 1, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_5 = 0;
  __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_Call(__pyx_t_1, __pyx_t_3, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1121; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(__pyx_v_xi);
  __pyx_v_xi = __pyx_t_2;
  __pyx_t_2 = 0;

  
  __pyx_t_2 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_2, __pyx_n_s__ascontiguousarray); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_GetAttr(__pyx_v_xi, __pyx_n_s__astype); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_5 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__double); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_5 = 0;
  __pyx_t_5 = PyObject_Call(__pyx_t_2, __pyx_t_1, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_5 = 0;
  __pyx_t_5 = PyObject_Call(__pyx_t_3, __pyx_t_1, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  if (!(likely(((__pyx_t_5) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_5, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_6 = ((PyArrayObject *)__pyx_t_5);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_x);
    __pyx_t_7 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_x, (PyObject*)__pyx_t_6, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_7 < 0)) {
      PyErr_Fetch(&__pyx_t_8, &__pyx_t_9, &__pyx_t_10);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_x, (PyObject*)__pyx_v_x, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_8); Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_10);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_8, __pyx_t_9, __pyx_t_10);
      }
    }
    __pyx_bstride_0_x = __pyx_bstruct_x.strides[0]; __pyx_bstride_1_x = __pyx_bstruct_x.strides[1];
    __pyx_bshape_0_x = __pyx_bstruct_x.shape[0]; __pyx_bshape_1_x = __pyx_bstruct_x.shape[1];
    if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1122; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_6 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_v_x));
  __pyx_v_x = ((PyArrayObject *)__pyx_t_5);
  __pyx_t_5 = 0;

  
  __pyx_f_5scipy_7spatial_5qhull__get_delaunay_info((&__pyx_v_info), __pyx_v_self, 0, 0);

  
  __pyx_t_5 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_1 = PyObject_GetAttr(__pyx_t_5, __pyx_n_s__zeros); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __pyx_t_5 = __Pyx_PyInt_to_py_Py_intptr_t((__pyx_v_x->dimensions[0])); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_3 = PyInt_FromLong(__pyx_v_info.nsimplex); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_2 = PyTuple_New(2); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  PyTuple_SET_ITEM(__pyx_t_2, 1, __pyx_t_3);
  __Pyx_GIVEREF(__pyx_t_3);
  __pyx_t_5 = 0;
  __pyx_t_3 = 0;
  __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  PyTuple_SET_ITEM(__pyx_t_3, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyDict_New(); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_2));
  __pyx_t_5 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_11 = PyObject_GetAttr(__pyx_t_5, __pyx_n_s__double); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_11);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  if (PyDict_SetItem(__pyx_t_2, ((PyObject *)__pyx_n_s__dtype), __pyx_t_11) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_11); __pyx_t_11 = 0;
  __pyx_t_11 = PyEval_CallObjectWithKeywords(__pyx_t_1, __pyx_t_3, ((PyObject *)__pyx_t_2)); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1126; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_11);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_2)); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_v_out);
  __pyx_v_out = __pyx_t_11;
  __pyx_t_11 = 0;

  
  if (!(likely(((__pyx_v_out) == Py_None) || likely(__Pyx_TypeTest(__pyx_v_out, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1127; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_12 = ((PyArrayObject *)__pyx_v_out);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_out_);
    __pyx_t_7 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_out_, (PyObject*)__pyx_t_12, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack);
    if (unlikely(__pyx_t_7 < 0)) {
      PyErr_Fetch(&__pyx_t_10, &__pyx_t_9, &__pyx_t_8);
      if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_out_, (PyObject*)__pyx_v_out_, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES| PyBUF_WRITABLE, 2, 0, __pyx_stack) == -1)) {
        Py_XDECREF(__pyx_t_10); Py_XDECREF(__pyx_t_9); Py_XDECREF(__pyx_t_8);
        __Pyx_RaiseBufferFallbackError();
      } else {
        PyErr_Restore(__pyx_t_10, __pyx_t_9, __pyx_t_8);
      }
    }
    __pyx_bstride_0_out_ = __pyx_bstruct_out_.strides[0]; __pyx_bstride_1_out_ = __pyx_bstruct_out_.strides[1];
    __pyx_bshape_0_out_ = __pyx_bstruct_out_.shape[0]; __pyx_bshape_1_out_ = __pyx_bstruct_out_.shape[1];
    if (unlikely(__pyx_t_7 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1127; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  __pyx_t_12 = 0;
  __Pyx_INCREF(__pyx_v_out);
  __Pyx_DECREF(((PyObject *)__pyx_v_out_));
  __pyx_v_out_ = ((PyArrayObject *)__pyx_v_out);

  
  __pyx_t_13 = (__pyx_v_x->dimensions[0]);
  for (__pyx_t_7 = 0; __pyx_t_7 < __pyx_t_13; __pyx_t_7+=1) {
    __pyx_v_i = __pyx_t_7;

    
    __pyx_t_14 = __pyx_v_info.nsimplex;
    for (__pyx_t_15 = 0; __pyx_t_15 < __pyx_t_14; __pyx_t_15+=1) {
      __pyx_v_j = __pyx_t_15;

      
      __pyx_f_5scipy_7spatial_5qhull__lift_point((&__pyx_v_info), (((double *)__pyx_v_x->data) + (__pyx_v_info.ndim * __pyx_v_i)), __pyx_v_z);

      
      __pyx_t_16 = __pyx_v_i;
      __pyx_t_17 = __pyx_v_j;
      __pyx_t_18 = -1;
      if (__pyx_t_16 < 0) {
        __pyx_t_16 += __pyx_bshape_0_out_;
        if (unlikely(__pyx_t_16 < 0)) __pyx_t_18 = 0;
      } else if (unlikely(__pyx_t_16 >= __pyx_bshape_0_out_)) __pyx_t_18 = 0;
      if (__pyx_t_17 < 0) {
        __pyx_t_17 += __pyx_bshape_1_out_;
        if (unlikely(__pyx_t_17 < 0)) __pyx_t_18 = 1;
      } else if (unlikely(__pyx_t_17 >= __pyx_bshape_1_out_)) __pyx_t_18 = 1;
      if (unlikely(__pyx_t_18 != -1)) {
        __Pyx_RaiseBufferIndexError(__pyx_t_18);
        {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1132; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      *__Pyx_BufPtrStrided2d(__pyx_t_5numpy_double_t *, __pyx_bstruct_out_.buf, __pyx_t_16, __pyx_bstride_0_out_, __pyx_t_17, __pyx_bstride_1_out_) = __pyx_f_5scipy_7spatial_5qhull__distplane((&__pyx_v_info), __pyx_v_j, __pyx_v_z);
    }
  }

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_11 = PyObject_GetAttr(__pyx_v_out, __pyx_n_s__reshape); if (unlikely(!__pyx_t_11)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1134; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_11);
  __pyx_t_2 = PySequence_GetSlice(__pyx_v_xi_shape, 0, -1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1134; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_3 = PyObject_GetAttr(__pyx_v_self, __pyx_n_s__nsimplex); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1134; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1134; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_t_3);
  __Pyx_GIVEREF(__pyx_t_3);
  __pyx_t_3 = 0;
  __pyx_t_3 = PyNumber_Add(__pyx_t_2, __pyx_t_1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1134; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyTuple_New(1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1134; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  PyTuple_SET_ITEM(__pyx_t_1, 0, __pyx_t_3);
  __Pyx_GIVEREF(__pyx_t_3);
  __pyx_t_3 = 0;
  __pyx_t_3 = PyObject_Call(__pyx_t_11, __pyx_t_1, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1134; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_11); __pyx_t_11 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_r = __pyx_t_3;
  __pyx_t_3 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_XDECREF(__pyx_t_11);
  { PyObject *__pyx_type, *__pyx_value, *__pyx_tb;
    __Pyx_ErrFetch(&__pyx_type, &__pyx_value, &__pyx_tb);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_out_);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_x);
  __Pyx_ErrRestore(__pyx_type, __pyx_value, __pyx_tb);}
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.plane_distance");
  __pyx_r = NULL;
  goto __pyx_L2;
  __pyx_L0:;
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_out_);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_x);
  __pyx_L2:;
  __Pyx_DECREF((PyObject *)__pyx_v_x);
  __Pyx_DECREF((PyObject *)__pyx_v_out_);
  __Pyx_DECREF(__pyx_v_xi_shape);
  __Pyx_DECREF(__pyx_v_out);
  __Pyx_DECREF(__pyx_v_xi);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_lift_points(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); 
static char __pyx_doc_5scipy_7spatial_5qhull_8Delaunay_lift_points[] = "\n        lift_points(tri, x)\n\n        Lift points to the Qhull paraboloid.\n\n        ";
static PyMethodDef __pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_lift_points = {__Pyx_NAMESTR("lift_points"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_lift_points, METH_VARARGS|METH_KEYWORDS, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull_8Delaunay_lift_points)};
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_8Delaunay_lift_points(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyObject *__pyx_v_tri = 0;
  PyObject *__pyx_v_x = 0;
  PyObject *__pyx_v_z;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  static PyObject **__pyx_pyargnames[] = {&__pyx_n_s__tri,&__pyx_n_s__x,0};
  __Pyx_RefNannySetupContext("lift_points");
  __pyx_self = __pyx_self;
  if (unlikely(__pyx_kwds)) {
    Py_ssize_t kw_args = PyDict_Size(__pyx_kwds);
    PyObject* values[2] = {0,0};
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  2: values[1] = PyTuple_GET_ITEM(__pyx_args, 1);
      case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
      case  0: break;
      default: goto __pyx_L5_argtuple_error;
    }
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  0:
      values[0] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__tri);
      if (likely(values[0])) kw_args--;
      else goto __pyx_L5_argtuple_error;
      case  1:
      values[1] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__x);
      if (likely(values[1])) kw_args--;
      else {
        __Pyx_RaiseArgtupleInvalid("lift_points", 1, 2, 2, 1); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1136; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
      }
    }
    if (unlikely(kw_args > 0)) {
      if (unlikely(__Pyx_ParseOptionalKeywords(__pyx_kwds, __pyx_pyargnames, 0, values, PyTuple_GET_SIZE(__pyx_args), "lift_points") < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1136; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    }
    __pyx_v_tri = values[0];
    __pyx_v_x = values[1];
  } else if (PyTuple_GET_SIZE(__pyx_args) != 2) {
    goto __pyx_L5_argtuple_error;
  } else {
    __pyx_v_tri = PyTuple_GET_ITEM(__pyx_args, 0);
    __pyx_v_x = PyTuple_GET_ITEM(__pyx_args, 1);
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L5_argtuple_error:;
  __Pyx_RaiseArgtupleInvalid("lift_points", 1, 2, 2, PyTuple_GET_SIZE(__pyx_args)); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1136; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  __pyx_L3_error:;
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.lift_points");
  __Pyx_RefNannyFinishContext();
  return NULL;
  __pyx_L4_argument_unpacking_done:;
  __pyx_v_z = Py_None; __Pyx_INCREF(Py_None);

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__zeros); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_x, __pyx_n_s__shape); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_3 = PySequence_GetSlice(__pyx_t_1, 0, -1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_x, __pyx_n_s__shape); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_4 = __Pyx_GetItemInt(__pyx_t_1, -1, sizeof(long), PyInt_FromLong); if (!__pyx_t_4) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyNumber_Add(__pyx_t_4, __pyx_int_1); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  PyTuple_SET_ITEM(__pyx_t_4, 0, __pyx_t_1);
  __Pyx_GIVEREF(__pyx_t_1);
  __pyx_t_1 = 0;
  __pyx_t_1 = PyNumber_Add(__pyx_t_3, __pyx_t_4); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  PyTuple_SET_ITEM(__pyx_t_4, 0, __pyx_t_1);
  __Pyx_GIVEREF(__pyx_t_1);
  __pyx_t_1 = 0;
  __pyx_t_1 = PyDict_New(); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_1));
  __pyx_t_3 = __Pyx_GetName(__pyx_m, __pyx_n_s__np); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_5 = PyObject_GetAttr(__pyx_t_3, __pyx_n_s__double); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_n_s__dtype), __pyx_t_5) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __pyx_t_5 = PyEval_CallObjectWithKeywords(__pyx_t_2, __pyx_t_4, ((PyObject *)__pyx_t_1)); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1143; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_1)); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_v_z);
  __pyx_v_z = __pyx_t_5;
  __pyx_t_5 = 0;

  
  __pyx_t_5 = PySlice_New(Py_None, __pyx_int_neg_1, Py_None); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1144; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_1 = PyTuple_New(2); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1144; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_INCREF(Py_Ellipsis);
  PyTuple_SET_ITEM(__pyx_t_1, 0, Py_Ellipsis);
  __Pyx_GIVEREF(Py_Ellipsis);
  PyTuple_SET_ITEM(__pyx_t_1, 1, __pyx_t_5);
  __Pyx_GIVEREF(__pyx_t_5);
  __pyx_t_5 = 0;
  if (PyObject_SetItem(__pyx_v_z, __pyx_t_1, __pyx_v_x) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1144; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_t_1 = PyNumber_Power(__pyx_v_x, __pyx_int_2, Py_None); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1145; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_5 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__sum); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1145; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyDict_New(); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1145; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_1));
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_n_s__axis), __pyx_int_neg_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1145; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_4 = PyEval_CallObjectWithKeywords(__pyx_t_5, ((PyObject *)__pyx_empty_tuple), ((PyObject *)__pyx_t_1)); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1145; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_1)); __pyx_t_1 = 0;
  __pyx_t_1 = PyTuple_New(2); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1145; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_INCREF(Py_Ellipsis);
  PyTuple_SET_ITEM(__pyx_t_1, 0, Py_Ellipsis);
  __Pyx_GIVEREF(Py_Ellipsis);
  __Pyx_INCREF(__pyx_int_neg_1);
  PyTuple_SET_ITEM(__pyx_t_1, 1, __pyx_int_neg_1);
  __Pyx_GIVEREF(__pyx_int_neg_1);
  if (PyObject_SetItem(__pyx_v_z, __pyx_t_1, __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1145; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;

  
  __pyx_t_4 = PyObject_GetAttr(__pyx_v_tri, __pyx_n_s__paraboloid_scale); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1146; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_1 = PyTuple_New(2); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1146; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_INCREF(Py_Ellipsis);
  PyTuple_SET_ITEM(__pyx_t_1, 0, Py_Ellipsis);
  __Pyx_GIVEREF(Py_Ellipsis);
  __Pyx_INCREF(__pyx_int_neg_1);
  PyTuple_SET_ITEM(__pyx_t_1, 1, __pyx_int_neg_1);
  __Pyx_GIVEREF(__pyx_int_neg_1);
  __pyx_t_5 = PyObject_GetItem(__pyx_v_z, __pyx_t_1); if (!__pyx_t_5) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1146; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_2 = PyNumber_InPlaceMultiply(__pyx_t_5, __pyx_t_4); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1146; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  if (PyObject_SetItem(__pyx_v_z, __pyx_t_1, __pyx_t_2) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1146; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_tri, __pyx_n_s__paraboloid_shift); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1147; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyTuple_New(2); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1147; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_INCREF(Py_Ellipsis);
  PyTuple_SET_ITEM(__pyx_t_2, 0, Py_Ellipsis);
  __Pyx_GIVEREF(Py_Ellipsis);
  __Pyx_INCREF(__pyx_int_neg_1);
  PyTuple_SET_ITEM(__pyx_t_2, 1, __pyx_int_neg_1);
  __Pyx_GIVEREF(__pyx_int_neg_1);
  __pyx_t_5 = PyObject_GetItem(__pyx_v_z, __pyx_t_2); if (!__pyx_t_5) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1147; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_5);
  __pyx_t_4 = PyNumber_InPlaceAdd(__pyx_t_5, __pyx_t_1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1147; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
  if (PyObject_SetItem(__pyx_v_z, __pyx_t_2, __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1147; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;

  
  __Pyx_XDECREF(__pyx_r);
  __Pyx_INCREF(__pyx_v_z);
  __pyx_r = __pyx_v_z;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_AddTraceback("scipy.spatial.qhull.Delaunay.lift_points");
  __pyx_r = NULL;
  __pyx_L0:;
  __Pyx_DECREF(__pyx_v_z);
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static PyObject *__pyx_pf_5scipy_7spatial_5qhull_tsearch(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds); 
static char __pyx_doc_5scipy_7spatial_5qhull_tsearch[] = "\n    tsearch(tri, xi)\n\n    Find simplices containing the given points. This function does the\n    same thing as Delaunay.find_simplex.\n\n    .. versionadded:: 0.9\n\n    See Also\n    --------\n    Delaunay.find_simplex\n\n    ";
static PyObject *__pyx_pf_5scipy_7spatial_5qhull_tsearch(PyObject *__pyx_self, PyObject *__pyx_args, PyObject *__pyx_kwds) {
  PyObject *__pyx_v_tri = 0;
  PyObject *__pyx_v_xi = 0;
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  static PyObject **__pyx_pyargnames[] = {&__pyx_n_s__tri,&__pyx_n_s__xi,0};
  __Pyx_RefNannySetupContext("tsearch");
  __pyx_self = __pyx_self;
  if (unlikely(__pyx_kwds)) {
    Py_ssize_t kw_args = PyDict_Size(__pyx_kwds);
    PyObject* values[2] = {0,0};
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  2: values[1] = PyTuple_GET_ITEM(__pyx_args, 1);
      case  1: values[0] = PyTuple_GET_ITEM(__pyx_args, 0);
      case  0: break;
      default: goto __pyx_L5_argtuple_error;
    }
    switch (PyTuple_GET_SIZE(__pyx_args)) {
      case  0:
      values[0] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__tri);
      if (likely(values[0])) kw_args--;
      else goto __pyx_L5_argtuple_error;
      case  1:
      values[1] = PyDict_GetItem(__pyx_kwds, __pyx_n_s__xi);
      if (likely(values[1])) kw_args--;
      else {
        __Pyx_RaiseArgtupleInvalid("tsearch", 1, 2, 2, 1); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1151; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
      }
    }
    if (unlikely(kw_args > 0)) {
      if (unlikely(__Pyx_ParseOptionalKeywords(__pyx_kwds, __pyx_pyargnames, 0, values, PyTuple_GET_SIZE(__pyx_args), "tsearch") < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1151; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
    }
    __pyx_v_tri = values[0];
    __pyx_v_xi = values[1];
  } else if (PyTuple_GET_SIZE(__pyx_args) != 2) {
    goto __pyx_L5_argtuple_error;
  } else {
    __pyx_v_tri = PyTuple_GET_ITEM(__pyx_args, 0);
    __pyx_v_xi = PyTuple_GET_ITEM(__pyx_args, 1);
  }
  goto __pyx_L4_argument_unpacking_done;
  __pyx_L5_argtuple_error:;
  __Pyx_RaiseArgtupleInvalid("tsearch", 1, 2, 2, PyTuple_GET_SIZE(__pyx_args)); {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1151; __pyx_clineno = __LINE__; goto __pyx_L3_error;}
  __pyx_L3_error:;
  __Pyx_AddTraceback("scipy.spatial.qhull.tsearch");
  __Pyx_RefNannyFinishContext();
  return NULL;
  __pyx_L4_argument_unpacking_done:;

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_tri, __pyx_n_s__find_simplex); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1165; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyTuple_New(1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1165; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_INCREF(__pyx_v_xi);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_v_xi);
  __Pyx_GIVEREF(__pyx_v_xi);
  __pyx_t_3 = PyObject_Call(__pyx_t_1, __pyx_t_2, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1165; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  __pyx_r = __pyx_t_3;
  __pyx_t_3 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_AddTraceback("scipy.spatial.qhull.tsearch");
  __pyx_r = NULL;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static  void __pyx_f_5scipy_7spatial_5qhull__get_delaunay_info(__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *__pyx_v_info, PyObject *__pyx_v_obj, int __pyx_v_compute_transform, int __pyx_v_compute_vertex_to_simplex) {
  PyArrayObject *__pyx_v_transform;
  PyArrayObject *__pyx_v_vertex_to_simplex;
  PyArrayObject *__pyx_v_points = 0;
  PyArrayObject *__pyx_v_vertices = 0;
  PyArrayObject *__pyx_v_neighbors = 0;
  PyArrayObject *__pyx_v_equations = 0;
  PyArrayObject *__pyx_v_min_bound = 0;
  PyArrayObject *__pyx_v_max_bound = 0;
  Py_buffer __pyx_bstruct_neighbors;
  Py_ssize_t __pyx_bstride_0_neighbors = 0;
  Py_ssize_t __pyx_bstride_1_neighbors = 0;
  Py_ssize_t __pyx_bshape_0_neighbors = 0;
  Py_ssize_t __pyx_bshape_1_neighbors = 0;
  Py_buffer __pyx_bstruct_transform;
  Py_ssize_t __pyx_bstride_0_transform = 0;
  Py_ssize_t __pyx_bstride_1_transform = 0;
  Py_ssize_t __pyx_bstride_2_transform = 0;
  Py_ssize_t __pyx_bshape_0_transform = 0;
  Py_ssize_t __pyx_bshape_1_transform = 0;
  Py_ssize_t __pyx_bshape_2_transform = 0;
  Py_buffer __pyx_bstruct_vertices;
  Py_ssize_t __pyx_bstride_0_vertices = 0;
  Py_ssize_t __pyx_bstride_1_vertices = 0;
  Py_ssize_t __pyx_bshape_0_vertices = 0;
  Py_ssize_t __pyx_bshape_1_vertices = 0;
  Py_buffer __pyx_bstruct_points;
  Py_ssize_t __pyx_bstride_0_points = 0;
  Py_ssize_t __pyx_bstride_1_points = 0;
  Py_ssize_t __pyx_bshape_0_points = 0;
  Py_ssize_t __pyx_bshape_1_points = 0;
  Py_buffer __pyx_bstruct_vertex_to_simplex;
  Py_ssize_t __pyx_bstride_0_vertex_to_simplex = 0;
  Py_ssize_t __pyx_bshape_0_vertex_to_simplex = 0;
  Py_buffer __pyx_bstruct_min_bound;
  Py_ssize_t __pyx_bstride_0_min_bound = 0;
  Py_ssize_t __pyx_bshape_0_min_bound = 0;
  Py_buffer __pyx_bstruct_max_bound;
  Py_ssize_t __pyx_bstride_0_max_bound = 0;
  Py_ssize_t __pyx_bshape_0_max_bound = 0;
  Py_buffer __pyx_bstruct_equations;
  Py_ssize_t __pyx_bstride_0_equations = 0;
  Py_ssize_t __pyx_bstride_1_equations = 0;
  Py_ssize_t __pyx_bshape_0_equations = 0;
  Py_ssize_t __pyx_bshape_1_equations = 0;
  PyObject *__pyx_t_1 = NULL;
  PyArrayObject *__pyx_t_2 = NULL;
  PyArrayObject *__pyx_t_3 = NULL;
  PyArrayObject *__pyx_t_4 = NULL;
  PyArrayObject *__pyx_t_5 = NULL;
  PyArrayObject *__pyx_t_6 = NULL;
  PyArrayObject *__pyx_t_7 = NULL;
  double __pyx_t_8;
  PyArrayObject *__pyx_t_9 = NULL;
  int __pyx_t_10;
  PyObject *__pyx_t_11 = NULL;
  PyObject *__pyx_t_12 = NULL;
  PyObject *__pyx_t_13 = NULL;
  PyArrayObject *__pyx_t_14 = NULL;
  __Pyx_RefNannySetupContext("_get_delaunay_info");
  __pyx_v_transform = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_vertex_to_simplex = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_bstruct_transform.buf = NULL;
  __pyx_bstruct_vertex_to_simplex.buf = NULL;
  __pyx_bstruct_points.buf = NULL;
  __pyx_bstruct_vertices.buf = NULL;
  __pyx_bstruct_neighbors.buf = NULL;
  __pyx_bstruct_equations.buf = NULL;
  __pyx_bstruct_min_bound.buf = NULL;
  __pyx_bstruct_max_bound.buf = NULL;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__points); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1178; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1178; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_2 = ((PyArrayObject *)__pyx_t_1);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_points, (PyObject*)__pyx_t_2, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
      __pyx_v_points = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None); __pyx_bstruct_points.buf = NULL;
      {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1178; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    } else {__pyx_bstride_0_points = __pyx_bstruct_points.strides[0]; __pyx_bstride_1_points = __pyx_bstruct_points.strides[1];
      __pyx_bshape_0_points = __pyx_bstruct_points.shape[0]; __pyx_bshape_1_points = __pyx_bstruct_points.shape[1];
    }
  }
  __pyx_t_2 = 0;
  __pyx_v_points = ((PyArrayObject *)__pyx_t_1);
  __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__vertices); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1179; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1179; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_3 = ((PyArrayObject *)__pyx_t_1);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_vertices, (PyObject*)__pyx_t_3, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
      __pyx_v_vertices = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None); __pyx_bstruct_vertices.buf = NULL;
      {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1179; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    } else {__pyx_bstride_0_vertices = __pyx_bstruct_vertices.strides[0]; __pyx_bstride_1_vertices = __pyx_bstruct_vertices.strides[1];
      __pyx_bshape_0_vertices = __pyx_bstruct_vertices.shape[0]; __pyx_bshape_1_vertices = __pyx_bstruct_vertices.shape[1];
    }
  }
  __pyx_t_3 = 0;
  __pyx_v_vertices = ((PyArrayObject *)__pyx_t_1);
  __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__neighbors); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1180; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1180; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_4 = ((PyArrayObject *)__pyx_t_1);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_neighbors, (PyObject*)__pyx_t_4, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
      __pyx_v_neighbors = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None); __pyx_bstruct_neighbors.buf = NULL;
      {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1180; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    } else {__pyx_bstride_0_neighbors = __pyx_bstruct_neighbors.strides[0]; __pyx_bstride_1_neighbors = __pyx_bstruct_neighbors.strides[1];
      __pyx_bshape_0_neighbors = __pyx_bstruct_neighbors.shape[0]; __pyx_bshape_1_neighbors = __pyx_bstruct_neighbors.shape[1];
    }
  }
  __pyx_t_4 = 0;
  __pyx_v_neighbors = ((PyArrayObject *)__pyx_t_1);
  __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__equations); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1181; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1181; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_5 = ((PyArrayObject *)__pyx_t_1);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_equations, (PyObject*)__pyx_t_5, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 2, 0, __pyx_stack) == -1)) {
      __pyx_v_equations = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None); __pyx_bstruct_equations.buf = NULL;
      {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1181; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    } else {__pyx_bstride_0_equations = __pyx_bstruct_equations.strides[0]; __pyx_bstride_1_equations = __pyx_bstruct_equations.strides[1];
      __pyx_bshape_0_equations = __pyx_bstruct_equations.shape[0]; __pyx_bshape_1_equations = __pyx_bstruct_equations.shape[1];
    }
  }
  __pyx_t_5 = 0;
  __pyx_v_equations = ((PyArrayObject *)__pyx_t_1);
  __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__min_bound); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1182; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1182; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_6 = ((PyArrayObject *)__pyx_t_1);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_min_bound, (PyObject*)__pyx_t_6, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 1, 0, __pyx_stack) == -1)) {
      __pyx_v_min_bound = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None); __pyx_bstruct_min_bound.buf = NULL;
      {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1182; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    } else {__pyx_bstride_0_min_bound = __pyx_bstruct_min_bound.strides[0];
      __pyx_bshape_0_min_bound = __pyx_bstruct_min_bound.shape[0];
    }
  }
  __pyx_t_6 = 0;
  __pyx_v_min_bound = ((PyArrayObject *)__pyx_t_1);
  __pyx_t_1 = 0;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__max_bound); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1183; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1183; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_7 = ((PyArrayObject *)__pyx_t_1);
  {
    __Pyx_BufFmt_StackElem __pyx_stack[1];
    if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_max_bound, (PyObject*)__pyx_t_7, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 1, 0, __pyx_stack) == -1)) {
      __pyx_v_max_bound = ((PyArrayObject *)Py_None); __Pyx_INCREF(Py_None); __pyx_bstruct_max_bound.buf = NULL;
      {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1183; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    } else {__pyx_bstride_0_max_bound = __pyx_bstruct_max_bound.strides[0];
      __pyx_bshape_0_max_bound = __pyx_bstruct_max_bound.shape[0];
    }
  }
  __pyx_t_7 = 0;
  __pyx_v_max_bound = ((PyArrayObject *)__pyx_t_1);
  __pyx_t_1 = 0;

  
  __pyx_v_info->ndim = (__pyx_v_points->dimensions[1]);

  
  __pyx_v_info->npoints = (__pyx_v_points->dimensions[0]);

  
  __pyx_v_info->nsimplex = (__pyx_v_vertices->dimensions[0]);

  
  __pyx_v_info->points = ((double *)__pyx_v_points->data);

  
  __pyx_v_info->vertices = ((int *)__pyx_v_vertices->data);

  
  __pyx_v_info->neighbors = ((int *)__pyx_v_neighbors->data);

  
  __pyx_v_info->equations = ((double *)__pyx_v_equations->data);

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__paraboloid_scale); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1192; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_8 = __pyx_PyFloat_AsDouble(__pyx_t_1); if (unlikely((__pyx_t_8 == (double)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1192; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_v_info->paraboloid_scale = __pyx_t_8;

  
  __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__paraboloid_shift); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1193; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_8 = __pyx_PyFloat_AsDouble(__pyx_t_1); if (unlikely((__pyx_t_8 == (double)-1) && PyErr_Occurred())) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1193; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_v_info->paraboloid_shift = __pyx_t_8;

  
  if (__pyx_v_compute_transform) {

    
    __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__transform); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1195; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1195; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __pyx_t_9 = ((PyArrayObject *)__pyx_t_1);
    {
      __Pyx_BufFmt_StackElem __pyx_stack[1];
      __Pyx_SafeReleaseBuffer(&__pyx_bstruct_transform);
      __pyx_t_10 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_transform, (PyObject*)__pyx_t_9, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 3, 0, __pyx_stack);
      if (unlikely(__pyx_t_10 < 0)) {
        PyErr_Fetch(&__pyx_t_11, &__pyx_t_12, &__pyx_t_13);
        if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_transform, (PyObject*)__pyx_v_transform, &__Pyx_TypeInfo_nn___pyx_t_5numpy_double_t, PyBUF_FORMAT| PyBUF_STRIDES, 3, 0, __pyx_stack) == -1)) {
          Py_XDECREF(__pyx_t_11); Py_XDECREF(__pyx_t_12); Py_XDECREF(__pyx_t_13);
          __Pyx_RaiseBufferFallbackError();
        } else {
          PyErr_Restore(__pyx_t_11, __pyx_t_12, __pyx_t_13);
        }
      }
      __pyx_bstride_0_transform = __pyx_bstruct_transform.strides[0]; __pyx_bstride_1_transform = __pyx_bstruct_transform.strides[1]; __pyx_bstride_2_transform = __pyx_bstruct_transform.strides[2];
      __pyx_bshape_0_transform = __pyx_bstruct_transform.shape[0]; __pyx_bshape_1_transform = __pyx_bstruct_transform.shape[1]; __pyx_bshape_2_transform = __pyx_bstruct_transform.shape[2];
      if (unlikely(__pyx_t_10 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1195; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    }
    __pyx_t_9 = 0;
    __Pyx_DECREF(((PyObject *)__pyx_v_transform));
    __pyx_v_transform = ((PyArrayObject *)__pyx_t_1);
    __pyx_t_1 = 0;

    
    __pyx_v_info->transform = ((double *)__pyx_v_transform->data);
    goto __pyx_L3;
  }
   {

    
    __pyx_v_info->transform = NULL;
  }
  __pyx_L3:;

  
  if (__pyx_v_compute_vertex_to_simplex) {

    
    __pyx_t_1 = PyObject_GetAttr(__pyx_v_obj, __pyx_n_s__vertex_to_simplex); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1200; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_1);
    if (!(likely(((__pyx_t_1) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_1, __pyx_ptype_5numpy_ndarray))))) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1200; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __pyx_t_14 = ((PyArrayObject *)__pyx_t_1);
    {
      __Pyx_BufFmt_StackElem __pyx_stack[1];
      __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertex_to_simplex);
      __pyx_t_10 = __Pyx_GetBufferAndValidate(&__pyx_bstruct_vertex_to_simplex, (PyObject*)__pyx_t_14, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 1, 0, __pyx_stack);
      if (unlikely(__pyx_t_10 < 0)) {
        PyErr_Fetch(&__pyx_t_13, &__pyx_t_12, &__pyx_t_11);
        if (unlikely(__Pyx_GetBufferAndValidate(&__pyx_bstruct_vertex_to_simplex, (PyObject*)__pyx_v_vertex_to_simplex, &__Pyx_TypeInfo_nn_npy_int, PyBUF_FORMAT| PyBUF_STRIDES, 1, 0, __pyx_stack) == -1)) {
          Py_XDECREF(__pyx_t_13); Py_XDECREF(__pyx_t_12); Py_XDECREF(__pyx_t_11);
          __Pyx_RaiseBufferFallbackError();
        } else {
          PyErr_Restore(__pyx_t_13, __pyx_t_12, __pyx_t_11);
        }
      }
      __pyx_bstride_0_vertex_to_simplex = __pyx_bstruct_vertex_to_simplex.strides[0];
      __pyx_bshape_0_vertex_to_simplex = __pyx_bstruct_vertex_to_simplex.shape[0];
      if (unlikely(__pyx_t_10 < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1200; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    }
    __pyx_t_14 = 0;
    __Pyx_DECREF(((PyObject *)__pyx_v_vertex_to_simplex));
    __pyx_v_vertex_to_simplex = ((PyArrayObject *)__pyx_t_1);
    __pyx_t_1 = 0;

    
    __pyx_v_info->vertex_to_simplex = ((int *)__pyx_v_vertex_to_simplex->data);
    goto __pyx_L4;
  }
   {

    
    __pyx_v_info->vertex_to_simplex = NULL;
  }
  __pyx_L4:;

  
  __pyx_v_info->min_bound = ((double *)__pyx_v_min_bound->data);

  
  __pyx_v_info->max_bound = ((double *)__pyx_v_max_bound->data);

  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  { PyObject *__pyx_type, *__pyx_value, *__pyx_tb;
    __Pyx_ErrFetch(&__pyx_type, &__pyx_value, &__pyx_tb);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_neighbors);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_transform);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_points);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertex_to_simplex);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_min_bound);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_max_bound);
    __Pyx_SafeReleaseBuffer(&__pyx_bstruct_equations);
  __Pyx_ErrRestore(__pyx_type, __pyx_value, __pyx_tb);}
  __Pyx_WriteUnraisable("scipy.spatial.qhull._get_delaunay_info");
  goto __pyx_L2;
  __pyx_L0:;
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_neighbors);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_transform);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertices);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_points);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_vertex_to_simplex);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_min_bound);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_max_bound);
  __Pyx_SafeReleaseBuffer(&__pyx_bstruct_equations);
  __pyx_L2:;
  __Pyx_DECREF((PyObject *)__pyx_v_transform);
  __Pyx_DECREF((PyObject *)__pyx_v_vertex_to_simplex);
  __Pyx_XDECREF((PyObject *)__pyx_v_points);
  __Pyx_XDECREF((PyObject *)__pyx_v_vertices);
  __Pyx_XDECREF((PyObject *)__pyx_v_neighbors);
  __Pyx_XDECREF((PyObject *)__pyx_v_equations);
  __Pyx_XDECREF((PyObject *)__pyx_v_min_bound);
  __Pyx_XDECREF((PyObject *)__pyx_v_max_bound);
  __Pyx_RefNannyFinishContext();
}



static CYTHON_UNUSED int __pyx_pf_5numpy_7ndarray___getbuffer__(PyObject *__pyx_v_self, Py_buffer *__pyx_v_info, int __pyx_v_flags); 
static CYTHON_UNUSED int __pyx_pf_5numpy_7ndarray___getbuffer__(PyObject *__pyx_v_self, Py_buffer *__pyx_v_info, int __pyx_v_flags) {
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
  int __pyx_t_1;
  int __pyx_t_2;
  int __pyx_t_3;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  int __pyx_t_6;
  int __pyx_t_7;
  int __pyx_t_8;
  char *__pyx_t_9;
  __Pyx_RefNannySetupContext("__getbuffer__");
  if (__pyx_v_info == NULL) return 0;
  __pyx_v_info->obj = Py_None; __Pyx_INCREF(Py_None);
  __Pyx_GIVEREF(__pyx_v_info->obj);

  
  __pyx_v_endian_detector = 1;

  
  __pyx_v_little_endian = ((((char *)(&__pyx_v_endian_detector))[0]) != 0);

  
  __pyx_v_ndim = PyArray_NDIM(((PyArrayObject *)__pyx_v_self));

  
  __pyx_t_1 = ((sizeof(npy_intp)) != (sizeof(Py_ssize_t)));
  if (__pyx_t_1) {

    
    __pyx_v_copy_shape = 1;
    goto __pyx_L5;
  }
   {

    
    __pyx_v_copy_shape = 0;
  }
  __pyx_L5:;

  
  __pyx_t_1 = ((__pyx_v_flags & PyBUF_C_CONTIGUOUS) == PyBUF_C_CONTIGUOUS);
  if (__pyx_t_1) {

    
    __pyx_t_2 = (!PyArray_CHKFLAGS(((PyArrayObject *)__pyx_v_self), NPY_C_CONTIGUOUS));
    __pyx_t_3 = __pyx_t_2;
  } else {
    __pyx_t_3 = __pyx_t_1;
  }
  if (__pyx_t_3) {

    
    __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 206; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __Pyx_INCREF(((PyObject *)__pyx_kp_u_13));
    PyTuple_SET_ITEM(__pyx_t_4, 0, ((PyObject *)__pyx_kp_u_13));
    __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_13));
    __pyx_t_5 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_4, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 206; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_5);
    __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
    __Pyx_Raise(__pyx_t_5, 0, 0);
    __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
    {__pyx_filename = __pyx_f[1]; __pyx_lineno = 206; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L6;
  }
  __pyx_L6:;

  
  __pyx_t_3 = ((__pyx_v_flags & PyBUF_F_CONTIGUOUS) == PyBUF_F_CONTIGUOUS);
  if (__pyx_t_3) {

    
    __pyx_t_1 = (!PyArray_CHKFLAGS(((PyArrayObject *)__pyx_v_self), NPY_F_CONTIGUOUS));
    __pyx_t_2 = __pyx_t_1;
  } else {
    __pyx_t_2 = __pyx_t_3;
  }
  if (__pyx_t_2) {

    
    __pyx_t_5 = PyTuple_New(1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 210; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_5);
    __Pyx_INCREF(((PyObject *)__pyx_kp_u_14));
    PyTuple_SET_ITEM(__pyx_t_5, 0, ((PyObject *)__pyx_kp_u_14));
    __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_14));
    __pyx_t_4 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_5, NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 210; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
    __Pyx_Raise(__pyx_t_4, 0, 0);
    __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
    {__pyx_filename = __pyx_f[1]; __pyx_lineno = 210; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    goto __pyx_L7;
  }
  __pyx_L7:;

  
  __pyx_v_info->buf = PyArray_DATA(((PyArrayObject *)__pyx_v_self));

  
  __pyx_v_info->ndim = __pyx_v_ndim;

  
  if (__pyx_v_copy_shape) {

    
    __pyx_v_info->strides = ((Py_ssize_t *)malloc((((sizeof(Py_ssize_t)) * __pyx_v_ndim) * 2)));

    
    __pyx_v_info->shape = (__pyx_v_info->strides + __pyx_v_ndim);

    
    __pyx_t_6 = __pyx_v_ndim;
    for (__pyx_t_7 = 0; __pyx_t_7 < __pyx_t_6; __pyx_t_7+=1) {
      __pyx_v_i = __pyx_t_7;

      
      (__pyx_v_info->strides[__pyx_v_i]) = (PyArray_STRIDES(((PyArrayObject *)__pyx_v_self))[__pyx_v_i]);

      
      (__pyx_v_info->shape[__pyx_v_i]) = (PyArray_DIMS(((PyArrayObject *)__pyx_v_self))[__pyx_v_i]);
    }
    goto __pyx_L8;
  }
   {

    
    __pyx_v_info->strides = ((Py_ssize_t *)PyArray_STRIDES(((PyArrayObject *)__pyx_v_self)));

    
    __pyx_v_info->shape = ((Py_ssize_t *)PyArray_DIMS(((PyArrayObject *)__pyx_v_self)));
  }
  __pyx_L8:;

  
  __pyx_v_info->suboffsets = NULL;

  
  __pyx_v_info->itemsize = PyArray_ITEMSIZE(((PyArrayObject *)__pyx_v_self));

  
  __pyx_v_info->readonly = (!PyArray_ISWRITEABLE(((PyArrayObject *)__pyx_v_self)));

  
  __pyx_v_f = NULL;

  
  __Pyx_INCREF(((PyObject *)((PyArrayObject *)__pyx_v_self)->descr));
  __pyx_v_descr = ((PyArrayObject *)__pyx_v_self)->descr;

  
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
    goto __pyx_L11;
  }
   {

    
    __Pyx_INCREF(__pyx_v_self);
    __Pyx_GIVEREF(__pyx_v_self);
    __Pyx_GOTREF(__pyx_v_info->obj);
    __Pyx_DECREF(__pyx_v_info->obj);
    __pyx_v_info->obj = __pyx_v_self;
  }
  __pyx_L11:;

  
  __pyx_t_1 = (!__pyx_v_hasfields);
  if (__pyx_t_1) {

    
    __pyx_v_t = __pyx_v_descr->type_num;

    
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
        __pyx_t_8 = __pyx_t_3;
      } else {
        __pyx_t_8 = __pyx_t_1;
      }
      __pyx_t_1 = __pyx_t_8;
    } else {
      __pyx_t_1 = __pyx_t_2;
    }
    if (__pyx_t_1) {

      
      __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 248; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_4);
      __Pyx_INCREF(((PyObject *)__pyx_kp_u_15));
      PyTuple_SET_ITEM(__pyx_t_4, 0, ((PyObject *)__pyx_kp_u_15));
      __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_15));
      __pyx_t_5 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_4, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 248; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
      __Pyx_Raise(__pyx_t_5, 0, 0);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      {__pyx_filename = __pyx_f[1]; __pyx_lineno = 248; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      goto __pyx_L13;
    }
    __pyx_L13:;

    
    __pyx_t_1 = (__pyx_v_t == NPY_BYTE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__b;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_UBYTE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__B;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_SHORT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__h;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_USHORT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__H;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_INT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__i;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_UINT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__I;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_LONG);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__l;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_ULONG);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__L;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_LONGLONG);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__q;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_ULONGLONG);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__Q;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_FLOAT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__f;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_DOUBLE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__d;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_LONGDOUBLE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__g;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_CFLOAT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__Zf;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_CDOUBLE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__Zd;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_CLONGDOUBLE);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__Zg;
      goto __pyx_L14;
    }

    
    __pyx_t_1 = (__pyx_v_t == NPY_OBJECT);
    if (__pyx_t_1) {
      __pyx_v_f = __pyx_k__O;
      goto __pyx_L14;
    }
     {

      
      __pyx_t_5 = PyInt_FromLong(__pyx_v_t); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 267; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_4 = PyNumber_Remainder(((PyObject *)__pyx_kp_u_16), __pyx_t_5); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 267; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(((PyObject *)__pyx_t_4));
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_5 = PyTuple_New(1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 267; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      PyTuple_SET_ITEM(__pyx_t_5, 0, ((PyObject *)__pyx_t_4));
      __Pyx_GIVEREF(((PyObject *)__pyx_t_4));
      __pyx_t_4 = 0;
      __pyx_t_4 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_5, NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 267; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_4);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __Pyx_Raise(__pyx_t_4, 0, 0);
      __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
      {__pyx_filename = __pyx_f[1]; __pyx_lineno = 267; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    }
    __pyx_L14:;

    
    __pyx_v_info->format = __pyx_v_f;

    
    __pyx_r = 0;
    goto __pyx_L0;
    goto __pyx_L12;
  }
   {

    
    __pyx_v_info->format = ((char *)malloc(255));

    
    (__pyx_v_info->format[0]) = '^';

    
    __pyx_v_offset = 0;

    
    __pyx_t_9 = __pyx_f_5numpy__util_dtypestring(__pyx_v_descr, (__pyx_v_info->format + 1), (__pyx_v_info->format + 255), (&__pyx_v_offset)); if (unlikely(__pyx_t_9 == NULL)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 274; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __pyx_v_f = __pyx_t_9;

    
    (__pyx_v_f[0]) = 0;
  }
  __pyx_L12:;

  __pyx_r = 0;
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_AddTraceback("numpy.ndarray.__getbuffer__");
  __pyx_r = -1;
  __Pyx_GOTREF(__pyx_v_info->obj);
  __Pyx_DECREF(__pyx_v_info->obj); __pyx_v_info->obj = NULL;
  goto __pyx_L2;
  __pyx_L0:;
  if (__pyx_v_info->obj == Py_None) {
    __Pyx_GOTREF(Py_None);
    __Pyx_DECREF(Py_None); __pyx_v_info->obj = NULL;
  }
  __pyx_L2:;
  __Pyx_XDECREF((PyObject *)__pyx_v_descr);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_UNUSED void __pyx_pf_5numpy_7ndarray___releasebuffer__(PyObject *__pyx_v_self, Py_buffer *__pyx_v_info); 
static CYTHON_UNUSED void __pyx_pf_5numpy_7ndarray___releasebuffer__(PyObject *__pyx_v_self, Py_buffer *__pyx_v_info) {
  int __pyx_t_1;
  __Pyx_RefNannySetupContext("__releasebuffer__");

  
  __pyx_t_1 = PyArray_HASFIELDS(((PyArrayObject *)__pyx_v_self));
  if (__pyx_t_1) {

    
    free(__pyx_v_info->format);
    goto __pyx_L5;
  }
  __pyx_L5:;

  
  __pyx_t_1 = ((sizeof(npy_intp)) != (sizeof(Py_ssize_t)));
  if (__pyx_t_1) {

    
    free(__pyx_v_info->strides);
    goto __pyx_L6;
  }
  __pyx_L6:;

  __Pyx_RefNannyFinishContext();
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew1(PyObject *__pyx_v_a) {
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew1");

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(1, ((void *)__pyx_v_a)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 757; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew1");
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew2(PyObject *__pyx_v_a, PyObject *__pyx_v_b) {
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew2");

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(2, ((void *)__pyx_v_a), ((void *)__pyx_v_b)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 760; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew2");
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew3(PyObject *__pyx_v_a, PyObject *__pyx_v_b, PyObject *__pyx_v_c) {
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew3");

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(3, ((void *)__pyx_v_a), ((void *)__pyx_v_b), ((void *)__pyx_v_c)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 763; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew3");
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew4(PyObject *__pyx_v_a, PyObject *__pyx_v_b, PyObject *__pyx_v_c, PyObject *__pyx_v_d) {
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew4");

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(4, ((void *)__pyx_v_a), ((void *)__pyx_v_b), ((void *)__pyx_v_c), ((void *)__pyx_v_d)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 766; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew4");
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE PyObject *__pyx_f_5numpy_PyArray_MultiIterNew5(PyObject *__pyx_v_a, PyObject *__pyx_v_b, PyObject *__pyx_v_c, PyObject *__pyx_v_d, PyObject *__pyx_v_e) {
  PyObject *__pyx_r = NULL;
  PyObject *__pyx_t_1 = NULL;
  __Pyx_RefNannySetupContext("PyArray_MultiIterNew5");

  
  __Pyx_XDECREF(__pyx_r);
  __pyx_t_1 = PyArray_MultiIterNew(5, ((void *)__pyx_v_a), ((void *)__pyx_v_b), ((void *)__pyx_v_c), ((void *)__pyx_v_d), ((void *)__pyx_v_e)); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 769; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_r = __pyx_t_1;
  __pyx_t_1 = 0;
  goto __pyx_L0;

  __pyx_r = Py_None; __Pyx_INCREF(Py_None);
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_AddTraceback("numpy.PyArray_MultiIterNew5");
  __pyx_r = 0;
  __pyx_L0:;
  __Pyx_XGIVEREF(__pyx_r);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE char *__pyx_f_5numpy__util_dtypestring(PyArray_Descr *__pyx_v_descr, char *__pyx_v_f, char *__pyx_v_end, int *__pyx_v_offset) {
  PyArray_Descr *__pyx_v_child;
  int __pyx_v_endian_detector;
  int __pyx_v_little_endian;
  PyObject *__pyx_v_fields;
  PyObject *__pyx_v_childname;
  PyObject *__pyx_v_new_offset;
  PyObject *__pyx_v_t;
  char *__pyx_r;
  Py_ssize_t __pyx_t_1;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  PyObject *__pyx_t_4 = NULL;
  PyObject *__pyx_t_5 = NULL;
  int __pyx_t_6;
  int __pyx_t_7;
  int __pyx_t_8;
  int __pyx_t_9;
  char *__pyx_t_10;
  __Pyx_RefNannySetupContext("_util_dtypestring");
  __pyx_v_child = ((PyArray_Descr *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_fields = ((PyObject *)Py_None); __Pyx_INCREF(Py_None);
  __pyx_v_childname = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_new_offset = Py_None; __Pyx_INCREF(Py_None);
  __pyx_v_t = Py_None; __Pyx_INCREF(Py_None);

  
  __pyx_v_endian_detector = 1;

  
  __pyx_v_little_endian = ((((char *)(&__pyx_v_endian_detector))[0]) != 0);

  
  if (likely(((PyObject *)__pyx_v_descr->names) != Py_None)) {
    __pyx_t_1 = 0; __pyx_t_2 = ((PyObject *)__pyx_v_descr->names); __Pyx_INCREF(__pyx_t_2);
  } else {
    PyErr_SetString(PyExc_TypeError, "'NoneType' object is not iterable"); {__pyx_filename = __pyx_f[1]; __pyx_lineno = 782; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  }
  for (;;) {
    if (__pyx_t_1 >= PyTuple_GET_SIZE(__pyx_t_2)) break;
    __pyx_t_3 = PyTuple_GET_ITEM(__pyx_t_2, __pyx_t_1); __Pyx_INCREF(__pyx_t_3); __pyx_t_1++;
    __Pyx_DECREF(__pyx_v_childname);
    __pyx_v_childname = __pyx_t_3;
    __pyx_t_3 = 0;

    
    __pyx_t_3 = PyObject_GetItem(__pyx_v_descr->fields, __pyx_v_childname); if (!__pyx_t_3) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 783; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    if (!(likely(PyTuple_CheckExact(__pyx_t_3))||((__pyx_t_3) == Py_None)||(PyErr_Format(PyExc_TypeError, "Expected tuple, got %.200s", Py_TYPE(__pyx_t_3)->tp_name), 0))) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 783; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(((PyObject *)__pyx_v_fields));
    __pyx_v_fields = ((PyObject *)__pyx_t_3);
    __pyx_t_3 = 0;

    
    if (likely(((PyObject *)__pyx_v_fields) != Py_None) && likely(PyTuple_GET_SIZE(((PyObject *)__pyx_v_fields)) == 2)) {
      PyObject* tuple = ((PyObject *)__pyx_v_fields);
      __pyx_t_3 = PyTuple_GET_ITEM(tuple, 0); __Pyx_INCREF(__pyx_t_3);
      if (!(likely(((__pyx_t_3) == Py_None) || likely(__Pyx_TypeTest(__pyx_t_3, __pyx_ptype_5numpy_dtype))))) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 784; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __pyx_t_4 = PyTuple_GET_ITEM(tuple, 1); __Pyx_INCREF(__pyx_t_4);
      __Pyx_DECREF(((PyObject *)__pyx_v_child));
      __pyx_v_child = ((PyArray_Descr *)__pyx_t_3);
      __pyx_t_3 = 0;
      __Pyx_DECREF(__pyx_v_new_offset);
      __pyx_v_new_offset = __pyx_t_4;
      __pyx_t_4 = 0;
    } else {
      __Pyx_UnpackTupleError(((PyObject *)__pyx_v_fields), 2);
      {__pyx_filename = __pyx_f[1]; __pyx_lineno = 784; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    }

    
    __pyx_t_4 = PyInt_FromLong((__pyx_v_end - __pyx_v_f)); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 786; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_4);
    __pyx_t_3 = PyInt_FromLong((__pyx_v_offset[0])); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 786; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __pyx_t_5 = PyNumber_Subtract(__pyx_v_new_offset, __pyx_t_3); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 786; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_5);
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    __pyx_t_3 = PyNumber_Subtract(__pyx_t_4, __pyx_t_5); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 786; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_3);
    __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
    __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
    __pyx_t_5 = PyObject_RichCompare(__pyx_t_3, __pyx_int_15, Py_LT); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 786; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_GOTREF(__pyx_t_5);
    __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
    __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 786; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
    __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
    if (__pyx_t_6) {

      
      __pyx_t_5 = PyTuple_New(1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 787; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_INCREF(((PyObject *)__pyx_kp_u_17));
      PyTuple_SET_ITEM(__pyx_t_5, 0, ((PyObject *)__pyx_kp_u_17));
      __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_17));
      __pyx_t_3 = PyObject_Call(__pyx_builtin_RuntimeError, __pyx_t_5, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 787; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __Pyx_Raise(__pyx_t_3, 0, 0);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      {__pyx_filename = __pyx_f[1]; __pyx_lineno = 787; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      goto __pyx_L5;
    }
    __pyx_L5:;

    
    __pyx_t_6 = (__pyx_v_child->byteorder == '>');
    if (__pyx_t_6) {
      __pyx_t_7 = __pyx_v_little_endian;
    } else {
      __pyx_t_7 = __pyx_t_6;
    }
    if (!__pyx_t_7) {

      
      __pyx_t_6 = (__pyx_v_child->byteorder == '<');
      if (__pyx_t_6) {
        __pyx_t_8 = (!__pyx_v_little_endian);
        __pyx_t_9 = __pyx_t_8;
      } else {
        __pyx_t_9 = __pyx_t_6;
      }
      __pyx_t_6 = __pyx_t_9;
    } else {
      __pyx_t_6 = __pyx_t_7;
    }
    if (__pyx_t_6) {

      
      __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 791; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_INCREF(((PyObject *)__pyx_kp_u_15));
      PyTuple_SET_ITEM(__pyx_t_3, 0, ((PyObject *)__pyx_kp_u_15));
      __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_15));
      __pyx_t_5 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_3, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 791; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __Pyx_Raise(__pyx_t_5, 0, 0);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      {__pyx_filename = __pyx_f[1]; __pyx_lineno = 791; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      goto __pyx_L6;
    }
    __pyx_L6:;

    
    while (1) {
      __pyx_t_5 = PyInt_FromLong((__pyx_v_offset[0])); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 801; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_t_5, __pyx_v_new_offset, Py_LT); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 801; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 801; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (!__pyx_t_6) break;

      
      (__pyx_v_f[0]) = 120;

      
      __pyx_v_f += 1;

      
      (__pyx_v_offset[0]) += 1;
    }

    
    (__pyx_v_offset[0]) += __pyx_v_child->elsize;

    
    __pyx_t_6 = (!PyDataType_HASFIELDS(__pyx_v_child));
    if (__pyx_t_6) {

      
      __pyx_t_3 = PyInt_FromLong(__pyx_v_child->type_num); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 809; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_v_t);
      __pyx_v_t = __pyx_t_3;
      __pyx_t_3 = 0;

      
      __pyx_t_6 = ((__pyx_v_end - __pyx_v_f) < 5);
      if (__pyx_t_6) {

        
        __pyx_t_3 = PyTuple_New(1); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 811; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        __Pyx_GOTREF(__pyx_t_3);
        __Pyx_INCREF(((PyObject *)__pyx_kp_u_18));
        PyTuple_SET_ITEM(__pyx_t_3, 0, ((PyObject *)__pyx_kp_u_18));
        __Pyx_GIVEREF(((PyObject *)__pyx_kp_u_18));
        __pyx_t_5 = PyObject_Call(__pyx_builtin_RuntimeError, __pyx_t_3, NULL); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 811; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        __Pyx_GOTREF(__pyx_t_5);
        __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
        __Pyx_Raise(__pyx_t_5, 0, 0);
        __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
        {__pyx_filename = __pyx_f[1]; __pyx_lineno = 811; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        goto __pyx_L10;
      }
      __pyx_L10:;

      
      __pyx_t_5 = PyInt_FromLong(NPY_BYTE); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 814; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 814; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 814; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 98;
        goto __pyx_L11;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_UBYTE); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 815; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 815; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 815; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 66;
        goto __pyx_L11;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_SHORT); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 816; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 816; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 816; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 104;
        goto __pyx_L11;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_USHORT); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 817; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 817; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 817; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 72;
        goto __pyx_L11;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_INT); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 818; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 818; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 818; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 105;
        goto __pyx_L11;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_UINT); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 819; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 819; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 819; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 73;
        goto __pyx_L11;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_LONG); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 820; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 820; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 820; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 108;
        goto __pyx_L11;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_ULONG); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 821; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 821; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 821; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 76;
        goto __pyx_L11;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_LONGLONG); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 822; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 822; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 822; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 113;
        goto __pyx_L11;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_ULONGLONG); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 823; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 823; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 823; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 81;
        goto __pyx_L11;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_FLOAT); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 824; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 824; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 824; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 102;
        goto __pyx_L11;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_DOUBLE); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 825; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 825; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 825; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 100;
        goto __pyx_L11;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_LONGDOUBLE); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 826; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 826; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 826; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 103;
        goto __pyx_L11;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_CFLOAT); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 827; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 827; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 827; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 90;
        (__pyx_v_f[1]) = 102;
        __pyx_v_f += 1;
        goto __pyx_L11;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_CDOUBLE); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 828; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 828; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 828; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 90;
        (__pyx_v_f[1]) = 100;
        __pyx_v_f += 1;
        goto __pyx_L11;
      }

      
      __pyx_t_3 = PyInt_FromLong(NPY_CLONGDOUBLE); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 829; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __pyx_t_5 = PyObject_RichCompare(__pyx_v_t, __pyx_t_3, Py_EQ); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 829; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_5); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 829; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 90;
        (__pyx_v_f[1]) = 103;
        __pyx_v_f += 1;
        goto __pyx_L11;
      }

      
      __pyx_t_5 = PyInt_FromLong(NPY_OBJECT); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 830; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_5);
      __pyx_t_3 = PyObject_RichCompare(__pyx_v_t, __pyx_t_5, Py_EQ); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 830; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_GOTREF(__pyx_t_3);
      __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
      __pyx_t_6 = __Pyx_PyObject_IsTrue(__pyx_t_3); if (unlikely(__pyx_t_6 < 0)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 830; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
      if (__pyx_t_6) {
        (__pyx_v_f[0]) = 79;
        goto __pyx_L11;
      }
       {

        
        __pyx_t_3 = PyNumber_Remainder(((PyObject *)__pyx_kp_u_16), __pyx_v_t); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 832; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        __Pyx_GOTREF(((PyObject *)__pyx_t_3));
        __pyx_t_5 = PyTuple_New(1); if (unlikely(!__pyx_t_5)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 832; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        __Pyx_GOTREF(__pyx_t_5);
        PyTuple_SET_ITEM(__pyx_t_5, 0, ((PyObject *)__pyx_t_3));
        __Pyx_GIVEREF(((PyObject *)__pyx_t_3));
        __pyx_t_3 = 0;
        __pyx_t_3 = PyObject_Call(__pyx_builtin_ValueError, __pyx_t_5, NULL); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 832; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
        __Pyx_GOTREF(__pyx_t_3);
        __Pyx_DECREF(__pyx_t_5); __pyx_t_5 = 0;
        __Pyx_Raise(__pyx_t_3, 0, 0);
        __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
        {__pyx_filename = __pyx_f[1]; __pyx_lineno = 832; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      }
      __pyx_L11:;

      
      __pyx_v_f += 1;
      goto __pyx_L9;
    }
     {

      
      __pyx_t_10 = __pyx_f_5numpy__util_dtypestring(__pyx_v_child, __pyx_v_f, __pyx_v_end, __pyx_v_offset); if (unlikely(__pyx_t_10 == NULL)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 837; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
      __pyx_v_f = __pyx_t_10;
    }
    __pyx_L9:;
  }
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;

  
  __pyx_r = __pyx_v_f;
  goto __pyx_L0;

  __pyx_r = 0;
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_4);
  __Pyx_XDECREF(__pyx_t_5);
  __Pyx_AddTraceback("numpy._util_dtypestring");
  __pyx_r = NULL;
  __pyx_L0:;
  __Pyx_DECREF((PyObject *)__pyx_v_child);
  __Pyx_DECREF(__pyx_v_fields);
  __Pyx_DECREF(__pyx_v_childname);
  __Pyx_DECREF(__pyx_v_new_offset);
  __Pyx_DECREF(__pyx_v_t);
  __Pyx_RefNannyFinishContext();
  return __pyx_r;
}



static CYTHON_INLINE void __pyx_f_5numpy_set_array_base(PyArrayObject *__pyx_v_arr, PyObject *__pyx_v_base) {
  PyObject *__pyx_v_baseptr;
  int __pyx_t_1;
  __Pyx_RefNannySetupContext("set_array_base");

  
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
  int __pyx_t_1;
  __Pyx_RefNannySetupContext("get_array_base");

  
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

static PyObject *__pyx_tp_new_5scipy_7spatial_5qhull_RidgeIter2D(PyTypeObject *t, PyObject *a, PyObject *k) {
  struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *p;
  PyObject *o = (*t->tp_alloc)(t, 0);
  if (!o) return 0;
  p = ((struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)o);
  p->delaunay = Py_None; Py_INCREF(Py_None);
  return o;
}

static void __pyx_tp_dealloc_5scipy_7spatial_5qhull_RidgeIter2D(PyObject *o) {
  struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *p = (struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)o;
  Py_XDECREF(p->delaunay);
  (*Py_TYPE(o)->tp_free)(o);
}

static int __pyx_tp_traverse_5scipy_7spatial_5qhull_RidgeIter2D(PyObject *o, visitproc v, void *a) {
  int e;
  struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *p = (struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)o;
  if (p->delaunay) {
    e = (*v)(p->delaunay, a); if (e) return e;
  }
  return 0;
}

static int __pyx_tp_clear_5scipy_7spatial_5qhull_RidgeIter2D(PyObject *o) {
  struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *p = (struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D *)o;
  PyObject* tmp;
  tmp = ((PyObject*)p->delaunay);
  p->delaunay = Py_None; Py_INCREF(Py_None);
  Py_XDECREF(tmp);
  return 0;
}

static PyMethodDef __pyx_methods_5scipy_7spatial_5qhull_RidgeIter2D[] = {
  {__Pyx_NAMESTR("__next__"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___next__, METH_NOARGS|METH_COEXIST, __Pyx_DOCSTR(0)},
  {0, 0, 0, 0}
};

static PyNumberMethods __pyx_tp_as_number_RidgeIter2D = {
  0, 
  0, 
  0, 
  #if PY_MAJOR_VERSION < 3
  0, 
  #endif
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  #if PY_MAJOR_VERSION < 3
  0, 
  #endif
  0, 
  #if PY_MAJOR_VERSION < 3
  0, 
  #else
  0, 
  #endif
  0, 
  #if PY_MAJOR_VERSION < 3
  0, 
  #endif
  #if PY_MAJOR_VERSION < 3
  0, 
  #endif
  0, 
  0, 
  0, 
  #if PY_MAJOR_VERSION < 3
  0, 
  #endif
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  #if PY_VERSION_HEX >= 0x02050000
  0, 
  #endif
};

static PySequenceMethods __pyx_tp_as_sequence_RidgeIter2D = {
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
};

static PyMappingMethods __pyx_tp_as_mapping_RidgeIter2D = {
  0, 
  0, 
  0, 
};

static PyBufferProcs __pyx_tp_as_buffer_RidgeIter2D = {
  #if PY_MAJOR_VERSION < 3
  0, 
  #endif
  #if PY_MAJOR_VERSION < 3
  0, 
  #endif
  #if PY_MAJOR_VERSION < 3
  0, 
  #endif
  #if PY_MAJOR_VERSION < 3
  0, 
  #endif
  #if PY_VERSION_HEX >= 0x02060000
  0, 
  #endif
  #if PY_VERSION_HEX >= 0x02060000
  0, 
  #endif
};

PyTypeObject __pyx_type_5scipy_7spatial_5qhull_RidgeIter2D = {
  PyVarObject_HEAD_INIT(0, 0)
  __Pyx_NAMESTR("scipy.spatial.qhull.RidgeIter2D"), 
  sizeof(struct __pyx_obj_5scipy_7spatial_5qhull_RidgeIter2D), 
  0, 
  __pyx_tp_dealloc_5scipy_7spatial_5qhull_RidgeIter2D, 
  0, 
  0, 
  0, 
  #if PY_MAJOR_VERSION < 3
  0, 
  #else
  0, 
  #endif
  0, 
  &__pyx_tp_as_number_RidgeIter2D, 
  &__pyx_tp_as_sequence_RidgeIter2D, 
  &__pyx_tp_as_mapping_RidgeIter2D, 
  0, 
  0, 
  0, 
  0, 
  0, 
  &__pyx_tp_as_buffer_RidgeIter2D, 
  Py_TPFLAGS_DEFAULT|Py_TPFLAGS_CHECKTYPES|Py_TPFLAGS_BASETYPE|Py_TPFLAGS_HAVE_NEWBUFFER|Py_TPFLAGS_HAVE_GC, 
  0, 
  __pyx_tp_traverse_5scipy_7spatial_5qhull_RidgeIter2D, 
  __pyx_tp_clear_5scipy_7spatial_5qhull_RidgeIter2D, 
  0, 
  0, 
  __pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___iter__, 
  __pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___next__, 
  __pyx_methods_5scipy_7spatial_5qhull_RidgeIter2D, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  __pyx_pf_5scipy_7spatial_5qhull_11RidgeIter2D___init__, 
  0, 
  __pyx_tp_new_5scipy_7spatial_5qhull_RidgeIter2D, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  0, 
  #if PY_VERSION_HEX >= 0x02060000
  0, 
  #endif
};

static PyMethodDef __pyx_methods[] = {
  {__Pyx_NAMESTR("_construct_delaunay"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull__construct_delaunay, METH_O, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull__construct_delaunay)},
  {__Pyx_NAMESTR("_qhull_get_facet_array"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull__qhull_get_facet_array, METH_VARARGS|METH_KEYWORDS, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull__qhull_get_facet_array)},
  {__Pyx_NAMESTR("_get_barycentric_transforms"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull__get_barycentric_transforms, METH_VARARGS|METH_KEYWORDS, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull__get_barycentric_transforms)},
  {__Pyx_NAMESTR("tsearch"), (PyCFunction)__pyx_pf_5scipy_7spatial_5qhull_tsearch, METH_VARARGS|METH_KEYWORDS, __Pyx_DOCSTR(__pyx_doc_5scipy_7spatial_5qhull_tsearch)},
  {0, 0, 0, 0}
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef __pyx_moduledef = {
    PyModuleDef_HEAD_INIT,
    __Pyx_NAMESTR("qhull"),
    __Pyx_DOCSTR(__pyx_k_19), 
    -1, 
    __pyx_methods ,
    NULL, 
    NULL, 
    NULL, 
    NULL 
};
#endif

static __Pyx_StringTabEntry __pyx_string_tab[] = {
  {&__pyx_kp_s_11, __pyx_k_11, sizeof(__pyx_k_11), 0, 0, 1, 0},
  {&__pyx_kp_s_12, __pyx_k_12, sizeof(__pyx_k_12), 0, 0, 1, 0},
  {&__pyx_kp_u_13, __pyx_k_13, sizeof(__pyx_k_13), 0, 1, 0, 0},
  {&__pyx_kp_u_14, __pyx_k_14, sizeof(__pyx_k_14), 0, 1, 0, 0},
  {&__pyx_kp_u_15, __pyx_k_15, sizeof(__pyx_k_15), 0, 1, 0, 0},
  {&__pyx_kp_u_16, __pyx_k_16, sizeof(__pyx_k_16), 0, 1, 0, 0},
  {&__pyx_kp_u_17, __pyx_k_17, sizeof(__pyx_k_17), 0, 1, 0, 0},
  {&__pyx_kp_u_18, __pyx_k_18, sizeof(__pyx_k_18), 0, 1, 0, 0},
  {&__pyx_kp_s_2, __pyx_k_2, sizeof(__pyx_k_2), 0, 0, 1, 0},
  {&__pyx_kp_s_20, __pyx_k_20, sizeof(__pyx_k_20), 0, 0, 1, 0},
  {&__pyx_kp_u_21, __pyx_k_21, sizeof(__pyx_k_21), 0, 1, 0, 0},
  {&__pyx_kp_u_22, __pyx_k_22, sizeof(__pyx_k_22), 0, 1, 0, 0},
  {&__pyx_kp_u_23, __pyx_k_23, sizeof(__pyx_k_23), 0, 1, 0, 0},
  {&__pyx_kp_u_24, __pyx_k_24, sizeof(__pyx_k_24), 0, 1, 0, 0},
  {&__pyx_kp_u_25, __pyx_k_25, sizeof(__pyx_k_25), 0, 1, 0, 0},
  {&__pyx_kp_u_26, __pyx_k_26, sizeof(__pyx_k_26), 0, 1, 0, 0},
  {&__pyx_kp_u_27, __pyx_k_27, sizeof(__pyx_k_27), 0, 1, 0, 0},
  {&__pyx_kp_u_28, __pyx_k_28, sizeof(__pyx_k_28), 0, 1, 0, 0},
  {&__pyx_kp_u_29, __pyx_k_29, sizeof(__pyx_k_29), 0, 1, 0, 0},
  {&__pyx_kp_s_3, __pyx_k_3, sizeof(__pyx_k_3), 0, 0, 1, 0},
  {&__pyx_kp_u_30, __pyx_k_30, sizeof(__pyx_k_30), 0, 1, 0, 0},
  {&__pyx_kp_s_4, __pyx_k_4, sizeof(__pyx_k_4), 0, 0, 1, 0},
  {&__pyx_n_s_5, __pyx_k_5, sizeof(__pyx_k_5), 0, 0, 1, 1},
  {&__pyx_kp_s_6, __pyx_k_6, sizeof(__pyx_k_6), 0, 0, 1, 0},
  {&__pyx_kp_s_7, __pyx_k_7, sizeof(__pyx_k_7), 0, 0, 1, 0},
  {&__pyx_kp_s_8, __pyx_k_8, sizeof(__pyx_k_8), 0, 0, 1, 0},
  {&__pyx_n_s_9, __pyx_k_9, sizeof(__pyx_k_9), 0, 0, 1, 1},
  {&__pyx_n_s__Delaunay, __pyx_k__Delaunay, sizeof(__pyx_k__Delaunay), 0, 0, 1, 1},
  {&__pyx_n_s__Lock, __pyx_k__Lock, sizeof(__pyx_k__Lock), 0, 0, 1, 1},
  {&__pyx_n_s__NOerrexit, __pyx_k__NOerrexit, sizeof(__pyx_k__NOerrexit), 0, 0, 1, 1},
  {&__pyx_n_s__RuntimeError, __pyx_k__RuntimeError, sizeof(__pyx_k__RuntimeError), 0, 0, 1, 1},
  {&__pyx_n_s__SCALElast, __pyx_k__SCALElast, sizeof(__pyx_k__SCALElast), 0, 0, 1, 1},
  {&__pyx_n_s__StopIteration, __pyx_k__StopIteration, sizeof(__pyx_k__StopIteration), 0, 0, 1, 1},
  {&__pyx_n_s__ValueError, __pyx_k__ValueError, sizeof(__pyx_k__ValueError), 0, 0, 1, 1},
  {&__pyx_n_s____all__, __pyx_k____all__, sizeof(__pyx_k____all__), 0, 0, 1, 1},
  {&__pyx_n_s____init__, __pyx_k____init__, sizeof(__pyx_k____init__), 0, 0, 1, 1},
  {&__pyx_n_s____main__, __pyx_k____main__, sizeof(__pyx_k____main__), 0, 0, 1, 1},
  {&__pyx_n_s____test__, __pyx_k____test__, sizeof(__pyx_k____test__), 0, 0, 1, 1},
  {&__pyx_n_s___construct_delaunay, __pyx_k___construct_delaunay, sizeof(__pyx_k___construct_delaunay), 0, 0, 1, 1},
  {&__pyx_n_s___qhull_lock, __pyx_k___qhull_lock, sizeof(__pyx_k___qhull_lock), 0, 0, 1, 1},
  {&__pyx_n_s___transform, __pyx_k___transform, sizeof(__pyx_k___transform), 0, 0, 1, 1},
  {&__pyx_n_s___vertex_to_simplex, __pyx_k___vertex_to_simplex, sizeof(__pyx_k___vertex_to_simplex), 0, 0, 1, 1},
  {&__pyx_n_s__acquire, __pyx_k__acquire, sizeof(__pyx_k__acquire), 0, 0, 1, 1},
  {&__pyx_n_s__asanyarray, __pyx_k__asanyarray, sizeof(__pyx_k__asanyarray), 0, 0, 1, 1},
  {&__pyx_n_s__ascontiguousarray, __pyx_k__ascontiguousarray, sizeof(__pyx_k__ascontiguousarray), 0, 0, 1, 1},
  {&__pyx_n_s__astype, __pyx_k__astype, sizeof(__pyx_k__astype), 0, 0, 1, 1},
  {&__pyx_n_s__axis, __pyx_k__axis, sizeof(__pyx_k__axis), 0, 0, 1, 1},
  {&__pyx_n_s__base, __pyx_k__base, sizeof(__pyx_k__base), 0, 0, 1, 1},
  {&__pyx_n_s__bruteforce, __pyx_k__bruteforce, sizeof(__pyx_k__bruteforce), 0, 0, 1, 1},
  {&__pyx_n_s__buf, __pyx_k__buf, sizeof(__pyx_k__buf), 0, 0, 1, 1},
  {&__pyx_n_s__byteorder, __pyx_k__byteorder, sizeof(__pyx_k__byteorder), 0, 0, 1, 1},
  {&__pyx_n_s__convex_hull, __pyx_k__convex_hull, sizeof(__pyx_k__convex_hull), 0, 0, 1, 1},
  {&__pyx_n_s__data, __pyx_k__data, sizeof(__pyx_k__data), 0, 0, 1, 1},
  {&__pyx_n_s__delaunay, __pyx_k__delaunay, sizeof(__pyx_k__delaunay), 0, 0, 1, 1},
  {&__pyx_n_s__descr, __pyx_k__descr, sizeof(__pyx_k__descr), 0, 0, 1, 1},
  {&__pyx_n_s__double, __pyx_k__double, sizeof(__pyx_k__double), 0, 0, 1, 1},
  {&__pyx_n_s__dtype, __pyx_k__dtype, sizeof(__pyx_k__dtype), 0, 0, 1, 1},
  {&__pyx_n_s__e, __pyx_k__e, sizeof(__pyx_k__e), 0, 0, 1, 1},
  {&__pyx_n_s__empty, __pyx_k__empty, sizeof(__pyx_k__empty), 0, 0, 1, 1},
  {&__pyx_n_s__eps, __pyx_k__eps, sizeof(__pyx_k__eps), 0, 0, 1, 1},
  {&__pyx_n_s__equations, __pyx_k__equations, sizeof(__pyx_k__equations), 0, 0, 1, 1},
  {&__pyx_n_s__facet_id, __pyx_k__facet_id, sizeof(__pyx_k__facet_id), 0, 0, 1, 1},
  {&__pyx_n_s__facet_list, __pyx_k__facet_list, sizeof(__pyx_k__facet_list), 0, 0, 1, 1},
  {&__pyx_n_s__fields, __pyx_k__fields, sizeof(__pyx_k__fields), 0, 0, 1, 1},
  {&__pyx_n_s__fill, __pyx_k__fill, sizeof(__pyx_k__fill), 0, 0, 1, 1},
  {&__pyx_n_s__find_simplex, __pyx_k__find_simplex, sizeof(__pyx_k__find_simplex), 0, 0, 1, 1},
  {&__pyx_n_s__finfo, __pyx_k__finfo, sizeof(__pyx_k__finfo), 0, 0, 1, 1},
  {&__pyx_n_s__format, __pyx_k__format, sizeof(__pyx_k__format), 0, 0, 1, 1},
  {&__pyx_n_s__id, __pyx_k__id, sizeof(__pyx_k__id), 0, 0, 1, 1},
  {&__pyx_n_s__index, __pyx_k__index, sizeof(__pyx_k__index), 0, 0, 1, 1},
  {&__pyx_n_s__info, __pyx_k__info, sizeof(__pyx_k__info), 0, 0, 1, 1},
  {&__pyx_n_s__intc, __pyx_k__intc, sizeof(__pyx_k__intc), 0, 0, 1, 1},
  {&__pyx_n_s__it, __pyx_k__it, sizeof(__pyx_k__it), 0, 0, 1, 1},
  {&__pyx_n_s__itemsize, __pyx_k__itemsize, sizeof(__pyx_k__itemsize), 0, 0, 1, 1},
  {&__pyx_n_s__ivertex, __pyx_k__ivertex, sizeof(__pyx_k__ivertex), 0, 0, 1, 1},
  {&__pyx_n_s__last_high, __pyx_k__last_high, sizeof(__pyx_k__last_high), 0, 0, 1, 1},
  {&__pyx_n_s__last_low, __pyx_k__last_low, sizeof(__pyx_k__last_low), 0, 0, 1, 1},
  {&__pyx_n_s__last_newhigh, __pyx_k__last_newhigh, sizeof(__pyx_k__last_newhigh), 0, 0, 1, 1},
  {&__pyx_n_s__lift_points, __pyx_k__lift_points, sizeof(__pyx_k__lift_points), 0, 0, 1, 1},
  {&__pyx_n_s__max, __pyx_k__max, sizeof(__pyx_k__max), 0, 0, 1, 1},
  {&__pyx_n_s__max_bound, __pyx_k__max_bound, sizeof(__pyx_k__max_bound), 0, 0, 1, 1},
  {&__pyx_n_s__min, __pyx_k__min, sizeof(__pyx_k__min), 0, 0, 1, 1},
  {&__pyx_n_s__min_bound, __pyx_k__min_bound, sizeof(__pyx_k__min_bound), 0, 0, 1, 1},
  {&__pyx_n_s__names, __pyx_k__names, sizeof(__pyx_k__names), 0, 0, 1, 1},
  {&__pyx_n_s__nan, __pyx_k__nan, sizeof(__pyx_k__nan), 0, 0, 1, 1},
  {&__pyx_n_s__ndim, __pyx_k__ndim, sizeof(__pyx_k__ndim), 0, 0, 1, 1},
  {&__pyx_n_s__neighbors, __pyx_k__neighbors, sizeof(__pyx_k__neighbors), 0, 0, 1, 1},
  {&__pyx_n_s__next, __pyx_k__next, sizeof(__pyx_k__next), 0, 0, 1, 1},
  {&__pyx_n_s__normal, __pyx_k__normal, sizeof(__pyx_k__normal), 0, 0, 1, 1},
  {&__pyx_n_s__np, __pyx_k__np, sizeof(__pyx_k__np), 0, 0, 1, 1},
  {&__pyx_n_s__npoints, __pyx_k__npoints, sizeof(__pyx_k__npoints), 0, 0, 1, 1},
  {&__pyx_n_s__nsimplex, __pyx_k__nsimplex, sizeof(__pyx_k__nsimplex), 0, 0, 1, 1},
  {&__pyx_n_s__numpoints, __pyx_k__numpoints, sizeof(__pyx_k__numpoints), 0, 0, 1, 1},
  {&__pyx_n_s__numpy, __pyx_k__numpy, sizeof(__pyx_k__numpy), 0, 0, 1, 1},
  {&__pyx_n_s__obj, __pyx_k__obj, sizeof(__pyx_k__obj), 0, 0, 1, 1},
  {&__pyx_n_s__object, __pyx_k__object, sizeof(__pyx_k__object), 0, 0, 1, 1},
  {&__pyx_n_s__offset, __pyx_k__offset, sizeof(__pyx_k__offset), 0, 0, 1, 1},
  {&__pyx_n_s__p, __pyx_k__p, sizeof(__pyx_k__p), 0, 0, 1, 1},
  {&__pyx_n_s__paraboloid_scale, __pyx_k__paraboloid_scale, sizeof(__pyx_k__paraboloid_scale), 0, 0, 1, 1},
  {&__pyx_n_s__paraboloid_shift, __pyx_k__paraboloid_shift, sizeof(__pyx_k__paraboloid_shift), 0, 0, 1, 1},
  {&__pyx_n_s__plane_distance, __pyx_k__plane_distance, sizeof(__pyx_k__plane_distance), 0, 0, 1, 1},
  {&__pyx_n_s__point, __pyx_k__point, sizeof(__pyx_k__point), 0, 0, 1, 1},
  {&__pyx_n_s__points, __pyx_k__points, sizeof(__pyx_k__points), 0, 0, 1, 1},
  {&__pyx_n_s__prod, __pyx_k__prod, sizeof(__pyx_k__prod), 0, 0, 1, 1},
  {&__pyx_n_s__property, __pyx_k__property, sizeof(__pyx_k__property), 0, 0, 1, 1},
  {&__pyx_n_s__range, __pyx_k__range, sizeof(__pyx_k__range), 0, 0, 1, 1},
  {&__pyx_n_s__readonly, __pyx_k__readonly, sizeof(__pyx_k__readonly), 0, 0, 1, 1},
  {&__pyx_n_s__release, __pyx_k__release, sizeof(__pyx_k__release), 0, 0, 1, 1},
  {&__pyx_n_s__reshape, __pyx_k__reshape, sizeof(__pyx_k__reshape), 0, 0, 1, 1},
  {&__pyx_n_s__resize, __pyx_k__resize, sizeof(__pyx_k__resize), 0, 0, 1, 1},
  {&__pyx_n_s__restart, __pyx_k__restart, sizeof(__pyx_k__restart), 0, 0, 1, 1},
  {&__pyx_n_s__self, __pyx_k__self, sizeof(__pyx_k__self), 0, 0, 1, 1},
  {&__pyx_n_s__shape, __pyx_k__shape, sizeof(__pyx_k__shape), 0, 0, 1, 1},
  {&__pyx_n_s__simplicial, __pyx_k__simplicial, sizeof(__pyx_k__simplicial), 0, 0, 1, 1},
  {&__pyx_n_s__start_index, __pyx_k__start_index, sizeof(__pyx_k__start_index), 0, 0, 1, 1},
  {&__pyx_n_s__start_triangle, __pyx_k__start_triangle, sizeof(__pyx_k__start_triangle), 0, 0, 1, 1},
  {&__pyx_n_s__strides, __pyx_k__strides, sizeof(__pyx_k__strides), 0, 0, 1, 1},
  {&__pyx_n_s__suboffsets, __pyx_k__suboffsets, sizeof(__pyx_k__suboffsets), 0, 0, 1, 1},
  {&__pyx_n_s__sum, __pyx_k__sum, sizeof(__pyx_k__sum), 0, 0, 1, 1},
  {&__pyx_n_s__threading, __pyx_k__threading, sizeof(__pyx_k__threading), 0, 0, 1, 1},
  {&__pyx_n_s__transform, __pyx_k__transform, sizeof(__pyx_k__transform), 0, 0, 1, 1},
  {&__pyx_n_s__tri, __pyx_k__tri, sizeof(__pyx_k__tri), 0, 0, 1, 1},
  {&__pyx_n_s__triangle, __pyx_k__triangle, sizeof(__pyx_k__triangle), 0, 0, 1, 1},
  {&__pyx_n_s__tsearch, __pyx_k__tsearch, sizeof(__pyx_k__tsearch), 0, 0, 1, 1},
  {&__pyx_n_s__type_num, __pyx_k__type_num, sizeof(__pyx_k__type_num), 0, 0, 1, 1},
  {&__pyx_n_s__upperdelaunay, __pyx_k__upperdelaunay, sizeof(__pyx_k__upperdelaunay), 0, 0, 1, 1},
  {&__pyx_n_s__vertex, __pyx_k__vertex, sizeof(__pyx_k__vertex), 0, 0, 1, 1},
  {&__pyx_n_s__vertex2, __pyx_k__vertex2, sizeof(__pyx_k__vertex2), 0, 0, 1, 1},
  {&__pyx_n_s__vertex_to_simplex, __pyx_k__vertex_to_simplex, sizeof(__pyx_k__vertex_to_simplex), 0, 0, 1, 1},
  {&__pyx_n_s__vertices, __pyx_k__vertices, sizeof(__pyx_k__vertices), 0, 0, 1, 1},
  {&__pyx_n_s__x, __pyx_k__x, sizeof(__pyx_k__x), 0, 0, 1, 1},
  {&__pyx_n_s__xi, __pyx_k__xi, sizeof(__pyx_k__xi), 0, 0, 1, 1},
  {&__pyx_n_s__xrange, __pyx_k__xrange, sizeof(__pyx_k__xrange), 0, 0, 1, 1},
  {&__pyx_n_s__zeros, __pyx_k__zeros, sizeof(__pyx_k__zeros), 0, 0, 1, 1},
  {0, 0, 0, 0, 0, 0, 0}
};
static int __Pyx_InitCachedBuiltins(void) {
  __pyx_builtin_object = __Pyx_GetName(__pyx_b, __pyx_n_s__object); if (!__pyx_builtin_object) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 855; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_builtin_property = __Pyx_GetName(__pyx_b, __pyx_n_s__property); if (!__pyx_builtin_property) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 933; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_builtin_ValueError = __Pyx_GetName(__pyx_b, __pyx_n_s__ValueError); if (!__pyx_builtin_ValueError) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 158; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_builtin_RuntimeError = __Pyx_GetName(__pyx_b, __pyx_n_s__RuntimeError); if (!__pyx_builtin_RuntimeError) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 170; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #if PY_MAJOR_VERSION >= 3
  __pyx_builtin_xrange = __Pyx_GetName(__pyx_b, __pyx_n_s__range); if (!__pyx_builtin_xrange) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 252; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #else
  __pyx_builtin_xrange = __Pyx_GetName(__pyx_b, __pyx_n_s__xrange); if (!__pyx_builtin_xrange) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 252; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #endif
  __pyx_builtin_StopIteration = __Pyx_GetName(__pyx_b, __pyx_n_s__StopIteration); if (!__pyx_builtin_StopIteration) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 610; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_builtin_range = __Pyx_GetName(__pyx_b, __pyx_n_s__range); if (!__pyx_builtin_range) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 219; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  return 0;
  __pyx_L1_error:;
  return -1;
}

static int __Pyx_InitGlobals(void) {
  if (__Pyx_InitStrings(__pyx_string_tab) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  __pyx_int_0 = PyInt_FromLong(0); if (unlikely(!__pyx_int_0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  __pyx_int_1 = PyInt_FromLong(1); if (unlikely(!__pyx_int_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  __pyx_int_2 = PyInt_FromLong(2); if (unlikely(!__pyx_int_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  __pyx_int_neg_1 = PyInt_FromLong(-1); if (unlikely(!__pyx_int_neg_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  __pyx_int_10 = PyInt_FromLong(10); if (unlikely(!__pyx_int_10)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  __pyx_int_15 = PyInt_FromLong(15); if (unlikely(!__pyx_int_15)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  return 0;
  __pyx_L1_error:;
  return -1;
}

#if PY_MAJOR_VERSION < 3
PyMODINIT_FUNC initqhull(void); 
PyMODINIT_FUNC initqhull(void)
#else
PyMODINIT_FUNC PyInit_qhull(void); 
PyMODINIT_FUNC PyInit_qhull(void)
#endif
{
  PyObject *__pyx_t_1 = NULL;
  PyObject *__pyx_t_2 = NULL;
  PyObject *__pyx_t_3 = NULL;
  PyObject *__pyx_t_4 = NULL;
  #if CYTHON_REFNANNY
  void* __pyx_refnanny = NULL;
  __Pyx_RefNanny = __Pyx_RefNannyImportAPI("refnanny");
  if (!__Pyx_RefNanny) {
      PyErr_Clear();
      __Pyx_RefNanny = __Pyx_RefNannyImportAPI("Cython.Runtime.refnanny");
      if (!__Pyx_RefNanny)
          Py_FatalError("failed to import 'refnanny' module");
  }
  __pyx_refnanny = __Pyx_RefNanny->SetupContext("PyMODINIT_FUNC PyInit_qhull(void)", __LINE__, __FILE__);
  #endif
  __pyx_empty_tuple = PyTuple_New(0); if (unlikely(!__pyx_empty_tuple)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_empty_bytes = PyBytes_FromStringAndSize("", 0); if (unlikely(!__pyx_empty_bytes)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #ifdef __pyx_binding_PyCFunctionType_USED
  if (__pyx_binding_PyCFunctionType_init() < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  #endif
  
  
  #if defined(__PYX_FORCE_INIT_THREADS) && __PYX_FORCE_INIT_THREADS
  #ifdef WITH_THREAD 
  PyEval_InitThreads();
  #endif
  #endif
  
  #if PY_MAJOR_VERSION < 3
  __pyx_m = Py_InitModule4(__Pyx_NAMESTR("qhull"), __pyx_methods, __Pyx_DOCSTR(__pyx_k_19), 0, PYTHON_API_VERSION);
  #else
  __pyx_m = PyModule_Create(&__pyx_moduledef);
  #endif
  if (!__pyx_m) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  #if PY_MAJOR_VERSION < 3
  Py_INCREF(__pyx_m);
  #endif
  __pyx_b = PyImport_AddModule(__Pyx_NAMESTR(__Pyx_BUILTIN_MODULE_NAME));
  if (!__pyx_b) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  if (__Pyx_SetAttrString(__pyx_m, "__builtins__", __pyx_b) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  
  if (unlikely(__Pyx_InitGlobals() < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__pyx_module_is_main_scipy__spatial__qhull) {
    if (__Pyx_SetAttrString(__pyx_m, "__name__", __pyx_n_s____main__) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;};
  }
  
  if (unlikely(__Pyx_InitCachedBuiltins() < 0)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  
  
  if (__Pyx_ExportFunction("_get_delaunay_info", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__get_delaunay_info, "void (__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, PyObject *, int, int)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_barycentric_inside", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__barycentric_inside, "int (int, double *, double *, double *, double)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_barycentric_coordinate_single", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__barycentric_coordinate_single, "void (int, double *, double *, double *, int)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_barycentric_coordinates", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__barycentric_coordinates, "void (int, double *, double *, double *)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_lift_point", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__lift_point, "void (__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double *)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_distplane", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__distplane, "double (__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, int, double *)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_is_point_fully_outside", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__is_point_fully_outside, "int (__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_find_simplex_bruteforce", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__find_simplex_bruteforce, "int (__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double *, double)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_find_simplex_directed", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__find_simplex_directed, "int (__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double *, int *, double)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_find_simplex", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__find_simplex, "int (__pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, double *, double *, int *, double)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_RidgeIter2D_init", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__RidgeIter2D_init, "void (__pyx_t_5scipy_7spatial_5qhull_RidgeIter2D_t *, __pyx_t_5scipy_7spatial_5qhull_DelaunayInfo_t *, int)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_ExportFunction("_RidgeIter2D_next", (void (*)(void))__pyx_f_5scipy_7spatial_5qhull__RidgeIter2D_next, "void (__pyx_t_5scipy_7spatial_5qhull_RidgeIter2D_t *)") < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  
  if (PyType_Ready(&__pyx_type_5scipy_7spatial_5qhull_RidgeIter2D) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 593; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  if (__Pyx_SetAttrString(__pyx_m, "RidgeIter2D", (PyObject *)&__pyx_type_5scipy_7spatial_5qhull_RidgeIter2D) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 593; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5scipy_7spatial_5qhull_RidgeIter2D = &__pyx_type_5scipy_7spatial_5qhull_RidgeIter2D;
  
  __pyx_ptype_5numpy_dtype = __Pyx_ImportType("numpy", "dtype", sizeof(PyArray_Descr), 0); if (unlikely(!__pyx_ptype_5numpy_dtype)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 149; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5numpy_flatiter = __Pyx_ImportType("numpy", "flatiter", sizeof(PyArrayIterObject), 0); if (unlikely(!__pyx_ptype_5numpy_flatiter)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 159; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5numpy_broadcast = __Pyx_ImportType("numpy", "broadcast", sizeof(PyArrayMultiIterObject), 0); if (unlikely(!__pyx_ptype_5numpy_broadcast)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 163; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5numpy_ndarray = __Pyx_ImportType("numpy", "ndarray", sizeof(PyArrayObject), 0); if (unlikely(!__pyx_ptype_5numpy_ndarray)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 172; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_ptype_5numpy_ufunc = __Pyx_ImportType("numpy", "ufunc", sizeof(PyUFuncObject), 0); if (unlikely(!__pyx_ptype_5numpy_ufunc)) {__pyx_filename = __pyx_f[1]; __pyx_lineno = 849; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  
  

  
  __pyx_t_1 = __Pyx_Import(((PyObject *)__pyx_n_s__threading), 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 13; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__threading, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 13; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_t_1 = __Pyx_Import(((PyObject *)__pyx_n_s__numpy), 0); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 14; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__np, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 14; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_t_1 = PyList_New(2); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 19; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_1));
  __Pyx_INCREF(((PyObject *)__pyx_n_s__Delaunay));
  PyList_SET_ITEM(__pyx_t_1, 0, ((PyObject *)__pyx_n_s__Delaunay));
  __Pyx_GIVEREF(((PyObject *)__pyx_n_s__Delaunay));
  __Pyx_INCREF(((PyObject *)__pyx_n_s__tsearch));
  PyList_SET_ITEM(__pyx_t_1, 1, ((PyObject *)__pyx_n_s__tsearch));
  __Pyx_GIVEREF(((PyObject *)__pyx_n_s__tsearch));
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s____all__, ((PyObject *)__pyx_t_1)) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 19; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(((PyObject *)__pyx_t_1)); __pyx_t_1 = 0;

  
  __pyx_t_1 = __Pyx_GetName(__pyx_m, __pyx_n_s__threading); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 118; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __pyx_t_2 = PyObject_GetAttr(__pyx_t_1, __pyx_n_s__Lock); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 118; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;
  __pyx_t_1 = PyObject_Call(__pyx_t_2, ((PyObject *)__pyx_empty_tuple), NULL); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 118; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_1);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s___qhull_lock, __pyx_t_1) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 118; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_1); __pyx_t_1 = 0;

  
  __pyx_t_1 = PyDict_New(); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 855; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_1));
  __pyx_t_2 = PyTuple_New(1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 855; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_INCREF(__pyx_builtin_object);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_builtin_object);
  __Pyx_GIVEREF(__pyx_builtin_object);
  if (PyDict_SetItemString(((PyObject *)__pyx_t_1), "__doc__", ((PyObject *)__pyx_kp_s_20)) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 855; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __pyx_t_3 = __Pyx_CreateClass(__pyx_t_2, ((PyObject *)__pyx_t_1), __pyx_n_s__Delaunay, "scipy.spatial.qhull"); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 855; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;

  
  __pyx_t_2 = PyCFunction_New(&__pyx_mdef_5scipy_7spatial_5qhull_8Delaunay___init__, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 914; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_4 = PyMethod_New(__pyx_t_2, 0, __pyx_t_3); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 914; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s____init__, __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 914; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;

  
  __pyx_t_4 = PyCFunction_New(&__pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_transform, NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 934; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_2 = PyMethod_New(__pyx_t_4, 0, __pyx_t_3); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 934; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s__transform, __pyx_t_2) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 934; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;

  
  __pyx_t_2 = __Pyx_GetName(__pyx_t_3, __pyx_n_s__transform); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 934; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 933; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  PyTuple_SET_ITEM(__pyx_t_4, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_Call(__pyx_builtin_property, __pyx_t_4, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 933; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s__transform, __pyx_t_2) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 934; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;

  
  __pyx_t_2 = PyCFunction_New(&__pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_vertex_to_simplex, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 957; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_4 = PyMethod_New(__pyx_t_2, 0, __pyx_t_3); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 957; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s__vertex_to_simplex, __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 957; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;

  
  __pyx_t_4 = __Pyx_GetName(__pyx_t_3, __pyx_n_s__vertex_to_simplex); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 957; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_2 = PyTuple_New(1); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 956; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  PyTuple_SET_ITEM(__pyx_t_2, 0, __pyx_t_4);
  __Pyx_GIVEREF(__pyx_t_4);
  __pyx_t_4 = 0;
  __pyx_t_4 = PyObject_Call(__pyx_builtin_property, __pyx_t_2, NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 956; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s__vertex_to_simplex, __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 957; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;

  
  __pyx_t_4 = PyCFunction_New(&__pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_convex_hull, NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 987; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_2 = PyMethod_New(__pyx_t_4, 0, __pyx_t_3); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 987; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s__convex_hull, __pyx_t_2) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 987; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;

  
  __pyx_t_2 = __Pyx_GetName(__pyx_t_3, __pyx_n_s__convex_hull); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 987; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_4 = PyTuple_New(1); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 985; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  PyTuple_SET_ITEM(__pyx_t_4, 0, __pyx_t_2);
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyObject_Call(__pyx_builtin_property, __pyx_t_4, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 985; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s__convex_hull, __pyx_t_2) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 987; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;

  
  __pyx_t_2 = __Pyx_PyBool_FromLong(0); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1033; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_k_10 = __pyx_t_2;
  __Pyx_GIVEREF(__pyx_t_2);
  __pyx_t_2 = 0;
  __pyx_t_2 = PyCFunction_New(&__pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_find_simplex, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1033; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_4 = PyMethod_New(__pyx_t_2, 0, __pyx_t_3); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1033; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s__find_simplex, __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1033; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;

  
  __pyx_t_4 = PyCFunction_New(&__pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_plane_distance, NULL); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1103; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_2 = PyMethod_New(__pyx_t_4, 0, __pyx_t_3); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1103; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s__plane_distance, __pyx_t_2) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1103; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;

  
  __pyx_t_2 = PyCFunction_New(&__pyx_mdef_5scipy_7spatial_5qhull_8Delaunay_lift_points, NULL); if (unlikely(!__pyx_t_2)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1136; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_2);
  __pyx_t_4 = PyMethod_New(__pyx_t_2, 0, __pyx_t_3); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1136; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_2); __pyx_t_2 = 0;
  if (PyObject_SetAttr(__pyx_t_3, __pyx_n_s__lift_points, __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1136; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s__Delaunay, __pyx_t_3) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 855; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __Pyx_DECREF(((PyObject *)__pyx_t_1)); __pyx_t_1 = 0;

  
  __pyx_t_1 = PyDict_New(); if (unlikely(!__pyx_t_1)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(((PyObject *)__pyx_t_1));
  __pyx_t_3 = PyObject_GetAttr(__pyx_m, __pyx_n_s___construct_delaunay); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_4 = __Pyx_GetAttrString(__pyx_t_3, "__doc__"); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_21), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyObject_GetAttr(__pyx_m, __pyx_n_s_5); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_3 = __Pyx_GetAttrString(__pyx_t_4, "__doc__"); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_22), __pyx_t_3) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  __pyx_t_3 = PyObject_GetAttr(__pyx_m, __pyx_n_s_9); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __pyx_t_4 = __Pyx_GetAttrString(__pyx_t_3, "__doc__"); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_23), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyObject_GetAttr(__pyx_m, __pyx_n_s__Delaunay); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__transform); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = __Pyx_GetAttrString(__pyx_t_3, "__doc__"); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_24), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyObject_GetAttr(__pyx_m, __pyx_n_s__Delaunay); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__vertex_to_simplex); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = __Pyx_GetAttrString(__pyx_t_3, "__doc__"); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_25), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyObject_GetAttr(__pyx_m, __pyx_n_s__Delaunay); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__convex_hull); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = __Pyx_GetAttrString(__pyx_t_3, "__doc__"); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_26), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyObject_GetAttr(__pyx_m, __pyx_n_s__Delaunay); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__find_simplex); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = __Pyx_GetAttrString(__pyx_t_3, "__doc__"); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_27), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyObject_GetAttr(__pyx_m, __pyx_n_s__Delaunay); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__plane_distance); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = __Pyx_GetAttrString(__pyx_t_3, "__doc__"); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_28), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyObject_GetAttr(__pyx_m, __pyx_n_s__Delaunay); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_3 = PyObject_GetAttr(__pyx_t_4, __pyx_n_s__lift_points); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = __Pyx_GetAttrString(__pyx_t_3, "__doc__"); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_29), __pyx_t_4) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  __pyx_t_4 = PyObject_GetAttr(__pyx_m, __pyx_n_s__tsearch); if (unlikely(!__pyx_t_4)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_4);
  __pyx_t_3 = __Pyx_GetAttrString(__pyx_t_4, "__doc__"); if (unlikely(!__pyx_t_3)) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_GOTREF(__pyx_t_3);
  __Pyx_DECREF(__pyx_t_4); __pyx_t_4 = 0;
  if (PyDict_SetItem(__pyx_t_1, ((PyObject *)__pyx_kp_u_30), __pyx_t_3) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(__pyx_t_3); __pyx_t_3 = 0;
  if (PyObject_SetAttr(__pyx_m, __pyx_n_s____test__, ((PyObject *)__pyx_t_1)) < 0) {__pyx_filename = __pyx_f[0]; __pyx_lineno = 1; __pyx_clineno = __LINE__; goto __pyx_L1_error;}
  __Pyx_DECREF(((PyObject *)__pyx_t_1)); __pyx_t_1 = 0;

  
  goto __pyx_L0;
  __pyx_L1_error:;
  __Pyx_XDECREF(__pyx_t_1);
  __Pyx_XDECREF(__pyx_t_2);
  __Pyx_XDECREF(__pyx_t_3);
  __Pyx_XDECREF(__pyx_t_4);
  if (__pyx_m) {
    __Pyx_AddTraceback("init scipy.spatial.qhull");
    Py_DECREF(__pyx_m); __pyx_m = 0;
  } else if (!PyErr_Occurred()) {
    PyErr_SetString(PyExc_ImportError, "init scipy.spatial.qhull");
  }
  __pyx_L0:;
  __Pyx_RefNannyFinishContext();
  #if PY_MAJOR_VERSION < 3
  return;
  #else
  return __pyx_m;
  #endif
}



static PyObject *__Pyx_GetName(PyObject *dict, PyObject *name) {
    PyObject *result;
    result = PyObject_GetAttr(dict, name);
    if (!result)
        PyErr_SetObject(PyExc_NameError, name);
    return result;
}

static int __Pyx_ArgTypeTest(PyObject *obj, PyTypeObject *type, int none_allowed,
    const char *name, int exact)
{
    if (!type) {
        PyErr_Format(PyExc_SystemError, "Missing type object");
        return 0;
    }
    if (none_allowed && obj == Py_None) return 1;
    else if (exact) {
        if (Py_TYPE(obj) == type) return 1;
    }
    else {
        if (PyObject_TypeCheck(obj, type)) return 1;
    }
    PyErr_Format(PyExc_TypeError,
        "Argument '%s' has incorrect type (expected %s, got %s)",
        name, type->tp_name, Py_TYPE(obj)->tp_name);
    return 0;
}

static CYTHON_INLINE int __Pyx_IsLittleEndian(void) {
  unsigned int n = 1;
  return *(unsigned char*)(&n) != 0;
}

typedef struct {
  __Pyx_StructField root;
  __Pyx_BufFmt_StackElem* head;
  size_t fmt_offset;
  int new_count, enc_count;
  int is_complex;
  char enc_type;
  char packmode;
} __Pyx_BufFmt_Context;

static void __Pyx_BufFmt_Init(__Pyx_BufFmt_Context* ctx,
                              __Pyx_BufFmt_StackElem* stack,
                              __Pyx_TypeInfo* type) {
  stack[0].field = &ctx->root;
  stack[0].parent_offset = 0;
  ctx->root.type = type;
  ctx->root.name = "buffer dtype";
  ctx->root.offset = 0;
  ctx->head = stack;
  ctx->head->field = &ctx->root;
  ctx->fmt_offset = 0;
  ctx->head->parent_offset = 0;
  ctx->packmode = '@';
  ctx->new_count = 1;
  ctx->enc_count = 0;
  ctx->enc_type = 0;
  ctx->is_complex = 0;
  while (type->typegroup == 'S') {
    ++ctx->head;
    ctx->head->field = type->fields;
    ctx->head->parent_offset = 0;
    type = type->fields->type;
  }
}

static int __Pyx_BufFmt_ParseNumber(const char** ts) {
    int count;
    const char* t = *ts;
    if (*t < '0' || *t > '9') {
      return -1;
    } else {
        count = *t++ - '0';
        while (*t >= '0' && *t < '9') {
            count *= 10;
            count += *t++ - '0';
        }
    }
    *ts = t;
    return count;
}

static void __Pyx_BufFmt_RaiseUnexpectedChar(char ch) {
  char msg[] = {ch, 0};
  PyErr_Format(PyExc_ValueError, "Unexpected format string character: '%s'", msg);
}

static const char* __Pyx_BufFmt_DescribeTypeChar(char ch, int is_complex) {
  switch (ch) {
    case 'b': return "'char'";
    case 'B': return "'unsigned char'";
    case 'h': return "'short'";
    case 'H': return "'unsigned short'";
    case 'i': return "'int'";
    case 'I': return "'unsigned int'";
    case 'l': return "'long'";
    case 'L': return "'unsigned long'";
    case 'q': return "'long long'";
    case 'Q': return "'unsigned long long'";
    case 'f': return (is_complex ? "'complex float'" : "'float'");
    case 'd': return (is_complex ? "'complex double'" : "'double'");
    case 'g': return (is_complex ? "'complex long double'" : "'long double'");
    case 'T': return "a struct";
    case 'O': return "Python object";
    case 'P': return "a pointer";
    case 0: return "end";
    default: return "unparseable format string";
  }
}

static size_t __Pyx_BufFmt_TypeCharToStandardSize(char ch, int is_complex) {
  switch (ch) {
    case '?': case 'c': case 'b': case 'B': return 1;
    case 'h': case 'H': return 2;
    case 'i': case 'I': case 'l': case 'L': return 4;
    case 'q': case 'Q': return 8;
    case 'f': return (is_complex ? 8 : 4);
    case 'd': return (is_complex ? 16 : 8);
    case 'g': {
      PyErr_SetString(PyExc_ValueError, "Python does not define a standard format string size for long double ('g')..");
      return 0;
    }
    case 'O': case 'P': return sizeof(void*);
    default:
      __Pyx_BufFmt_RaiseUnexpectedChar(ch);
      return 0;
    }
}

static size_t __Pyx_BufFmt_TypeCharToNativeSize(char ch, int is_complex) {
  switch (ch) {
    case 'c': case 'b': case 'B': return 1;
    case 'h': case 'H': return sizeof(short);
    case 'i': case 'I': return sizeof(int);
    case 'l': case 'L': return sizeof(long);
    #ifdef HAVE_LONG_LONG
    case 'q': case 'Q': return sizeof(PY_LONG_LONG);
    #endif
    case 'f': return sizeof(float) * (is_complex ? 2 : 1);
    case 'd': return sizeof(double) * (is_complex ? 2 : 1);
    case 'g': return sizeof(long double) * (is_complex ? 2 : 1);
    case 'O': case 'P': return sizeof(void*);
    default: {
      __Pyx_BufFmt_RaiseUnexpectedChar(ch);
      return 0;
    }    
  }
}

typedef struct { char c; short x; } __Pyx_st_short;
typedef struct { char c; int x; } __Pyx_st_int;
typedef struct { char c; long x; } __Pyx_st_long;
typedef struct { char c; float x; } __Pyx_st_float;
typedef struct { char c; double x; } __Pyx_st_double;
typedef struct { char c; long double x; } __Pyx_st_longdouble;
typedef struct { char c; void *x; } __Pyx_st_void_p;
#ifdef HAVE_LONG_LONG
typedef struct { char c; PY_LONG_LONG x; } __Pyx_s_long_long;
#endif

static size_t __Pyx_BufFmt_TypeCharToAlignment(char ch, int is_complex) {
  switch (ch) {
    case '?': case 'c': case 'b': case 'B': return 1;
    case 'h': case 'H': return sizeof(__Pyx_st_short) - sizeof(short);
    case 'i': case 'I': return sizeof(__Pyx_st_int) - sizeof(int);
    case 'l': case 'L': return sizeof(__Pyx_st_long) - sizeof(long);
#ifdef HAVE_LONG_LONG
    case 'q': case 'Q': return sizeof(__Pyx_s_long_long) - sizeof(PY_LONG_LONG);
#endif
    case 'f': return sizeof(__Pyx_st_float) - sizeof(float);
    case 'd': return sizeof(__Pyx_st_double) - sizeof(double);
    case 'g': return sizeof(__Pyx_st_longdouble) - sizeof(long double);
    case 'P': case 'O': return sizeof(__Pyx_st_void_p) - sizeof(void*);
    default:
      __Pyx_BufFmt_RaiseUnexpectedChar(ch);
      return 0;
    }
}

static size_t __Pyx_BufFmt_TypeCharToGroup(char ch, int is_complex) {
  switch (ch) {
    case 'c': case 'b': case 'h': case 'i': case 'l': case 'q': return 'I';
    case 'B': case 'H': case 'I': case 'L': case 'Q': return 'U';
    case 'f': case 'd': case 'g': return (is_complex ? 'C' : 'R');
    case 'O': return 'O';
    case 'P': return 'P';
    default: {
      __Pyx_BufFmt_RaiseUnexpectedChar(ch);
      return 0;
    }    
  }
}

static void __Pyx_BufFmt_RaiseExpected(__Pyx_BufFmt_Context* ctx) {
  if (ctx->head == NULL || ctx->head->field == &ctx->root) {
    const char* expected;
    const char* quote;
    if (ctx->head == NULL) {
      expected = "end";
      quote = "";
    } else {
      expected = ctx->head->field->type->name;
      quote = "'";
    }
    PyErr_Format(PyExc_ValueError,
                 "Buffer dtype mismatch, expected %s%s%s but got %s",
                 quote, expected, quote,
                 __Pyx_BufFmt_DescribeTypeChar(ctx->enc_type, ctx->is_complex));
  } else {
    __Pyx_StructField* field = ctx->head->field;
    __Pyx_StructField* parent = (ctx->head - 1)->field;
    PyErr_Format(PyExc_ValueError,
                 "Buffer dtype mismatch, expected '%s' but got %s in '%s.%s'",
                 field->type->name, __Pyx_BufFmt_DescribeTypeChar(ctx->enc_type, ctx->is_complex),
                 parent->type->name, field->name);
  }
}

static int __Pyx_BufFmt_ProcessTypeChunk(__Pyx_BufFmt_Context* ctx) {
  char group;
  size_t size, offset;
  if (ctx->enc_type == 0) return 0;
  group = __Pyx_BufFmt_TypeCharToGroup(ctx->enc_type, ctx->is_complex);
  do {
    __Pyx_StructField* field = ctx->head->field;
    __Pyx_TypeInfo* type = field->type;
  
    if (ctx->packmode == '@' || ctx->packmode == '^') {
      size = __Pyx_BufFmt_TypeCharToNativeSize(ctx->enc_type, ctx->is_complex);
    } else {
      size = __Pyx_BufFmt_TypeCharToStandardSize(ctx->enc_type, ctx->is_complex);
    }
    if (ctx->packmode == '@') {
      int align_at = __Pyx_BufFmt_TypeCharToAlignment(ctx->enc_type, ctx->is_complex);
      int align_mod_offset;
      if (align_at == 0) return -1;
      align_mod_offset = ctx->fmt_offset % align_at;
      if (align_mod_offset > 0) ctx->fmt_offset += align_at - align_mod_offset;
    }

    if (type->size != size || type->typegroup != group) {
      if (type->typegroup == 'C' && type->fields != NULL) {
        
        size_t parent_offset = ctx->head->parent_offset + field->offset;
        ++ctx->head;
        ctx->head->field = type->fields;
        ctx->head->parent_offset = parent_offset;
        continue;
      }
    
      __Pyx_BufFmt_RaiseExpected(ctx);
      return -1;
    }

    offset = ctx->head->parent_offset + field->offset;
    if (ctx->fmt_offset != offset) {
      PyErr_Format(PyExc_ValueError,
                   "Buffer dtype mismatch; next field is at offset %"PY_FORMAT_SIZE_T"d "
                   "but %"PY_FORMAT_SIZE_T"d expected", ctx->fmt_offset, offset);
      return -1;
    }

    ctx->fmt_offset += size;
  
    --ctx->enc_count; 

    
    while (1) {
      if (field == &ctx->root) {
        ctx->head = NULL;
        if (ctx->enc_count != 0) {
          __Pyx_BufFmt_RaiseExpected(ctx);
          return -1;
        }
        break; 
      }
      ctx->head->field = ++field;
      if (field->type == NULL) {
        --ctx->head;
        field = ctx->head->field;
        continue;
      } else if (field->type->typegroup == 'S') {
        size_t parent_offset = ctx->head->parent_offset + field->offset;
        if (field->type->fields->type == NULL) continue; 
        field = field->type->fields;
        ++ctx->head;
        ctx->head->field = field;
        ctx->head->parent_offset = parent_offset;
        break;
      } else {
        break;
      }
    }
  } while (ctx->enc_count);
  ctx->enc_type = 0;
  ctx->is_complex = 0;
  return 0;    
}

static int __Pyx_BufFmt_FirstPack(__Pyx_BufFmt_Context* ctx) {
  if (ctx->enc_type != 0 || ctx->packmode != '@') {
    PyErr_SetString(PyExc_ValueError, "Buffer packing mode currently only allowed at beginning of format string (this is a defect)");
    return -1;
  }
  return 0;
}

static const char* __Pyx_BufFmt_CheckString(__Pyx_BufFmt_Context* ctx, const char* ts) {
  int got_Z = 0;
  while (1) {
    switch(*ts) {
      case 0:
        if (ctx->enc_type != 0 && ctx->head == NULL) {
          __Pyx_BufFmt_RaiseExpected(ctx);
          return NULL;
        }
        if (__Pyx_BufFmt_ProcessTypeChunk(ctx) == -1) return NULL;
        if (ctx->head != NULL) {
          __Pyx_BufFmt_RaiseExpected(ctx);
          return NULL;
        }
        return ts;
      case ' ':
      case 10:
      case 13:
        ++ts;
        break;
      case '<':
        if (!__Pyx_IsLittleEndian()) {
          PyErr_SetString(PyExc_ValueError, "Little-endian buffer not supported on big-endian compiler");
          return NULL;
        }
        if (__Pyx_BufFmt_FirstPack(ctx) == -1) return NULL;
        ctx->packmode = '=';
        ++ts;
        break;
      case '>':
      case '!':
        if (__Pyx_IsLittleEndian()) {
          PyErr_SetString(PyExc_ValueError, "Big-endian buffer not supported on little-endian compiler");
          return NULL;
        }
        if (__Pyx_BufFmt_FirstPack(ctx) == -1) return NULL;
        ctx->packmode = '=';
        ++ts;
        break;
      case '=':
      case '@':
      case '^':
        if (__Pyx_BufFmt_FirstPack(ctx) == -1) return NULL;
        ctx->packmode = *ts++;
        break;
      case 'T': 
        {
          int i;
          const char* ts_after_sub;
          int struct_count = ctx->new_count;
          ctx->new_count = 1;
          ++ts;
          if (*ts != '{') {
            PyErr_SetString(PyExc_ValueError, "Buffer acquisition: Expected '{' after 'T'");
            return NULL;
          }
          ++ts;
          ts_after_sub = ts;
          for (i = 0; i != struct_count; ++i) {
            ts_after_sub = __Pyx_BufFmt_CheckString(ctx, ts);
            if (!ts_after_sub) return NULL;
          }
          ts = ts_after_sub;
        }
        break;
      case '}': 
        ++ts;
        return ts;
      case 'x':
        if (__Pyx_BufFmt_ProcessTypeChunk(ctx) == -1) return NULL;
        ctx->fmt_offset += ctx->new_count;
        ctx->new_count = 1;
        ctx->enc_count = 0;
        ctx->enc_type = 0;
        ++ts;
        break;
      case 'Z':
        got_Z = 1;
        ++ts;
        if (*ts != 'f' && *ts != 'd' && *ts != 'g') {
          __Pyx_BufFmt_RaiseUnexpectedChar('Z');
          return NULL;
        }        
      case 'c': case 'b': case 'B': case 'h': case 'H': case 'i': case 'I':
      case 'l': case 'L': case 'q': case 'Q':
      case 'f': case 'd': case 'g':
      case 'O':
        if (ctx->enc_type == *ts && got_Z == ctx->is_complex) {
          
          ctx->enc_count += ctx->new_count;
        } else {
          
          if (__Pyx_BufFmt_ProcessTypeChunk(ctx) == -1) return NULL;
          ctx->enc_count = ctx->new_count;
          ctx->enc_type = *ts;
          ctx->is_complex = got_Z;
        }
        ++ts;
        ctx->new_count = 1;
        got_Z = 0;
        break;
      default:
        {
          ctx->new_count = __Pyx_BufFmt_ParseNumber(&ts);
          if (ctx->new_count == -1) { 
            char msg[2] = { *ts, 0 };
            PyErr_Format(PyExc_ValueError,
                         "Does not understand character buffer dtype format string ('%s')", msg);
            return NULL;
          }
        }
      
    }
  }
}

static CYTHON_INLINE void __Pyx_ZeroBuffer(Py_buffer* buf) {
  buf->buf = NULL;
  buf->obj = NULL;
  buf->strides = __Pyx_zeros;
  buf->shape = __Pyx_zeros;
  buf->suboffsets = __Pyx_minusones;
}

static CYTHON_INLINE int __Pyx_GetBufferAndValidate(Py_buffer* buf, PyObject* obj, __Pyx_TypeInfo* dtype, int flags, int nd, int cast, __Pyx_BufFmt_StackElem* stack) {
  if (obj == Py_None) {
    __Pyx_ZeroBuffer(buf);
    return 0;
  }
  buf->buf = NULL;
  if (__Pyx_GetBuffer(obj, buf, flags) == -1) goto fail;
  if (buf->ndim != nd) {
    PyErr_Format(PyExc_ValueError,
                 "Buffer has wrong number of dimensions (expected %d, got %d)",
                 nd, buf->ndim);
    goto fail;
  }
  if (!cast) {
    __Pyx_BufFmt_Context ctx;
    __Pyx_BufFmt_Init(&ctx, stack, dtype);
    if (!__Pyx_BufFmt_CheckString(&ctx, buf->format)) goto fail;
  }
  if ((unsigned)buf->itemsize != dtype->size) {
    PyErr_Format(PyExc_ValueError,
      "Item size of buffer (%"PY_FORMAT_SIZE_T"d byte%s) does not match size of '%s' (%"PY_FORMAT_SIZE_T"d byte%s)",
      buf->itemsize, (buf->itemsize > 1) ? "s" : "",
      dtype->name,
      dtype->size, (dtype->size > 1) ? "s" : "");
    goto fail;
  }
  if (buf->suboffsets == NULL) buf->suboffsets = __Pyx_minusones;
  return 0;
fail:;
  __Pyx_ZeroBuffer(buf);
  return -1;
}

static CYTHON_INLINE void __Pyx_SafeReleaseBuffer(Py_buffer* info) {
  if (info->buf == NULL) return;
  if (info->suboffsets == __Pyx_minusones) info->suboffsets = NULL;
  __Pyx_ReleaseBuffer(info);
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

static void __Pyx_RaiseBufferFallbackError(void) {
  PyErr_Format(PyExc_ValueError,
     "Buffer acquisition failed on assignment; and then reacquiring the old buffer failed too!");
}


static CYTHON_INLINE void __Pyx_RaiseNeedMoreValuesError(Py_ssize_t index) {
    PyErr_Format(PyExc_ValueError,
        #if PY_VERSION_HEX < 0x02050000
                 "need more than %d value%s to unpack", (int)index,
        #else
                 "need more than %zd value%s to unpack", index,
        #endif
                 (index == 1) ? "" : "s");
}

static CYTHON_INLINE void __Pyx_RaiseTooManyValuesError(Py_ssize_t expected) {
    PyErr_Format(PyExc_ValueError,
        #if PY_VERSION_HEX < 0x02050000
            "too many values to unpack (expected %d)", (int)expected);
        #else
            "too many values to unpack (expected %zd)", expected);
        #endif
}

static PyObject *__Pyx_UnpackItem(PyObject *iter, Py_ssize_t index) {
    PyObject *item;
    if (!(item = PyIter_Next(iter))) {
        if (!PyErr_Occurred()) {
            __Pyx_RaiseNeedMoreValuesError(index);
        }
    }
    return item;
}

static int __Pyx_EndUnpack(PyObject *iter, Py_ssize_t expected) {
    PyObject *item;
    if ((item = PyIter_Next(iter))) {
        Py_DECREF(item);
        __Pyx_RaiseTooManyValuesError(expected);
        return -1;
    }
    else if (!PyErr_Occurred())
        return 0;
    else
        return -1;
}

static CYTHON_INLINE void __Pyx_ErrRestore(PyObject *type, PyObject *value, PyObject *tb) {
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
}

static CYTHON_INLINE void __Pyx_ErrFetch(PyObject **type, PyObject **value, PyObject **tb) {
    PyThreadState *tstate = PyThreadState_GET();
    *type = tstate->curexc_type;
    *value = tstate->curexc_value;
    *tb = tstate->curexc_traceback;

    tstate->curexc_type = 0;
    tstate->curexc_value = 0;
    tstate->curexc_traceback = 0;
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
        PyString_AS_STRING(kw_name));
        #endif
}

static void __Pyx_RaiseArgtupleInvalid(
    const char* func_name,
    int exact,
    Py_ssize_t num_min,
    Py_ssize_t num_max,
    Py_ssize_t num_found)
{
    Py_ssize_t num_expected;
    const char *number, *more_or_less;

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
    number = (num_expected == 1) ? "" : "s";
    PyErr_Format(PyExc_TypeError,
        #if PY_VERSION_HEX < 0x02050000
            "%s() takes %s %d positional argument%s (%d given)",
        #else
            "%s() takes %s %zd positional argument%s (%zd given)",
        #endif
        func_name, more_or_less, num_expected, number, num_found);
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
        } else {
            #if PY_MAJOR_VERSION < 3
            if (unlikely(!PyString_CheckExact(key)) && unlikely(!PyString_Check(key))) {
            #else
            if (unlikely(!PyUnicode_CheckExact(key)) && unlikely(!PyUnicode_Check(key))) {
            #endif
                goto invalid_keyword_type;
            } else {
                for (name = first_kw_arg; *name; name++) {
                    #if PY_MAJOR_VERSION >= 3
                    if (PyUnicode_GET_SIZE(**name) == PyUnicode_GET_SIZE(key) &&
                        PyUnicode_Compare(**name, key) == 0) break;
                    #else
                    if (PyString_GET_SIZE(**name) == PyString_GET_SIZE(key) &&
                        _PyString_Eq(**name, key)) break;
                    #endif
                }
                if (*name) {
                    values[name-argnames] = value;
                } else {
                    
                    for (name=argnames; name != first_kw_arg; name++) {
                        if (**name == key) goto arg_passed_twice;
                        #if PY_MAJOR_VERSION >= 3
                        if (PyUnicode_GET_SIZE(**name) == PyUnicode_GET_SIZE(key) &&
                            PyUnicode_Compare(**name, key) == 0) goto arg_passed_twice;
                        #else
                        if (PyString_GET_SIZE(**name) == PyString_GET_SIZE(key) &&
                            _PyString_Eq(**name, key)) goto arg_passed_twice;
                        #endif
                    }
                    if (kwds2) {
                        if (unlikely(PyDict_SetItem(kwds2, key, value))) goto bad;
                    } else {
                        goto invalid_keyword;
                    }
                }
            }
        }
    }
    return 0;
arg_passed_twice:
    __Pyx_RaiseDoubleKeywordsError(function_name, **name);
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
static void __Pyx_RaiseBufferIndexError(int axis) {
  PyErr_Format(PyExc_IndexError,
     "Out of bounds on buffer access (axis %d)", axis);
}



static CYTHON_INLINE void __Pyx_RaiseNoneNotIterableError(void) {
    PyErr_SetString(PyExc_TypeError, "'NoneType' object is not iterable");
}

static void __Pyx_UnpackTupleError(PyObject *t, Py_ssize_t index) {
    if (t == Py_None) {
      __Pyx_RaiseNoneNotIterableError();
    } else if (PyTuple_GET_SIZE(t) < index) {
      __Pyx_RaiseNeedMoreValuesError(PyTuple_GET_SIZE(t));
    } else {
      __Pyx_RaiseTooManyValuesError(index);
    }
}

#if PY_MAJOR_VERSION < 3
static int __Pyx_GetBuffer(PyObject *obj, Py_buffer *view, int flags) {
  #if PY_VERSION_HEX >= 0x02060000
  if (Py_TYPE(obj)->tp_flags & Py_TPFLAGS_HAVE_NEWBUFFER)
      return PyObject_GetBuffer(obj, view, flags);
  #endif
  if (PyObject_TypeCheck(obj, __pyx_ptype_5numpy_ndarray)) return __pyx_pf_5numpy_7ndarray___getbuffer__(obj, view, flags);
  else {
  PyErr_Format(PyExc_TypeError, "'%100s' does not have the buffer interface", Py_TYPE(obj)->tp_name);
  return -1;
    }
}

static void __Pyx_ReleaseBuffer(Py_buffer *view) {
  PyObject* obj = view->obj;
  if (obj) {
if (PyObject_TypeCheck(obj, __pyx_ptype_5numpy_ndarray)) __pyx_pf_5numpy_7ndarray___releasebuffer__(obj, view);
    Py_DECREF(obj);
    view->obj = NULL;
  }
}

#endif

static PyObject *__Pyx_Import(PyObject *name, PyObject *from_list) {
    PyObject *py_import = 0;
    PyObject *empty_list = 0;
    PyObject *module = 0;
    PyObject *global_dict = 0;
    PyObject *empty_dict = 0;
    PyObject *list;
    py_import = __Pyx_GetAttrString(__pyx_b, "__import__");
    if (!py_import)
        goto bad;
    if (from_list)
        list = from_list;
    else {
        empty_list = PyList_New(0);
        if (!empty_list)
            goto bad;
        list = empty_list;
    }
    global_dict = PyModule_GetDict(__pyx_m);
    if (!global_dict)
        goto bad;
    empty_dict = PyDict_New();
    if (!empty_dict)
        goto bad;
    module = PyObject_CallFunctionObjArgs(py_import,
        name, global_dict, empty_dict, list, NULL);
bad:
    Py_XDECREF(empty_list);
    Py_XDECREF(py_import);
    Py_XDECREF(empty_dict);
    return module;
}

static PyObject *__Pyx_CreateClass(
    PyObject *bases, PyObject *dict, PyObject *name, const char *modname)
{
    PyObject *py_modname;
    PyObject *result = 0;

    #if PY_MAJOR_VERSION < 3
    py_modname = PyString_FromString(modname);
    #else
    py_modname = PyUnicode_FromString(modname);
    #endif
    if (!py_modname)
        goto bad;
    if (PyDict_SetItemString(dict, "__module__", py_modname) < 0)
        goto bad;
    #if PY_MAJOR_VERSION < 3
    result = PyClass_New(bases, dict, name);
    #else
    result = PyObject_CallFunctionObjArgs((PyObject *)&PyType_Type, name, bases, dict, NULL);
    #endif
bad:
    Py_XDECREF(py_modname);
    return result;
}

#if PY_MAJOR_VERSION < 3
static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb) {
    Py_XINCREF(type);
    Py_XINCREF(value);
    Py_XINCREF(tb);
    
    if (tb == Py_None) {
        Py_DECREF(tb);
        tb = 0;
    }
    else if (tb != NULL && !PyTraceBack_Check(tb)) {
        PyErr_SetString(PyExc_TypeError,
            "raise: arg 3 must be a traceback or None");
        goto raise_error;
    }
    
    if (value == NULL) {
        value = Py_None;
        Py_INCREF(value);
    }
    #if PY_VERSION_HEX < 0x02050000
    if (!PyClass_Check(type))
    #else
    if (!PyType_Check(type))
    #endif
    {
        
        if (value != Py_None) {
            PyErr_SetString(PyExc_TypeError,
                "instance exception may not have a separate value");
            goto raise_error;
        }
        
        Py_DECREF(value);
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

static void __Pyx_Raise(PyObject *type, PyObject *value, PyObject *tb) {
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
    } else if (!PyExceptionClass_Check(type)) {
        PyErr_SetString(PyExc_TypeError,
            "raise: exception class must be a subclass of BaseException");
        goto bad;
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
    return;
}
#endif

static CYTHON_INLINE PyObject *__Pyx_PyInt_to_py_flagT(flagT val) {
    const flagT neg_one = (flagT)-1, const_zero = (flagT)0;
    const int is_unsigned = const_zero < neg_one;
    if ((sizeof(flagT) == sizeof(char))  ||
        (sizeof(flagT) == sizeof(short))) {
        return PyInt_FromLong((long)val);
    } else if ((sizeof(flagT) == sizeof(int)) ||
               (sizeof(flagT) == sizeof(long))) {
        if (is_unsigned)
            return PyLong_FromUnsignedLong((unsigned long)val);
        else
            return PyInt_FromLong((long)val);
    } else if (sizeof(flagT) == sizeof(PY_LONG_LONG)) {
        if (is_unsigned)
            return PyLong_FromUnsignedLongLong((unsigned PY_LONG_LONG)val);
        else
            return PyLong_FromLongLong((PY_LONG_LONG)val);
    } else {
        int one = 1; int little = (int)*(unsigned char *)&one;
        unsigned char *bytes = (unsigned char *)&val;
        return _PyLong_FromByteArray(bytes, sizeof(flagT), 
                                     little, !is_unsigned);
    }
}

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

#endif

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

#endif

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
            return PyLong_AsUnsignedLong(x);
        } else {
            return PyLong_AsLong(x);
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
            return PyLong_AsUnsignedLongLong(x);
        } else {
            return PyLong_AsLongLong(x);
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
            return PyLong_AsUnsignedLong(x);
        } else {
            return PyLong_AsLong(x);
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
            return PyLong_AsUnsignedLongLong(x);
        } else {
            return PyLong_AsLongLong(x);
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
            return PyLong_AsUnsignedLong(x);
        } else {
            return PyLong_AsLong(x);
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
            return PyLong_AsUnsignedLongLong(x);
        } else {
            return PyLong_AsLongLong(x);
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

static void __Pyx_WriteUnraisable(const char *name) {
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

static int __Pyx_ExportFunction(const char *name, void (*f)(void), const char *sig) {
    PyObject *d = 0;
    PyObject *cobj = 0;
    union {
        void (*fp)(void);
        void *p;
    } tmp;

    d = PyObject_GetAttrString(__pyx_m, (char *)"__pyx_capi__");
    if (!d) {
        PyErr_Clear();
        d = PyDict_New();
        if (!d)
            goto bad;
        Py_INCREF(d);
        if (PyModule_AddObject(__pyx_m, (char *)"__pyx_capi__", d) < 0)
            goto bad;
    }
    tmp.fp = f;
#if PY_VERSION_HEX >= 0x02070000 && !(PY_MAJOR_VERSION==3&&PY_MINOR_VERSION==0)
    cobj = PyCapsule_New(tmp.p, sig, 0);
#else
    cobj = PyCObject_FromVoidPtrAndDesc(tmp.p, (void *)sig, 0);
#endif
    if (!cobj)
        goto bad;
    if (PyDict_SetItemString(d, name, cobj) < 0)
        goto bad;
    Py_DECREF(cobj);
    Py_DECREF(d);
    return 0;
bad:
    Py_XDECREF(cobj);
    Py_XDECREF(d);
    return -1;
}

#ifndef __PYX_HAVE_RT_ImportType
#define __PYX_HAVE_RT_ImportType
static PyTypeObject *__Pyx_ImportType(const char *module_name, const char *class_name,
    long size, int strict)
{
    PyObject *py_module = 0;
    PyObject *result = 0;
    PyObject *py_name = 0;
    char warning[200];

    py_module = __Pyx_ImportModule(module_name);
    if (!py_module)
        goto bad;
    #if PY_MAJOR_VERSION < 3
    py_name = PyString_FromString(class_name);
    #else
    py_name = PyUnicode_FromString(class_name);
    #endif
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
    if (!strict && ((PyTypeObject *)result)->tp_basicsize > size) {
        PyOS_snprintf(warning, sizeof(warning), 
            "%s.%s size changed, may indicate binary incompatibility",
            module_name, class_name);
        #if PY_VERSION_HEX < 0x02050000
        PyErr_Warn(NULL, warning);
        #else
        PyErr_WarnEx(NULL, warning, 0);
        #endif
    }
    else if (((PyTypeObject *)result)->tp_basicsize != size) {
        PyErr_Format(PyExc_ValueError, 
            "%s.%s has the wrong size, try recompiling",
            module_name, class_name);
        goto bad;
    }
    return (PyTypeObject *)result;
bad:
    Py_XDECREF(py_module);
    Py_XDECREF(result);
    return 0;
}
#endif

#ifndef __PYX_HAVE_RT_ImportModule
#define __PYX_HAVE_RT_ImportModule
static PyObject *__Pyx_ImportModule(const char *name) {
    PyObject *py_name = 0;
    PyObject *py_module = 0;

    #if PY_MAJOR_VERSION < 3
    py_name = PyString_FromString(name);
    #else
    py_name = PyUnicode_FromString(name);
    #endif
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

#include "compile.h"
#include "frameobject.h"
#include "traceback.h"

static void __Pyx_AddTraceback(const char *funcname) {
    PyObject *py_srcfile = 0;
    PyObject *py_funcname = 0;
    PyObject *py_globals = 0;
    PyCodeObject *py_code = 0;
    PyFrameObject *py_frame = 0;

    #if PY_MAJOR_VERSION < 3
    py_srcfile = PyString_FromString(__pyx_filename);
    #else
    py_srcfile = PyUnicode_FromString(__pyx_filename);
    #endif
    if (!py_srcfile) goto bad;
    if (__pyx_clineno) {
        #if PY_MAJOR_VERSION < 3
        py_funcname = PyString_FromFormat( "%s (%s:%d)", funcname, __pyx_cfilenm, __pyx_clineno);
        #else
        py_funcname = PyUnicode_FromFormat( "%s (%s:%d)", funcname, __pyx_cfilenm, __pyx_clineno);
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
    py_globals = PyModule_GetDict(__pyx_m);
    if (!py_globals) goto bad;
    py_code = PyCode_New(
        0,            
        #if PY_MAJOR_VERSION >= 3
        0,            
        #endif
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
        __pyx_lineno,   
        __pyx_empty_bytes  
    );
    if (!py_code) goto bad;
    py_frame = PyFrame_New(
        PyThreadState_GET(), 
        py_code,             
        py_globals,          
        0                    
    );
    if (!py_frame) goto bad;
    py_frame->f_lineno = __pyx_lineno;
    PyTraceBack_Here(py_frame);
bad:
    Py_XDECREF(py_srcfile);
    Py_XDECREF(py_funcname);
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
