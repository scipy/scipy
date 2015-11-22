/*  
 *  Utility code for creating a NumPy array from the memory
 *  in a std::vector and setting the array's base attribute
 *  to the owner of the std::vector.
 *
 */


#ifndef CKDTREE_CONVERTER
#define CKDTREE_CONVERTER

#include <Python.h>
#include "numpy/arrayobject.h"

static PyObject *
numpyarray_with_base(void *data, PyObject *base_obj, npy_intp n, 
    PyObject *dtype, int aligned)
{
    PyObject *x = NULL;
    
    int flags = NPY_ARRAY_C_CONTIGUOUS | NPY_ARRAY_WRITABLE;
    if (aligned) flags |= NPY_ARRAY_ALIGNED;
    
    Py_XINCREF(base_obj);
    Py_INCREF(dtype);
    
    if (!PyArray_DescrCheck(dtype)) {    
        PyErr_SetString(PyExc_ValueError, 
            "dtype is not an instance of numpy.dtype.");
        goto error;   
    } 
    
    x = PyArray_NewFromDescr(&PyArray_Type, (PyArray_Descr *)dtype, 
            1, &n, NULL, data, flags, NULL);
    if (x == NULL) goto error;

#if NPY_API_VERSION >= 0x00000007
    if (PyArray_SetBaseObject((PyArrayObject *)x, base_obj) == -1) {
        PyErr_SetString(PyExc_ValueError, 
            "Failed to set base object on newly formed array.");
        goto error;
    }
#else            
    if (PyArray_BASE(x) != NULL) {
        PyErr_SetString(PyExc_ValueError, 
            "The newly formed array already had a base object.");
        goto error;  
    }
    PyArray_BASE(x) = base_obj;
#endif
    return x;
error:
    Py_DECREF(dtype);
    Py_XDECREF(base_obj);
    Py_XDECREF(x);
    return NULL;
}


static PyObject *
numpyarray_with_base(void *data, PyObject *base_obj, npy_intp n, 
    int typenum, int aligned)
{
    PyObject *x = NULL;
    
    Py_XINCREF(base_obj);
    
    x = PyArray_SimpleNewFromData(1, &n, typenum, data);
    if (x == NULL) goto error;

#if NPY_API_VERSION >= 0x00000007
    if (PyArray_SetBaseObject((PyArrayObject*)x, base_obj) == -1) {
        PyErr_SetString(PyExc_ValueError, 
            "Failed to set base object on newly formed array.");
        goto error;
    }
#else            
    if (PyArray_BASE(x) != NULL) {
        PyErr_SetString(PyExc_ValueError, 
            "The newly formed array already had a base object.");
        goto error;  
    }
    PyArray_BASE(x) = base_obj;
#endif    
    return x;
error:
    Py_XDECREF(base_obj);
    Py_XDECREF(x);
    return NULL;
}

inline PyObject *
intp_array_with_base(std::vector<npy_intp> *v, PyObject *base_obj)
{
    std::vector<npy_intp> &tmp = *v;
    void *data = (void*)(&tmp[0]);
    return numpyarray_with_base(data, base_obj, tmp.size(), NPY_INTP, 1);
}

inline PyObject *
float64_array_with_base(std::vector<npy_float64> *v, PyObject *base_obj)
{
    std::vector<npy_float64> &tmp = *v;
    void *data = (void*)(&tmp[0]);
    return numpyarray_with_base(data, base_obj, tmp.size(), NPY_FLOAT64, 1);
}

inline PyObject *
record_array_with_base(std::vector<T> *v, PyObject *b, PyObject *dt, int aligned)
{
    std::vector<T> &tmp = *v;
    void *data = (void*)(&tmp[0]);
    return numpyarray_with_base(data, b, tmp.size(), dt, aligned);
}

    
#endif

