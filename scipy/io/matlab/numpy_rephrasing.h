#include <numpy/arrayobject.h>
#define PyArray_Set_BASE(arr, obj) PyArray_BASE(arr) = obj
#define PyArray_PyANewFromDescr(descr, nd, dims, data, parent)                 \
        PyArray_NewFromDescr(&PyArray_Type, descr, nd, dims,                  \
                             NULL, data, 0, parent)
