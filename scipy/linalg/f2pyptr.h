/* Header used to extract function pointers from f2py wrappers. */
 
#ifndef F2PYPTR_H_
#define F2PYPTR_H_

#include <Python.h>

void *f2py_ptr(PyObject *obj) 
{
#if PY_VERSION_HEX < 0x03000000
    if (PyCObject_Check(obj)) {
        return PyCObject_AsVoidPtr(obj);}
#endif

#if PY_VERSION_HEX >= 0x02070000
    if (PyCapsule_CheckExact(obj)) {
        return PyCapsule_GetPointer(obj, NULL);}
#endif
    PyErr_SetString(PyExc_ValueError, "Not an object containing a void ptr");
    return NULL;
}
#endif
