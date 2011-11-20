#include <Python.h>
#define PY_ARRAY_UNIQUE_SYMBOL _scipy_signal_ARRAY_API
#define NO_IMPORT_ARRAY
#include <numpy/noprefix.h>

/*
 *  int check_all_numbers(PyObject *x, char *name)
 *
 *  Determine if all the objects in the object array x are numbers.
 *  x must be a pointer to an ndarray of object dtype.
 *  name must be a null-terminated string that holds the name of the
 *. variable being checked.  It is used in the exception message if any
 *  non-numeric objects are found in x.
 *
 *  Return values:
 *    0:  all objects pass PyNumber_Check()
 *    1:  at least one object does not satisfy PyNumber_Check()
 *   -1:  internal error (PyArray_IterNew(x) failed).
 *
 */

int check_all_numbers(PyObject *x, const char *name)
{
    PyArrayIterObject *iter;
    int result;

    iter = (PyArrayIterObject *) PyArray_IterNew(x);
    if (iter == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "internal error, possibly out of memory.");
        return -1;
    }
    result = 0;
    while (iter->index < iter->size) {
        if (!PyNumber_Check(*((PyObject **) (iter->dataptr)))) {
            PyErr_Format(PyExc_ValueError,
                         "%s is an object array containing objects that are not numbers.",
                         name);
            result = 1;
            break;
        }
        PyArray_ITER_NEXT(iter);
    }
    Py_DECREF(iter);
    return result;
}
