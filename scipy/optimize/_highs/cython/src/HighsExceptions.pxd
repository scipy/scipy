# distutils: language=c++
# cython: language_level=3

cdef extern from "Python.h" nogil:
    ctypedef struct PyObject

cdef extern from *:
    """
    #include <Python.h>
    #include <stdexcept>
    #include "HighsExceptions.h"
    #include <ios>

    PyObject *PresolveExceptionError;

    void create_highs_exceptions() {
        PresolveExceptionError = PyErr_NewException("highs.PresolveExceptionError", NULL, NULL);
    }

    void highs_exception_handler() {
        try {
            if (PyErr_Occurred()) {
                ; // let the latest Python exn pass through and ignore the current one
            } else {
                throw;
            }
        }  catch (const PresolveTooLarge& exn) {
            // Add mapping of PresolveTooLarge -> PresolveException
            PyErr_SetString(PresolveExceptionError, exn.what());
        } catch (...) {
            PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
        }
    }
    """
    cdef PyObject* PresolveExceptionError
    cdef void create_highs_exceptions()
    cdef void highs_exception_handler()
