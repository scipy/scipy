/*
 * Catch C++ standard exceptions and raise corresponding
 * Python exceptions. Use simplified GIL API to ensure 
 * we have the GIL.
 */

#include <new>
#include <typeinfo>
#include <stdexcept>
#include <ios>

#include <Python.h>

extern "C" void
translate_cpp_exception() 
{
    try {
        if (PyErr_Occurred())
            return; 
        else
            throw;
    } catch (const std::bad_alloc &exn) {
        PyErr_NoMemory();
    } catch (const std::bad_cast &exn) {
        PyErr_SetString(PyExc_TypeError, exn.what());
    } catch (const std::domain_error &exn) {
        PyErr_SetString(PyExc_ValueError, exn.what());
    } catch (const std::invalid_argument &exn) {
        PyErr_SetString(PyExc_ValueError, exn.what());
    } catch (const std::ios_base::failure &exn) {
        PyErr_SetString(PyExc_IOError, exn.what());
    } catch (const std::out_of_range &exn) {
        PyErr_SetString(PyExc_IndexError, exn.what());
    } catch (const std::overflow_error &exn) {
        PyErr_SetString(PyExc_OverflowError, exn.what());
    } catch (const std::range_error &exn) {
        PyErr_SetString(PyExc_ArithmeticError, exn.what());
    } catch (const std::underflow_error &exn) {
        PyErr_SetString(PyExc_ArithmeticError, exn.what());
    } catch (const std::logic_error &exn) {
        PyErr_SetString(PyExc_RuntimeError, exn.what());
    } catch (const std::exception& exn) {
        PyErr_SetString(PyExc_RuntimeError, exn.what());
    } catch (...) {
        PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
    }
}

extern "C" void
translate_cpp_exception_with_gil() 
{
    PyGILState_STATE state = PyGILState_Ensure();
    translate_cpp_exception();
    PyGILState_Release(state);
}
