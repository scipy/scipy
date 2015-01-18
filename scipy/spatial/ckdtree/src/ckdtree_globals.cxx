/*
 * This would break SciPy with NumPy 1.6 so just accept the compiler
 * warning for now.
 * #define NPY_NO_DEPRECATED_API NPY_1_9_API_VERSION 
 *
 */
 
#include <Python.h>
#include "numpy/arrayobject.h"

npy_float64 infinity;
int number_of_processors;

