/*
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following
   disclaimer in the documentation and/or other materials provided
   with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "numpy/arrayobject.h"
#include "numpy/ndarraytypes.h"
#include "rectangular_lsap/rectangular_lsap.h"

static PyObject*
calculate_assignment(PyObject* self, PyObject* args)
{
    PyObject* x = NULL;
    PyObject* obj_cost = NULL;
    if (!PyArg_ParseTuple(args, "O", &obj_cost))
        return NULL;

    PyArrayObject* obj_cont =
      (PyArrayObject*)PyArray_ContiguousFromAny(obj_cost, NPY_DOUBLE, 2, 2);
    if (obj_cont == NULL) {
        PyErr_SetString(PyExc_TypeError, "invalid cost matrix object");
        return NULL;
    }

    double* cost_matrix = (double*)PyArray_DATA(obj_cont);
    if (cost_matrix == NULL) {
        PyErr_SetString(PyExc_TypeError, "invalid cost matrix object");
        goto cleanup;
    }

    int num_rows = PyArray_DIM(obj_cont, 0);
    int num_cols = PyArray_DIM(obj_cont, 1);

    npy_intp dim[1] = { num_rows };
    x = PyArray_SimpleNew(1, dim, NPY_INT64);
    int result = solve_rectangular_linear_sum_assignment(
      num_rows, num_cols, cost_matrix, PyArray_DATA((PyArrayObject*)x));
    if (result != 0) {
        PyErr_SetString(PyExc_ValueError, "cost matrix is infeasible");
        Py_DECREF(x);
        x = NULL;
    }

cleanup:
    Py_DECREF((PyObject*)obj_cont);
    return x;
}

static PyMethodDef lsap_module_methods[] = {
    { "calculate_assignment", calculate_assignment, METH_VARARGS,
      "Solves the rectangular linear sum assignment problem." },
    { NULL, NULL, 0, NULL }
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_lsap_module",
    "Solves the rectangular linear sum assignment.",
    -1,
    lsap_module_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyObject*
PyInit__lsap_module(void)
{
    PyObject* m;
    m = PyModule_Create(&moduledef);
    import_array();
    return m;
}
