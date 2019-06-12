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

#include "jonkervolgenant/jonkervolgenant.h"
#include "numpy/arrayobject.h"
#include "numpy/ndarraytypes.h"

static PyObject *calculate_assignment(PyObject *self, PyObject *args) {
    PyObject *obj_cost = NULL;
    if (!PyArg_ParseTuple(args, "O", &obj_cost))
        return NULL;

    PyObject *obj_cont = PyArray_ContiguousFromAny(obj_cost, NPY_DOUBLE, 2, 2);
    if (obj_cont == NULL)
        return NULL;

    int num_rows = PyArray_DIM(obj_cont, 0);
    int num_cols = PyArray_DIM(obj_cont, 1);

    double *cost_matrix = (double *)PyArray_DATA(obj_cont);
    if (cost_matrix == NULL)
        return NULL;

    npy_intp dim[1] = {num_rows};
    PyObject *x = PyArray_SimpleNew(1, dim, NPY_INT64);
    int result = solve_jonker_volgenant(num_rows, num_cols, cost_matrix,
                                        PyArray_DATA((PyArrayObject *)x));
    if (result != 0) {
        PyErr_SetString(PyExc_ValueError, "cost matrix is infeasible");
        Py_DECREF(x);
        return NULL;
    } else {
        return x;
    }
}

static PyMethodDef jonkervolgenant_methods[] = {
    {"calculate_assignment", calculate_assignment, METH_VARARGS,
     "Solves rectangular linear sum assignment problem."},
    {NULL, NULL, 0, NULL}};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_jonkervolgenant",
    "Jonker-Volgenant algorithm for rectangular linear sum assignment.",
    -1,
    jonkervolgenant_methods,
    NULL,
    NULL,
    NULL,
    NULL,
};

PyObject *PyInit__jonkervolgenant(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    import_array();
    return m;
}
