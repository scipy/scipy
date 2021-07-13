#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "_direct/direct.h"

static PyObject *
direct(PyObject *self, PyObject *args)
{
    const char *command;
    int sts;

    if (!PyArg_ParseTuple(args, "s", &command))
        return NULL;
    sts = system(command);
    return PyLong_FromLong(sts);
}