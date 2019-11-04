#include <Python.h>

#include <stdio.h>


#ifdef _MSC_VER
#pragma float_control(precise, on)
#pragma fenv_access (on)
#endif


static char get_fpu_mode_doc[] = (
    "get_fpu_mode()\n"
    "\n"
    "Get the current FPU control word, in a platform-dependent format.\n"
    "Returns None if not implemented on current platform.");


static PyObject *
get_fpu_mode(PyObject *self, PyObject *args)
{
    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

#if defined(_MSC_VER)
    {
        unsigned int result = 0;
        result = _controlfp(0, 0);
        return PyLong_FromLongLong(result);
    }
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
    {
        unsigned short cw = 0;
        __asm__("fstcw %w0" : "=m" (cw));
        return PyLong_FromLongLong(cw);
    }
#else
    Py_RETURN_NONE;
#endif
}


static struct PyMethodDef methods[] = {
    {"get_fpu_mode", get_fpu_mode, METH_VARARGS, get_fpu_mode_doc},
    {NULL, NULL, 0, NULL}
};


#if PY_MAJOR_VERSION >= 3

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_fpumode",
    NULL,
    -1,
    methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *PyInit__fpumode(void)
{
    return PyModule_Create(&moduledef);
}

#else

PyMODINIT_FUNC init_fpumode(void)
{
    Py_InitModule("_fpumode", methods);
}

#endif
