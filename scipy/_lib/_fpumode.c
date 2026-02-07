#include <Python.h>

#include <stdio.h>


#ifdef _MSC_VER
#include <float.h>
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

static struct PyModuleDef_Slot _fpumode_slots[] = {
#if PY_VERSION_HEX >= 0x030c00f0  // Python 3.12+
    // signal that this module can be imported in isolated subinterpreters
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#endif
#if PY_VERSION_HEX >= 0x030d00f0  // Python 3.13+
    // signal that this module supports running without an active GIL
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};

static struct PyModuleDef moduledef = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "_fpumode",
    .m_size = 0,
    .m_methods = methods,
    .m_slots = _fpumode_slots,
};

PyMODINIT_FUNC
PyInit__fpumode(void)
{
    return PyModuleDef_Init(&moduledef);
}
