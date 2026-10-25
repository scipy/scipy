// Raw C callbacks from xsf for use with scipy.LowLevelCallable. Each is exported as a PyCapsule
// whose name is the C signature.

#include <xsf/numpy.h>
#include <xsf/cpu/stats.h>


// user_data points to {p, a, b}.
static double
_genhyperbolic_pdf_callback(double x, void *user_data)
{
    const double *params = static_cast<const double *>(user_data);
    return xsf::cpu::genhyperbolic_pdf(x, params[0], params[1], params[2]);
}


static int
_special_callbacks_module_exec(PyObject *module)
{
    // Needed by wrap_PyUFunc_getfperr, which sf_error.cc calls.
    if (_import_array() < 0) { return -1; }
    if (_import_umath() < 0) { return -1; }

    PyObject *capsule = PyCapsule_New(
        (void *)_genhyperbolic_pdf_callback, "double (double, void *)", NULL);
    if (capsule == NULL) {
        return -1;
    }
    int ret = PyModule_AddObjectRef(module, "_genhyperbolic_pdf_callback", capsule);
    Py_DECREF(capsule);
    return ret;
}


static PyModuleDef_Slot _special_callbacks_slots[] = {
    {Py_mod_exec, (void *)_special_callbacks_module_exec},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  /* Python 3.13+ */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, NULL},
};


static struct PyModuleDef _special_callbacks_def = {
    /* m_base     */ PyModuleDef_HEAD_INIT,
    /* m_name     */ "_special_callbacks",
    /* m_doc      */ NULL,
    /* m_size     */ 0,
    /* m_methods  */ NULL,
    /* m_slots    */ _special_callbacks_slots,
    /* m_traverse */ NULL,
    /* m_clear    */ NULL,
    /* m_free     */ NULL
};


PyMODINIT_FUNC PyInit__special_callbacks()
{
    return PyModuleDef_Init(&_special_callbacks_def);
}
