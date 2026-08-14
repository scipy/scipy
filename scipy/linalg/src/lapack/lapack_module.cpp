/**
 * @file
 * @brief Module definition for the C++ LAPACK wrappers (`_flapack_cpp` / `_flapack_cpp_64`).
 *
 * This is the only file here that calls `_import_array()`; every other .cpp defines
 * `NO_IMPORT_ARRAY` so the extension shares one numpy C-API table.
 */
#define PY_ARRAY_UNIQUE_SYMBOL scipy_lapack_ARRAY_API
#include "lapack_helpers.hpp"

/* The ILP64 build differs from the LP64 one only in the module's name. */
#ifdef HAVE_BLAS_ILP64
#define FLAPACK_CPP_MODULE_STRING "_flapack_cpp_64"
#define FLAPACK_CPP_PYINIT        PyInit__flapack_cpp_64
#else
#define FLAPACK_CPP_MODULE_STRING "_flapack_cpp"
#define FLAPACK_CPP_PYINIT        PyInit__flapack_cpp
#endif

/* One table per `.pyf.src` group, each contributed by the matching wrapper .cpp.  Adding a
 * group means one more `extern` here and one more line in the exec slot below. */
namespace lapack {
namespace capi {
    extern PyMethodDef gen_methods[];
}
}

static int lapack_module_exec(PyObject *module)
{
    if (_import_array() < 0) { return -1; }
    return PyModule_AddFunctions(module, lapack::capi::gen_methods);
}

/* The wrappers keep no module or process state; the one callback slot is thread-local. */
static PyModuleDef_Slot lapack_module_slots[] = {
    {Py_mod_exec, reinterpret_cast<void *>(lapack_module_exec)},
    {Py_mod_multiple_interpreters, Py_MOD_PER_INTERPRETER_GIL_SUPPORTED},
#if PY_VERSION_HEX >= 0x030d00f0  /* Python 3.13+ */
    {Py_mod_gil, Py_MOD_GIL_NOT_USED},
#endif
    {0, nullptr},
};

static PyModuleDef lapack_moduledef = {
    PyModuleDef_HEAD_INIT,
    FLAPACK_CPP_MODULE_STRING,  /* m_name    */
    nullptr,                    /* m_doc     */
    0,                          /* m_size    */
    nullptr,                    /* m_methods, added in the exec slot */
    lapack_module_slots,        /* m_slots   */
    nullptr,                    /* m_traverse */
    nullptr,                    /* m_clear   */
    nullptr,                    /* m_free    */
};

PyMODINIT_FUNC FLAPACK_CPP_PYINIT(void)
{
    return PyModuleDef_Init(&lapack_moduledef);
}
