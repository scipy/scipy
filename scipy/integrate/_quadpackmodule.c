/*
  From Multipack project
 */
#include "__quadpack.h"

static struct PyMethodDef quadpack_module_methods[] = {
{"_qagse", quadpack_qagse, METH_VARARGS, doc_qagse},
{"_qagie", quadpack_qagie, METH_VARARGS, doc_qagie},
{"_qagpe", quadpack_qagpe, METH_VARARGS, doc_qagpe},
{"_qawoe", quadpack_qawoe, METH_VARARGS, doc_qawoe},
{"_qawfe", quadpack_qawfe, METH_VARARGS, doc_qawfe},
{"_qawse", quadpack_qawse, METH_VARARGS, doc_qawse},
{"_qawce", quadpack_qawce, METH_VARARGS, doc_qawce},
{NULL,		NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_quadpack",
    NULL,
    -1,
    quadpack_module_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyMODINIT_FUNC
PyInit__quadpack(void)
{
    PyObject *module, *mdict;

    import_array();

    module = PyModule_Create(&moduledef);
    if (module == NULL) {
        return NULL;
    }

    mdict = PyModule_GetDict(module);
    if (mdict == NULL) {
        return NULL;
    }

    quadpack_error = PyErr_NewException ("_quadpack.error", NULL, NULL);
    if (quadpack_error == NULL) {
        return NULL;
    }
    if (PyDict_SetItemString(mdict, "error", quadpack_error)) {
        return NULL;
    }

    return module;
}
