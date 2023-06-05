extern "C" {

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <stdbool.h>
} // extern "C"

#include "HConst.h"
#include "SimplexConst.h"

extern "C" {
static PyMethodDef highsConstantsMethods[] = {
    {nullptr, nullptr, 0, nullptr} /* Sentinel */
};


static struct PyModuleDef _highs_constants_ext = {
    PyModuleDef_HEAD_INIT, "_highs_constants", /* name of module */
    nullptr,                              /* module documentation, may be NULL */
    -1, /* size of per-interpreter state of the module, or -1 if the module
   keeps state in global variables. */
    highsConstantsMethods};


PyMODINIT_FUNC PyInit__highs_constants(void) {
    PyObject *mod = PyModule_Create(&_highs_constants_ext);

    // constants that should come from highs_bindings but currently don't
    // FIXME: error handling
    PyModule_AddIntConstant(mod, "kHighsDebugLevelNone", HighsDebugLevel::kHighsDebugLevelNone);
    PyModule_AddIntConstant(mod, "kSimplexEdgeWeightStrategyDantzig", SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategyDantzig);
    PyModule_AddIntConstant(mod, "kSimplexEdgeWeightStrategyDevex", SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategyDevex);
    PyModule_AddIntConstant(mod, "kSimplexEdgeWeightStrategyChoose", SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategyChoose);
    PyModule_AddIntConstant(mod, "kSimplexEdgeWeightStrategySteepestEdge", SimplexEdgeWeightStrategy::kSimplexEdgeWeightStrategySteepestEdge);
    PyModule_AddIntConstant(mod, "kSimplexStrategyDual", SimplexStrategy::kSimplexStrategyDual);
    PyModule_AddIntConstant(mod, "kSimplexCrashStrategyOff", SimplexCrashStrategy::kSimplexCrashStrategyOff);

    return mod;
}

} // extern "C"
