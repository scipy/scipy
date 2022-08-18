#include <Python.h>
#include "ccallback.h"
#include "unuran.h"

#define UNURAN_THUNK(CAST_FUNC, FUNCNAME, LEN)                                              \
    PyGILState_STATE gstate = PyGILState_Ensure();                                          \
    /* If an error has occured, return INFINITY. */                                         \
    if (PyErr_Occurred()) return UNUR_INFINITY;                                             \
    ccallback_t *callback = ccallback_obtain();                                             \
                                                                                            \
    PyObject *arg1 = NULL, *argobj = NULL, *res = NULL, *funcname = NULL;                   \
    double result = 0.;                                                                     \
    int error = 0;                                                                          \
                                                                                            \
    argobj = CAST_FUNC(x);                                                                  \
    if (argobj == NULL) {                                                                   \
        error = 1;                                                                          \
        goto done;                                                                          \
    }                                                                                       \
                                                                                            \
    funcname = Py_BuildValue("s#", FUNCNAME, LEN);                                          \
    if (funcname == NULL) {                                                                 \
        error = 1;                                                                          \
        goto done;                                                                          \
    }                                                                                       \
                                                                                            \
    arg1 = PyTuple_New(2);                                                                  \
    if (arg1 == NULL) {                                                                     \
        error = 1;                                                                          \
        goto done;                                                                          \
    }                                                                                       \
                                                                                            \
    PyTuple_SET_ITEM(arg1, 0, argobj);                                                      \
    PyTuple_SET_ITEM(arg1, 1, funcname);                                                    \
    argobj = NULL; funcname = NULL;                                                         \
                                                                                            \
    res = PyObject_CallObject(callback->py_function, arg1);                                 \
    if (res == NULL) {                                                                      \
        error = 1;                                                                          \
        goto done;                                                                          \
    }                                                                                       \
                                                                                            \
    result = PyFloat_AsDouble(res);                                                         \
                                                                                            \
    if (PyErr_Occurred()) {                                                                 \
        error = 1;                                                                          \
        goto done;                                                                          \
    }                                                                                       \
                                                                                            \
done:                                                                                       \
    PyGILState_Release(gstate);                                                             \
    Py_XDECREF(arg1);                                                                       \
    Py_XDECREF(argobj);                                                                     \
    Py_XDECREF(funcname);                                                                   \
    Py_XDECREF(res);                                                                        \
                                                                                            \
    if (error) {                                                                            \
        /* nonlocal return causes memory leaks. So, if the Python error variable has been   \
           set, just return INFINITY. This will cause an error in UNU.RAN and it            \
           will return an errorcode or NULL value. We can then raise the error from Cython  \
           by checking if the Python error variable variable has been set when an           \
           errorcode/NULL value is returned. */                                             \
        return UNUR_INFINITY;                                                               \
    }                                                                                       \
                                                                                            \
    return result

void error_handler(const char *objid, const char *file, int line, const char *errortype, int unur_errno, const char *reason)
{
    if ( unur_errno != UNUR_SUCCESS ) {
        if (PyErr_Occurred()) {
            return;
        }
        FILE *LOG = unur_get_stream();
        char objid_[256], reason_[256];
        (objid == NULL || strcmp(objid, "") == 0) ? strcpy(objid_, "unknown") : strcpy(objid_, objid);
        (reason == NULL || strcmp(reason, "") == 0) ? strcpy(reason_, "unknown error!") : strcpy(reason_, reason);
        const char *errno_msg = unur_get_strerror(unur_errno);
        if ( strcmp(errortype, "error") == 0 ) {
            fprintf(LOG, "[objid: %s] %d : %s => %s", objid_, unur_errno, reason_, errno_msg);
        }
        else {
            PyErr_WarnFormat(PyExc_RuntimeWarning, 1, "[objid: %s] %d : %s => %s", objid_, unur_errno, reason_, errno_msg);
        }
    }
}

static ccallback_signature_t unuran_call_signatures[] = {
    {NULL}
};

int init_unuran_callback(ccallback_t *callback, PyObject *fcn)
{
    int ret;
    int flags = CCALLBACK_OBTAIN;

    ret = ccallback_prepare(callback, unuran_call_signatures, fcn, flags);
    if (ret == -1) {
        return -1;
    }

    callback->info_p = NULL;

    return 0;
}

int release_unuran_callback(ccallback_t *callback) {
    int ret = ccallback_release(callback);
    return ret;
}


/* ********************************************************************************** */
/* ********************************* UNU.RAN Thunks ********************************* */
/* ********************************************************************************** */

double pmf_thunk(int x, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyLong_FromLong, "pmf", 3);
}

double pdf_thunk(double x, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyFloat_FromDouble, "pdf", 3);
}

double dpdf_thunk(double x, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyFloat_FromDouble, "dpdf", 4);
}

double cont_cdf_thunk(double x, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyFloat_FromDouble, "cdf", 3);
}

double discr_cdf_thunk(int x, const struct unur_distr *distr)
{
    UNURAN_THUNK(PyLong_FromLong, "cdf", 3);
}
