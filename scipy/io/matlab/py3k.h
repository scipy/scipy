#include <Python.h>
#include <stdio.h>

#if PY_VERSION_HEX < 0x03000000

#include "cStringIO.h"

#define npy_PyFile_Dup(file, mode) PyFile_AsFile(file)
#define npy_PyFile_DupClose(file, handle) (0)
#define npy_PyFile_Check PyFile_Check

#else

/*
 * No-op implementation -- always fall back to the generic one.
 */

static struct PycStringIO_CAPI {
    int(*cread)(PyObject *, char **, Py_ssize_t);
    int(*creadline)(PyObject *, char **);
    int(*cwrite)(PyObject *, const char *, Py_ssize_t);
    PyObject *(*cgetvalue)(PyObject *);
    PyObject *(*NewOutput)(int);
    PyObject *(*NewInput)(PyObject *);
    PyTypeObject *InputType, *OutputType;
} *PycStringIO;

static void PycString_IMPORT() {}

#define PycStringIO_InputCheck(O) 0
#define PycStringIO_OutputCheck(O) 0

/*
 * PyFile_* compatibility
 */

/*
 * Get a FILE* handle to the file represented by the Python object
 */
static FILE*
npy_PyFile_Dup(PyObject *file, char *mode)
{
    int fd, fd2;
    PyObject *ret, *os;
    /* Flush first to ensure things end up in the file in the correct order */
    ret = PyObject_CallMethod(file, "flush", "");
    if (ret == NULL) {
        return NULL;
    }
    Py_DECREF(ret);
    fd = PyObject_AsFileDescriptor(file);
    if (fd == -1) {
        return NULL;
    }
    os = PyImport_ImportModule("os");
    if (os == NULL) {
        return NULL;
    }
    ret = PyObject_CallMethod(os, "dup", "i", fd);
    Py_DECREF(os);
    if (ret == NULL) {
        return NULL;
    }
    fd2 = PyNumber_AsSsize_t(ret, NULL);
    Py_DECREF(ret);
#ifdef _WIN32
    return _fdopen(fd2, mode);
#else
    return fdopen(fd2, mode);
#endif
}

/*
 * Close the dup-ed file handle, and seek the Python one to the current position
 */
static int
npy_PyFile_DupClose(PyObject *file, FILE* handle)
{
    PyObject *ret;
    long position;
    position = ftell(handle);
    fclose(handle);

    ret = PyObject_CallMethod(file, "seek", "li", position, 0);
    if (ret == NULL) {
        return -1;
    }
    Py_DECREF(ret);
    return 0;
}

static int
npy_PyFile_Check(PyObject *file)
{
    int fd;
    fd = PyObject_AsFileDescriptor(file);
    if (fd == -1) {
        PyErr_Clear();
        return 0;
    }
    return 1;
}

#endif
