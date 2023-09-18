/*
 * This file implements the API for the voidptr type.
 *
 * Copyright (c) 2022 Riverbank Computing Limited <info@riverbankcomputing.com>
 *
 * This file is part of SIP.
 *
 * This copy of SIP is licensed for use under the terms of the SIP License
 * Agreement.  See the file LICENSE for more details.
 *
 * This copy of SIP may also used under the terms of the GNU General Public
 * License v2 or v3 as published by the Free Software Foundation which can be
 * found in the files LICENSE-GPL2 and LICENSE-GPL3 included in this package.
 *
 * SIP is supplied WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 */


#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stddef.h>
#include <string.h>

#include "sipint.h"
#include "sip_array.h"


/* The object data structure. */
typedef struct {
    PyObject_HEAD
    void *voidptr;
    Py_ssize_t size;
    int rw;
} sipVoidPtrObject;


/* The structure used to hold the results of a voidptr conversion. */
struct vp_values {
    void *voidptr;
    Py_ssize_t size;
    int rw;
};


static int check_size(PyObject *self);
static int check_rw(PyObject *self);
static int check_index(PyObject *self, Py_ssize_t idx);
static void bad_key(PyObject *key);
static int check_slice_size(Py_ssize_t size, Py_ssize_t value_size);
static PyObject *make_voidptr(void *voidptr, Py_ssize_t size, int rw);
static int vp_convertor(PyObject *arg, struct vp_values *vp);
static Py_ssize_t get_size_from_arg(sipVoidPtrObject *v, Py_ssize_t size);


/*
 * Implement ascapsule() for the type.
 */
static PyObject *sipVoidPtr_ascapsule(sipVoidPtrObject *v, PyObject *arg)
{
    (void)arg;

    return PyCapsule_New(v->voidptr, NULL, NULL);
}


/*
 * Implement asarray() for the type.
 */
static PyObject *sipVoidPtr_asarray(sipVoidPtrObject *v, PyObject *args,
        PyObject *kw)
{
    static char *kwlist[] = {"size", NULL};

    Py_ssize_t size = -1;

    if (!PyArg_ParseTupleAndKeywords(args, kw, "|n:asarray", kwlist, &size))
        return NULL;

    if ((size = get_size_from_arg(v, size)) < 0)
        return NULL;

    return sip_api_convert_to_array(v->voidptr, "B", size,
            (v->rw ? 0 : SIP_READ_ONLY));
}


/*
 * Implement asstring() for the type.
 */
static PyObject *sipVoidPtr_asstring(sipVoidPtrObject *v, PyObject *args,
        PyObject *kw)
{
    static char *kwlist[] = {"size", NULL};

    Py_ssize_t size = -1;

    if (!PyArg_ParseTupleAndKeywords(args, kw, "|n:asstring", kwlist, &size))
        return NULL;

    if ((size = get_size_from_arg(v, size)) < 0)
        return NULL;

    return PyBytes_FromStringAndSize(v->voidptr, size);
}


/*
 * Implement getsize() for the type.
 */
static PyObject *sipVoidPtr_getsize(sipVoidPtrObject *v, PyObject *arg)
{
    (void)arg;

    return PyLong_FromSsize_t(v->size);
}


/*
 * Implement setsize() for the type.
 */
static PyObject *sipVoidPtr_setsize(sipVoidPtrObject *v, PyObject *arg)
{
    Py_ssize_t size = PyLong_AsSsize_t(arg);

    if (PyErr_Occurred())
        return NULL;

    v->size = size;

    Py_INCREF(Py_None);
    return Py_None;
}


/*
 * Implement getwriteable() for the type.
 */
static PyObject *sipVoidPtr_getwriteable(sipVoidPtrObject *v, PyObject *arg)
{
    (void)arg;

    return PyBool_FromLong(v->rw);
}


/*
 * Implement setwriteable() for the type.
 */
static PyObject *sipVoidPtr_setwriteable(sipVoidPtrObject *v, PyObject *arg)
{
    int rw;

    if ((rw = PyObject_IsTrue(arg)) < 0)
        return NULL;

    v->rw = rw;

    Py_INCREF(Py_None);
    return Py_None;
}


/* The methods data structure. */
static PyMethodDef sipVoidPtr_Methods[] = {
    {"asarray", (PyCFunction)sipVoidPtr_asarray, METH_VARARGS|METH_KEYWORDS, NULL},
    {"ascapsule", (PyCFunction)sipVoidPtr_ascapsule, METH_NOARGS, NULL},
    {"asstring", (PyCFunction)sipVoidPtr_asstring, METH_VARARGS|METH_KEYWORDS, NULL},
    {"getsize", (PyCFunction)sipVoidPtr_getsize, METH_NOARGS, NULL},
    {"setsize", (PyCFunction)sipVoidPtr_setsize, METH_O, NULL},
    {"getwriteable", (PyCFunction)sipVoidPtr_getwriteable, METH_NOARGS, NULL},
    {"setwriteable", (PyCFunction)sipVoidPtr_setwriteable, METH_O, NULL},
    {NULL, NULL, 0, NULL}
};


/*
 * Implement bool() for the type.
 */
static int sipVoidPtr_bool(PyObject *self)
{
    return (((sipVoidPtrObject *)self)->voidptr != NULL);
}


/*
 * Implement int() for the type.
 */
static PyObject *sipVoidPtr_int(PyObject *self)
{
    return PyLong_FromVoidPtr(((sipVoidPtrObject *)self)->voidptr);
}


/* The number methods data structure. */
static PyNumberMethods sipVoidPtr_NumberMethods = {
    0,                      /* nb_add */
    0,                      /* nb_subtract */
    0,                      /* nb_multiply */
    0,                      /* nb_remainder */
    0,                      /* nb_divmod */
    0,                      /* nb_power */
    0,                      /* nb_negative */
    0,                      /* nb_positive */
    0,                      /* nb_absolute */
    sipVoidPtr_bool,        /* nb_bool */
    0,                      /* nb_invert */
    0,                      /* nb_lshift */
    0,                      /* nb_rshift */
    0,                      /* nb_and */
    0,                      /* nb_xor */
    0,                      /* nb_or */
    sipVoidPtr_int,         /* nb_int */
    0,                      /* nb_reserved */
    0,                      /* nb_float */
    0,                      /* nb_inplace_add */
    0,                      /* nb_inplace_subtract */
    0,                      /* nb_inplace_multiply */
    0,                      /* nb_inplace_remainder */
    0,                      /* nb_inplace_power */
    0,                      /* nb_inplace_lshift */
    0,                      /* nb_inplace_rshift */
    0,                      /* nb_inplace_and */
    0,                      /* nb_inplace_xor */
    0,                      /* nb_inplace_or */
    0,                      /* nb_floor_divide */
    0,                      /* nb_true_divide */
    0,                      /* nb_inplace_floor_divide */
    0,                      /* nb_inplace_true_divide */
    0,                      /* nb_index */
    0,                      /* nb_matrix_multiply */
    0,                      /* nb_inplace_matrix_multiply */
};


/*
 * Implement len() for the type.
 */
static Py_ssize_t sipVoidPtr_length(PyObject *self)
{
    if (check_size(self) < 0)
        return -1;

    return ((sipVoidPtrObject *)self)->size;
}


/*
 * Implement sequence item sub-script for the type.
 */
static PyObject *sipVoidPtr_item(PyObject *self, Py_ssize_t idx)
{
    if (check_size(self) < 0 || check_index(self, idx) < 0)
        return NULL;

    return PyBytes_FromStringAndSize(
            (char *)((sipVoidPtrObject *)self)->voidptr + idx, 1);
}


/* The sequence methods data structure. */
static PySequenceMethods sipVoidPtr_SequenceMethods = {
    sipVoidPtr_length,      /* sq_length */
    0,                      /* sq_concat */
    0,                      /* sq_repeat */
    sipVoidPtr_item,        /* sq_item */
    0,                      /* sq_slice */
    0,                      /* sq_ass_item */
    0,                      /* sq_ass_slice */
    0,                      /* sq_contains */
    0,                      /* sq_inplace_concat */
    0,                      /* sq_inplace_repeat */
};


/*
 * Implement mapping sub-script for the type.
 */
static PyObject *sipVoidPtr_subscript(PyObject *self, PyObject *key)
{
    sipVoidPtrObject *v;

    if (check_size(self) < 0)
        return NULL;

    v = (sipVoidPtrObject *)self;

    if (PyIndex_Check(key))
    {
        Py_ssize_t idx = PyNumber_AsSsize_t(key, PyExc_IndexError);

        if (idx == -1 && PyErr_Occurred())
            return NULL;

        if (idx < 0)
            idx += v->size;

        return sipVoidPtr_item(self, idx);
    }

    if (PySlice_Check(key))
    {
        Py_ssize_t start, stop, step, slicelength;

        if (sip_api_convert_from_slice_object(key, v->size, &start, &stop, &step, &slicelength) < 0)
            return NULL;

        if (step != 1)
        {
            PyErr_SetNone(PyExc_NotImplementedError);
            return NULL;
        }

        return make_voidptr((char *)v->voidptr + start, slicelength, v->rw);
    }

    bad_key(key);

    return NULL;
}


/*
 * Implement mapping assignment sub-script for the type.
 */
static int sipVoidPtr_ass_subscript(PyObject *self, PyObject *key,
        PyObject *value)
{
    sipVoidPtrObject *v;
    Py_ssize_t start, size;
    Py_buffer value_view;

    if (check_rw(self) < 0 || check_size(self) < 0)
        return -1;

    v = (sipVoidPtrObject *)self;

    if (PyIndex_Check(key))
    {
        start = PyNumber_AsSsize_t(key, PyExc_IndexError);

        if (start == -1 && PyErr_Occurred())
            return -1;

        if (start < 0)
            start += v->size;

        if (check_index(self, start) < 0)
            return -1;

        size = 1;
    }
    else if (PySlice_Check(key))
    {
        Py_ssize_t stop, step;

        if (sip_api_convert_from_slice_object(key, v->size, &start, &stop, &step, &size) < 0)
            return -1;

        if (step != 1)
        {
            PyErr_SetNone(PyExc_NotImplementedError);
            return -1;
        }
    }
    else
    {
        bad_key(key);

        return -1;
    }

    if (PyObject_GetBuffer(value, &value_view, PyBUF_CONTIG_RO) < 0)
        return -1;

    /* We could allow any item size... */
    if (value_view.itemsize != 1)
    {
        PyErr_Format(PyExc_TypeError, "'%s' must have an item size of 1",
                Py_TYPE(value_view.obj)->tp_name);

        PyBuffer_Release(&value_view);
        return -1;
    }

    if (check_slice_size(size, value_view.len) < 0)
    {
        PyBuffer_Release(&value_view);
        return -1;
    }

    memmove((char *)v->voidptr + start, value_view.buf, size);

    PyBuffer_Release(&value_view);

    return 0;
}


/* The mapping methods data structure. */
static PyMappingMethods sipVoidPtr_MappingMethods = {
    sipVoidPtr_length,      /* mp_length */
    sipVoidPtr_subscript,   /* mp_subscript */
    sipVoidPtr_ass_subscript,   /* mp_ass_subscript */
};


/*
 * The buffer implementation for Python v2.6.3 and later.
 */
static int sipVoidPtr_getbuffer(PyObject *self, Py_buffer *buf, int flags)
{
    sipVoidPtrObject *v;

    if (check_size(self) < 0)
        return -1;

    v = (sipVoidPtrObject *)self;

    return PyBuffer_FillInfo(buf, self, v->voidptr, v->size, !v->rw, flags);
}


/* The buffer methods data structure. */
static PyBufferProcs sipVoidPtr_BufferProcs = {
    sipVoidPtr_getbuffer,   /* bf_getbuffer */
    0                       /* bf_releasebuffer */
};


/*
 * Implement __new__ for the type.
 */
static PyObject *sipVoidPtr_new(PyTypeObject *subtype, PyObject *args,
        PyObject *kw)
{
    static char *kwlist[] = {"address", "size", "writeable", NULL};

    struct vp_values vp_conversion;
    Py_ssize_t size = -1;
    int rw = -1;
    PyObject *obj;

    if (!PyArg_ParseTupleAndKeywords(args, kw, "O&|ni:voidptr", kwlist, vp_convertor, &vp_conversion, &size, &rw))
        return NULL;

    /* Use the explicit size if one was given. */
    if (size >= 0)
        vp_conversion.size = size;

    /* Use the explicit writeable flag if one was given. */
    if (rw >= 0)
        vp_conversion.rw = rw;

    /* Create the instance. */
    if ((obj = subtype->tp_alloc(subtype, 0)) == NULL)
        return NULL;

    /* Save the values. */
    ((sipVoidPtrObject *)obj)->voidptr = vp_conversion.voidptr;
    ((sipVoidPtrObject *)obj)->size = vp_conversion.size;
    ((sipVoidPtrObject *)obj)->rw = vp_conversion.rw;

    return obj;
}


/* The type data structure. */
PyTypeObject sipVoidPtr_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "sip.voidptr",          /* tp_name */
    sizeof (sipVoidPtrObject),  /* tp_basicsize */
    0,                      /* tp_itemsize */
    0,                      /* tp_dealloc */
    0,                      /* tp_print */
    0,                      /* tp_getattr */
    0,                      /* tp_setattr */
    0,                      /* tp_reserved (Python v3), tp_compare (Python v2) */
    0,                      /* tp_repr */
    &sipVoidPtr_NumberMethods,  /* tp_as_number */
    &sipVoidPtr_SequenceMethods,    /* tp_as_sequence */
    &sipVoidPtr_MappingMethods, /* tp_as_mapping */
    0,                      /* tp_hash */
    0,                      /* tp_call */
    0,                      /* tp_str */
    0,                      /* tp_getattro */
    0,                      /* tp_setattro */
    &sipVoidPtr_BufferProcs,    /* tp_as_buffer */
#if defined(Py_TPFLAGS_HAVE_NEWBUFFER)
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_NEWBUFFER,   /* tp_flags */
#else
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
#endif
    0,                      /* tp_doc */
    0,                      /* tp_traverse */
    0,                      /* tp_clear */
    0,                      /* tp_richcompare */
    0,                      /* tp_weaklistoffset */
    0,                      /* tp_iter */
    0,                      /* tp_iternext */
    sipVoidPtr_Methods,     /* tp_methods */
    0,                      /* tp_members */
    0,                      /* tp_getset */
    0,                      /* tp_base */
    0,                      /* tp_dict */
    0,                      /* tp_descr_get */
    0,                      /* tp_descr_set */
    0,                      /* tp_dictoffset */
    0,                      /* tp_init */
    0,                      /* tp_alloc */
    sipVoidPtr_new,         /* tp_new */
    0,                      /* tp_free */
    0,                      /* tp_is_gc */
    0,                      /* tp_bases */
    0,                      /* tp_mro */
    0,                      /* tp_cache */
    0,                      /* tp_subclasses */
    0,                      /* tp_weaklist */
    0,                      /* tp_del */
    0,                      /* tp_version_tag */
    0,                      /* tp_finalize */
#if PY_VERSION_HEX >= 0x03080000
    0,                      /* tp_vectorcall */
#endif
};


/*
 * A convenience function to convert a C/C++ void pointer from a Python object.
 */
void *sip_api_convert_to_void_ptr(PyObject *obj)
{
    struct vp_values vp;

    if (obj == NULL)
    {
        PyErr_SetString(PyExc_TypeError, "sip.voidptr is NULL");
        return NULL;
    }

    if (vp_convertor(obj, &vp))
        return vp.voidptr;

    return PyLong_AsVoidPtr(obj);
}


/*
 * Convert a C/C++ void pointer to a sip.voidptr object.
 */
PyObject *sip_api_convert_from_void_ptr(void *val)
{
    return make_voidptr(val, -1, TRUE);
}


/*
 * Convert a C/C++ void pointer to a sip.voidptr object.
 */
PyObject *sip_api_convert_from_const_void_ptr(const void *val)
{
    return make_voidptr((void *)val, -1, FALSE);
}


/*
 * Convert a sized C/C++ void pointer to a sip.voidptr object.
 */
PyObject *sip_api_convert_from_void_ptr_and_size(void *val, Py_ssize_t size)
{
    return make_voidptr(val, size, TRUE);
}


/*
 * Convert a sized C/C++ const void pointer to a sip.voidptr object.
 */
PyObject *sip_api_convert_from_const_void_ptr_and_size(const void *val,
        Py_ssize_t size)
{
    return make_voidptr((void *)val, size, FALSE);
}


/*
 * Check that a void pointer has an explicit size and raise an exception if it
 * hasn't.
 */
static int check_size(PyObject *self)
{
    if (((sipVoidPtrObject *)self)->size >= 0)
        return 0;

    PyErr_SetString(PyExc_IndexError,
            "sip.voidptr object has an unknown size");

    return -1;
}


/*
 * Check that a void pointer is writable.
 */
static int check_rw(PyObject *self)
{
    if (((sipVoidPtrObject *)self)->rw)
        return 0;

    PyErr_SetString(PyExc_TypeError,
            "cannot modify a read-only sip.voidptr object");

    return -1;
}


/*
 * Check that an index is valid for a void pointer.
 */
static int check_index(PyObject *self, Py_ssize_t idx)
{
    if (idx >= 0 && idx < ((sipVoidPtrObject *)self)->size)
        return 0;

    PyErr_SetString(PyExc_IndexError, "index out of bounds");

    return -1;
}


/*
 * Raise an exception about a bad sub-script key.
 */
static void bad_key(PyObject *key)
{
    PyErr_Format(PyExc_TypeError,
            "cannot index a sip.voidptr object using '%s'",
            Py_TYPE(key)->tp_name);
}


/*
 * Check that the size of a value is the same as the size of the slice it is
 * replacing.
 */
static int check_slice_size(Py_ssize_t size, Py_ssize_t value_size)
{
    if (value_size == size)
        return 0;

    PyErr_SetString(PyExc_ValueError,
            "cannot modify the size of a sip.voidptr object");

    return -1;
}


/*
 * Do the work of converting a void pointer.
 */
static PyObject *make_voidptr(void *voidptr, Py_ssize_t size, int rw)
{
    sipVoidPtrObject *self;

    if (voidptr == NULL)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }

    if ((self = PyObject_NEW(sipVoidPtrObject, &sipVoidPtr_Type)) == NULL)
        return NULL;

    self->voidptr = voidptr;
    self->size = size;
    self->rw = rw;

    return (PyObject *)self;
}


/*
 * Convert a Python object to the values needed to create a voidptr.
 */
static int vp_convertor(PyObject *arg, struct vp_values *vp)
{
    void *ptr;
    Py_ssize_t size = -1;
    int rw = TRUE;

    if (arg == Py_None)
    {
        ptr = NULL;
    }
    else if (PyCapsule_CheckExact(arg))
    {
        ptr = PyCapsule_GetPointer(arg, NULL);
    }
    else if (PyObject_TypeCheck(arg, &sipVoidPtr_Type))
    {
        ptr = ((sipVoidPtrObject *)arg)->voidptr;
        size = ((sipVoidPtrObject *)arg)->size;
        rw = ((sipVoidPtrObject *)arg)->rw;
    }
    else if (PyObject_CheckBuffer(arg))
    {
        Py_buffer view;

        if (PyObject_GetBuffer(arg, &view, PyBUF_SIMPLE) < 0)
            return 0;

        ptr = view.buf;
        size = view.len;
        rw = !view.readonly;

        PyBuffer_Release(&view);
    }
    else
    {
        PyErr_Clear();
        ptr = PyLong_AsVoidPtr(arg);

        if (PyErr_Occurred())
        {
            PyErr_SetString(PyExc_TypeError, "a single integer, Capsule, None, bytes-like object or another sip.voidptr object is required");
            return 0;
        }
    }

    vp->voidptr = ptr;
    vp->size = size;
    vp->rw = rw;

    return 1;
}


/*
 * Get a size possibly supplied as an argument, otherwise get it from the
 * object.  Raise an exception if there was no size specified.
 */
static Py_ssize_t get_size_from_arg(sipVoidPtrObject *v, Py_ssize_t size)
{
    /* Use the current size if one wasn't explicitly given. */
    if (size < 0)
        size = v->size;

    if (size < 0)
        PyErr_SetString(PyExc_ValueError,
                "a size must be given or the sip.voidptr object must have a size");

    return size;
}
