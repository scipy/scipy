/*
 * This file implements the API for the array type.
 *
 * Copyright (c) 2023 Riverbank Computing Limited <info@riverbankcomputing.com>
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

#include "sip_core.h"

#include "sip_array.h"


/* The object data structure. */
typedef struct {
    PyObject_HEAD
    void *data;
    const sipTypeDef *td;
    const char *format;
    size_t stride;
    Py_ssize_t len;
    int flags;
    PyObject *owner;
} sipArrayObject;


static void bad_key(PyObject *key);
static int check_index(sipArrayObject *array, Py_ssize_t idx);
static int check_writable(sipArrayObject *array);
static PyObject *create_array(void *data, const sipTypeDef *td,
        const char *format, size_t stride, Py_ssize_t len, int flags,
        PyObject *owner);
static void *element(sipArrayObject *array, Py_ssize_t idx);
static void *get_slice(sipArrayObject *array, PyObject *value, Py_ssize_t len);
static const char *get_type_name(sipArrayObject *array);
static void *get_value(sipArrayObject *array, PyObject *value);
static void init_array(sipArrayObject *array, void *data, const sipTypeDef *td,
        const char *format, size_t stride, Py_ssize_t len, int flags,
        PyObject *owner);


/*
 * Implement len() for the type.
 */
static Py_ssize_t sipArray_length(PyObject *self)
{
    return ((sipArrayObject *)self)->len;
}


/*
 * Implement sequence item sub-script for the type.
 */
static PyObject *sipArray_item(PyObject *self, Py_ssize_t idx)
{
    sipArrayObject *array = (sipArrayObject *)self;
    PyObject *py_item;
    void *data;

    if (check_index(array, idx) < 0)
        return NULL;

    data = element(array, idx);

    if (array->td != NULL)
    {
        py_item = sip_api_convert_from_type(data, array->td, NULL);
    }
    else
    {
        switch (*array->format)
        {
        case 'b':
            py_item = PyLong_FromLong(*(char *)data);
            break;

        case 'B':
            py_item = PyLong_FromUnsignedLong(*(unsigned char *)data);
            break;

        case 'h':
            py_item = PyLong_FromLong(*(short *)data);
            break;

        case 'H':
            py_item = PyLong_FromUnsignedLong(*(unsigned short *)data);
            break;

        case 'i':
            py_item = PyLong_FromLong(*(int *)data);
            break;

        case 'I':
            py_item = PyLong_FromUnsignedLong(*(unsigned int *)data);
            break;

        case 'f':
            py_item = PyFloat_FromDouble(*(float *)data);
            break;

        case 'd':
            py_item = PyFloat_FromDouble(*(double *)data);
            break;

        default:
            py_item = NULL;
        }
    }

    return py_item;
}


/* The sequence methods data structure. */
static PySequenceMethods sipArray_SequenceMethods = {
    sipArray_length,        /* sq_length */
    0,                      /* sq_concat */
    0,                      /* sq_repeat */
    sipArray_item,          /* sq_item */
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
static PyObject *sipArray_subscript(PyObject *self, PyObject *key)
{
    sipArrayObject *array = (sipArrayObject *)self;

    if (PyIndex_Check(key))
    {
        Py_ssize_t idx = PyNumber_AsSsize_t(key, PyExc_IndexError);

        if (idx == -1 && PyErr_Occurred())
            return NULL;

        if (idx < 0)
            idx += array->len;

        return sipArray_item(self, idx);
    }

    if (PySlice_Check(key))
    {
        Py_ssize_t start, stop, step, slicelength;

        if (sip_api_convert_from_slice_object(key, array->len, &start, &stop, &step, &slicelength) < 0)
            return NULL;

        if (step != 1)
        {
            PyErr_SetNone(PyExc_NotImplementedError);
            return NULL;
        }

        return create_array(element(array, start), array->td, array->format,
                array->stride, slicelength, (array->flags & ~SIP_OWNS_MEMORY),
                array->owner);
    }

    bad_key(key);

    return NULL;
}


/*
 * Implement mapping assignment sub-script for the type.
 */
static int sipArray_ass_subscript(PyObject *self, PyObject *key,
        PyObject *value)
{
    sipArrayObject *array = (sipArrayObject *)self;
    Py_ssize_t start, len;
    void *value_data;

    if (check_writable(array) < 0)
        return -1;

    if (PyIndex_Check(key))
    {
        start = PyNumber_AsSsize_t(key, PyExc_IndexError);

        if (start == -1 && PyErr_Occurred())
            return -1;

        if (start < 0)
            start += array->len;

        if (check_index(array, start) < 0)
            return -1;

        if ((value_data = get_value(array, value)) == NULL)
            return -1;

        len = 1;
    }
    else if (PySlice_Check(key))
    {
        Py_ssize_t stop, step;

        if (sip_api_convert_from_slice_object(key, array->len, &start, &stop, &step, &len) < 0)
            return -1;

        if (step != 1)
        {
            PyErr_SetNone(PyExc_NotImplementedError);
            return -1;
        }

        if ((value_data = get_slice(array, value, len)) == NULL)
            return -1;
    }
    else
    {
        bad_key(key);

        return -1;
    }

    if (array->td != NULL)
    {
        const sipClassTypeDef *ctd = (const sipClassTypeDef *)(array->td);
        sipAssignFunc assign;
        Py_ssize_t i;

        if ((assign = ctd->ctd_assign) == NULL)
        {
            PyErr_Format(PyExc_TypeError,
                    "a " _SIP_MODULE_FQ_NAME ".array cannot copy '%s'",
                Py_TYPE(self)->tp_name);
            return -1;
        }

        for (i = 0; i < len; ++i)
        {
            assign(array->data, start + i, value_data);
            value_data = (char *)value_data + array->stride;
        }
    }
    else
    {
        memmove(element(array, start), value_data, len * array->stride);
    }

    return 0;
}


/* The mapping methods data structure. */
static PyMappingMethods sipArray_MappingMethods = {
    sipArray_length,        /* mp_length */
    sipArray_subscript,     /* mp_subscript */
    sipArray_ass_subscript, /* mp_ass_subscript */
};


/*
 * The buffer implementation.
 */
static int sipArray_getbuffer(PyObject *self, Py_buffer *view, int flags)
{
    sipArrayObject *array = (sipArrayObject *)self;
    const char *format;
    Py_ssize_t itemsize;

    if (view == NULL)
        return 0;

    if (((flags & PyBUF_WRITABLE) == PyBUF_WRITABLE) && (array->flags & SIP_READ_ONLY))
    {
        PyErr_SetString(PyExc_BufferError, "object is not writable");
        return -1;
    }

    view->obj = self;
    Py_INCREF(self);

    /*
     * If there is no format, ie. it is a wrapped type, then present it as
     * bytes.
     */
    if ((format = array->format) == NULL)
    {
        format = "B";
        itemsize = sizeof (unsigned char);
    }
    else
    {
        itemsize = array->stride;
    }

    view->buf = array->data;
    view->len = array->len * array->stride;
    view->readonly = (array->flags & SIP_READ_ONLY);
    view->itemsize = itemsize;

    view->format = NULL;
    if ((flags & PyBUF_FORMAT) == PyBUF_FORMAT)
        view->format = format;

    view->ndim = 1;

    view->shape = NULL;
    if ((flags & PyBUF_ND) == PyBUF_ND)
        view->shape = &view->len;

    view->strides = NULL;
    if ((flags & PyBUF_STRIDES) == PyBUF_STRIDES)
        view->strides = &view->itemsize;

    view->suboffsets = NULL;
    view->internal = NULL;

    return 0;
}


/* The buffer methods data structure. */
static PyBufferProcs sipArray_BufferProcs = {
    sipArray_getbuffer,     /* bf_getbuffer */
    0                       /* bf_releasebuffer */
};


/*
 * The instance deallocation function.
 */
static void sipArray_dealloc(PyObject *self)
{
    sipArrayObject *array = (sipArrayObject *)self;

    if (array->flags & SIP_OWNS_MEMORY)
    {
        if (array->td != NULL)
            ((const sipClassTypeDef *)(array->td))->ctd_array_delete(array->data);
        else
            PyMem_Free(array->data);
    }
    else
    {
        Py_XDECREF(array->owner);
    }
}


/*
 * Implement __repr__ for the type.
 */
static PyObject *sipArray_repr(PyObject *self)
{
    sipArrayObject *array = (sipArrayObject *)self;

    return PyUnicode_FromFormat(_SIP_MODULE_FQ_NAME ".array(%s, %zd)",
            get_type_name(array), array->len);
}


/*
 * Implement __new__ for the type.
 */
static PyObject *sipArray_new(PyTypeObject *cls, PyObject *args, PyObject *kw)
{
    static char *kwlist[] = {"", "", NULL};

    Py_ssize_t length;
    PyObject *array, *type;
    const sipClassTypeDef *ctd;

    if (!PyArg_ParseTupleAndKeywords(args, kw, "O!n:array", kwlist, &sipWrapperType_Type, &type, &length))
        return NULL;

    ctd = (const sipClassTypeDef *)((sipWrapperType *)type)->wt_td;

    /* We require the array delete helper which was added in ABI v13.4. */
    if (ctd->ctd_base.td_module->em_api_minor < 4)
    {
        PyErr_SetString(PyExc_TypeError,
                "a " _SIP_MODULE_FQ_NAME ".array can only be created for types using ABI v13.4 or later");
        return NULL;
    }

    if (ctd->ctd_array == NULL || ctd->ctd_sizeof == 0)
    {
        PyErr_Format(PyExc_TypeError,
                "a " _SIP_MODULE_FQ_NAME ".array cannot be created for '%s'",
                Py_TYPE(type)->tp_name);
        return NULL;
    }

    if (length < 0)
    {
        PyErr_SetString(PyExc_ValueError,
                "a " _SIP_MODULE_FQ_NAME ".array length cannot be negative");
        return NULL;
    }

    /* Create the instance. */
    if ((array = cls->tp_alloc(cls, 0)) == NULL)
        return NULL;

    init_array((sipArrayObject *)array, ctd->ctd_array(length), &ctd->ctd_base,
            NULL, ctd->ctd_sizeof, length, SIP_OWNS_MEMORY, NULL);

    return array;
}


/* The type data structure. */
PyTypeObject sipArray_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    _SIP_MODULE_FQ_NAME ".array",   /* tp_name */
    sizeof (sipArrayObject),    /* tp_basicsize */
    0,                      /* tp_itemsize */
    sipArray_dealloc,       /* tp_dealloc */
    0,                      /* tp_print */
    0,                      /* tp_getattr */
    0,                      /* tp_setattr */
    0,                      /* tp_reserved */
    sipArray_repr,          /* tp_repr */
    0,                      /* tp_as_number */
    &sipArray_SequenceMethods,  /* tp_as_sequence */
    &sipArray_MappingMethods,   /* tp_as_mapping */
    0,                      /* tp_hash */
    0,                      /* tp_call */
    0,                      /* tp_str */
    0,                      /* tp_getattro */
    0,                      /* tp_setattro */
    &sipArray_BufferProcs,  /* tp_as_buffer */
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
    0,                      /* tp_methods */
    0,                      /* tp_members */
    0,                      /* tp_getset */
    0,                      /* tp_base */
    0,                      /* tp_dict */
    0,                      /* tp_descr_get */
    0,                      /* tp_descr_set */
    0,                      /* tp_dictoffset */
    0,                      /* tp_init */
    0,                      /* tp_alloc */
    sipArray_new,           /* tp_new */
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
 * Return TRUE if an object is a sip.array with elements of a given type.
 */
int sip_array_can_convert(PyObject *obj, const sipTypeDef *td)
{
    if (!PyObject_TypeCheck(obj, &sipArray_Type))
        return FALSE;

    return (((sipArrayObject *)obj)->td == td);
}


/*
 * Return the address and number of elements of a sip.array for which
 * sip_array_can_convert has already returned TRUE.
 */
void sip_array_convert(PyObject *obj, void **data, Py_ssize_t *size)
{
    sipArrayObject *array = (sipArrayObject *)obj;

    *data = array->data;
    *size = array->len;
}


/*
 * Check that an array is writable.
 */
static int check_writable(sipArrayObject *array)
{
    if (array->flags & SIP_READ_ONLY)
    {
        PyErr_SetString(PyExc_TypeError,
                _SIP_MODULE_FQ_NAME ".array object is read-only");
        return -1;
    }

    return 0;
}


/*
 * Check that an index is valid for an array.
 */
static int check_index(sipArrayObject *array, Py_ssize_t idx)
{
    if (idx >= 0 && idx < array->len)
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
            "cannot index a " _SIP_MODULE_FQ_NAME ".array object using '%s'",
            Py_TYPE(key)->tp_name);
}


/*
 * Get the address of an element of an array.
 */
static void *element(sipArrayObject *array, Py_ssize_t idx)
{
    return (unsigned char *)(array->data) + idx * array->stride;
}


/*
 * Get the address of a value that will be copied to an array.
 */
static void *get_value(sipArrayObject *array, PyObject *value)
{
    static union {
        signed char s_char_t;
        unsigned char u_char_t;
        signed short s_short_t;
        unsigned short u_short_t;
        signed int s_int_t;
        unsigned int u_int_t;
        float float_t;
        double double_t;
    } static_data;

    void *data;

    if (array->td != NULL)
    {
        int iserr = FALSE;

        data = sip_api_force_convert_to_type_us(value, array->td, NULL,
                SIP_NOT_NONE|SIP_NO_CONVERTORS, NULL, NULL, &iserr);
    }
    else
    {
        PyErr_Clear();

        switch (*array->format)
        {
        case 'b':
            static_data.s_char_t = sip_api_long_as_char(value);
            data = &static_data.s_char_t;
            break;

        case 'B':
            static_data.u_char_t = sip_api_long_as_unsigned_char(value);
            data = &static_data.u_char_t;
            break;

        case 'h':
            static_data.s_short_t = sip_api_long_as_short(value);
            data = &static_data.s_short_t;
            break;

        case 'H':
            static_data.u_short_t = sip_api_long_as_unsigned_short(value);
            data = &static_data.u_short_t;
            break;

        case 'i':
            static_data.s_int_t = sip_api_long_as_int(value);
            data = &static_data.s_int_t;
            break;

        case 'I':
            static_data.u_int_t = sip_api_long_as_unsigned_int(value);
            data = &static_data.u_int_t;
            break;

        case 'f':
            static_data.float_t = (float)PyFloat_AsDouble(value);
            data = &static_data.float_t;
            break;

        case 'd':
            static_data.double_t = PyFloat_AsDouble(value);
            data = &static_data.double_t;
            break;

        default:
            data = NULL;
        }

        if (PyErr_Occurred())
            data = NULL;
    }

    return data;
}


/*
 * Get the address of an value that will be copied to an array slice.
 */
static void *get_slice(sipArrayObject *array, PyObject *value, Py_ssize_t len)
{
    sipArrayObject *other = (sipArrayObject *)value;

    if (!PyObject_IsInstance(value, (PyObject *)&sipArray_Type) || array->td != other->td || strcmp(array->format, other->format) != 0)
    {
        PyErr_Format(PyExc_TypeError,
                "can only assign another array of %s to the slice",
                get_type_name(array));
        return NULL;
    }

    if (other->len != len)
    {
        PyErr_Format(PyExc_TypeError,
                "the array being assigned must have length %zd", len);
        return NULL;
    }

    if (other->stride == array->stride)
    {
        PyErr_Format(PyExc_TypeError,
                "the array being assigned must have stride %zu",
                array->stride);
        return NULL;
    }

    return other->data;
}


/*
 * Get the name of the type of an element of an array.
 */
static const char *get_type_name(sipArrayObject *array)
{
    const char *type_name;

    if (array->td != NULL)
    {
        type_name = sipTypeName(array->td);
    }
    else
    {
        switch (*array->format)
        {
        case 'b':
            type_name = "char";
            break;

        case 'B':
            type_name = "unsigned char";
            break;

        case 'h':
            type_name = "short";
            break;

        case 'H':
            type_name = "unsigned short";
            break;

        case 'i':
            type_name = "int";
            break;

        case 'I':
            type_name = "unsigned int";
            break;

        case 'f':
            type_name = "float";
            break;

        case 'd':
            type_name = "double";
            break;

        default:
            type_name = "";
        }
    }

    return type_name;
}



/*
 * Create an array for the C API.
 */
static PyObject *create_array(void *data, const sipTypeDef *td,
        const char *format, size_t stride, Py_ssize_t len, int flags,
        PyObject *owner)
{
    sipArrayObject *array;

    if ((array = PyObject_NEW(sipArrayObject, &sipArray_Type)) == NULL)
        return NULL;

    init_array(array, data, td, format, stride, len, flags, owner);

    return (PyObject *)array;
}


/*
 * Initialise an array.
 */
static void init_array(sipArrayObject *array, void *data, const sipTypeDef *td,
        const char *format, size_t stride, Py_ssize_t len, int flags,
        PyObject *owner)
{
    array->data = data;
    array->td = td;
    array->format = format;
    array->stride = stride;
    array->len = len;
    array->flags = flags;

    if (flags & SIP_OWNS_MEMORY)
    {
        /* This is a borrowed reference to itself. */
        array->owner = (PyObject *)array;
    }
    else
    {
        Py_XINCREF(owner);
        array->owner = owner;
    }
}


/*
 * Wrap an array of instances of a fundamental type.  At the moment format must
 * be either "b" (char), "B" (unsigned char), "h" (short), "H" (unsigned
 * short), "i" (int), "I" (unsigned int), "f" (float) or "d" (double).
 */
PyObject *sip_api_convert_to_array(void *data, const char *format,
        Py_ssize_t len, int flags)
{
    size_t stride;

    assert(len >= 0);

    if (data == NULL)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }

    switch (*format)
    {
    case 'b':
        stride = sizeof (char);
        break;

    case 'B':
        stride = sizeof (unsigned char);
        break;

    case 'h':
        stride = sizeof (short);
        break;

    case 'H':
        stride = sizeof (unsigned short);
        break;

    case 'i':
        stride = sizeof (int);
        break;

    case 'I':
        stride = sizeof (unsigned int);
        break;

    case 'f':
        stride = sizeof (float);
        break;

    case 'd':
        stride = sizeof (double);
        break;

    default:
        PyErr_Format(PyExc_ValueError, "'%c' is not a supported format",
                format);
        return NULL;
    }

    return create_array(data, NULL, format, stride, len, flags, NULL);
}


/*
 * Wrap an array of instances of a defined type.
 */
PyObject *sip_api_convert_to_typed_array(void *data, const sipTypeDef *td,
        const char *format, size_t stride, Py_ssize_t len, int flags)
{
    if (data == NULL)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }

    assert(stride > 0);
    assert(len >= 0);

    return create_array(data, td, format, stride, len, flags, NULL);
}
