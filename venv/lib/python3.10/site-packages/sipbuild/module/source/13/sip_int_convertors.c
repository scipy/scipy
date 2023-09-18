/*
 * The implementation of the Python object to C/C++ integer convertors.
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

#include <limits.h>

#include "sip_core.h"


static long long long_as_long_long(PyObject *o, long long min, long long max);
static unsigned long long_as_unsigned_long(PyObject *o, unsigned long max);
static void raise_signed_overflow(long long min, long long max);
static void raise_unsigned_overflow(unsigned long long max);


/*
 * Convert a Python object to a C++ bool (returned as an int).
 */
int sip_api_convert_to_bool(PyObject *o)
{
    int v;

    /* Convert the object to an int while checking for overflow. */
    v = sip_api_long_as_int(o);

    if (PyErr_Occurred())
    {
        if (PyErr_ExceptionMatches(PyExc_OverflowError))
        {
            PyErr_Clear();

            /* The value must have been non-zero. */
            v = 1;
        }
        else
        {
            PyErr_Format(PyExc_TypeError, "a 'bool' is expected not '%s'",
                    Py_TYPE(o)->tp_name);

            v = -1;
        }
    }
    else if (v != 0)
    {
        v = 1;
    }

    return v;
}


/*
 * Convert a Python object to a C char.
 */
char sip_api_long_as_char(PyObject *o)
{
    return (char)long_as_long_long(o, CHAR_MIN, CHAR_MAX);
}


/*
 * Convert a Python object to a C signed char.
 */
signed char sip_api_long_as_signed_char(PyObject *o)
{
    return (signed char)long_as_long_long(o, SCHAR_MIN, SCHAR_MAX);
}


/*
 * Convert a Python object to a C unsigned char.
 */
unsigned char sip_api_long_as_unsigned_char(PyObject *o)
{
    return (unsigned char)long_as_unsigned_long(o, UCHAR_MAX);
}


/*
 * Convert a Python object to a C short.
 */
short sip_api_long_as_short(PyObject *o)
{
    return (short)long_as_long_long(o, SHRT_MIN, SHRT_MAX);
}


/*
 * Convert a Python object to a C unsigned short.
 */
unsigned short sip_api_long_as_unsigned_short(PyObject *o)
{
    return (unsigned short)long_as_unsigned_long(o, USHRT_MAX);
}


/*
 * Convert a Python object to a C int.
 */
int sip_api_long_as_int(PyObject *o)
{
    return (int)long_as_long_long(o, INT_MIN, INT_MAX);
}


/*
 * Convert a Python object to a C unsigned int.
 */
unsigned sip_api_long_as_unsigned_int(PyObject *o)
{
    return (unsigned)long_as_unsigned_long(o, UINT_MAX);
}


/*
 * Convert a Python object to a C size_t.
 */
size_t sip_api_long_as_size_t(PyObject *o)
{
    return (size_t)long_as_unsigned_long(o, SIZE_MAX);
}


/*
 * Convert a Python object to a C long.
 */
long sip_api_long_as_long(PyObject *o)
{
    return (long)long_as_long_long(o, LONG_MIN, LONG_MAX);
}


/*
 * Convert a Python object to a C unsigned long.
 */
unsigned long sip_api_long_as_unsigned_long(PyObject *o)
{
    return long_as_unsigned_long(o, ULONG_MAX);
}


/*
 * Convert a Python object to a C long long.
 */
long long sip_api_long_as_long_long(PyObject *o)
{
    return long_as_long_long(o, LLONG_MIN, LLONG_MAX);
}


/*
 * Convert a Python object to a C unsigned long long.
 */
unsigned long long sip_api_long_as_unsigned_long_long(PyObject *o)
{
    unsigned long long value;

    PyErr_Clear();

    value = PyLong_AsUnsignedLongLong(o);

    if (PyErr_Occurred())
    {
        /* Provide a better exception message. */
        if (PyErr_ExceptionMatches(PyExc_OverflowError))
            raise_unsigned_overflow(ULLONG_MAX);
    }

    return value;
}


/*
 * Convert a Python object to a long long checking that the value is within a
 * range if overflow checking is enabled.
 */
static long long long_as_long_long(PyObject *o, long long min, long long max)
{
    long long value;

    PyErr_Clear();

    value = PyLong_AsLongLong(o);

    if (PyErr_Occurred())
    {
        /* Provide a better exception message. */
        if (PyErr_ExceptionMatches(PyExc_OverflowError))
            raise_signed_overflow(min, max);
    }
    else if (value < min || value > max)
    {
        raise_signed_overflow(min, max);
    }

    return value;
}


/*
 * Convert a Python object to an unsigned long checking that the value is
 * within a range if overflow checking is enabled.
 */
static unsigned long long_as_unsigned_long(PyObject *o, unsigned long max)
{
    unsigned long value;

    PyErr_Clear();

    value = PyLong_AsUnsignedLong(o);

    if (PyErr_Occurred())
    {
        /* Provide a better exception message. */
        if (PyErr_ExceptionMatches(PyExc_OverflowError))
            raise_unsigned_overflow(max);
    }
    else if (value > max)
    {
        raise_unsigned_overflow(max);
    }

    return value;
}


/*
 * Raise an overflow exception if a signed value is out of range.
 */
static void raise_signed_overflow(long long min, long long max)
{
    PyErr_Format(PyExc_OverflowError,
            "value must be in the range %lld to %lld", min, max);
}


/*
 * Raise an overflow exception if an unsigned value is out of range.
 */
static void raise_unsigned_overflow(unsigned long long max)
{
    PyErr_Format(PyExc_OverflowError, "value must be in the range 0 to %llu",
            max);
}
