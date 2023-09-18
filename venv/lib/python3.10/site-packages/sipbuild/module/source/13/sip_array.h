/*
 * This file defines the API for the array type.
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


#ifndef _SIP_ARRAY_H
#define _SIP_ARRAY_H


#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "sip.h"


#ifdef __cplusplus
extern "C" {
#endif


extern PyTypeObject sipArray_Type;

PyObject *sip_api_convert_to_array(void *data, const char *format,
        Py_ssize_t len, int flags);
PyObject *sip_api_convert_to_typed_array(void *data, const sipTypeDef *td,
        const char *format, size_t stride, Py_ssize_t len, int flags);

int sip_array_can_convert(PyObject *obj, const sipTypeDef *td);
void sip_array_convert(PyObject *obj, void **data, Py_ssize_t *size);


#ifdef __cplusplus
}
#endif

#endif
