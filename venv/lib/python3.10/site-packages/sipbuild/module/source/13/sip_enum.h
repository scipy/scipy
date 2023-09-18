/*
 * This file defines the API for the enum support.
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


#ifndef _SIP_ENUM_H
#define _SIP_ENUM_H


#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "sip.h"


#ifdef __cplusplus
extern "C" {
#endif


PyObject *sip_api_convert_from_enum(int member, const sipTypeDef *td);
int sip_api_convert_to_enum(PyObject *obj, const sipTypeDef *td);
int sip_api_is_enum_flag(PyObject *obj);

int sip_enum_create(sipExportedModuleDef *client, sipEnumTypeDef *etd,
        sipIntInstanceDef **next_int_p, PyObject *dict);
const sipTypeDef *sip_enum_get_generated_type(PyObject *obj);
int sip_enum_init(void);
int sip_enum_is_enum(PyObject *obj);


#ifdef __cplusplus
}
#endif

#endif
