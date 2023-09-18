/*
 * This file defines the SIP library internal interfaces.
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


#ifndef _SIPINT_H
#define _SIPINT_H


#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include "sip.h"


#ifdef __cplusplus
extern "C" {
#endif

#undef  TRUE
#define TRUE        1

#undef  FALSE
#define FALSE       0


/*
 * This defines a single entry in an object map's hash table.
 */
typedef struct
{
    void *key;                  /* The C/C++ address. */
    sipSimpleWrapper *first;    /* The first object at this address. */
} sipHashEntry;


/*
 * This defines the interface to a hash table class for mapping C/C++ addresses
 * to the corresponding wrapped Python object.
 */
typedef struct
{
    int primeIdx;               /* Index into table sizes. */
    unsigned long size;         /* Size of hash table. */
    unsigned long unused;       /* Nr. unused in hash table. */
    unsigned long stale;        /* Nr. stale in hash table. */
    sipHashEntry *hash_array;   /* Current hash table. */
} sipObjectMap;


/*
 * Support for the descriptors.
 */
extern PyTypeObject sipMethodDescr_Type;
PyObject *sipMethodDescr_New(PyMethodDef *pmd);
PyObject *sipMethodDescr_Copy(PyObject *orig, PyObject *mixin_name);

extern PyTypeObject sipVariableDescr_Type;
PyObject *sipVariableDescr_New(sipVariableDef *vd, const sipTypeDef *td,
    const sipContainerDef *cod);
PyObject *sipVariableDescr_Copy(PyObject *orig, PyObject *mixin_name);


/*
 * Support for API versions.
 */
PyObject *sipGetAPI(PyObject *self, PyObject *args);
PyObject *sipSetAPI(PyObject *self, PyObject *args);
int sip_api_is_api_enabled(const char *name, int from, int to);
int sipIsRangeEnabled(sipExportedModuleDef *em, int range_index);
int sipInitAPI(sipExportedModuleDef *em, PyObject *mod_dict);


/*
 * Support for void pointers.
 */
extern PyTypeObject sipVoidPtr_Type;
void *sip_api_convert_to_void_ptr(PyObject *obj);
PyObject *sip_api_convert_from_void_ptr(void *val);
PyObject *sip_api_convert_from_const_void_ptr(const void *val);
PyObject *sip_api_convert_from_void_ptr_and_size(void *val, Py_ssize_t size);
PyObject *sip_api_convert_from_const_void_ptr_and_size(const void *val,
        Py_ssize_t size);


/*
 * Support for int convertors.
 */
PyObject *sipEnableOverflowChecking(PyObject *self, PyObject *args);
int sip_api_enable_overflow_checking(int enable);
int sip_api_convert_to_bool(PyObject *o);
char sip_api_long_as_char(PyObject *o);
signed char sip_api_long_as_signed_char(PyObject *o);
unsigned char sip_api_long_as_unsigned_char(PyObject *o);
short sip_api_long_as_short(PyObject *o);
unsigned short sip_api_long_as_unsigned_short(PyObject *o);
int sip_api_long_as_int(PyObject *o);
unsigned int sip_api_long_as_unsigned_int(PyObject *o);
long sip_api_long_as_long(PyObject *o);
unsigned long sip_api_long_as_unsigned_long(PyObject *o);
#if defined(HAVE_LONG_LONG)
PY_LONG_LONG sip_api_long_as_long_long(PyObject *o);
unsigned PY_LONG_LONG sip_api_long_as_unsigned_long_long(PyObject *o);
#endif
size_t sip_api_long_as_size_t(PyObject *o);


extern sipQtAPI *sipQtSupport;  /* The Qt support API. */
extern PyTypeObject sipWrapperType_Type;        /* The wrapper type type. */
extern sipWrapperType sipSimpleWrapper_Type;    /* The simple wrapper type. */
extern sipTypeDef *sipQObjectType;  /* The QObject type. */

void *sipGetRx(sipSimpleWrapper *txSelf, const char *sigargs, PyObject *rxObj,
        const char *slot, const char **memberp);
PyObject *sip_api_connect_rx(PyObject *txObj, const char *sig, PyObject *rxObj,
        const char *slot, int type);
PyObject *sip_api_disconnect_rx(PyObject *txObj, const char *sig,
        PyObject *rxObj,const char *slot);


/*
 * These are part of the SIP API but are also used within the SIP module.
 */
void *sip_api_malloc(size_t nbytes);
void sip_api_free(void *mem);
void *sip_api_get_address(sipSimpleWrapper *w);
void *sip_api_get_cpp_ptr(sipSimpleWrapper *w, const sipTypeDef *td);
PyObject *sip_api_convert_from_type(void *cppPtr, const sipTypeDef *td,
        PyObject *transferObj);
void sip_api_instance_destroyed(sipSimpleWrapper *sipSelf);
void sip_api_end_thread(void);
void *sip_api_force_convert_to_type(PyObject *pyObj, const sipTypeDef *td,
        PyObject *transferObj, int flags, int *statep, int *iserrp);
void sip_api_free_sipslot(sipSlot *slot);
int sip_api_same_slot(const sipSlot *sp, PyObject *rxObj, const char *slot);
PyObject *sip_api_invoke_slot(const sipSlot *slot, PyObject *sigargs);
PyObject *sip_api_invoke_slot_ex(const sipSlot *slot, PyObject *sigargs,
        int no_receiver_check);
void *sip_api_convert_rx(sipWrapper *txSelf, const char *sigargs,
        PyObject *rxObj, const char *slot, const char **memberp, int flags);
int sip_api_save_slot(sipSlot *sp, PyObject *rxObj, const char *slot);
int sip_api_convert_from_slice_object(PyObject *slice, Py_ssize_t length,
        Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step,
        Py_ssize_t *slicelength);
int sip_api_deprecated(const char *classname, const char *method);


/*
 * These are not part of the SIP API but are used within the SIP module.
 */
sipClassTypeDef *sipGetGeneratedClassType(const sipEncodedTypeDef *enc,
        const sipClassTypeDef *ctd);
void sipSaveMethod(sipPyMethod *pm,PyObject *meth);
int sipGetPending(void **pp, sipWrapper **op, int *fp);
int sipIsPending(void);
PyObject *sipWrapInstance(void *cpp,  PyTypeObject *py_type, PyObject *args,
        sipWrapper *owner, int flags);
void *sipConvertRxEx(sipWrapper *txSelf, const char *sigargs,
        PyObject *rxObj, const char *slot, const char **memberp, int flags);

void sipOMInit(sipObjectMap *om);
void sipOMFinalise(sipObjectMap *om);
sipSimpleWrapper *sipOMFindObject(sipObjectMap *om, void *key,
        const sipTypeDef *td);
void sipOMAddObject(sipObjectMap *om, sipSimpleWrapper *val);
int sipOMRemoveObject(sipObjectMap *om, sipSimpleWrapper *val);

#if defined(__STDC_VERSION__) && __STDC_VERSION__ >= 199901L
#define sipSetBool(p, v)    (*(_Bool *)(p) = (v))
#endif


#ifdef __cplusplus
}
#endif

#endif
