/*
 * SIP library code.
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
#include <datetime.h>
#include <frameobject.h>

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stddef.h>
#include <string.h>

#include "sip.h"
#include "sipint.h"
#include "sip_array.h"


/* There doesn't seem to be a standard way of checking for C99 support. */
#if !defined(va_copy)
#define va_copy(d, s)   ((d) = (s))
#endif


/*
 * The Python metatype for a C++ wrapper type.  We inherit everything from the
 * standard Python metatype except the init and getattro methods and the size
 * of the type object created is increased to accomodate the extra information
 * we associate with a wrapped type.
 */

static PyObject *sipWrapperType_alloc(PyTypeObject *self, Py_ssize_t nitems);
static PyObject *sipWrapperType_getattro(PyObject *self, PyObject *name);
static int sipWrapperType_init(sipWrapperType *self, PyObject *args,
        PyObject *kwds);
static int sipWrapperType_setattro(PyObject *self, PyObject *name,
        PyObject *value);

PyTypeObject sipWrapperType_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "sip.wrappertype",      /* tp_name */
    sizeof (sipWrapperType),    /* tp_basicsize */
    0,                      /* tp_itemsize */
    0,                      /* tp_dealloc */
    0,                      /* tp_print */
    0,                      /* tp_getattr */
    0,                      /* tp_setattr */
    0,                      /* tp_as_async (Python v3.5), tp_compare (Python v2) */
    0,                      /* tp_repr */
    0,                      /* tp_as_number */
    0,                      /* tp_as_sequence */
    0,                      /* tp_as_mapping */
    0,                      /* tp_hash */
    0,                      /* tp_call */
    0,                      /* tp_str */
    sipWrapperType_getattro,    /* tp_getattro */
    sipWrapperType_setattro,    /* tp_setattro */
    0,                      /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
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
    (initproc)sipWrapperType_init,  /* tp_init */
    sipWrapperType_alloc,   /* tp_alloc */
    0,                      /* tp_new */
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
 * The Python type that is the super-type for all C++ wrapper types that
 * support parent/child relationships.
 */

static int sipWrapper_clear(sipWrapper *self);
static void sipWrapper_dealloc(sipWrapper *self);
static int sipWrapper_traverse(sipWrapper *self, visitproc visit, void *arg);

static sipWrapperType sipWrapper_Type = {
#if !defined(STACKLESS)
    {
#endif
        {
            PyVarObject_HEAD_INIT(&sipWrapperType_Type, 0)
            "sip.wrapper",  /* tp_name */
            sizeof (sipWrapper),    /* tp_basicsize */
            0,              /* tp_itemsize */
            (destructor)sipWrapper_dealloc, /* tp_dealloc */
            0,              /* tp_print */
            0,              /* tp_getattr */
            0,              /* tp_setattr */
            0,              /* tp_as_async (Python v3.5), tp_compare (Python v2) */
            0,              /* tp_repr */
            0,              /* tp_as_number */
            0,              /* tp_as_sequence */
            0,              /* tp_as_mapping */
            0,              /* tp_hash */
            0,              /* tp_call */
            0,              /* tp_str */
            0,              /* tp_getattro */
            0,              /* tp_setattro */
            0,              /* tp_as_buffer */
            Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,  /* tp_flags */
            0,              /* tp_doc */
            (traverseproc)sipWrapper_traverse,  /* tp_traverse */
            (inquiry)sipWrapper_clear,  /* tp_clear */
            0,              /* tp_richcompare */
            0,              /* tp_weaklistoffset */
            0,              /* tp_iter */
            0,              /* tp_iternext */
            0,              /* tp_methods */
            0,              /* tp_members */
            0,              /* tp_getset */
            0,              /* tp_base */
            0,              /* tp_dict */
            0,              /* tp_descr_get */
            0,              /* tp_descr_set */
            0,              /* tp_dictoffset */
            0,              /* tp_init */
            0,              /* tp_alloc */
            0,              /* tp_new */
            0,              /* tp_free */
            0,              /* tp_is_gc */
            0,              /* tp_bases */
            0,              /* tp_mro */
            0,              /* tp_cache */
            0,              /* tp_subclasses */
            0,              /* tp_weaklist */
            0,              /* tp_del */
            0,              /* tp_version_tag */
            0,              /* tp_finalize */
#if PY_VERSION_HEX >= 0x03080000
            0,              /* tp_vectorcall */
#endif
        },
        {
            0,              /* am_await */
            0,              /* am_aiter */
            0,              /* am_anext */
        },
        {
            0,              /* nb_add */
            0,              /* nb_subtract */
            0,              /* nb_multiply */
            0,              /* nb_remainder */
            0,              /* nb_divmod */
            0,              /* nb_power */
            0,              /* nb_negative */
            0,              /* nb_positive */
            0,              /* nb_absolute */
            0,              /* nb_bool */
            0,              /* nb_invert */
            0,              /* nb_lshift */
            0,              /* nb_rshift */
            0,              /* nb_and */
            0,              /* nb_xor */
            0,              /* nb_or */
            0,              /* nb_int */
            0,              /* nb_reserved */
            0,              /* nb_float */
            0,              /* nb_inplace_add */
            0,              /* nb_inplace_subtract */
            0,              /* nb_inplace_multiply */
            0,              /* nb_inplace_remainder */
            0,              /* nb_inplace_power */
            0,              /* nb_inplace_lshift */
            0,              /* nb_inplace_rshift */
            0,              /* nb_inplace_and */
            0,              /* nb_inplace_xor */
            0,              /* nb_inplace_or */
            0,              /* nb_floor_divide */
            0,              /* nb_true_divide */
            0,              /* nb_inplace_floor_divide */
            0,              /* nb_inplace_true_divide */
            0,              /* nb_index */
            0,              /* nb_matrix_multiply */
            0,              /* nb_inplace_matrix_multiply */
        },
        {
            0,              /* mp_length */
            0,              /* mp_subscript */
            0,              /* mp_ass_subscript */
        },
        {
            0,              /* sq_length */
            0,              /* sq_concat */
            0,              /* sq_repeat */
            0,              /* sq_item */
            0,              /* was_sq_slice */
            0,              /* sq_ass_item */
            0,              /* was_sq_ass_slice */
            0,              /* sq_contains */
            0,              /* sq_inplace_concat */
            0,              /* sq_inplace_repeat */
        },
        {
            0,              /* bf_getbuffer */
            0,              /* bf_releasebuffer */
        },
        0,                  /* ht_name */
        0,                  /* ht_slots */
        0,                  /* ht_qualname */
        0,                  /* ht_cached_keys */
#if PY_VERSION_HEX >= 0x03090000
        0,                  /* ht_module */
#endif
#if !defined(STACKLESS)
    },
#endif
    0,                      /* wt_user_type */
    0,                      /* wt_dict_complete */
    0,                      /* wt_unused */
    0,                      /* wt_td */
    0,                      /* wt_iextend */
    0,                      /* wt_new_user_type_handler */
    0,                      /* wt_user_data */
};


static void sip_api_bad_catcher_result(PyObject *method);
static void sip_api_bad_length_for_slice(Py_ssize_t seqlen,
        Py_ssize_t slicelen);
static PyObject *sip_api_build_result(int *isErr, const char *fmt, ...);
static PyObject *sip_api_call_method(int *isErr, PyObject *method,
        const char *fmt, ...);
static void sip_api_call_procedure_method(sip_gilstate_t gil_state,
        sipVirtErrorHandlerFunc error_handler, sipSimpleWrapper *py_self,
        PyObject *method, const char *fmt, ...);
static Py_ssize_t sip_api_convert_from_sequence_index(Py_ssize_t idx,
        Py_ssize_t len);
static int sip_api_can_convert_to_type(PyObject *pyObj, const sipTypeDef *td,
        int flags);
static void *sip_api_convert_to_type(PyObject *pyObj, const sipTypeDef *td,
        PyObject *transferObj, int flags, int *statep, int *iserrp);
static int sip_api_can_convert_to_enum(PyObject *pyObj, const sipTypeDef *td);
static int sip_api_convert_to_enum(PyObject *pyObj, const sipTypeDef *td);
static void sip_api_release_type(void *cpp, const sipTypeDef *td, int state);
static PyObject *sip_api_convert_from_new_type(void *cpp, const sipTypeDef *td,
        PyObject *transferObj);
static PyObject *sip_api_convert_from_new_pytype(void *cpp,
        PyTypeObject *py_type, sipWrapper *owner, sipSimpleWrapper **selfp,
        const char *fmt, ...);
static int sip_api_get_state(PyObject *transferObj);
static PyObject *sip_api_get_pyobject(void *cppPtr, const sipTypeDef *td);
static sipWrapperType *sip_api_map_int_to_class(int typeInt,
        const sipIntTypeClassMap *map, int maplen);
static sipWrapperType *sip_api_map_string_to_class(const char *typeString,
        const sipStringTypeClassMap *map, int maplen);
static int sip_api_parse_result_ex(sip_gilstate_t gil_state,
        sipVirtErrorHandlerFunc error_handler, sipSimpleWrapper *py_self,
        PyObject *method, PyObject *res, const char *fmt, ...);
static int sip_api_parse_result(int *isErr, PyObject *method, PyObject *res,
        const char *fmt, ...);
static void sip_api_call_error_handler(sipVirtErrorHandlerFunc error_handler,
        sipSimpleWrapper *py_self, sip_gilstate_t gil_state);
static void sip_api_trace(unsigned mask,const char *fmt,...);
static void sip_api_transfer_back(PyObject *self);
static void sip_api_transfer_to(PyObject *self, PyObject *owner);
static int sip_api_export_module(sipExportedModuleDef *client,
        unsigned api_major, unsigned api_minor, void *unused);
static int sip_api_init_module(sipExportedModuleDef *client,
        PyObject *mod_dict);
static int sip_api_parse_args(PyObject **parseErrp, PyObject *sipArgs,
        const char *fmt, ...);
static int sip_api_parse_kwd_args(PyObject **parseErrp, PyObject *sipArgs,
        PyObject *sipKwdArgs, const char **kwdlist, PyObject **unused,
        const char *fmt, ...);
static int sip_api_parse_pair(PyObject **parseErrp, PyObject *sipArg0,
        PyObject *sipArg1, const char *fmt, ...);
static void sip_api_no_function(PyObject *parseErr, const char *func,
        const char *doc);
static void sip_api_no_method(PyObject *parseErr, const char *scope,
        const char *method, const char *doc);
static void sip_api_abstract_method(const char *classname, const char *method);
static void sip_api_bad_class(const char *classname);
static void *sip_api_get_complex_cpp_ptr(sipSimpleWrapper *sw);
static PyObject *sip_api_is_py_method(sip_gilstate_t *gil, char *pymc,
        sipSimpleWrapper *sipSelf, const char *cname, const char *mname);
static PyObject *sip_api_is_py_method_12_8(sip_gilstate_t *gil, char *pymc,
        sipSimpleWrapper **sipSelfp, const char *cname, const char *mname);
static void sip_api_call_hook(const char *hookname);
static void sip_api_raise_unknown_exception(void);
static void sip_api_raise_type_exception(const sipTypeDef *td, void *ptr);
static int sip_api_add_type_instance(PyObject *dict, const char *name,
        void *cppPtr, const sipTypeDef *td);
static sipErrorState sip_api_bad_callable_arg(int arg_nr, PyObject *arg);
static void sip_api_bad_operator_arg(PyObject *self, PyObject *arg,
        sipPySlotType st);
static PyObject *sip_api_pyslot_extend(sipExportedModuleDef *mod,
        sipPySlotType st, const sipTypeDef *td, PyObject *arg0,
        PyObject *arg1);
static void sip_api_add_delayed_dtor(sipSimpleWrapper *w);
static int sip_api_export_symbol(const char *name, void *sym);
static void *sip_api_import_symbol(const char *name);
static const sipTypeDef *sip_api_find_type(const char *type);
static sipWrapperType *sip_api_find_class(const char *type);
static const sipMappedType *sip_api_find_mapped_type(const char *type);
static PyTypeObject *sip_api_find_named_enum(const char *type);
static char sip_api_bytes_as_char(PyObject *obj);
static const char *sip_api_bytes_as_string(PyObject *obj);
static char sip_api_string_as_ascii_char(PyObject *obj);
static const char *sip_api_string_as_ascii_string(PyObject **obj);
static char sip_api_string_as_latin1_char(PyObject *obj);
static const char *sip_api_string_as_latin1_string(PyObject **obj);
static char sip_api_string_as_utf8_char(PyObject *obj);
static const char *sip_api_string_as_utf8_string(PyObject **obj);
#if defined(HAVE_WCHAR_H)
static wchar_t sip_api_unicode_as_wchar(PyObject *obj);
static wchar_t *sip_api_unicode_as_wstring(PyObject *obj);
#else
static int sip_api_unicode_as_wchar(PyObject *obj);
static int *sip_api_unicode_as_wstring(PyObject *obj);
#endif
static void sip_api_transfer_break(PyObject *self);
static int sip_api_register_py_type(PyTypeObject *supertype);
static PyObject *sip_api_convert_from_enum(int eval, const sipTypeDef *td);
static const sipTypeDef *sip_api_type_from_py_type_object(PyTypeObject *py_type);
static const sipTypeDef *sip_api_type_scope(const sipTypeDef *td);
static const char *sip_api_resolve_typedef(const char *name);
static int sip_api_register_attribute_getter(const sipTypeDef *td,
        sipAttrGetterFunc getter);
static void sip_api_clear_any_slot_reference(sipSlot *slot);
static int sip_api_visit_slot(sipSlot *slot, visitproc visit, void *arg);
static void sip_api_keep_reference(PyObject *self, int key, PyObject *obj);
static PyObject *sip_api_get_reference(PyObject *self, int key);
static int sip_api_is_owned_by_python(sipSimpleWrapper *sw);
static int sip_api_is_derived_class(sipSimpleWrapper *sw);
static void sip_api_add_exception(sipErrorState es, PyObject **parseErrp);
static void sip_api_set_destroy_on_exit(int value);
static int sip_api_enable_autoconversion(const sipTypeDef *td, int enable);
static int sip_api_init_mixin(PyObject *self, PyObject *args, PyObject *kwds,
        const sipClassTypeDef *ctd);
static void *sip_api_get_mixin_address(sipSimpleWrapper *w,
        const sipTypeDef *td);
static int sip_api_register_proxy_resolver(const sipTypeDef *td,
        sipProxyResolverFunc resolver);
static PyInterpreterState *sip_api_get_interpreter(void);
static sipNewUserTypeFunc sip_api_set_new_user_type_handler(
        const sipTypeDef *td, sipNewUserTypeFunc handler);
static void sip_api_set_type_user_data(sipWrapperType *wt, void *data);
static void *sip_api_get_type_user_data(const sipWrapperType *wt);
static PyObject *sip_api_py_type_dict(const PyTypeObject *py_type);
static const char *sip_api_py_type_name(const PyTypeObject *py_type);
static int sip_api_get_method(PyObject *obj, sipMethodDef *method);
static PyObject *sip_api_from_method(const sipMethodDef *method);
static int sip_api_get_c_function(PyObject *obj, sipCFunctionDef *c_function);
static int sip_api_get_date(PyObject *obj, sipDateDef *date);
static PyObject *sip_api_from_date(const sipDateDef *date);
static int sip_api_get_datetime(PyObject *obj, sipDateDef *date,
        sipTimeDef *time);
static PyObject *sip_api_from_datetime(const sipDateDef *date,
        const sipTimeDef *time);
static int sip_api_get_time(PyObject *obj, sipTimeDef *time);
static PyObject *sip_api_from_time(const sipTimeDef *time);
static int sip_api_is_user_type(const sipWrapperType *wt);
static struct _frame *sip_api_get_frame(int);
static int sip_api_check_plugin_for_type(const sipTypeDef *td,
        const char *name);
static PyObject *sip_api_unicode_new(Py_ssize_t len, unsigned maxchar,
        int *kind, void **data);
static void sip_api_unicode_write(int kind, void *data, int index,
        unsigned value);
static void *sip_api_unicode_data(PyObject *obj, int *char_size,
        Py_ssize_t *len);
static int sip_api_get_buffer_info(PyObject *obj, sipBufferInfoDef *bi);
static void sip_api_release_buffer_info(sipBufferInfoDef *bi);
static PyObject *sip_api_get_user_object(const sipSimpleWrapper *sw);
static void sip_api_set_user_object(sipSimpleWrapper *sw, PyObject *user);
static int sip_api_enable_gc(int enable);
static void sip_api_print_object(PyObject *o);
static int sip_api_register_event_handler(sipEventType type,
        const sipTypeDef *td, void *handler);
static void sip_api_instance_destroyed_ex(sipSimpleWrapper **sipSelfp);
static void sip_api_visit_wrappers(sipWrapperVisitorFunc visitor,
        void *closure);
static int sip_api_register_exit_notifier(PyMethodDef *md);
static sipExceptionHandler sip_api_next_exception_handler(void **statep);


/*
 * The data structure that represents the SIP API.
 */
static const sipAPIDef sip_api = {
    /* This must be first. */
    sip_api_export_module,
    /*
     * The following are part of the public API.
     */
    (PyTypeObject *)&sipSimpleWrapper_Type,
    (PyTypeObject *)&sipWrapper_Type,
    &sipWrapperType_Type,
    &sipVoidPtr_Type,

    sip_api_bad_catcher_result,
    sip_api_bad_length_for_slice,
    sip_api_build_result,
    sip_api_call_method,
    sip_api_call_procedure_method,
    sip_api_connect_rx,
    sip_api_convert_from_sequence_index,
    sip_api_can_convert_to_type,
    sip_api_convert_to_type,
    sip_api_force_convert_to_type,
    sip_api_can_convert_to_enum,
    sip_api_release_type,
    sip_api_convert_from_type,
    sip_api_convert_from_new_type,
    sip_api_convert_from_enum,
    sip_api_get_state,
    sip_api_disconnect_rx,
    sip_api_free,
    sip_api_get_pyobject,
    sip_api_malloc,
    sip_api_parse_result,
    sip_api_trace,
    sip_api_transfer_back,
    sip_api_transfer_to,
    sip_api_transfer_break,
    sip_api_long_as_unsigned_long,
    sip_api_convert_from_void_ptr,
    sip_api_convert_from_const_void_ptr,
    sip_api_convert_from_void_ptr_and_size,
    sip_api_convert_from_const_void_ptr_and_size,
    sip_api_convert_to_void_ptr,
    sip_api_export_symbol,
    sip_api_import_symbol,
    sip_api_find_type,
    sip_api_register_py_type,
    sip_api_type_from_py_type_object,
    sip_api_type_scope,
    sip_api_resolve_typedef,
    sip_api_register_attribute_getter,
    sip_api_is_api_enabled,
    sip_api_bad_callable_arg,
    sip_api_get_address,
    sip_api_set_destroy_on_exit,
    sip_api_enable_autoconversion,
    sip_api_get_mixin_address,
    sip_api_convert_from_new_pytype,
    sip_api_convert_to_typed_array,
    sip_api_convert_to_array,
    sip_api_register_proxy_resolver,
    sip_api_get_interpreter,
    sip_api_set_new_user_type_handler,
    sip_api_set_type_user_data,
    sip_api_get_type_user_data,
    sip_api_py_type_dict,
    sip_api_py_type_name,
    sip_api_get_method,
    sip_api_from_method,
    sip_api_get_c_function,
    sip_api_get_date,
    sip_api_from_date,
    sip_api_get_datetime,
    sip_api_from_datetime,
    sip_api_get_time,
    sip_api_from_time,
    sip_api_is_user_type,
    sip_api_get_frame,
    sip_api_check_plugin_for_type,
    sip_api_unicode_new,
    sip_api_unicode_write,
    sip_api_unicode_data,
    sip_api_get_buffer_info,
    sip_api_release_buffer_info,
    sip_api_get_user_object,
    sip_api_set_user_object,
    /*
     * The following are not part of the public API.
     */
    sip_api_init_module,
    sip_api_parse_args,
    sip_api_parse_pair,
    /*
     * The following are part of the public API.
     */
    sip_api_instance_destroyed,
    /*
     * The following are not part of the public API.
     */
    sip_api_no_function,
    sip_api_no_method,
    sip_api_abstract_method,
    sip_api_bad_class,
    sip_api_get_cpp_ptr,
    sip_api_get_complex_cpp_ptr,
    sip_api_is_py_method,
    sip_api_call_hook,
    sip_api_end_thread,
    sip_api_raise_unknown_exception,
    sip_api_raise_type_exception,
    sip_api_add_type_instance,
    sip_api_bad_operator_arg,
    sip_api_pyslot_extend,
    sip_api_add_delayed_dtor,
    sip_api_bytes_as_char,
    sip_api_bytes_as_string,
    sip_api_string_as_ascii_char,
    sip_api_string_as_ascii_string,
    sip_api_string_as_latin1_char,
    sip_api_string_as_latin1_string,
    sip_api_string_as_utf8_char,
    sip_api_string_as_utf8_string,
    sip_api_unicode_as_wchar,
    sip_api_unicode_as_wstring,
    sip_api_deprecated,
    sip_api_keep_reference,
    sip_api_parse_kwd_args,
    sip_api_add_exception,
    sip_api_parse_result_ex,
    sip_api_call_error_handler,
    sip_api_init_mixin,
    sip_api_get_reference,
    /*
     * The following are part of the public API.
     */
    sip_api_is_owned_by_python,
    /*
     * The following are not part of the public API.
     */
    sip_api_is_derived_class,
    /*
     * The following may be used by Qt support code but by no other handwritten
     * code.
     */
    sip_api_free_sipslot,
    sip_api_same_slot,
    sip_api_convert_rx,
    sip_api_invoke_slot,
    sip_api_invoke_slot_ex,
    sip_api_save_slot,
    sip_api_clear_any_slot_reference,
    sip_api_visit_slot,
    /*
     * The following are deprecated parts of the public API.
     */
    sip_api_find_named_enum,
    sip_api_find_mapped_type,
    sip_api_find_class,
    sip_api_map_int_to_class,
    sip_api_map_string_to_class,
    /*
     * The following are part of the public API.
     */
    sip_api_enable_gc,
    sip_api_print_object,
    sip_api_register_event_handler,
    sip_api_convert_to_enum,
    sip_api_convert_to_bool,
    sip_api_enable_overflow_checking,
    sip_api_long_as_char,
    sip_api_long_as_signed_char,
    sip_api_long_as_unsigned_char,
    sip_api_long_as_short,
    sip_api_long_as_unsigned_short,
    sip_api_long_as_int,
    sip_api_long_as_unsigned_int,
    sip_api_long_as_long,
#if defined(HAVE_LONG_LONG)
    sip_api_long_as_long_long,
    sip_api_long_as_unsigned_long_long,
#else
    0,
    0,
#endif
    /*
     * The following are not part of the public API.
     */
    sip_api_instance_destroyed_ex,
    /*
     * The following are part of the public API.
     */
    sip_api_convert_from_slice_object,
    sip_api_long_as_size_t,
    sip_api_visit_wrappers,
    sip_api_register_exit_notifier,
    /*
     * The following are not part of the public API.
     */
    sip_api_is_py_method_12_8,
    sip_api_next_exception_handler,
};


#define AUTO_DOCSTRING          '\1'    /* Marks an auto class docstring. */


/*
 * These are the format flags supported by argument parsers.
 */
#define FMT_AP_DEREF            0x01    /* The pointer will be dereferenced. */
#define FMT_AP_TRANSFER         0x02    /* Implement /Transfer/. */
#define FMT_AP_TRANSFER_BACK    0x04    /* Implement /TransferBack/. */
#define FMT_AP_NO_CONVERTORS    0x08    /* Suppress any convertors. */
#define FMT_AP_TRANSFER_THIS    0x10    /* Support for /TransferThis/. */


/*
 * These are the format flags supported by result parsers.  Deprecated values
 * have a _DEPR suffix.
 */
#define FMT_RP_DEREF            0x01    /* The pointer will be dereferenced. */
#define FMT_RP_FACTORY          0x02    /* /Factory/ or /TransferBack/. */
#define FMT_RP_MAKE_COPY        0x04    /* Return a copy of the value. */
#define FMT_RP_NO_STATE_DEPR    0x04    /* Don't return the C/C++ state. */


/*
 * The different reasons for failing to parse an overload.  These include
 * internal (i.e. non-user) errors.
 */
typedef enum {
    Ok, Unbound, TooFew, TooMany, UnknownKeyword, Duplicate, WrongType, Raised,
    KeywordNotString, Exception, Overflow
} sipParseFailureReason;


/*
 * The description of a failure to parse an overload because of a user error.
 */
typedef struct _sipParseFailure {
    sipParseFailureReason reason;   /* The reason for the failure. */
    const char *detail_str;         /* The detail if a string. */
    PyObject *detail_obj;           /* The detail if a Python object. */
    int arg_nr;                     /* The wrong positional argument. */
    const char *arg_name;           /* The wrong keyword argument. */
    int overflow_arg_nr;            /* The overflowed positional argument. */
    const char *overflow_arg_name;  /* The overflowed keyword argument. */
} sipParseFailure;


/*
 * An entry in a linked list of name/symbol pairs.
 */
typedef struct _sipSymbol {
    const char *name;           /* The name. */
    void *symbol;               /* The symbol. */
    struct _sipSymbol *next;    /* The next in the list. */
} sipSymbol;


/*
 * An entry in a linked list of Python objects.
 */
typedef struct _sipPyObject {
    PyObject *object;           /* The Python object. */
    struct _sipPyObject *next;  /* The next in the list. */
} sipPyObject;


/*
 * An entry in the linked list of attribute getters.
 */
typedef struct _sipAttrGetter {
    PyTypeObject *type;             /* The Python type being handled. */
    sipAttrGetterFunc getter;       /* The getter. */
    struct _sipAttrGetter *next;    /* The next in the list. */
} sipAttrGetter;


/*
 * An entry in the linked list of proxy resolvers.
 */
typedef struct _sipProxyResolver {
    const sipTypeDef *td;           /* The type the resolver handles. */
    sipProxyResolverFunc resolver;  /* The resolver. */
    struct _sipProxyResolver *next; /* The next in the list. */
} sipProxyResolver;


/*
 * An entry in the linked list of event handlers.
 */
typedef struct _sipEventHandler {
    const sipClassTypeDef *ctd;     /* The type the handler handles. */
    void *handler;                  /* The handler. */
    struct _sipEventHandler *next;  /* The next in the list. */
} sipEventHandler;


/*****************************************************************************
 * The structures to support a Python type to hold a named enum.
 *****************************************************************************/

static PyObject *sipEnumType_alloc(PyTypeObject *self, Py_ssize_t nitems);
static PyObject *sipEnumType_getattro(PyObject *self, PyObject *name);


/*
 * The type data structure.  We inherit everything from the standard Python
 * metatype and the size of the type object created is increased to accomodate
 * the extra information we associate with a named enum type.
 */
static PyTypeObject sipEnumType_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "sip.enumtype",         /* tp_name */
    sizeof (sipEnumTypeObject), /* tp_basicsize */
    0,                      /* tp_itemsize */
    0,                      /* tp_dealloc */
    0,                      /* tp_print */
    0,                      /* tp_getattr */
    0,                      /* tp_setattr */
    0,                      /* tp_as_async (Python v3.5), tp_compare (Python v2) */
    0,                      /* tp_repr */
    0,                      /* tp_as_number */
    0,                      /* tp_as_sequence */
    0,                      /* tp_as_mapping */
    0,                      /* tp_hash */
    0,                      /* tp_call */
    0,                      /* tp_str */
    sipEnumType_getattro,   /* tp_getattro */
    0,                      /* tp_setattro */
    0,                      /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
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
    sipEnumType_alloc,      /* tp_alloc */
    0,                      /* tp_new */
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


sipQtAPI *sipQtSupport = NULL;
sipTypeDef *sipQObjectType;

static int got_kw_handler = FALSE;
static int (*kw_handler)(PyObject *, void *, PyObject *);


/*
 * Various strings as Python objects created as and when needed.
 */
static PyObject *licenseName = NULL;
static PyObject *licenseeName = NULL;
static PyObject *typeName = NULL;
static PyObject *timestampName = NULL;
static PyObject *signatureName = NULL;

static sipObjectMap cppPyMap;           /* The C/C++ to Python map. */
static sipExportedModuleDef *moduleList = NULL; /* List of registered modules. */
static unsigned traceMask = 0;          /* The current trace mask. */

static sipTypeDef *currentType = NULL;  /* The type being created. */
static PyObject **unused_backdoor = NULL;   /* For passing dict of unused arguments. */

static PyObject *init_name = NULL;      /* '__init__'. */
static PyObject *empty_tuple;           /* The empty tuple. */
static PyObject *type_unpickler;        /* The type unpickler function. */
static PyObject *enum_unpickler;        /* The enum unpickler function. */
static sipSymbol *sipSymbolList = NULL; /* The list of published symbols. */
static sipAttrGetter *sipAttrGetters = NULL;  /* The list of attribute getters. */
static sipProxyResolver *proxyResolvers = NULL; /* The list of proxy resolvers. */
static sipPyObject *sipRegisteredPyTypes = NULL;    /* Registered Python types. */
static sipPyObject *sipDisabledAutoconversions = NULL;  /* Python types whose auto-conversion is disabled. */
static PyInterpreterState *sipInterpreter = NULL;   /* The interpreter. */
static int destroy_on_exit = TRUE;      /* Destroy owned objects on exit. */
static sipEventHandler *event_handlers[sipEventNrEvents];   /* The event handler lists. */

static void addClassSlots(sipWrapperType *wt, const sipClassTypeDef *ctd);
static void addTypeSlots(PyHeapTypeObject *heap_to, sipPySlotDef *slots);
static void *findSlot(PyObject *self, sipPySlotType st);
static void *findSlotInClass(const sipClassTypeDef *psd, sipPySlotType st);
static void *findSlotInSlotList(sipPySlotDef *psd, sipPySlotType st);
static int objobjargprocSlot(PyObject *self, PyObject *arg1, PyObject *arg2,
        sipPySlotType st);
static int ssizeobjargprocSlot(PyObject *self, Py_ssize_t arg1,
        PyObject *arg2, sipPySlotType st);
static PyObject *buildObject(PyObject *tup, const char *fmt, va_list va);
static int parseKwdArgs(PyObject **parseErrp, PyObject *sipArgs,
        PyObject *sipKwdArgs, const char **kwdlist, PyObject **unused,
        const char *fmt, va_list va_orig);
static int parsePass1(PyObject **parseErrp, sipSimpleWrapper **selfp,
        int *selfargp, PyObject *sipArgs, PyObject *sipKwdArgs,
        const char **kwdlist, PyObject **unused, const char *fmt, va_list va);
static int parsePass2(sipSimpleWrapper *self, int selfarg, PyObject *sipArgs,
        PyObject *sipKwdArgs, const char **kwdlist, const char *fmt,
        va_list va);
static int parseResult(PyObject *method, PyObject *res,
        sipSimpleWrapper *py_self, const char *fmt, va_list va);
static PyObject *signature_FromDocstring(const char *doc, Py_ssize_t line);
static PyObject *detail_FromFailure(PyObject *failure_obj);
static int isQObject(PyObject *obj);
static int canConvertFromSequence(PyObject *seq, const sipTypeDef *td);
static int convertFromSequence(PyObject *seq, const sipTypeDef *td,
        void **array, Py_ssize_t *nr_elem);
static PyObject *convertToSequence(void *array, Py_ssize_t nr_elem,
        const sipTypeDef *td);
static int getSelfFromArgs(sipTypeDef *td, PyObject *args, int argnr,
        sipSimpleWrapper **selfp);
static int compareTypedefName(const void *key, const void *el);
static int checkPointer(void *ptr, sipSimpleWrapper *sw);
static void *cast_cpp_ptr(void *ptr, PyTypeObject *src_type,
        const sipTypeDef *dst_type);
static void finalise(void);
static PyObject *getDefaultBase(void);
static PyObject *getDefaultSimpleBase(void);
static PyObject *getScopeDict(sipTypeDef *td, PyObject *mod_dict,
        sipExportedModuleDef *client);
static PyObject *createContainerType(sipContainerDef *cod, sipTypeDef *td,
        PyObject *bases, PyObject *metatype, PyObject *mod_dict,
        PyObject *type_dict, sipExportedModuleDef *client);
static int createClassType(sipExportedModuleDef *client, sipClassTypeDef *ctd,
        PyObject *mod_dict);
static int createMappedType(sipExportedModuleDef *client,
        sipMappedTypeDef *mtd, PyObject *mod_dict);
static sipExportedModuleDef *getModule(PyObject *mname_obj);
static PyObject *pickle_type(PyObject *obj, PyObject *args);
static PyObject *unpickle_type(PyObject *obj, PyObject *args);
static PyObject *pickle_enum(PyObject *obj, PyObject *args);
static PyObject *unpickle_enum(PyObject *obj, PyObject *args);
static int setReduce(PyTypeObject *type, PyMethodDef *pickler);
static int createEnum(sipExportedModuleDef *client, sipEnumTypeDef *etd,
        int enum_nr, PyObject *mod_dict);
static PyObject *createUnscopedEnum(sipExportedModuleDef *client,
        sipEnumTypeDef *etd, PyObject *name);
static PyObject *createScopedEnum(sipExportedModuleDef *client,
        sipEnumTypeDef *etd, int enum_nr, PyObject *name);
static PyObject *createTypeDict(sipExportedModuleDef *em);
static sipTypeDef *getGeneratedType(const sipEncodedTypeDef *enc,
        sipExportedModuleDef *em);
static const sipTypeDef *convertSubClass(const sipTypeDef *td, void **cppPtr);
static int convertPass(const sipTypeDef **tdp, void **cppPtr);
static void *getPtrTypeDef(sipSimpleWrapper *self,
        const sipClassTypeDef **ctd);
static int addInstances(PyObject *dict, sipInstancesDef *id);
static int addVoidPtrInstances(PyObject *dict, sipVoidPtrInstanceDef *vi);
static int addCharInstances(PyObject *dict, sipCharInstanceDef *ci);
static int addStringInstances(PyObject *dict, sipStringInstanceDef *si);
static int addIntInstances(PyObject *dict, sipIntInstanceDef *ii);
static int addLongInstances(PyObject *dict, sipLongInstanceDef *li);
static int addUnsignedLongInstances(PyObject *dict,
        sipUnsignedLongInstanceDef *uli);
static int addLongLongInstances(PyObject *dict, sipLongLongInstanceDef *lli);
static int addUnsignedLongLongInstances(PyObject *dict,
        sipUnsignedLongLongInstanceDef *ulli);
static int addDoubleInstances(PyObject *dict, sipDoubleInstanceDef *di);
static int addTypeInstances(PyObject *dict, sipTypeInstanceDef *ti);
static int addSingleTypeInstance(PyObject *dict, const char *name,
        void *cppPtr, const sipTypeDef *td, int initflags);
static int addLicense(PyObject *dict, sipLicenseDef *lc);
static PyObject *assign(PyObject *self, PyObject *args);
static PyObject *cast(PyObject *self, PyObject *args);
static PyObject *callDtor(PyObject *self, PyObject *args);
static PyObject *dumpWrapper(PyObject *self, PyObject *arg);
static PyObject *enableAutoconversion(PyObject *self, PyObject *args);
static PyObject *isDeleted(PyObject *self, PyObject *args);
static PyObject *isPyCreated(PyObject *self, PyObject *args);
static PyObject *isPyOwned(PyObject *self, PyObject *args);
static PyObject *setDeleted(PyObject *self, PyObject *args);
static PyObject *setTraceMask(PyObject *self, PyObject *args);
static PyObject *wrapInstance(PyObject *self, PyObject *args);
static PyObject *unwrapInstance(PyObject *self, PyObject *args);
static PyObject *transferBack(PyObject *self, PyObject *args);
static PyObject *transferTo(PyObject *self, PyObject *args);
static PyObject *setDestroyOnExit(PyObject *self, PyObject *args);
static void clear_wrapper(sipSimpleWrapper *sw);
static void print_object(const char *label, PyObject *obj);
static void addToParent(sipWrapper *self, sipWrapper *owner);
static void removeFromParent(sipWrapper *self);
static void detachChildren(sipWrapper *self);
static void release(void *addr, const sipTypeDef *td, int state);
static void callPyDtor(sipSimpleWrapper *self);
static int parseBytes_AsCharArray(PyObject *obj, const char **ap,
        Py_ssize_t *aszp);
static int parseBytes_AsChar(PyObject *obj, char *ap);
static int parseBytes_AsString(PyObject *obj, const char **ap);
static int parseString_AsASCIIChar(PyObject *obj, char *ap);
static PyObject *parseString_AsASCIIString(PyObject *obj, const char **ap);
static int parseString_AsLatin1Char(PyObject *obj, char *ap);
static PyObject *parseString_AsLatin1String(PyObject *obj, const char **ap);
static int parseString_AsUTF8Char(PyObject *obj, char *ap);
static PyObject *parseString_AsUTF8String(PyObject *obj, const char **ap);
static int parseString_AsEncodedChar(PyObject *bytes, PyObject *obj, char *ap);
static PyObject *parseString_AsEncodedString(PyObject *bytes, PyObject *obj,
        const char **ap);
#if defined(HAVE_WCHAR_H)
static int parseWCharArray(PyObject *obj, wchar_t **ap, Py_ssize_t *aszp);
static int convertToWCharArray(PyObject *obj, wchar_t **ap, Py_ssize_t *aszp);
static int parseWChar(PyObject *obj, wchar_t *ap);
static int convertToWChar(PyObject *obj, wchar_t *ap);
static int parseWCharString(PyObject *obj, wchar_t **ap);
static int convertToWCharString(PyObject *obj, wchar_t **ap);
#else
static void raiseNoWChar();
#endif
static void *getComplexCppPtr(sipSimpleWrapper *w, const sipTypeDef *td);
static PyObject *findPyType(const char *name);
static int addPyObjectToList(sipPyObject **head, PyObject *object);
static PyObject *getDictFromObject(PyObject *obj);
static void forgetObject(sipSimpleWrapper *sw);
static int add_lazy_container_attrs(sipTypeDef *td, sipContainerDef *cod,
        PyObject *dict);
static int add_lazy_attrs(sipTypeDef *td);
static int add_all_lazy_attrs(sipTypeDef *td);
static int objectify(const char *s, PyObject **objp);
static void add_failure(PyObject **parseErrp, sipParseFailure *failure);
static PyObject *bad_type_str(int arg_nr, PyObject *arg);
static void *explicit_access_func(sipSimpleWrapper *sw, AccessFuncOp op);
static void *indirect_access_func(sipSimpleWrapper *sw, AccessFuncOp op);
static void clear_access_func(sipSimpleWrapper *sw);
static int check_encoded_string(PyObject *obj);
static int isNonlazyMethod(PyMethodDef *pmd);
static int addMethod(PyObject *dict, PyMethodDef *pmd);
static PyObject *create_property(sipVariableDef *vd);
static PyObject *create_function(PyMethodDef *ml);
static PyObject *sip_exit(PyObject *self, PyObject *args);
static sipConvertFromFunc get_from_convertor(const sipTypeDef *td);
static sipPyObject **autoconversion_disabled(const sipTypeDef *td);
static void fix_slots(PyTypeObject *py_type, sipPySlotDef *psd);
static sipFinalFunc find_finalisation(sipClassTypeDef *ctd);
static sipNewUserTypeFunc find_new_user_type_handler(sipWrapperType *wt);
static PyObject *next_in_mro(PyObject *self, PyObject *after);
static int super_init(PyObject *self, PyObject *args, PyObject *kwds,
        PyObject *type);
static sipSimpleWrapper *deref_mixin(sipSimpleWrapper *w);
static PyObject *wrap_simple_instance(void *cpp, const sipTypeDef *td,
        sipWrapper *owner, int flags);
static void *resolve_proxy(const sipTypeDef *td, void *proxy);
static PyObject *call_method(PyObject *method, const char *fmt, va_list va);
static int importTypes(sipExportedModuleDef *client, sipImportedModuleDef *im,
        sipExportedModuleDef *em);
static int importErrorHandlers(sipExportedModuleDef *client,
        sipImportedModuleDef *im, sipExportedModuleDef *em);
static int importExceptions(sipExportedModuleDef *client,
        sipImportedModuleDef *im, sipExportedModuleDef *em);
static int is_subtype(const sipClassTypeDef *ctd,
        const sipClassTypeDef *base_ctd);
static PyObject *import_module_attr(const char *module, const char *attr);
static const sipContainerDef *get_container(const sipTypeDef *td);
static PyObject *get_qualname(const sipTypeDef *td, PyObject *name);
static int convert_to_enum(PyObject *obj, const sipTypeDef *td, int allow_int);
static void handle_failed_int_conversion(sipParseFailure *pf, PyObject *arg);
static void enum_expected(PyObject *obj, const sipTypeDef *td);
static int long_as_nonoverflow_int(PyObject *val_obj);
static int dict_set_and_discard(PyObject *dict, const char *name,
        PyObject *obj);


/*
 * Initialise the module as a library.
 */
const sipAPIDef *sip_init_library(PyObject *mod_dict)
{
    static PyMethodDef methods[] = {
        /* This must be first. */
        {"_unpickle_enum", unpickle_enum, METH_VARARGS, NULL},
        /* This must be second. */
        {"_unpickle_type", unpickle_type, METH_VARARGS, NULL},
        {"assign", assign, METH_VARARGS, NULL},
        {"cast", cast, METH_VARARGS, NULL},
        {"delete", callDtor, METH_VARARGS, NULL},
        {"dump", dumpWrapper, METH_O, NULL},
        {"enableautoconversion", enableAutoconversion, METH_VARARGS, NULL},
        {"enableoverflowchecking", sipEnableOverflowChecking, METH_VARARGS, NULL},
        {"getapi", sipGetAPI, METH_VARARGS, NULL},
        {"isdeleted", isDeleted, METH_VARARGS, NULL},
        {"ispycreated", isPyCreated, METH_VARARGS, NULL},
        {"ispyowned", isPyOwned, METH_VARARGS, NULL},
        {"setapi", sipSetAPI, METH_VARARGS, NULL},
        {"setdeleted", setDeleted, METH_VARARGS, NULL},
        {"setdestroyonexit", setDestroyOnExit, METH_VARARGS, NULL},
        {"settracemask", setTraceMask, METH_VARARGS, NULL},
        {"transferback", transferBack, METH_VARARGS, NULL},
        {"transferto", transferTo, METH_VARARGS, NULL},
        {"wrapinstance", wrapInstance, METH_VARARGS, NULL},
        {"unwrapinstance", unwrapInstance, METH_VARARGS, NULL},
        {NULL, NULL, 0, NULL}
    };

    static PyMethodDef sip_exit_md = {
        "_sip_exit", sip_exit, METH_NOARGS, NULL
    };

    PyObject *obj;
    PyMethodDef *md;

    /* Add the SIP version number. */
    obj = PyLong_FromLong(SIP_VERSION);

    if (dict_set_and_discard(mod_dict, "SIP_VERSION", obj) < 0)
        return NULL;

    obj = PyUnicode_FromString(SIP_VERSION_STR);

    if (dict_set_and_discard(mod_dict, "SIP_VERSION_STR", obj) < 0)
        return NULL;

    /* Add the methods. */
    for (md = methods; md->ml_name != NULL; ++md)
    {
        PyObject *meth = PyCFunction_New(md, NULL);

        if (dict_set_and_discard(mod_dict, md->ml_name, meth) < 0)
            return NULL;

        if (md == &methods[0])
        {
            Py_INCREF(meth);
            enum_unpickler = meth;
        }
        else if (md == &methods[1])
        {
            Py_INCREF(meth);
            type_unpickler = meth;
        }
    }

    /* Initialise the types. */
    sipWrapperType_Type.tp_base = &PyType_Type;

    if (PyType_Ready(&sipWrapperType_Type) < 0)
        return NULL;

    if (PyType_Ready((PyTypeObject *)&sipSimpleWrapper_Type) < 0)
        return NULL;

    if (sip_api_register_py_type((PyTypeObject *)&sipSimpleWrapper_Type) < 0)
        return NULL;

#if defined(STACKLESS)
    sipWrapper_Type.super.tp_base = (PyTypeObject *)&sipSimpleWrapper_Type;
#else
    sipWrapper_Type.super.ht_type.tp_base = (PyTypeObject *)&sipSimpleWrapper_Type;
#endif

    if (PyType_Ready((PyTypeObject *)&sipWrapper_Type) < 0)
        return NULL;

    if (PyType_Ready(&sipMethodDescr_Type) < 0)
        return NULL;

    if (PyType_Ready(&sipVariableDescr_Type) < 0)
        return NULL;

    sipEnumType_Type.tp_base = &PyType_Type;

    if (PyType_Ready(&sipEnumType_Type) < 0)
        return NULL;

    if (PyType_Ready(&sipVoidPtr_Type) < 0)
        return NULL;

    if (PyType_Ready(&sipArray_Type) < 0)
        return NULL;

    /* Add the public types. */
    if (PyDict_SetItemString(mod_dict, "wrappertype", (PyObject *)&sipWrapperType_Type) < 0)
        return NULL;

    if (PyDict_SetItemString(mod_dict, "simplewrapper", (PyObject *)&sipSimpleWrapper_Type) < 0)
        return NULL;

    if (PyDict_SetItemString(mod_dict, "wrapper", (PyObject *)&sipWrapper_Type) < 0)
        return NULL;

    if (PyDict_SetItemString(mod_dict, "voidptr", (PyObject *)&sipVoidPtr_Type) < 0)
        return NULL;

    if (PyDict_SetItemString(mod_dict, "array", (PyObject *)&sipArray_Type) < 0)
        return NULL;

    /* These will always be needed. */
    if (objectify("__init__", &init_name) < 0)
        return NULL;

    if ((empty_tuple = PyTuple_New(0)) == NULL)
        return NULL;

    /* Initialise the object map. */
    sipOMInit(&cppPyMap);

    /* Make sure we are notified at the end of the exit process. */
    if (Py_AtExit(finalise) < 0)
        return NULL;

    /* Make sure we are notified at the start of the exit process. */
    if (sip_api_register_exit_notifier(&sip_exit_md) < 0)
        return NULL;

    /*
     * Get the current interpreter.  This will be shared between all threads.
     */
    sipInterpreter = PyThreadState_Get()->interp;

    return &sip_api;
}


/*
 * Set a dictionary item and discard the reference to the item even if there
 * was an error.
 */
static int dict_set_and_discard(PyObject *dict, const char *name, PyObject *obj)
{
    int rc;

    if (obj == NULL)
        return -1;

    rc = PyDict_SetItemString(dict, name, obj);

    Py_DECREF(obj);

    return rc;
}


#if _SIP_MODULE_SHARED
/*
 * The Python module initialisation function.
 */
#if defined(SIP_STATIC_MODULE)
PyObject *_SIP_MODULE_ENTRY(void)
#else
PyMODINIT_FUNC _SIP_MODULE_ENTRY(void)
#endif
{
    static PyModuleDef module_def = {
        PyModuleDef_HEAD_INIT,
        _SIP_MODULE_FQ_NAME,    /* m_name */
        NULL,                   /* m_doc */
        -1,                     /* m_size */
        NULL,                   /* m_methods */
        NULL,                   /* m_reload */
        NULL,                   /* m_traverse */
        NULL,                   /* m_clear */
        NULL,                   /* m_free */
    };

    const sipAPIDef *api;
    PyObject *mod, *mod_dict, *api_obj;

    /* Create the module. */
    if ((mod = PyModule_Create(&module_def)) == NULL)
        return NULL;

    mod_dict = PyModule_GetDict(mod);

    /* Initialise the module dictionary and static variables. */
    if ((api = sip_init_library(mod_dict)) == NULL)
        return NULL;

    /* Publish the SIP API. */
    api_obj = PyCapsule_New((void *)api, _SIP_MODULE_FQ_NAME "._C_API", NULL);

    if (dict_set_and_discard(mod_dict, "_C_API", api_obj) < 0)
    {
        Py_DECREF(mod);
        return NULL;
    }

#if _SIP_MODULE_LEGACY
    {
        /*
         * Also install the package-specific module at the top level for
         * backwards compatibility.
         */
        PyObject *modules = PySys_GetObject("modules");

        if (modules != NULL)
            PyDict_SetItemString(modules, "sip", mod);
    }
#endif

    return mod;
}
#endif


/*
 * Return the current interpreter, if there is one.
 */
static PyInterpreterState *sip_api_get_interpreter(void)
{
    return sipInterpreter;
}


/*
 * Display a printf() style message to stderr according to the current trace
 * mask.
 */
static void sip_api_trace(unsigned mask, const char *fmt, ...)
{
    va_list ap;

    va_start(ap,fmt);

    if (mask & traceMask)
        vfprintf(stderr, fmt, ap);

    va_end(ap);
}


/*
 * Set the trace mask.
 */
static PyObject *setTraceMask(PyObject *self, PyObject *args)
{
    unsigned new_mask;

    (void)self;

    if (PyArg_ParseTuple(args, "I:settracemask", &new_mask))
    {
        traceMask = new_mask;

        Py_INCREF(Py_None);
        return Py_None;
    }

    return NULL;
}


/*
 * Dump various bits of potentially useful information to stdout.  Note that we
 * use the same calling convention as sys.getrefcount() so that it has the
 * same caveat regarding the reference count.
 */
static PyObject *dumpWrapper(PyObject *self, PyObject *arg)
{
    sipSimpleWrapper *sw;

    (void)self;

    if (!PyObject_TypeCheck(arg, (PyTypeObject *)&sipSimpleWrapper_Type))
    {
        PyErr_Format(PyExc_TypeError,
                "dump() argument 1 must be sip.simplewrapper, not %s",
                Py_TYPE(arg)->tp_name);
        return NULL;
    }

    sw = (sipSimpleWrapper *)arg;

    print_object(NULL, (PyObject *)sw);

    printf("    Reference count: %" PY_FORMAT_SIZE_T "d\n", Py_REFCNT(sw));
    printf("    Address of wrapped object: %p\n", sip_api_get_address(sw));
    printf("    Created by: %s\n", (sipIsDerived(sw) ? "Python" : "C/C++"));
    printf("    To be destroyed by: %s\n", (sipIsPyOwned(sw) ? "Python" : "C/C++"));

    if (PyObject_TypeCheck((PyObject *)sw, (PyTypeObject *)&sipWrapper_Type))
    {
        sipWrapper *w = (sipWrapper *)sw;

        print_object("Parent wrapper", (PyObject *)w->parent);
        print_object("Next sibling wrapper", (PyObject *)w->sibling_next);
        print_object("Previous sibling wrapper", (PyObject *)w->sibling_prev);
        print_object("First child wrapper", (PyObject *)w->first_child);
    }

    Py_INCREF(Py_None);
    return Py_None;
}


/*
 * Write a reference to a wrapper to stdout.
 */
static void print_object(const char *label, PyObject *obj)
{
    if (label != NULL)
        printf("    %s: ", label);

    if (obj != NULL)
        PyObject_Print(obj, stdout, 0);
    else
        printf("NULL");

    printf("\n");
}


/*
 * Transfer the ownership of an instance to C/C++.
 */
static PyObject *transferTo(PyObject *self, PyObject *args)
{
    PyObject *w, *owner;

    (void)self;

    if (PyArg_ParseTuple(args, "O!O:transferto", &sipWrapper_Type, &w, &owner))
    {
        if (owner == Py_None)
        {
            /*
             * Note that the Python API is different to the C API when the
             * owner is None.
             */
            owner = NULL;
        }
        else if (!PyObject_TypeCheck(owner, (PyTypeObject *)&sipWrapper_Type))
        {
            PyErr_Format(PyExc_TypeError, "transferto() argument 2 must be sip.wrapper, not %s", Py_TYPE(owner)->tp_name);
            return NULL;
        }

        sip_api_transfer_to(w, owner);

        Py_INCREF(Py_None);
        return Py_None;
    }

    return NULL;
}


/*
 * Transfer the ownership of an instance to Python.
 */
static PyObject *transferBack(PyObject *self, PyObject *args)
{
    PyObject *w;

    (void)self;

    if (PyArg_ParseTuple(args, "O!:transferback", &sipWrapper_Type, &w))
    {
        sip_api_transfer_back(w);

        Py_INCREF(Py_None);
        return Py_None;
    }

    return NULL;
}


/*
 * Invoke the assignment operator for a C++ instance.
 */
static PyObject *assign(PyObject *self, PyObject *args)
{
    sipSimpleWrapper *dst, *src;
    PyTypeObject *dst_type, *src_type;
    const sipTypeDef *td, *super_td;
    sipAssignFunc assign_helper;
    void *dst_addr, *src_addr;

    (void)self;

    if (!PyArg_ParseTuple(args, "O!O!:assign", &sipSimpleWrapper_Type, &dst, &sipSimpleWrapper_Type, &src))
        return NULL;

    /* Get the assignment helper. */
    dst_type = Py_TYPE(dst);
    td = ((sipWrapperType *)dst_type)->wt_td;

    if (sipTypeIsMapped(td))
        assign_helper = ((const sipMappedTypeDef *)td)->mtd_assign;
    else
        assign_helper = ((const sipClassTypeDef *)td)->ctd_assign;

    if (assign_helper == NULL)
    {
        PyErr_SetString(PyExc_TypeError,
                "argument 1 of assign() does not support assignment");
        return NULL;
    }

    /* Check the types are compatible. */
    src_type = Py_TYPE(src);

    if (src_type == dst_type)
    {
        super_td = NULL;
    }
    else if (PyType_IsSubtype(src_type, dst_type))
    {
        super_td = td;
    }
    else
    {
        PyErr_SetString(PyExc_TypeError,
                "type of argument 1 of assign() must be a super-type of type of argument 2");
        return NULL;
    }

    /* Get the addresses. */
    if ((dst_addr = sip_api_get_cpp_ptr(dst, NULL)) == NULL)
        return NULL;

    if ((src_addr = sip_api_get_cpp_ptr(src, super_td)) == NULL)
        return NULL;

    /* Do the assignment. */
    assign_helper(dst_addr, 0, src_addr);

    Py_INCREF(Py_None);
    return Py_None;
}


/*
 * Cast an instance to one of it's sub or super-classes by returning a new
 * Python object with the superclass type wrapping the same C++ instance.
 */
static PyObject *cast(PyObject *self, PyObject *args)
{
    sipSimpleWrapper *sw;
    sipWrapperType *wt;
    const sipTypeDef *td;
    void *addr;
    PyTypeObject *ft, *tt;

    (void)self;

    if (!PyArg_ParseTuple(args, "O!O!:cast", &sipSimpleWrapper_Type, &sw, &sipWrapperType_Type, &wt))
        return NULL;

    ft = Py_TYPE(sw);
    tt = (PyTypeObject *)wt;

    if (ft == tt || PyType_IsSubtype(tt, ft))
        td = NULL;
    else if (PyType_IsSubtype(ft, tt))
        td = wt->wt_td;
    else
    {
        PyErr_SetString(PyExc_TypeError, "argument 1 of cast() must be an instance of a sub or super-type of argument 2");
        return NULL;
    }

    if ((addr = sip_api_get_cpp_ptr(sw, td)) == NULL)
        return NULL;

    /*
     * We don't put this new object into the map so that the original object is
     * always found.  It would also totally confuse the map logic.
     */
    return wrap_simple_instance(addr, wt->wt_td, NULL,
            (sw->sw_flags | SIP_NOT_IN_MAP) & ~SIP_PY_OWNED);
}


/*
 * Call an instance's dtor.
 */
static PyObject *callDtor(PyObject *self, PyObject *args)
{
    sipSimpleWrapper *sw;
    void *addr;
    const sipClassTypeDef *ctd;

    (void)self;

    if (!PyArg_ParseTuple(args, "O!:delete", &sipSimpleWrapper_Type, &sw))
        return NULL;

    addr = getPtrTypeDef(sw, &ctd);

    if (checkPointer(addr, sw) < 0)
        return NULL;

    clear_wrapper(sw);

    release(addr, (const sipTypeDef *)ctd, sw->sw_flags);

    Py_INCREF(Py_None);
    return Py_None;
}


/*
 * Check if an instance still exists without raising an exception.
 */
static PyObject *isDeleted(PyObject *self, PyObject *args)
{
    sipSimpleWrapper *sw;
    PyObject *res;

    (void)self;

    if (!PyArg_ParseTuple(args, "O!:isdeleted", &sipSimpleWrapper_Type, &sw))
        return NULL;

    res = (sip_api_get_address(sw) == NULL ? Py_True : Py_False);

    Py_INCREF(res);
    return res;
}


/*
 * Check if an instance was created by Python.
 */
static PyObject *isPyCreated(PyObject *self, PyObject *args)
{
    sipSimpleWrapper *sw;
    PyObject *res;

    (void)self;

    if (!PyArg_ParseTuple(args, "O!:ispycreated", &sipSimpleWrapper_Type, &sw))
        return NULL;

    /* sipIsDerived() is a misnomer. */
    res = (sipIsDerived(sw) ? Py_True : Py_False);

    Py_INCREF(res);
    return res;
}


/*
 * Check if an instance is owned by Python or C/C++.
 */
static PyObject *isPyOwned(PyObject *self, PyObject *args)
{
    sipSimpleWrapper *sw;
    PyObject *res;

    (void)self;

    if (!PyArg_ParseTuple(args, "O!:ispyowned", &sipSimpleWrapper_Type, &sw))
        return NULL;

    res = (sipIsPyOwned(sw) ? Py_True : Py_False);

    Py_INCREF(res);
    return res;
}


/*
 * Mark an instance as having been deleted.
 */
static PyObject *setDeleted(PyObject *self, PyObject *args)
{
    sipSimpleWrapper *sw;

    (void)self;

    if (!PyArg_ParseTuple(args, "O!:setdeleted", &sipSimpleWrapper_Type, &sw))
        return NULL;

    clear_wrapper(sw);

    Py_INCREF(Py_None);
    return Py_None;
}


/*
 * Unwrap an instance.
 */
static PyObject *unwrapInstance(PyObject *self, PyObject *args)
{
    sipSimpleWrapper *sw;

    (void)self;

    if (PyArg_ParseTuple(args, "O!:unwrapinstance", &sipSimpleWrapper_Type, &sw))
    {
        void *addr;

        /*
         * We just get the pointer but don't try and cast it (which isn't
         * needed and wouldn't work with the way casts are currently
         * implemented if we are unwrapping something derived from a wrapped
         * class).
         */
        if ((addr = sip_api_get_cpp_ptr(sw, NULL)) == NULL)
            return NULL;

        return PyLong_FromVoidPtr(addr);
    }

    return NULL;
}


/*
 * Wrap an instance.
 */
static PyObject *wrapInstance(PyObject *self, PyObject *args)
{
    unsigned PY_LONG_LONG addr;
    sipWrapperType *wt;

    (void)self;

    if (PyArg_ParseTuple(args, "KO!:wrapinstance", &addr, &sipWrapperType_Type, &wt))
        return sip_api_convert_from_type((void *)addr, wt->wt_td, NULL);

    return NULL;
}


/*
 * Set the destroy on exit flag from Python code.
 */
static PyObject *setDestroyOnExit(PyObject *self, PyObject *args)
{
    (void)self;

    if (PyArg_ParseTuple(args, "i:setdestroyonexit", &destroy_on_exit))
    {
        Py_INCREF(Py_None);
        return Py_None;
    }

    return NULL;
}


/*
 * Set the destroy on exit flag from C++ code.
 */
static void sip_api_set_destroy_on_exit(int value)
{
    destroy_on_exit = value;
}


/*
 * Register a client module.  A negative value is returned and an exception
 * raised if there was an error.
 */
static int sip_api_export_module(sipExportedModuleDef *client,
        unsigned api_major, unsigned api_minor, void *unused)
{
    sipExportedModuleDef *em;
    const char *full_name = sipNameOfModule(client);

    (void)unused;

    /* Check that we can support it. */

    if (api_major != SIP_API_MAJOR_NR || api_minor > SIP_API_MINOR_NR)
    {
#if SIP_API_MINOR_NR > 0
        PyErr_Format(PyExc_RuntimeError,
                "the sip module implements API v%d.0 to v%d.%d but the %s module requires API v%d.%d",
                SIP_API_MAJOR_NR, SIP_API_MAJOR_NR, SIP_API_MINOR_NR,
                full_name, api_major, api_minor);
#else
        PyErr_Format(PyExc_RuntimeError,
                "the sip module implements API v%d.0 but the %s module requires API v%d.%d",
                SIP_API_MAJOR_NR, full_name, api_major, api_minor);
#endif

        return -1;
    }

    /* Import any required modules. */
    if (client->em_imports != NULL)
    {
        sipImportedModuleDef *im = client->em_imports;

        while (im->im_name != NULL)
        {
            PyObject *mod;

            if ((mod = PyImport_ImportModule(im->im_name)) == NULL)
                return -1;

            for (em = moduleList; em != NULL; em = em->em_next)
                if (strcmp(sipNameOfModule(em), im->im_name) == 0)
                    break;

            if (em == NULL)
            {
                PyErr_Format(PyExc_RuntimeError,
                        "the %s module failed to register with the sip module",
                        im->im_name);

                return -1;
            }

            if (im->im_imported_types != NULL && importTypes(client, im, em) < 0)
                return -1;

            if (im->im_imported_veh != NULL && importErrorHandlers(client, im, em) < 0)
                return -1;

            if (im->im_imported_exceptions != NULL && importExceptions(client, im, em) < 0)
                return -1;

            ++im;
        }
    }

    for (em = moduleList; em != NULL; em = em->em_next)
    {
        /* SIP clients must have unique names. */
        if (strcmp(sipNameOfModule(em), full_name) == 0)
        {
            PyErr_Format(PyExc_RuntimeError,
                    "the sip module has already registered a module called %s",
                    full_name);

            return -1;
        }

        /* Only one module can claim to wrap QObject. */
        if (em->em_qt_api != NULL && client->em_qt_api != NULL)
        {
            PyErr_Format(PyExc_RuntimeError,
                    "the %s and %s modules both wrap the QObject class",
                    full_name, sipNameOfModule(em));

            return -1;
        }
    }

    /* Convert the module name to an object. */
    if ((client->em_nameobj = PyUnicode_FromString(full_name)) == NULL)
        return -1;

    /* Add it to the list of client modules. */
    client->em_next = moduleList;
    moduleList = client;

    /* Get any keyword handler. */
    if (!got_kw_handler)
    {
        kw_handler = sip_api_import_symbol("pyqt_kw_handler");
        got_kw_handler = TRUE;
    }

    return 0;
}


/*
 * Initialise the contents of a client module.  By this time anything that
 * this depends on should have been initialised.  A negative value is returned
 * and an exception raised if there was an error.
 */
static int sip_api_init_module(sipExportedModuleDef *client,
        PyObject *mod_dict)
{
    sipExportedModuleDef *em;
    sipEnumMemberDef *emd;
    int i;

    /* Handle any API. */
    if (sipInitAPI(client, mod_dict) < 0)
        return -1;

    /* Create the module's types. */
    for (i = 0; i < client->em_nrtypes; ++i)
    {
        sipTypeDef *td = client->em_types[i];

        /* Skip external classes. */
        if (td == NULL)
             continue;

        /* Skip if already initialised. */
        if (td->td_module != NULL)
            continue;

        /* If it is a stub then just set the module so we can get its name. */
        if (sipTypeIsStub(td))
        {
            td->td_module = client;
            continue;
        }

        if (sipTypeIsEnum(td) || sipTypeIsScopedEnum(td))
        {
            sipEnumTypeDef *etd = (sipEnumTypeDef *)td;

            if (td->td_version < 0 || sipIsRangeEnabled(client, td->td_version))
                if (createEnum(client, etd, i, mod_dict) < 0)
                    return -1;

            /*
             * Register the enum pickler for nested enums (non-nested enums
             * don't need special treatment).
             */
            if (sipTypeIsEnum(td) && etd->etd_scope >= 0)
            {
                static PyMethodDef md = {
                    "_pickle_enum", pickle_enum, METH_NOARGS, NULL
                };

                if (setReduce(sipTypeAsPyTypeObject(td), &md) < 0)
                    return -1;
            }
        }
        else if (sipTypeIsMapped(td))
        {
            sipMappedTypeDef *mtd = (sipMappedTypeDef *)td;

            /* If there is a name then we need a namespace. */
            if (mtd->mtd_container.cod_name >= 0)
            {
                if (createMappedType(client, mtd, mod_dict) < 0)
                    return -1;
            }
            else
            {
                td->td_module = client;
            }
        }
        else
        {
            sipClassTypeDef *ctd = (sipClassTypeDef *)td;

            /* See if this is a namespace extender. */
            if (ctd->ctd_container.cod_name < 0)
            {
                sipTypeDef *real_nspace;
                sipClassTypeDef **last;

                ctd->ctd_base.td_module = client;

                real_nspace = getGeneratedType(&ctd->ctd_container.cod_scope,
                        client);

                /* Append this type to the real one. */
                last = &((sipClassTypeDef *)real_nspace)->ctd_nsextender;

                while (*last != NULL)
                    last = &(*last)->ctd_nsextender;

                *last = ctd;

                /*
                 * Save the real namespace type so that it is the correct scope
                 * for any enums or classes defined in this module.
                 */
                client->em_types[i] = real_nspace;
            }
            else if (createClassType(client, ctd, mod_dict) < 0)
                return -1;
        }
    }

    /* Set any Qt support API. */
    if (client->em_qt_api != NULL)
    {
        sipQtSupport = client->em_qt_api;
        sipQObjectType = *sipQtSupport->qt_qobject;
    }

    /* Append any initialiser extenders to the relevant classes. */
    if (client->em_initextend != NULL)
    {
        sipInitExtenderDef *ie = client->em_initextend;

        while (ie->ie_extender != NULL)
        {
            sipTypeDef *td = getGeneratedType(&ie->ie_class, client);
            int enabled;

            if (ie->ie_api_range < 0)
                enabled = TRUE;
            else
                enabled = sipIsRangeEnabled(td->td_module, ie->ie_api_range);

            if (enabled)
            {
                sipWrapperType *wt = (sipWrapperType *)sipTypeAsPyTypeObject(td);

                ie->ie_next = wt->wt_iextend;
                wt->wt_iextend = ie;
            }

            ++ie;
        }
    }

    /* Set the base class object for any sub-class convertors. */
    if (client->em_convertors != NULL)
    {
        sipSubClassConvertorDef *scc = client->em_convertors;

        while (scc->scc_convertor != NULL)
        {
            scc->scc_basetype = getGeneratedType(&scc->scc_base, client);

            ++scc;
        }
    }

    /* Create the module's enum members. */
    for (emd = client->em_enummembers, i = 0; i < client->em_nrenummembers; ++i, ++emd)
    {
        sipTypeDef *etd = client->em_types[emd->em_enum];
        PyObject *mo;

        if (sipTypeIsScopedEnum(etd))
            continue;

        mo = sip_api_convert_from_enum(emd->em_val, etd);

        if (dict_set_and_discard(mod_dict, emd->em_name, mo) < 0)
            return -1;
    }


    /*
     * Add any class static instances.  We need to do this once all types are
     * fully formed because of potential interdependencies.
     */
    for (i = 0; i < client->em_nrtypes; ++i)
    {
        sipTypeDef *td = client->em_types[i];

        if (td != NULL && !sipTypeIsStub(td) && sipTypeIsClass(td))
            if (addInstances((sipTypeAsPyTypeObject(td))->tp_dict, &((sipClassTypeDef *)td)->ctd_container.cod_instances) < 0)
                return -1;
    }

    /* Add any global static instances. */
    if (addInstances(mod_dict, &client->em_instances) < 0)
        return -1;

    /* Add any license. */
    if (client->em_license != NULL && addLicense(mod_dict, client->em_license) < 0)
        return -1;

    /* See if the new module satisfies any outstanding external types. */
    for (em = moduleList; em != NULL; em = em->em_next)
    {
        sipExternalTypeDef *etd;

        if (em == client || em->em_external == NULL)
            continue;

        for (etd = em->em_external; etd->et_nr >= 0; ++etd)
        {
            if (etd->et_name == NULL)
                continue;

            for (i = 0; i < client->em_nrtypes; ++i)
            {
                sipTypeDef *td = client->em_types[i];

                if (td != NULL && !sipTypeIsStub(td) && sipTypeIsClass(td))
                {
                    const char *pyname = sipPyNameOfContainer(
                            &((sipClassTypeDef *)td)->ctd_container, td);

                    if (strcmp(etd->et_name, pyname) == 0)
                    {
                        em->em_types[etd->et_nr] = td;
                        etd->et_name = NULL;

                        break;
                    }
                }
            }
        }
    }

    return 0;
}


/*
 * Called by the interpreter to do any final clearing up, just in case the
 * interpreter will re-start.
 */
static void finalise(void)
{
    sipExportedModuleDef *em;

    /*
     * Mark the Python API as unavailable.  This should already have been done,
     * but just in case...
     */
    sipInterpreter = NULL;

    /* Handle any delayed dtors. */
    for (em = moduleList; em != NULL; em = em->em_next)
        if (em->em_ddlist != NULL)
        {
            em->em_delayeddtors(em->em_ddlist);

            /* Free the list. */
            do
            {
                sipDelayedDtor *dd = em->em_ddlist;

                em->em_ddlist = dd->dd_next;
                sip_api_free(dd);
            }
            while (em->em_ddlist != NULL);
        }

    licenseName = NULL;
    licenseeName = NULL;
    typeName = NULL;
    timestampName = NULL;
    signatureName = NULL;

    /* Release all memory we've allocated directly. */
    sipOMFinalise(&cppPyMap);

    /* Re-initialise those globals that (might) need it. */
    moduleList = NULL;
}


/*
 * Register the given Python type.
 */
static int sip_api_register_py_type(PyTypeObject *type)
{
    return addPyObjectToList(&sipRegisteredPyTypes, (PyObject *)type);
}


/*
 * Find the registered type with the given name.  Raise an exception if it
 * couldn't be found.
 */
static PyObject *findPyType(const char *name)
{
    sipPyObject *po;

    for (po = sipRegisteredPyTypes; po != NULL; po = po->next)
    {
        PyObject *type = po->object;

        if (strcmp(((PyTypeObject *)type)->tp_name, name) == 0)
            return type;
    }

    PyErr_Format(PyExc_RuntimeError, "%s is not a registered type", name);

    return NULL;
}


/*
 * Add a wrapped C/C++ pointer to the list of delayed dtors.
 */
static void sip_api_add_delayed_dtor(sipSimpleWrapper *sw)
{
    void *ptr;
    const sipClassTypeDef *ctd;
    sipExportedModuleDef *em;

    if ((ptr = getPtrTypeDef(sw, &ctd)) == NULL)
        return;

    /* Find the defining module. */
    for (em = moduleList; em != NULL; em = em->em_next)
    {
        int i;

        for (i = 0; i < em->em_nrtypes; ++i)
            if (em->em_types[i] == (const sipTypeDef *)ctd)
            {
                sipDelayedDtor *dd;

                if ((dd = sip_api_malloc(sizeof (sipDelayedDtor))) == NULL)
                    return;

                /* Add to the list. */
                dd->dd_ptr = ptr;
                dd->dd_name = sipPyNameOfContainer(&ctd->ctd_container,
                        (sipTypeDef *)ctd);
                dd->dd_isderived = sipIsDerived(sw);
                dd->dd_next = em->em_ddlist;

                em->em_ddlist = dd;

                return;
            }
    }
}


/*
 * A wrapper around the Python memory allocater that will raise an exception if
 * if the allocation fails.
 */
void *sip_api_malloc(size_t nbytes)
{
    void *mem;

    if ((mem = PyMem_RawMalloc(nbytes)) == NULL)
        PyErr_NoMemory();

    return mem;
}


/*
 * A wrapper around the Python memory de-allocater.
 */
void sip_api_free(void *mem)
{
    PyMem_RawFree(mem);
}


/*
 * Extend a Python slot by looking in other modules to see if there is an
 * extender function that can handle the arguments.
 */
static PyObject *sip_api_pyslot_extend(sipExportedModuleDef *mod,
        sipPySlotType st, const sipTypeDef *td, PyObject *arg0,
        PyObject *arg1)
{
    sipExportedModuleDef *em;

    /* Go through each module. */
    for (em = moduleList; em != NULL; em = em->em_next)
    {
        sipPySlotExtenderDef *ex;

        /* Skip the module that couldn't handle the arguments. */
        if (em == mod)
            continue;

        /* Skip if the module doesn't have any extenders. */
        if (em->em_slotextend == NULL)
            continue;

        /* Go through each extender. */
        for (ex = em->em_slotextend; ex->pse_func != NULL; ++ex)
        {
            PyObject *res;

            /* Skip if not the right slot type. */
            if (ex->pse_type != st)
                continue;

            /* Check against the type if one was given. */
            if (td != NULL && td != getGeneratedType(&ex->pse_class, NULL))
                continue;

            PyErr_Clear();

            res = ((binaryfunc)ex->pse_func)(arg0, arg1);

            if (res != Py_NotImplemented)
                return res;
        }
    }

    /* The arguments couldn't handled anywhere. */
    PyErr_Clear();

    Py_INCREF(Py_NotImplemented);
    return Py_NotImplemented;
}


/*
 * Convert a new C/C++ instance to a Python instance of a specific Python type..
 */
static PyObject *sip_api_convert_from_new_pytype(void *cpp,
        PyTypeObject *py_type, sipWrapper *owner, sipSimpleWrapper **selfp,
        const char *fmt, ...)
{
    PyObject *args, *res;
    va_list va;

    va_start(va, fmt);

    if ((args = PyTuple_New(strlen(fmt))) != NULL && buildObject(args, fmt, va) != NULL)
    {
        res = sipWrapInstance(cpp, py_type, args, owner,
                (selfp != NULL ? SIP_DERIVED_CLASS : 0));

        /* Initialise the rest of an instance of a derived class. */
        if (selfp != NULL)
            *selfp = (sipSimpleWrapper *)res;
    }
    else
    {
        res = NULL;
    }

    Py_XDECREF(args);

    va_end(va);

    return res;
}


/*
 * Call a method and return the result.
 */
static PyObject *call_method(PyObject *method, const char *fmt, va_list va)
{
    PyObject *args, *res;

    if ((args = PyTuple_New(strlen(fmt))) == NULL)
        return NULL;

    if (buildObject(args, fmt, va) != NULL)
        res = PyObject_CallObject(method, args);
    else
        res = NULL;

    Py_DECREF(args);

    return res;
}


/*
 * Call the Python re-implementation of a C++ virtual that does not return a
 * value and handle the result..
 */
static void sip_api_call_procedure_method(sip_gilstate_t gil_state,
        sipVirtErrorHandlerFunc error_handler, sipSimpleWrapper *py_self,
        PyObject *method, const char *fmt, ...)
{
    PyObject *res;
    va_list va;

    va_start(va, fmt);
    res = call_method(method, fmt, va);
    va_end(va);

    if (res != NULL)
    {
        Py_DECREF(res);

        if (res != Py_None)
        {
            sip_api_bad_catcher_result(method);
            res = NULL;
        }
    }

    Py_DECREF(method);

    if (res == NULL)
        sip_api_call_error_handler(error_handler, py_self, gil_state);

    SIP_RELEASE_GIL(gil_state);
}


/*
 * Call the Python re-implementation of a C++ virtual.
 */
static PyObject *sip_api_call_method(int *isErr, PyObject *method,
        const char *fmt, ...)
{
    PyObject *res;
    va_list va;

    va_start(va, fmt);
    res = call_method(method, fmt, va);
    va_end(va);

    if (res == NULL && isErr != NULL)
        *isErr = TRUE;

    return res;
}


/*
 * Build a result object based on a format string.
 */
static PyObject *sip_api_build_result(int *isErr, const char *fmt, ...)
{
    PyObject *res = NULL;
    int badfmt, tupsz;
    va_list va;

    va_start(va,fmt);

    /* Basic validation of the format string. */

    badfmt = FALSE;

    if (*fmt == '(')
    {
        char *ep;

        if ((ep = strchr(fmt,')')) == NULL || ep[1] != '\0')
            badfmt = TRUE;
        else
            tupsz = (int)(ep - fmt - 1);
    }
    else if (strlen(fmt) == 1)
        tupsz = -1;
    else
        badfmt = TRUE;

    if (badfmt)
        PyErr_Format(PyExc_SystemError,"sipBuildResult(): invalid format string \"%s\"",fmt);
    else if (tupsz < 0 || (res = PyTuple_New(tupsz)) != NULL)
        res = buildObject(res,fmt,va);

    va_end(va);

    if (res == NULL && isErr != NULL)
        *isErr = TRUE;

    return res;
}


/*
 * Get the values off the stack and put them into an object.
 */
static PyObject *buildObject(PyObject *obj, const char *fmt, va_list va)
{
    char ch, termch;
    int i;

    /*
     * The format string has already been checked that it is properly formed if
     * it is enclosed in parenthesis.
     */
    if (*fmt == '(')
    {
        termch = ')';
        ++fmt;
    }
    else
        termch = '\0';

    i = 0;

    while ((ch = *fmt++) != termch)
    {
        PyObject *el;

        switch (ch)
        {
        case 'g':
            {
                char *s;
                Py_ssize_t l;

                s = va_arg(va, char *);
                l = va_arg(va, Py_ssize_t);

                if (s != NULL)
                {
                    el = PyBytes_FromStringAndSize(s, l);
                }
                else
                {
                    Py_INCREF(Py_None);
                    el = Py_None;
                }
            }

            break;

        case 'G':
#if defined(HAVE_WCHAR_H)
            {
                wchar_t *s;
                Py_ssize_t l;

                s = va_arg(va, wchar_t *);
                l = va_arg(va, Py_ssize_t);

                if (s != NULL)
                    el = PyUnicode_FromWideChar(s, l);
                else
                {
                    Py_INCREF(Py_None);
                    el = Py_None;
                }
            }
#else
            raiseNoWChar();
            el = NULL;
#endif

            break;

        case 'b':
            el = PyBool_FromLong(va_arg(va,int));
            break;

        case 'c':
            {
                char c = va_arg(va, int);

                el = PyBytes_FromStringAndSize(&c, 1);
            }

            break;

        case 'a':
            {
                char c = va_arg(va, int);

                el = PyUnicode_FromStringAndSize(&c, 1);
            }

            break;

        case 'w':
#if defined(HAVE_WCHAR_H)
            {
                wchar_t c = va_arg(va, int);

                el = PyUnicode_FromWideChar(&c, 1);
            }
#else
            raiseNoWChar();
            el = NULL;
#endif

            break;

        case 'E':
            {
                int ev = va_arg(va, int);
                PyTypeObject *et = va_arg(va, PyTypeObject *);

                el = sip_api_convert_from_enum(ev,
                        ((const sipEnumTypeObject *)et)->type);
            }

            break;

        case 'F':
            {
                int ev = va_arg(va, int);
                const sipTypeDef *td = va_arg(va, const sipTypeDef *);

                el = sip_api_convert_from_enum(ev, td);
            }

            break;

        case 'd':
        case 'f':
            el = PyFloat_FromDouble(va_arg(va, double));
            break;

        case 'e':
        case 'h':
        case 'i':
        case 'L':
            el = PyLong_FromLong(va_arg(va, int));
            break;

        case 'l':
            el = PyLong_FromLong(va_arg(va, long));
            break;

        case 'm':
            el = PyLong_FromUnsignedLong(va_arg(va, unsigned long));
            break;

        case 'n':
#if defined(HAVE_LONG_LONG)
            el = PyLong_FromLongLong(va_arg(va, PY_LONG_LONG));
#else
            el = PyLong_FromLong(va_arg(va, long));
#endif
            break;

        case 'o':
#if defined(HAVE_LONG_LONG)
            el = PyLong_FromUnsignedLongLong(va_arg(va, unsigned PY_LONG_LONG));
#else
            el = PyLong_FromUnsignedLong(va_arg(va, unsigned long));
#endif
            break;

        case 's':
            {
                char *s = va_arg(va, char *);

                if (s != NULL)
                {
                    el = PyBytes_FromString(s);
                }
                else
                {
                    Py_INCREF(Py_None);
                    el = Py_None;
                }
            }

            break;

        case 'A':
            {
                char *s = va_arg(va, char *);

                if (s != NULL)
                {
                    el = PyUnicode_FromString(s);
                }
                else
                {
                    Py_INCREF(Py_None);
                    el = Py_None;
                }
            }

            break;

        case 'x':
#if defined(HAVE_WCHAR_H)
            {
                wchar_t *s = va_arg(va, wchar_t *);

                if (s != NULL)
                    el = PyUnicode_FromWideChar(s, (Py_ssize_t)wcslen(s));
                else
                {
                    Py_INCREF(Py_None);
                    el = Py_None;
                }
            }
#else
            raiseNoWChar();
            el = NULL;
#endif

            break;

        case 't':
        case 'u':
        case 'M':
            el = PyLong_FromUnsignedLong(va_arg(va, unsigned));
            break;

        case '=':
            el = PyLong_FromUnsignedLong(va_arg(va, size_t));
            break;

        case 'B':
            {
                void *p = va_arg(va,void *);
                sipWrapperType *wt = va_arg(va, sipWrapperType *);
                PyObject *xfer = va_arg(va, PyObject *);

                el = sip_api_convert_from_new_type(p, wt->wt_td, xfer);
            }

            break;

        case 'N':
            {
                void *p = va_arg(va, void *);
                const sipTypeDef *td = va_arg(va, const sipTypeDef *);
                PyObject *xfer = va_arg(va, PyObject *);

                el = sip_api_convert_from_new_type(p, td, xfer);
            }

            break;

        case 'C':
            {
                void *p = va_arg(va,void *);
                sipWrapperType *wt = va_arg(va, sipWrapperType *);
                PyObject *xfer = va_arg(va, PyObject *);

                el = sip_api_convert_from_type(p, wt->wt_td, xfer);
            }

            break;

        case 'D':
            {
                void *p = va_arg(va, void *);
                const sipTypeDef *td = va_arg(va, const sipTypeDef *);
                PyObject *xfer = va_arg(va, PyObject *);

                el = sip_api_convert_from_type(p, td, xfer);
            }

            break;

        case 'r':
            {
                void *p = va_arg(va, void *);
                Py_ssize_t l = va_arg(va, Py_ssize_t);
                const sipTypeDef *td = va_arg(va, const sipTypeDef *);

                el = convertToSequence(p, l, td);
            }

            break;

        case 'R':
            el = va_arg(va,PyObject *);
            break;

        case 'S':
            el = va_arg(va,PyObject *);
            Py_INCREF(el);
            break;

        case 'V':
            el = sip_api_convert_from_void_ptr(va_arg(va, void *));
            break;

        case 'z':
            {
                const char *name = va_arg(va, const char *);
                void *p = va_arg(va, void *);

                if (p == NULL)
                {
                    el = Py_None;
                    Py_INCREF(el);
                }
                else
                {
                    el = PyCapsule_New(p, name, NULL);
                }
            }

            break;

        default:
            PyErr_Format(PyExc_SystemError,"buildObject(): invalid format character '%c'",ch);
            el = NULL;
        }

        if (el == NULL)
        {
            Py_XDECREF(obj);
            return NULL;
        }

        if (obj == NULL)
            return el;

        PyTuple_SET_ITEM(obj,i,el);
        ++i;
    }

    return obj;
}


/*
 * Parse a result object based on a format string.  As of v9.0 of the API this
 * is only ever called by handwritten code.
 */
static int sip_api_parse_result(int *isErr, PyObject *method, PyObject *res,
        const char *fmt, ...)
{
    int rc;
    va_list va;

    va_start(va, fmt);
    rc = parseResult(method, res, NULL, fmt, va);
    va_end(va);

    if (isErr != NULL && rc < 0)
        *isErr = TRUE;

    return rc;
}


/*
 * Parse a result object based on a format string.
 */
static int sip_api_parse_result_ex(sip_gilstate_t gil_state,
        sipVirtErrorHandlerFunc error_handler, sipSimpleWrapper *py_self,
        PyObject *method, PyObject *res, const char *fmt, ...)
{
    int rc;

    if (res != NULL)
    {
        va_list va;

        va_start(va, fmt);
        rc = parseResult(method, res, deref_mixin(py_self), fmt, va);
        va_end(va);

        Py_DECREF(res);
    }
    else
    {
        rc = -1;
    }

    Py_DECREF(method);

    if (rc < 0)
        sip_api_call_error_handler(error_handler, py_self, gil_state);

    SIP_RELEASE_GIL(gil_state);

    return rc;
}


/*
 * Call a virtual error handler.  This is called with the GIL and from the
 * thread that raised the error.
 */
static void sip_api_call_error_handler(sipVirtErrorHandlerFunc error_handler,
        sipSimpleWrapper *py_self, sip_gilstate_t sipGILState)
{
    if (error_handler != NULL)
        error_handler(deref_mixin(py_self), sipGILState);
    else
        PyErr_Print();
}


/*
 * Do the main work of parsing a result object based on a format string.
 */
static int parseResult(PyObject *method, PyObject *res,
        sipSimpleWrapper *py_self, const char *fmt, va_list va)
{
    int tupsz, rc = 0;

    /* We rely on PyErr_Occurred(). */
    PyErr_Clear();

    /* Get self if it is provided as an argument. */
    if (*fmt == 'S')
    {
        py_self = va_arg(va, sipSimpleWrapper *);
        ++fmt;
    }

    /* Basic validation of the format string. */
    if (*fmt == '(')
    {
        char ch;
        const char *cp = ++fmt;
        int sub_format = FALSE;

        tupsz = 0;

        while ((ch = *cp++) != ')')
        {
            if (ch == '\0')
            {
                PyErr_Format(PyExc_SystemError, "sipParseResult(): invalid format string \"%s\"", fmt - 1);
                rc = -1;

                break;
            }

            if (sub_format)
            {
                sub_format = FALSE;
            }
            else
            {
                ++tupsz;

                /* Some format characters have a sub-format. */
                if (strchr("aAHDC", ch) != NULL)
                    sub_format = TRUE;
            }
        }

        if (rc == 0)
            if (!PyTuple_Check(res) || PyTuple_GET_SIZE(res) != tupsz)
            {
                sip_api_bad_catcher_result(method);
                rc = -1;
            }
    }
    else
        tupsz = -1;

    if (rc == 0)
    {
        char ch;
        int i = 0;

        while ((ch = *fmt++) != '\0' && ch != ')' && rc == 0)
        {
            PyObject *arg;
            int invalid = FALSE;

            if (tupsz > 0)
            {
                arg = PyTuple_GET_ITEM(res,i);
                ++i;
            }
            else
                arg = res;

            switch (ch)
            {
            case 'g':
                {
                    const char **p = va_arg(va, const char **);
                    Py_ssize_t *szp = va_arg(va, Py_ssize_t *);

                    if (parseBytes_AsCharArray(arg, p, szp) < 0)
                        invalid = TRUE;
                }

                break;

            case 'G':
#if defined(HAVE_WCHAR_H)
                {
                    wchar_t **p = va_arg(va, wchar_t **);
                    Py_ssize_t *szp = va_arg(va, Py_ssize_t *);

                    if (parseWCharArray(arg, p, szp) < 0)
                        invalid = TRUE;
                }
#else
                raiseNoWChar();
                invalid = TRUE;
#endif

                break;

            case 'b':
                {
                    char *p = va_arg(va, void *);
                    int v = sip_api_convert_to_bool(arg);

                    if (v < 0)
                        invalid = TRUE;
                    else if (p != NULL)
                        sipSetBool(p, v);
                }

                break;

            case 'c':
                {
                    char *p = va_arg(va, char *);

                    if (parseBytes_AsChar(arg, p) < 0)
                        invalid = TRUE;
                }

                break;

            case 'a':
                {
                    char *p = va_arg(va, char *);
                    int enc;

                    switch (*fmt++)
                    {
                    case 'A':
                        enc = parseString_AsASCIIChar(arg, p);
                        break;

                    case 'L':
                        enc = parseString_AsLatin1Char(arg, p);
                        break;

                    case '8':
                        enc = parseString_AsUTF8Char(arg, p);
                        break;

                    default:
                        enc = -1;
                    }

                    if (enc < 0)
                        invalid = TRUE;
                }

                break;

            case 'w':
#if defined(HAVE_WCHAR_H)
                {
                    wchar_t *p = va_arg(va, wchar_t *);

                    if (parseWChar(arg, p) < 0)
                        invalid = TRUE;
                }
#else
                raiseNoWChar();
                invalid = TRUE;
#endif

                break;

            case 'd':
                {
                    double *p = va_arg(va, double *);
                    double v = PyFloat_AsDouble(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'E':
                {
                    PyTypeObject *et = va_arg(va, PyTypeObject *);
                    int *p = va_arg(va, int *);
                    int v = sip_api_convert_to_enum(arg, ((sipEnumTypeObject *)et)->type);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'F':
                {
                    sipTypeDef *td = va_arg(va, sipTypeDef *);
                    int *p = va_arg(va, int *);
                    int v = sip_api_convert_to_enum(arg, td);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'f':
                {
                    float *p = va_arg(va, float *);
                    float v = (float)PyFloat_AsDouble(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'L':
                {
                    signed char *p = va_arg(va, signed char *);
                    signed char v = sip_api_long_as_signed_char(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'M':
                {
                    unsigned char *p = va_arg(va, unsigned char *);
                    unsigned char v = sip_api_long_as_unsigned_char(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'h':
                {
                    signed short *p = va_arg(va, signed short *);
                    signed short v = sip_api_long_as_short(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 't':
                {
                    unsigned short *p = va_arg(va, unsigned short *);
                    unsigned short v = sip_api_long_as_unsigned_short(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'e':
                {
                    int *p = va_arg(va, int *);
                    int v = long_as_nonoverflow_int(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'i':
                {
                    int *p = va_arg(va, int *);
                    int v = sip_api_long_as_int(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'u':
                {
                    unsigned *p = va_arg(va, unsigned *);
                    unsigned v = sip_api_long_as_unsigned_int(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case '=':
                {
                    size_t *p = va_arg(va, size_t *);
                    size_t v = sip_api_long_as_size_t(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'l':
                {
                    long *p = va_arg(va, long *);
                    long v = sip_api_long_as_long(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'm':
                {
                    unsigned long *p = va_arg(va, unsigned long *);
                    unsigned long v = sip_api_long_as_unsigned_long(arg);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'n':
                {
#if defined(HAVE_LONG_LONG)
                    PY_LONG_LONG *p = va_arg(va, PY_LONG_LONG *);
                    PY_LONG_LONG v = sip_api_long_as_long_long(arg);
#else
                    long *p = va_arg(va, long *);
                    long v = sip_api_long_as_long(arg);
#endif

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'o':
                {
#if defined(HAVE_LONG_LONG)
                    unsigned PY_LONG_LONG *p = va_arg(va, unsigned PY_LONG_LONG *);
                    unsigned PY_LONG_LONG v = sip_api_long_as_unsigned_long_long(arg);
#else
                    unsigned long *p = va_arg(va, unsigned long *);
                    unsigned long v = sip_api_long_as_unsigned_long(arg);
#endif

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 's':
                {
                    const char **p = va_arg(va, const char **);

                    if (parseBytes_AsString(arg, p) < 0)
                        invalid = TRUE;
                }

                break;

            case 'A':
                {
                    int key = va_arg(va, int);
                    const char **p = va_arg(va, const char **);
                    PyObject *keep;

                    switch (*fmt++)
                    {
                    case 'A':
                        keep = parseString_AsASCIIString(arg, p);
                        break;

                    case 'L':
                        keep = parseString_AsLatin1String(arg, p);
                        break;

                    case '8':
                        keep = parseString_AsUTF8String(arg, p);
                        break;

                    default:
                        keep = NULL;
                    }

                    if (keep == NULL)
                        invalid = TRUE;
                    else
                        sip_api_keep_reference((PyObject *)py_self, key, keep);
                }

                break;

            case 'B':
                {
                    int key = va_arg(va, int);
                    const char **p = va_arg(va, const char **);

                    if (parseBytes_AsString(arg, p) < 0)
                        invalid = TRUE;
                    else
                        sip_api_keep_reference((PyObject *)py_self, key, arg);
                }

                break;

            case 'x':
#if defined(HAVE_WCHAR_H)
                {
                    wchar_t **p = va_arg(va, wchar_t **);

                    if (parseWCharString(arg, p) < 0)
                        invalid = TRUE;
                }
#else
                raiseNoWChar();
                invalid = TRUE;
#endif

                break;

            case 'C':
                {
                    if (*fmt == '\0')
                    {
                        invalid = TRUE;
                    }
                    else
                    {
                        int flags = *fmt++ - '0';
                        int iserr = FALSE;
                        sipWrapperType *type;
                        void **cpp;
                        int *state;

                        type = va_arg(va, sipWrapperType *);

                        if (flags & FMT_RP_NO_STATE_DEPR)
                            state = NULL;
                        else
                            state = va_arg(va, int *);

                        cpp = va_arg(va, void **);

                        *cpp = sip_api_force_convert_to_type(arg, type->wt_td, (flags & FMT_RP_FACTORY ? arg : NULL), (flags & FMT_RP_DEREF ? SIP_NOT_NONE : 0), state, &iserr);

                        if (iserr)
                            invalid = TRUE;
                    }
                }

                break;

            case 'D':
                {
                    if (*fmt == '\0')
                    {
                        invalid = TRUE;
                    }
                    else
                    {
                        int flags = *fmt++ - '0';
                        int iserr = FALSE;
                        const sipTypeDef *td;
                        void **cpp;
                        int *state;

                        td = va_arg(va, const sipTypeDef *);

                        if (flags & FMT_RP_NO_STATE_DEPR)
                            state = NULL;
                        else
                            state = va_arg(va, int *);

                        cpp = va_arg(va, void **);

                        *cpp = sip_api_force_convert_to_type(arg, td, (flags & FMT_RP_FACTORY ? arg : NULL), (flags & FMT_RP_DEREF ? SIP_NOT_NONE : 0), state, &iserr);

                        if (iserr)
                            invalid = TRUE;
                    }
                }

                break;

            case 'H':
                {
                    if (*fmt == '\0')
                    {
                        invalid = TRUE;
                    }
                    else
                    {
                        int flags = *fmt++ - '0';
                        int iserr = FALSE, state;
                        const sipTypeDef *td;
                        void *cpp, *val;

                        td = va_arg(va, const sipTypeDef *);
                        cpp = va_arg(va, void **);

                        val = sip_api_force_convert_to_type(arg, td, (flags & FMT_RP_FACTORY ? arg : NULL), (flags & FMT_RP_DEREF ? SIP_NOT_NONE : 0), &state, &iserr);

                        if (iserr)
                        {
                            invalid = TRUE;
                        }
                        else if (flags & FMT_RP_MAKE_COPY)
                        {
                            sipAssignFunc assign_helper;

                            if (sipTypeIsMapped(td))
                                assign_helper = ((const sipMappedTypeDef *)td)->mtd_assign;
                            else
                                assign_helper = ((const sipClassTypeDef *)td)->ctd_assign;

                            assert(assign_helper != NULL);

                            if (cpp != NULL)
                                assign_helper(cpp, 0, val);

                            sip_api_release_type(val, td, state);
                        }
                        else if (cpp != NULL)
                        {
                            *(void **)cpp = val;
                        }
                    }
                }

                break;

            case 'N':
                {
                    PyTypeObject *type = va_arg(va, PyTypeObject *);
                    PyObject **p = va_arg(va, PyObject **);

                    if (arg == Py_None || PyObject_TypeCheck(arg, type))
                    {
                        if (p != NULL)
                        {
                            Py_INCREF(arg);
                            *p = arg;
                        }
                    }
                    else
                    {
                        invalid = TRUE;
                    }
                }

                break;

            case 'O':
                {
                    PyObject **p = va_arg(va, PyObject **);

                    if (p != NULL)
                    {
                        Py_INCREF(arg);
                        *p = arg;
                    }
                }

                break;

            case 'T':
                {
                    PyTypeObject *type = va_arg(va, PyTypeObject *);
                    PyObject **p = va_arg(va, PyObject **);

                    if (PyObject_TypeCheck(arg, type))
                    {
                        if (p != NULL)
                        {
                            Py_INCREF(arg);
                            *p = arg;
                        }
                    }
                    else
                    {
                        invalid = TRUE;
                    }
                }

                break;

            case 'V':
                {
                    void *v = sip_api_convert_to_void_ptr(arg);
                    void **p = va_arg(va, void **);

                    if (PyErr_Occurred())
                        invalid = TRUE;
                    else if (p != NULL)
                        *p = v;
                }

                break;

            case 'z':
                {
                    const char *name = va_arg(va, const char *);
                    void **p = va_arg(va, void **);

                    if (arg == Py_None)
                    {
                        if (p != NULL)
                            *p = NULL;
                    }
                    else
                    {
#if defined(SIP_USE_CAPSULE)
                        void *v = PyCapsule_GetPointer(arg, name);
#else
                        void *v = sip_api_convert_to_void_ptr(arg);

                        (void)name;
#endif

                        if (PyErr_Occurred())
                            invalid = TRUE;
                        else if (p != NULL)
                            *p = v;
                    }
                }

                break;

            case 'Z':
                if (arg != Py_None)
                    invalid = TRUE;

                break;

            case '!':
                {
                    PyObject **p = va_arg(va, PyObject **);

                    if (PyObject_CheckBuffer(arg))
                    {
                        if (p != NULL)
                        {
                            Py_INCREF(arg);
                            *p = arg;
                        }
                    }
                    else
                    {
                        invalid = TRUE;
                    }
                }

                break;

            case '$':
                {
                    PyObject **p = va_arg(va, PyObject **);

                    if (arg == Py_None || PyObject_CheckBuffer(arg))
                    {
                        if (p != NULL)
                        {
                            Py_INCREF(arg);
                            *p = arg;
                        }
                    }
                    else
                    {
                        invalid = TRUE;
                    }
                }

                break;

            default:
                PyErr_Format(PyExc_SystemError,"sipParseResult(): invalid format character '%c'",ch);
                rc = -1;
            }

            if (invalid)
            {
                sip_api_bad_catcher_result(method);
                rc = -1;
                break;
            }
        }
    }

    return rc;
}


/*
 * Parse the arguments to a C/C++ function without any side effects.
 */
static int sip_api_parse_args(PyObject **parseErrp, PyObject *sipArgs,
        const char *fmt, ...)
{
    int ok;
    va_list va;

    va_start(va, fmt);
    ok = parseKwdArgs(parseErrp, sipArgs, NULL, NULL, NULL, fmt, va);
    va_end(va);

    return ok;
}


/*
 * Parse the positional and/or keyword arguments to a C/C++ function without
 * any side effects.
 */
static int sip_api_parse_kwd_args(PyObject **parseErrp, PyObject *sipArgs,
        PyObject *sipKwdArgs, const char **kwdlist, PyObject **unused,
        const char *fmt, ...)
{
    int ok;
    va_list va;

    if (unused != NULL)
    {
        /*
         * Initialise the return of any unused keyword arguments.  This is
         * used by any ctor overload.
         */
        *unused = NULL;
    }

    va_start(va, fmt);
    ok = parseKwdArgs(parseErrp, sipArgs, sipKwdArgs, kwdlist, unused, fmt,
            va);
    va_end(va);

    /* Release any unused arguments if the parse failed. */
    if (!ok && unused != NULL)
    {
        Py_XDECREF(*unused);
    }

    return ok;
}


/*
 * Parse the arguments to a C/C++ function without any side effects.
 */
static int parseKwdArgs(PyObject **parseErrp, PyObject *sipArgs,
        PyObject *sipKwdArgs, const char **kwdlist, PyObject **unused,
        const char *fmt, va_list va_orig)
{
    int no_tmp_tuple, ok, selfarg;
    sipSimpleWrapper *self;
    PyObject *single_arg;
    va_list va;

    /* Previous second pass errors stop subsequent parses. */
    if (*parseErrp != NULL && !PyList_Check(*parseErrp))
        return FALSE;

    /*
     * See if we are parsing a single argument.  In current versions we are
     * told explicitly by the first character of the format string.  In earlier
     * versions we guessed (sometimes wrongly).
     */
    if (*fmt == '1')
    {
        ++fmt;
        no_tmp_tuple = FALSE;
    }
    else
        no_tmp_tuple = PyTuple_Check(sipArgs);

    if (no_tmp_tuple)
    {
        Py_INCREF(sipArgs);
    }
    else if ((single_arg = PyTuple_New(1)) != NULL)
    {
        Py_INCREF(sipArgs);
        PyTuple_SET_ITEM(single_arg, 0, sipArgs);

        sipArgs = single_arg;
    }
    else
    {
        /* Stop all parsing and indicate an exception has been raised. */
        Py_XDECREF(*parseErrp);
        *parseErrp = Py_None;
        Py_INCREF(Py_None);

        return FALSE;
    }

    /*
     * The first pass checks all the types and does conversions that are cheap
     * and have no side effects.
     */
    va_copy(va, va_orig);
    ok = parsePass1(parseErrp, &self, &selfarg, sipArgs, sipKwdArgs, kwdlist,
            unused, fmt, va);
    va_end(va);

    if (ok)
    {
        /*
         * The second pass does any remaining conversions now that we know we
         * have the right signature.
         */
        va_copy(va, va_orig);
        ok = parsePass2(self, selfarg, sipArgs, sipKwdArgs, kwdlist, fmt, va);
        va_end(va);

        /* Remove any previous failed parses. */
        Py_XDECREF(*parseErrp);

        if (ok)
        {
            *parseErrp = NULL;
        }
        else
        {
            /* Indicate that an exception has been raised. */
            *parseErrp = Py_None;
            Py_INCREF(Py_None);
        }
    }

    Py_DECREF(sipArgs);

    return ok;
}


/*
 * Return a string as a Python object that describes an argument with an
 * unexpected type.
 */
static PyObject *bad_type_str(int arg_nr, PyObject *arg)
{
    return PyUnicode_FromFormat("argument %d has unexpected type '%s'", arg_nr,
            Py_TYPE(arg)->tp_name);
}


/*
 * Adds a failure about an argument with an incorrect type to the current list
 * of exceptions.
 */
static sipErrorState sip_api_bad_callable_arg(int arg_nr, PyObject *arg)
{
    PyObject *detail = bad_type_str(arg_nr + 1, arg);

    if (detail == NULL)
        return sipErrorFail;

    PyErr_SetObject(PyExc_TypeError, detail);
    Py_DECREF(detail);

    return sipErrorContinue;
}


/*
 * Adds the current exception to the current list of exceptions (if it is a
 * user exception) or replace the current list of exceptions.
 */
static void sip_api_add_exception(sipErrorState es, PyObject **parseErrp)
{
    assert(*parseErrp == NULL);

    if (es == sipErrorContinue)
    {
        sipParseFailure failure;
        PyObject *e_type, *e_traceback;

        /* Get the value of the exception. */
        PyErr_Fetch(&e_type, &failure.detail_obj, &e_traceback);
        Py_XDECREF(e_type);
        Py_XDECREF(e_traceback);

        failure.reason = Exception;

        add_failure(parseErrp, &failure);

        if (failure.reason == Raised)
        {
            Py_XDECREF(failure.detail_obj);
            es = sipErrorFail;
        }
    }

    if (es == sipErrorFail)
    {
        Py_XDECREF(*parseErrp);
        *parseErrp = Py_None;
        Py_INCREF(Py_None);
    }
}


/*
 * The dtor for parse failure wrapped in a Python object.
 */
static void failure_dtor(PyObject *capsule)
{
    sipParseFailure *failure = (sipParseFailure *)PyCapsule_GetPointer(capsule, NULL);

    Py_XDECREF(failure->detail_obj);

    sip_api_free(failure);
}


/*
 * Add a parse failure to the current list of exceptions.
 */
static void add_failure(PyObject **parseErrp, sipParseFailure *failure)
{
    sipParseFailure *failure_copy;
    PyObject *failure_obj;

    /* Create the list if necessary. */
    if (*parseErrp == NULL && (*parseErrp = PyList_New(0)) == NULL)
    {
        failure->reason = Raised;
        return;
    }

    /*
     * Make a copy of the failure, convert it to a Python object and add it to
     * the list.  We do it this way to make it as lightweight as possible.
     */
    if ((failure_copy = sip_api_malloc(sizeof (sipParseFailure))) == NULL)
    {
        failure->reason = Raised;
        return;
    }

    *failure_copy = *failure;

    if ((failure_obj = PyCapsule_New(failure_copy, NULL, failure_dtor)) == NULL)
    {
        sip_api_free(failure_copy);
        failure->reason = Raised;
        return;
    }

    /* Ownership of any detail object is now with the wrapped failure. */
    failure->detail_obj = NULL;

    if (PyList_Append(*parseErrp, failure_obj) < 0)
    {
        Py_DECREF(failure_obj);
        failure->reason = Raised;
        return;
    }

    Py_DECREF(failure_obj);
}


/*
 * Parse one or a pair of arguments to a C/C++ function without any side
 * effects.
 */
static int sip_api_parse_pair(PyObject **parseErrp, PyObject *sipArg0,
        PyObject *sipArg1, const char *fmt, ...)
{
    int ok, selfarg;
    sipSimpleWrapper *self;
    PyObject *args;
    va_list va;

    /* Previous second pass errors stop subsequent parses. */
    if (*parseErrp != NULL && !PyList_Check(*parseErrp))
        return FALSE;

    if ((args = PyTuple_New(sipArg1 != NULL ? 2 : 1)) == NULL)
    {
        /* Stop all parsing and indicate an exception has been raised. */
        Py_XDECREF(*parseErrp);
        *parseErrp = Py_None;
        Py_INCREF(Py_None);

        return FALSE;
    }

    Py_INCREF(sipArg0);
    PyTuple_SET_ITEM(args, 0, sipArg0);

    if (sipArg1 != NULL)
    {
        Py_INCREF(sipArg1);
        PyTuple_SET_ITEM(args, 1, sipArg1);
    }

    /*
     * The first pass checks all the types and does conversions that are cheap
     * and have no side effects.
     */
    va_start(va, fmt);
    ok = parsePass1(parseErrp, &self, &selfarg, args, NULL, NULL, NULL, fmt,
            va);
    va_end(va);

    if (ok)
    {
        /*
         * The second pass does any remaining conversions now that we know we
         * have the right signature.
         */
        va_start(va, fmt);
        ok = parsePass2(self, selfarg, args, NULL, NULL, fmt, va);
        va_end(va);

        /* Remove any previous failed parses. */
        Py_XDECREF(*parseErrp);

        if (ok)
        {
            *parseErrp = NULL;
        }
        else
        {
            /* Indicate that an exception has been raised. */
            *parseErrp = Py_None;
            Py_INCREF(Py_None);
        }
    }

    Py_DECREF(args);

    return ok;
}


/*
 * First pass of the argument parse, converting those that can be done so
 * without any side effects.  Return TRUE if the arguments matched.
 */
static int parsePass1(PyObject **parseErrp, sipSimpleWrapper **selfp,
        int *selfargp, PyObject *sipArgs, PyObject *sipKwdArgs,
        const char **kwdlist, PyObject **unused, const char *fmt, va_list va)
{
    int compulsory, argnr, nr_args;
    Py_ssize_t nr_pos_args, nr_kwd_args, nr_kwd_args_used;
    sipParseFailure failure;

    failure.reason = Ok;
    failure.detail_obj = NULL;
    compulsory = TRUE;
    argnr = 0;
    nr_args = 0;
    nr_pos_args = PyTuple_GET_SIZE(sipArgs);
    nr_kwd_args = nr_kwd_args_used = 0;

    if (sipKwdArgs != NULL)
    {
        assert(PyDict_Check(sipKwdArgs));

        nr_kwd_args = PyDict_Size(sipKwdArgs);
    }

    /*
     * Handle those format characters that deal with the "self" argument.  They
     * will always be the first one.
     */
    *selfp = NULL;
    *selfargp = FALSE;

    switch (*fmt++)
    {
    case '#':
            /* A ctor has an argument with the /Transfer/ annotation. */
            *selfp = va_arg(va, PyObject *);
            break;

    case 'B':
    case 'p':
        {
            PyObject *self;
            sipTypeDef *td;

            self = *va_arg(va, PyObject **);
            td = va_arg(va, sipTypeDef *);
            va_arg(va, void **);

            if (self == NULL)
            {
                if (!getSelfFromArgs(td, sipArgs, argnr, selfp))
                {
                    failure.reason = Unbound;
                    failure.detail_str = sipPyNameOfContainer(
                            &((sipClassTypeDef *)td)->ctd_container, td);
                    break;
                }

                *selfargp = TRUE;
                ++argnr;
            }
            else
                *selfp = (sipSimpleWrapper *)self;

            break;
        }

    case 'C':
        *selfp = (sipSimpleWrapper *)va_arg(va,PyObject *);
        break;

    default:
        --fmt;
    }

    /*
     * Now handle the remaining arguments.  We continue to parse if we get an
     * overflow because that is, strictly speaking, a second pass error.
     */
    while (failure.reason == Ok || failure.reason == Overflow)
    {
        char ch;
        PyObject *arg;

        PyErr_Clear();

        /* See if the following arguments are optional. */
        if ((ch = *fmt++) == '|')
        {
            compulsory = FALSE;
            ch = *fmt++;
        }

        /* See if we don't expect anything else. */

        if (ch == '\0')
        {
            if (argnr < nr_pos_args)
            {
                /* There are still positional arguments. */
                failure.reason = TooMany;
            }
            else if (nr_kwd_args_used != nr_kwd_args)
            {
                /*
                 * Take a shortcut if no keyword arguments were used and we are
                 * interested in them.
                 */
                if (nr_kwd_args_used == 0 && unused != NULL)
                {
                    Py_INCREF(sipKwdArgs);
                    *unused = sipKwdArgs;
                }
                else
                {
                    PyObject *key, *value, *unused_dict = NULL;
                    Py_ssize_t pos = 0;

                    /*
                     * Go through the keyword arguments to find any that were
                     * duplicates of positional arguments.  For the remaining
                     * ones remember the unused ones if we are interested.
                     */
                    while (PyDict_Next(sipKwdArgs, &pos, &key, &value))
                    {
                        int a;

                        if (!PyUnicode_Check(key))
                        {
                            failure.reason = KeywordNotString;
                            failure.detail_obj = key;
                            Py_INCREF(key);
                            break;
                        }

                        if (kwdlist != NULL)
                        {
                            /* Get the argument's index if it is one. */
                            for (a = 0; a < nr_args; ++a)
                            {
                                const char *name = kwdlist[a];

                                if (name == NULL)
                                    continue;

                                if (PyUnicode_CompareWithASCIIString(key, name) == 0)
                                    break;
                            }
                        }
                        else
                        {
                            a = nr_args;
                        }

                        if (a == nr_args)
                        {
                            /*
                             * The name doesn't correspond to a keyword
                             * argument.
                             */
                            if (unused == NULL)
                            {
                                /*
                                 * It may correspond to a keyword argument of a
                                 * different overload.
                                 */
                                failure.reason = UnknownKeyword;
                                failure.detail_obj = key;
                                Py_INCREF(key);

                                break;
                            }

                            /*
                             * Add it to the dictionary of unused arguments
                             * creating it if necessary.  Note that if the
                             * unused arguments are actually used by a later
                             * overload then the parse will incorrectly
                             * succeed.  This should be picked up (perhaps with
                             * a misleading exception) so long as the code that
                             * handles the unused arguments checks that it can
                             * handle them all.
                             */
                            if (unused_dict == NULL && (*unused = unused_dict = PyDict_New()) == NULL)
                            {
                                failure.reason = Raised;
                                break;
                            }

                            if (PyDict_SetItem(unused_dict, key, value) < 0)
                            {
                                failure.reason = Raised;
                                break;
                            }
                        }
                        else if (a < nr_pos_args - *selfargp)
                        {
                            /*
                             * The argument has been given positionally and as
                             * a keyword.
                             */
                            failure.reason = Duplicate;
                            failure.detail_obj = key;
                            Py_INCREF(key);
                            break;
                        }
                    }
                }
            }

            break;
        }

        /* Get the next argument. */
        arg = NULL;
        failure.arg_nr = -1;
        failure.arg_name = NULL;

        if (argnr < nr_pos_args)
        {
            arg = PyTuple_GET_ITEM(sipArgs, argnr);
            failure.arg_nr = argnr + 1;
        }
        else if (nr_kwd_args != 0 && kwdlist != NULL)
        {
            const char *name = kwdlist[argnr - *selfargp];

            if (name != NULL)
            {
                arg = PyDict_GetItemString(sipKwdArgs, name);

                if (arg != NULL)
                    ++nr_kwd_args_used;

                failure.arg_name = name;
            }
        }

        ++argnr;
        ++nr_args;

        if (arg == NULL && compulsory)
        {
            if (ch == 'W')
            {
                /*
                 * A variable number of arguments was allowed but none were
                 * given.
                 */
                break;
            }

            /* An argument was required. */
            failure.reason = TooFew;

            /*
             * Check if there were any unused keyword arguments so that we give
             * a (possibly) more accurate diagnostic in the case that a keyword
             * argument has been mis-spelled.
             */
            if (unused == NULL && sipKwdArgs != NULL && nr_kwd_args_used != nr_kwd_args)
            {
                PyObject *key, *value;
                Py_ssize_t pos = 0;

                while (PyDict_Next(sipKwdArgs, &pos, &key, &value))
                {
                    int a;

                    if (!PyUnicode_Check(key))
                    {
                        failure.reason = KeywordNotString;
                        failure.detail_obj = key;
                        Py_INCREF(key);
                        break;
                    }

                    if (kwdlist != NULL)
                    {
                        /* Get the argument's index if it is one. */
                        for (a = 0; a < nr_args; ++a)
                        {
                            const char *name = kwdlist[a];

                            if (name == NULL)
                                continue;

                            if (PyUnicode_CompareWithASCIIString(key, name) == 0)
                                break;
                        }
                    }
                    else
                    {
                        a = nr_args;
                    }

                    if (a == nr_args)
                    {
                        failure.reason = UnknownKeyword;
                        failure.detail_obj = key;
                        Py_INCREF(key);

                        break;
                    }
                }
            }

            break;
        }

        /*
         * Handle the format character even if we don't have an argument so
         * that we skip the right number of arguments.
         */
        switch (ch)
        {
        case 'W':
            /* Ellipsis. */
            break;

        case '@':
            {
                /* Implement /GetWrapper/. */

                PyObject **p = va_arg(va, PyObject **);

                if (arg != NULL)
                    *p = arg;

                /* Process the same argument next time round. */
                --argnr;
                --nr_args;

                break;
            }

        case 's':
            {
                /* String from a Python bytes or None. */

                const char **p = va_arg(va, const char **);

                if (arg != NULL && parseBytes_AsString(arg, p) < 0)
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'A':
            {
                /* String from a Python string or None. */

                va_arg(va, PyObject **);
                va_arg(va, const char **);
                fmt++;

                if (arg != NULL && check_encoded_string(arg) < 0)
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'a':
            {
                /* Character from a Python string. */

                va_arg(va, char *);
                fmt++;

                if (arg != NULL && check_encoded_string(arg) < 0)
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'x':
#if defined(HAVE_WCHAR_H)
            {
                /* Wide string or None. */

                wchar_t **p = va_arg(va, wchar_t **);

                if (arg != NULL && parseWCharString(arg, p) < 0)
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }
#else
            raiseNoWChar();
            failure.reason = Raised;
            break;
#endif

        case 'U':
            {
                /* Slot name or callable, return the name or callable. */

                char **sname = va_arg(va, char **);
                PyObject **scall = va_arg(va, PyObject **);

                if (arg != NULL)
                {
                    *sname = NULL;
                    *scall = NULL;

                    if (PyBytes_Check(arg))
                    {
                        char *s = PyBytes_AS_STRING(arg);

                        if (*s == '1' || *s == '2' || *s == '9')
                        {
                            *sname = s;
                        }
                        else
                        {
                            failure.reason = WrongType;
                            failure.detail_obj = arg;
                            Py_INCREF(arg);
                        }
                    }
                    else if (PyCallable_Check(arg))
                    {
                         *scall = arg;
                    }
                    else if (arg != Py_None)
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }

                break;
            }

        case 'S':
            {
                /* Slot name, return the name. */

                char **p = va_arg(va, char **);

                if (arg != NULL)
                {
                    if (PyBytes_Check(arg))
                    {
                        char *s = PyBytes_AS_STRING(arg);

                        if (*s == '1' || *s == '2' || *s == '9')
                        {
                            *p = s;
                        }
                        else
                        {
                            failure.reason = WrongType;
                            failure.detail_obj = arg;
                            Py_INCREF(arg);
                        }
                    }
                    else
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }

                break;
            }

        case 'G':
            {
                /* Signal name, return the name. */

                char **p = va_arg(va, char **);

                if (arg != NULL)
                {
                    if (PyBytes_Check(arg))
                    {
                        char *s = PyBytes_AS_STRING(arg);

                        if (*s == '2' || *s == '9')
                        {
                            *p = s;
                        }
                        else
                        {
                            failure.reason = WrongType;
                            failure.detail_obj = arg;
                            Py_INCREF(arg);
                        }
                    }
                    else
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }

                break;
            }

        case 'r':
            {
                /*
                 * Sequence of mapped type instances.  For ABI v12.10 and
                 * earlier this is also used for class instances.
                 */

                const sipTypeDef *td;

                td = va_arg(va, const sipTypeDef *);
                va_arg(va, void **);
                va_arg(va, Py_ssize_t *);

                if (arg != NULL && !canConvertFromSequence(arg, td))
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case '>':
            {
                /*
                 * Sequence or sip.array of class instances.  This is only used
                 * by ABI v12.11 and later.
                 */

                const sipTypeDef *td;

                td = va_arg(va, const sipTypeDef *);
                va_arg(va, void **);
                va_arg(va, Py_ssize_t *);
                va_arg(va, int *);

                if (arg != NULL && !sip_array_can_convert(arg, td) && !canConvertFromSequence(arg, td))
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'J':
            {
                /* Class or mapped type instance. */

                char sub_fmt = *fmt++;
                const sipTypeDef *td;
                int flags = sub_fmt - '0';
                int iflgs = 0;

                td = va_arg(va, const sipTypeDef *);
                va_arg(va, void **);

                if (flags & FMT_AP_DEREF)
                    iflgs |= SIP_NOT_NONE;

                if (flags & FMT_AP_TRANSFER_THIS)
                    va_arg(va, PyObject **);

                if (flags & FMT_AP_NO_CONVERTORS)
                    iflgs |= SIP_NO_CONVERTORS;
                else
                    va_arg(va, int *);

                if (arg != NULL && !sip_api_can_convert_to_type(arg, td, iflgs))
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'N':
            {
                /* Python object of given type or None. */

                PyTypeObject *type = va_arg(va,PyTypeObject *);
                PyObject **p = va_arg(va,PyObject **);

                if (arg != NULL)
                {
                    if (arg == Py_None || PyObject_TypeCheck(arg,type))
                    {
                        *p = arg;
                    }
                    else
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }

                break;
            }

        case 'P':
            {
                /* Python object of any type with a sub-format. */

                va_arg(va, PyObject **);

                /* Skip the sub-format. */
                ++fmt;

                break;
            }

        case 'T':
            {
                /* Python object of given type. */

                PyTypeObject *type = va_arg(va, PyTypeObject *);
                PyObject **p = va_arg(va, PyObject **);

                if (arg != NULL)
                {
                    if (PyObject_TypeCheck(arg,type))
                    {
                        *p = arg;
                    }
                    else
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }

                break;
            }

        case 'R':
            {
                /* Sub-class of QObject. */

                PyObject **p = va_arg(va, PyObject **);

                if (arg != NULL)
                {
                    if (isQObject(arg))
                    {
                        *p = arg;
                    }
                    else
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }

                break;
            }

        case 'F':
            {
                /* Python callable object. */
 
                PyObject **p = va_arg(va, PyObject **);

                if (arg != NULL)
                {
                    if (PyCallable_Check(arg))
                    {
                        *p = arg;
                    }
                    else
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }
 
                break;
            }

        case 'H':
            {
                /* Python callable object or None. */
 
                PyObject **p = va_arg(va, PyObject **);

                if (arg != NULL)
                {
                    if (arg == Py_None || PyCallable_Check(arg))
                    {
                        *p = arg;
                    }
                    else
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }
 
                break;
            }

        case '!':
            {
                /* Python object that implements the buffer protocol. */
 
                PyObject **p = va_arg(va, PyObject **);

                if (arg != NULL)
                {
                    if (PyObject_CheckBuffer(arg))
                    {
                        *p = arg;
                    }
                    else
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }
 
                break;
            }

        case '$':
            {
                /*
                 * Python object that implements the buffer protocol or None.
                 */
 
                PyObject **p = va_arg(va, PyObject **);

                if (arg != NULL)
                {
                    if (arg == Py_None || PyObject_CheckBuffer(arg))
                    {
                        *p = arg;
                    }
                    else
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                }
 
                break;
            }

        case 'q':
            {
                /* Qt receiver to connect. */

                va_arg(va, char *);
                va_arg(va, void **);
                va_arg(va, const char **);

                if (arg != NULL && !isQObject(arg))
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'Q':
            {
                /* Qt receiver to disconnect. */

                va_arg(va, char *);
                va_arg(va, void **);
                va_arg(va, const char **);

                if (arg != NULL && !isQObject(arg))
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'g':
        case 'y':
            {
                /* Python slot to connect. */

                va_arg(va, char *);
                va_arg(va, void **);
                va_arg(va, const char **);

                if (arg != NULL && (sipQtSupport == NULL || !PyCallable_Check(arg)))
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'Y':
            {
                /* Python slot to disconnect. */

                va_arg(va, char *);
                va_arg(va, void **);
                va_arg(va, const char **);

                if (arg != NULL && (sipQtSupport == NULL || !PyCallable_Check(arg)))
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'k':
            {
                /* Char array or None. */

                const char **p = va_arg(va, const char **);
                Py_ssize_t *szp = va_arg(va, Py_ssize_t *);

                if (arg != NULL && parseBytes_AsCharArray(arg, p, szp) < 0)
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'K':
#if defined(HAVE_WCHAR_H)
            {
                /* Wide char array or None. */

                wchar_t **p = va_arg(va, wchar_t **);
                Py_ssize_t *szp = va_arg(va, Py_ssize_t *);

                if (arg != NULL && parseWCharArray(arg, p, szp) < 0)
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }
#else
            raiseNoWChar();
            failure.reason = Raised;
            break
#endif

        case 'c':
            {
                /* Character from a Python bytes. */

                char *p = va_arg(va, char *);

                if (arg != NULL && parseBytes_AsChar(arg, p) < 0)
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }

        case 'w':
#if defined(HAVE_WCHAR_H)
            {
                /* Wide character. */

                wchar_t *p = va_arg(va, wchar_t *);

                if (arg != NULL && parseWChar(arg, p) < 0)
                {
                    failure.reason = WrongType;
                    failure.detail_obj = arg;
                    Py_INCREF(arg);
                }

                break;
            }
#else
            raiseNoWChar();
            failure.reason = Raised;
            break
#endif

        case 'b':
            {
                /* Bool. */

                void *p = va_arg(va, void *);

                if (arg != NULL)
                {
                    int v = sip_api_convert_to_bool(arg);

                    if (v < 0)
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                    else
                    {
                        sipSetBool(p, v);
                    }
                }

                break;
            }

        case 'E':
            {
                /* Named enum or integer. */

                sipTypeDef *td = va_arg(va, sipTypeDef *);
                int *p = va_arg(va, int *);

                if (arg != NULL)
                {
                    int v = sip_api_convert_to_enum(arg, td);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }
            }

            break;

        case 'e':
            {
                /* Anonymous enum. */

                int *p = va_arg(va, int *);

                if (arg != NULL)
                {
                    int v = long_as_nonoverflow_int(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }
            }

            break;

        case 'i':
            {
                /* Integer. */

                int *p = va_arg(va, int *);

                if (arg != NULL)
                {
                    int v = sip_api_long_as_int(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 'u':
            {
                /* Unsigned integer. */

                unsigned *p = va_arg(va, unsigned *);

                if (arg != NULL)
                {
                    unsigned v = sip_api_long_as_unsigned_int(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case '=':
            {
                /* size_t integer. */

                size_t *p = va_arg(va, size_t *);

                if (arg != NULL)
                {
                    size_t v = sip_api_long_as_size_t(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 'L':
            {
                /* Signed char. */

                signed char *p = va_arg(va, signed char *);

                if (arg != NULL)
                {
                    signed char v = sip_api_long_as_signed_char(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 'M':
            {
                /* Unsigned char. */

                unsigned char *p = va_arg(va, unsigned char *);

                if (arg != NULL)
                {
                    unsigned char v = sip_api_long_as_unsigned_char(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 'h':
            {
                /* Short integer. */

                signed short *p = va_arg(va, signed short *);

                if (arg != NULL)
                {
                    signed short v = sip_api_long_as_short(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 't':
            {
                /* Unsigned short integer. */

                unsigned short *p = va_arg(va, unsigned short *);

                if (arg != NULL)
                {
                    unsigned short v = sip_api_long_as_unsigned_short(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 'l':
            {
                /* Long integer. */

                long *p = va_arg(va, long *);

                if (arg != NULL)
                {
                    long v = sip_api_long_as_long(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 'm':
            {
                /* Unsigned long integer. */

                unsigned long *p = va_arg(va, unsigned long *);

                if (arg != NULL)
                {
                    unsigned long v = sip_api_long_as_unsigned_long(arg);

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 'n':
            {
                /* Long long integer. */

#if defined(HAVE_LONG_LONG)
                PY_LONG_LONG *p = va_arg(va, PY_LONG_LONG *);
#else
                long *p = va_arg(va, long *);
#endif

                if (arg != NULL)
                {
#if defined(HAVE_LONG_LONG)
                    PY_LONG_LONG v = sip_api_long_as_long_long(arg);
#else
                    long v = sip_api_long_as_long(arg);
#endif

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 'o':
            {
                /* Unsigned long long integer. */

#if defined(HAVE_LONG_LONG)
                unsigned PY_LONG_LONG *p = va_arg(va, unsigned PY_LONG_LONG *);
#else
                unsigned long *p = va_arg(va, unsigned long *);
#endif

                if (arg != NULL)
                {
#if defined(HAVE_LONG_LONG)
                    unsigned PY_LONG_LONG v = sip_api_long_as_unsigned_long_long(arg);
#else
                    unsigned long v = sip_api_long_as_unsigned_long(arg);
#endif

                    if (PyErr_Occurred())
                        handle_failed_int_conversion(&failure, arg);
                    else
                        *p = v;
                }

                break;
            }

        case 'f':
            {
                /* Float. */

                float *p = va_arg(va, float *);

                if (arg != NULL)
                {
                    double v = PyFloat_AsDouble(arg);

                    if (PyErr_Occurred())
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                    else
                    {
                        *p = (float)v;
                    }
                }

                break;
            }

        case 'X':
            {
                /* Constrained types. */

                char sub_fmt = *fmt++;

                if (sub_fmt == 'E')
                {
                    /* Named enum. */

                    sipTypeDef *td = va_arg(va, sipTypeDef *);
                    int *p = va_arg(va, int *);

                    if (arg != NULL)
                    {
                        *p = convert_to_enum(arg, td, FALSE);

                        if (PyErr_Occurred())
                            handle_failed_int_conversion(&failure, arg);
                    }
                }
                else
                {
                    void *p = va_arg(va, void *);

                    if (arg != NULL)
                    {
                        switch (sub_fmt)
                        {
                        case 'b':
                            {
                                /* Boolean. */

                                if (PyBool_Check(arg))
                                {
                                    sipSetBool(p, (arg == Py_True));
                                }
                                else
                                {
                                    failure.reason = WrongType;
                                    failure.detail_obj = arg;
                                    Py_INCREF(arg);
                                }

                                break;
                            }

                        case 'd':
                            {
                                /* Double float. */

                                if (PyFloat_Check(arg))
                                {
                                    *(double *)p = PyFloat_AS_DOUBLE(arg);
                                }
                                else
                                {
                                    failure.reason = WrongType;
                                    failure.detail_obj = arg;
                                    Py_INCREF(arg);
                                }

                                break;
                            }

                        case 'f':
                            {
                                /* Float. */

                                if (PyFloat_Check(arg))
                                {
                                    *(float *)p = (float)PyFloat_AS_DOUBLE(arg);
                                }
                                else
                                {
                                    failure.reason = WrongType;
                                    failure.detail_obj = arg;
                                    Py_INCREF(arg);
                                }

                                break;
                            }

                        case 'i':
                            {
                                /* Integer. */

                                if (PyLong_Check(arg))
                                {
                                    *(int *)p = sip_api_long_as_int(arg);

                                    if (PyErr_Occurred())
                                        handle_failed_int_conversion(&failure,
                                                arg);
                                }
                                else
                                {
                                    failure.reason = WrongType;
                                    failure.detail_obj = arg;
                                    Py_INCREF(arg);
                                }

                                break;
                            }
                        }
                    }
                }

                break;
            }

        case 'd':
            {
                /* Double float. */

                double *p = va_arg(va,double *);

                if (arg != NULL)
                {
                    double v = PyFloat_AsDouble(arg);

                    if (PyErr_Occurred())
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                    else
                    {
                        *p = v;
                    }
                }

                break;
            }

        case 'v':
            {
                /* Void pointer. */

                void **p = va_arg(va, void **);

                if (arg != NULL)
                {
                    void *v = sip_api_convert_to_void_ptr(arg);

                    if (PyErr_Occurred())
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                    else
                    {
                        *p = v;
                    }
                }

                break;
            }

        case 'z':
            {
                /* Void pointer as a capsule. */

                const char *name = va_arg(va, const char *);
                void **p = va_arg(va, void **);

                if (arg == Py_None)
                {
                    *p = NULL;
                }
                else if (arg != NULL)
                {
                    void *v = PyCapsule_GetPointer(arg, name);

                    if (PyErr_Occurred())
                    {
                        failure.reason = WrongType;
                        failure.detail_obj = arg;
                        Py_INCREF(arg);
                    }
                    else
                    {
                        *p = v;
                    }
                }

                break;
            }
        }

        if ((failure.reason == Ok || failure.reason == Overflow) && ch == 'W')
        {
            /* An ellipsis matches everything and ends the parse. */
            break;
        }
    }

    /* Handle parse failures appropriately. */

    if (failure.reason == Ok)
        return TRUE;

    if (failure.reason == Overflow)
    {
        /*
         * We have successfully parsed the signature but one of the arguments
         * has been found to overflow.  Raise an appropriate exception and
         * ensure we don't parse any subsequent overloads.
         */
        if (failure.overflow_arg_nr >= 0)
        {
            PyErr_Format(PyExc_OverflowError, "argument %d overflowed: %S",
                    failure.overflow_arg_nr, failure.detail_obj);
        }
        else
        {
            PyErr_Format(PyExc_OverflowError, "argument '%s' overflowed: %S",
                    failure.overflow_arg_name, failure.detail_obj);
        }

        /* The overflow exception has now been raised. */
        failure.reason = Raised;
    }

    if (failure.reason != Raised)
        add_failure(parseErrp, &failure);

    if (failure.reason == Raised)
    {
        Py_XDECREF(failure.detail_obj);

        /*
         * Discard any previous errors and flag that the exception we want the
         * user to see has been raised.
         */
        Py_XDECREF(*parseErrp);
        *parseErrp = Py_None;
        Py_INCREF(Py_None);
    }

    return FALSE;
}


/*
 * Called after a failed conversion of an integer.
 */
static void handle_failed_int_conversion(sipParseFailure *pf, PyObject *arg)
{
    PyObject *xtype, *xvalue, *xtb;

    assert(pf->reason == Ok || pf->reason == Overflow);

    PyErr_Fetch(&xtype, &xvalue, &xtb);

    if (PyErr_GivenExceptionMatches(xtype, PyExc_OverflowError) && xvalue != NULL)
    {
        /* Remove any previous overflow exception. */
        Py_XDECREF(pf->detail_obj);

        pf->reason = Overflow;
        pf->overflow_arg_nr = pf->arg_nr;
        pf->overflow_arg_name = pf->arg_name;
        pf->detail_obj = xvalue;
        Py_INCREF(xvalue);
    }
    else
    {
        pf->reason = WrongType;
        pf->detail_obj = arg;
        Py_INCREF(arg);
    }

    PyErr_Restore(xtype, xvalue, xtb);
}


/*
 * Second pass of the argument parse, converting the remaining ones that might
 * have side effects.  Return TRUE if there was no error.
 */
static int parsePass2(sipSimpleWrapper *self, int selfarg, PyObject *sipArgs,
        PyObject *sipKwdArgs, const char **kwdlist, const char *fmt,
        va_list va)
{
    int a, ok;
    Py_ssize_t nr_pos_args;

    /* Handle the converions of "self" first. */
    switch (*fmt++)
    {
    case '#':
        va_arg(va, PyObject *);
        break;

    case 'B':
        {
            /*
             * The address of a C++ instance when calling one of its public
             * methods.
             */

            const sipTypeDef *td;
            void **p;

            *va_arg(va, PyObject **) = (PyObject *)self;
            td = va_arg(va, const sipTypeDef *);
            p = va_arg(va, void **);

            if ((*p = sip_api_get_cpp_ptr(self, td)) == NULL)
                return FALSE;

            break;
        }

    case 'p':
        {
            /*
             * The address of a C++ instance when calling one of its protected
             * methods.
             */

            const sipTypeDef *td;
            void **p;

            *va_arg(va, PyObject **) = (PyObject *)self;
            td = va_arg(va, const sipTypeDef *);
            p = va_arg(va, void **);

            if ((*p = getComplexCppPtr(self, td)) == NULL)
                return FALSE;

            break;
        }

    case 'C':
        va_arg(va, PyObject *);
        break;

    default:
        --fmt;
    }

    ok = TRUE;
    nr_pos_args = PyTuple_GET_SIZE(sipArgs);

    for (a = (selfarg ? 1 : 0); *fmt != '\0' && *fmt != 'W' && ok; ++a)
    {
        char ch;
        PyObject *arg;

        /* Skip the optional character. */
        if ((ch = *fmt++) == '|')
            ch = *fmt++;

        /* Get the next argument. */
        arg = NULL;

        if (a < nr_pos_args)
        {
            arg = PyTuple_GET_ITEM(sipArgs, a);
        }
        else if (sipKwdArgs != NULL)
        {
            const char *name = kwdlist[a - selfarg];

            if (name != NULL)
                arg = PyDict_GetItemString(sipKwdArgs, name);
        }

        /*
         * Do the outstanding conversions.  For most types it has already been
         * done, so we are just skipping the parameters.
         */
        switch (ch)
        {
        case '@':
            /* Implement /GetWrapper/. */
            va_arg(va, PyObject **);

            /* Process the same argument next time round. */
            --a;

            break;

        case 'q':
            {
                /* Qt receiver to connect. */

                char *sig = va_arg(va, char *);
                void **rx = va_arg(va, void **);
                const char **slot = va_arg(va, const char **);

                if (arg != NULL)
                {
                    *rx = sip_api_convert_rx((sipWrapper *)self, sig, arg,
                            *slot, slot, 0);

                    if (*rx == NULL)
                        return FALSE;
                }

                break;
            }

        case 'Q':
            {
                /* Qt receiver to disconnect. */

                char *sig = va_arg(va, char *);
                void **rx = va_arg(va, void **);
                const char **slot = va_arg(va, const char **);

                if (arg != NULL)
                    *rx = sipGetRx(self, sig, arg, *slot, slot);

                break;
            }

        case 'g':
            {
                /* Python single shot slot to connect. */

                char *sig = va_arg(va, char *);
                void **rx = va_arg(va, void **);
                const char **slot = va_arg(va, const char **);

                if (arg != NULL)
                {
                    *rx = sip_api_convert_rx((sipWrapper *)self, sig, arg,
                            NULL, slot, SIP_SINGLE_SHOT);

                    if (*rx == NULL)
                        return FALSE;
                }

                break;
            }

        case 'y':
            {
                /* Python slot to connect. */

                char *sig = va_arg(va, char *);
                void **rx = va_arg(va, void **);
                const char **slot = va_arg(va, const char **);

                if (arg != NULL)
                {
                    *rx = sip_api_convert_rx((sipWrapper *)self, sig, arg,
                            NULL, slot, 0);

                    if (*rx == NULL)
                        return FALSE;
                }

                break;
            }

        case 'Y':
            {
                /* Python slot to disconnect. */

                char *sig = va_arg(va, char *);
                void **rx = va_arg(va, void **);
                const char **slot = va_arg(va, const char **);

                if (arg != NULL)
                    *rx = sipGetRx(self, sig, arg, NULL, slot);

                break;
            }

        case 'r':
            {
                /* Sequence of mapped type instances. */

                const sipTypeDef *td;
                void **array;
                Py_ssize_t *nr_elem;

                td = va_arg(va, const sipTypeDef *);
                array = va_arg(va, void **);
                nr_elem = va_arg(va, Py_ssize_t *);

                if (arg != NULL && !convertFromSequence(arg, td, array, nr_elem))
                    return FALSE;

                break;
            }

        case '>':
            {
                /* Sequence or sip.array of class instances. */

                const sipTypeDef *td;
                void **array;
                Py_ssize_t *nr_elem;
                int *is_temp;

                td = va_arg(va, const sipTypeDef *);
                array = va_arg(va, void **);
                nr_elem = va_arg(va, Py_ssize_t *);
                is_temp = va_arg(va, int *);

                if (arg != NULL)
                {
                    if (sip_array_can_convert(arg, td))
                    {
                        sip_array_convert(arg, array, nr_elem);
                        *is_temp = FALSE;
                    }
                    else if (convertFromSequence(arg, td, array, nr_elem))
                    {
                        /*
                         * Note that this will leak if there is a subsequent
                         * error.
                         */
                        *is_temp = TRUE;
                    }
                    else
                    {
                        return FALSE;
                    }
                }

                break;
            }

        case 'J':
            {
                /* Class or mapped type instance. */

                int flags = *fmt++ - '0';
                const sipTypeDef *td;
                void **p;
                int iflgs = 0;
                int *state;
                PyObject *xfer, **owner;

                td = va_arg(va, const sipTypeDef *);
                p = va_arg(va, void **);

                if (flags & FMT_AP_TRANSFER)
                    xfer = (self ? (PyObject *)self : arg);
                else if (flags & FMT_AP_TRANSFER_BACK)
                    xfer = Py_None;
                else
                    xfer = NULL;

                if (flags & FMT_AP_DEREF)
                    iflgs |= SIP_NOT_NONE;

                if (flags & FMT_AP_TRANSFER_THIS)
                    owner = va_arg(va, PyObject **);

                if (flags & FMT_AP_NO_CONVERTORS)
                {
                    iflgs |= SIP_NO_CONVERTORS;
                    state = NULL;
                }
                else
                {
                    state = va_arg(va, int *);
                }

                if (arg != NULL)
                {
                    int iserr = FALSE;

                    *p = sip_api_convert_to_type(arg, td, xfer, iflgs, state,
                            &iserr);

                    if (iserr)
                        return FALSE;

                    if (flags & FMT_AP_TRANSFER_THIS && *p != NULL)
                        *owner = arg;
                }

                break;
            }

        case 'P':
            {
                /* Python object of any type with a sub-format. */

                PyObject **p = va_arg(va, PyObject **);
                int flags = *fmt++ - '0';

                if (arg != NULL)
                {
                    if (flags & FMT_AP_TRANSFER)
                    {
                        Py_XINCREF(arg);
                    }
                    else if (flags & FMT_AP_TRANSFER_BACK)
                    {
                        Py_XDECREF(arg);
                    }

                    *p = arg;
                }

                break;
            }

        case 'X':
            {
                /* Constrained types. */

                if (*fmt++ == 'E')
                    va_arg(va, void *);

                va_arg(va, void *);

                break;
            }

        case 'A':
            {
                /* String from a Python string or None. */

                PyObject **keep = va_arg(va, PyObject **);
                const char **p = va_arg(va, const char **);
                char sub_fmt = *fmt++;

                if (arg != NULL)
                {
                    PyObject *s = NULL;

                    switch (sub_fmt)
                    {
                    case 'A':
                        s = parseString_AsASCIIString(arg, p);
                        break;

                    case 'L':
                        s = parseString_AsLatin1String(arg, p);
                        break;

                    case '8':
                        s = parseString_AsUTF8String(arg, p);
                        break;
                    }

                    if (s == NULL)
                        return FALSE;

                    *keep = s;
                }

                break;
            }

        case 'a':
            {
                /* Character from a Python string. */

                char *p = va_arg(va, char *);
                char sub_fmt = *fmt++;

                if (arg != NULL)
                {
                    int enc;

                    switch (sub_fmt)
                    {
                    case 'A':
                        enc = parseString_AsASCIIChar(arg, p);
                        break;

                    case 'L':
                        enc = parseString_AsLatin1Char(arg, p);
                        break;

                    case '8':
                        enc = parseString_AsUTF8Char(arg, p);
                        break;
                    }

                    if (enc < 0)
                        return FALSE;
                }

                break;
            }

        /*
         * Every other argument is a pointer and only differ in how many there
         * are.
         */
        case 'N':
        case 'T':
        case 'k':
        case 'K':
        case 'U':
        case 'E':
            va_arg(va, void *);

            /* Drop through. */

        default:
            va_arg(va, void *);
        }
    }

    /* Handle any ellipsis argument. */
    if (*fmt == 'W')
    {
        PyObject *al;
        int da = 0;

        /* Create a tuple for any remaining arguments. */
        if ((al = PyTuple_New(nr_pos_args - a)) == NULL)
            return FALSE;

        while (a < nr_pos_args)
        {
            PyObject *arg = PyTuple_GET_ITEM(sipArgs, a);

            /* Add the remaining argument to the tuple. */
            Py_INCREF(arg);
            PyTuple_SET_ITEM(al, da, arg);

            ++a;
            ++da;
        }

        /* Return the tuple. */
        *va_arg(va, PyObject **) = al;
    }

    return TRUE;
}


/*
 * Return TRUE if an object is a QObject.
 */
static int isQObject(PyObject *obj)
{
    return (sipQtSupport != NULL && PyObject_TypeCheck(obj, sipTypeAsPyTypeObject(sipQObjectType)));
}


/*
 * See if a Python object is a sequence of a particular type.
 */
static int canConvertFromSequence(PyObject *seq, const sipTypeDef *td)
{
    Py_ssize_t i, size = PySequence_Size(seq);

    if (size < 0)
        return FALSE;

    /*
     * Check the type of each element.  Note that this is inconsistent with how
     * similiar situations are handled elsewhere.  We should instead just check
     * we have an iterator and assume (until the second pass) that the type is
     * correct.
     */
    for (i = 0; i < size; ++i)
    {
        int ok;
        PyObject *val_obj;

        if ((val_obj = PySequence_GetItem(seq, i)) == NULL)
            return FALSE;

        ok = sip_api_can_convert_to_type(val_obj, td,
                SIP_NO_CONVERTORS|SIP_NOT_NONE);

        Py_DECREF(val_obj);

        if (!ok)
            return FALSE;
    }

    return TRUE;
}


/*
 * Convert a Python sequence to an array that has already "passed"
 * canConvertFromSequence().  Return TRUE if the conversion was successful.
 */
static int convertFromSequence(PyObject *seq, const sipTypeDef *td,
        void **array, Py_ssize_t *nr_elem)
{
    int iserr = 0;
    Py_ssize_t i, size = PySequence_Size(seq);
    sipArrayFunc array_helper;
    sipAssignFunc assign_helper;
    void *array_mem;

    /* Get the type's helpers. */
    if (sipTypeIsMapped(td))
    {
        array_helper = ((const sipMappedTypeDef *)td)->mtd_array;
        assign_helper = ((const sipMappedTypeDef *)td)->mtd_assign;
    }
    else
    {
        array_helper = ((const sipClassTypeDef *)td)->ctd_array;
        assign_helper = ((const sipClassTypeDef *)td)->ctd_assign;
    }

    assert(array_helper != NULL);
    assert(assign_helper != NULL);

    /*
     * Create the memory for the array of values.  Note that this will leak if
     * there is an error.
     */
    array_mem = array_helper(size);

    for (i = 0; i < size; ++i)
    {
        PyObject *val_obj;
        void *val;

        if ((val_obj = PySequence_GetItem(seq, i)) == NULL)
            return FALSE;

        val = sip_api_convert_to_type(val_obj, td, NULL,
                SIP_NO_CONVERTORS|SIP_NOT_NONE, NULL, &iserr);

        Py_DECREF(val_obj);

        if (iserr)
            return FALSE;

        assign_helper(array_mem, i, val);
    }

    *array = array_mem;
    *nr_elem = size;

    return TRUE;
}


/*
 * Convert an array of a type to a Python sequence.
 */
static PyObject *convertToSequence(void *array, Py_ssize_t nr_elem,
        const sipTypeDef *td)
{
    Py_ssize_t i;
    PyObject *seq;
    sipCopyFunc copy_helper;

    /* Get the type's copy helper. */
    if (sipTypeIsMapped(td))
        copy_helper = ((const sipMappedTypeDef *)td)->mtd_copy;
    else
        copy_helper = ((const sipClassTypeDef *)td)->ctd_copy;

    assert(copy_helper != NULL);

    if ((seq = PyTuple_New(nr_elem)) == NULL)
        return NULL;

    for (i = 0; i < nr_elem; ++i)
    {
        void *el = copy_helper(array, i);
        PyObject *el_obj = sip_api_convert_from_new_type(el, td, NULL);

        if (el_obj == NULL)
        {
            release(el, td, 0);
            Py_DECREF(seq);
        }

        PyTuple_SET_ITEM(seq, i, el_obj);
    }

    return seq;
}


/*
 * Perform housekeeping after a C++ instance has been destroyed.
 */
void sip_api_instance_destroyed(sipSimpleWrapper *sw)
{
    sip_api_instance_destroyed_ex(&sw);
}


/*
 * Carry out actions common to all dtors.
 */
static void sip_api_instance_destroyed_ex(sipSimpleWrapper **sipSelfp)
{
    /* If there is no interpreter just to the minimum and get out. */
    if (sipInterpreter == NULL)
    {
        *sipSelfp = NULL;
        return;
    }

    SIP_BLOCK_THREADS

    sipSimpleWrapper *sipSelf = *sipSelfp;

    if (sipSelf != NULL)
    {
        PyObject *xtype, *xvalue, *xtb;

        /* We may be tidying up after an exception so preserve it. */
        PyErr_Fetch(&xtype, &xvalue, &xtb);
        callPyDtor(sipSelf);
        PyErr_Restore(xtype, xvalue, xtb);

        sipOMRemoveObject(&cppPyMap, sipSelf);

        /*
         * This no longer points to anything useful.  Actually it might do as
         * the partialy destroyed C++ instance may still be trying to invoke
         * reimplemented virtuals.
         */
        clear_access_func(sipSelf);

        /*
         * If C/C++ has a reference (and therefore no parent) then remove it.
         * Otherwise remove the object from any parent.
         */
        if (sipCppHasRef(sipSelf))
        {
            sipResetCppHasRef(sipSelf);
            Py_DECREF(sipSelf);
        }
        else if (PyObject_TypeCheck((PyObject *)sipSelf, (PyTypeObject *)&sipWrapper_Type))
        {
            removeFromParent((sipWrapper *)sipSelf);
        }

        /*
         * Normally this is done in the generated dealloc function.  However
         * this is only called if the pointer/access function has not been
         * reset (which it has).  It acts as a guard to prevent any further
         * invocations of reimplemented virtuals.
         */
        *sipSelfp = NULL;
    }

    SIP_UNBLOCK_THREADS
}


/*
 * Clear any access function so that sip_api_get_address() will always return a
 * NULL pointer.
 */
static void clear_access_func(sipSimpleWrapper *sw)
{
    if (sw->access_func != NULL)
    {
        sw->access_func(sw, ReleaseGuard);
        sw->access_func = NULL;
    }

    sw->data = NULL;
}


/*
 * Call self.__dtor__() if it is implemented.
 */
static void callPyDtor(sipSimpleWrapper *self)
{
    sip_gilstate_t sipGILState;
    char pymc = 0;
    PyObject *meth;

    meth = sip_api_is_py_method_12_8(&sipGILState, &pymc, &self, NULL,
            "__dtor__");

    if (meth != NULL)
    {
        PyObject *res = sip_api_call_method(0, meth, "", NULL);

        Py_DECREF(meth);

        /* Discard any result. */
        Py_XDECREF(res);

        /* Handle any error the best we can. */
        if (PyErr_Occurred())
            PyErr_Print();

        SIP_RELEASE_GIL(sipGILState);
    }
}


/*
 * Add a wrapper to it's parent owner.  The wrapper must not currently have a
 * parent and, therefore, no siblings.
 */
static void addToParent(sipWrapper *self, sipWrapper *owner)
{
    if (owner->first_child != NULL)
    {
        self->sibling_next = owner->first_child;
        owner->first_child->sibling_prev = self;
    }

    owner->first_child = self;
    self->parent = owner;

    /*
     * The owner holds a real reference so that the cyclic garbage collector
     * works properly.
     */
    Py_INCREF((sipSimpleWrapper *)self);
}


/*
 * Remove a wrapper from it's parent if it has one.
 */
static void removeFromParent(sipWrapper *self)
{
    if (self->parent != NULL)
    {
        if (self->parent->first_child == self)
            self->parent->first_child = self->sibling_next;

        if (self->sibling_next != NULL)
            self->sibling_next->sibling_prev = self->sibling_prev;

        if (self->sibling_prev != NULL)
            self->sibling_prev->sibling_next = self->sibling_next;

        self->parent = NULL;
        self->sibling_next = NULL;
        self->sibling_prev = NULL;

        /*
         * We must do this last, after all the pointers are correct, because
         * this is used by the clear slot.
         */
        Py_DECREF((sipSimpleWrapper *)self);
    }
}


/*
 * Detach and children of a parent.
 */
static void detachChildren(sipWrapper *self)
{
    while (self->first_child != NULL)
        removeFromParent(self->first_child);
}


/*
 * Convert a sequence index.  Return the index or a negative value if there was
 * an error.
 */
static Py_ssize_t sip_api_convert_from_sequence_index(Py_ssize_t idx,
        Py_ssize_t len)
{
    /* Negative indices start from the other end. */
    if (idx < 0)
        idx = len + idx;

    if (idx < 0 || idx >= len)
    {
        PyErr_Format(PyExc_IndexError, "sequence index out of range");
        return -1;
    }

    return idx;
}


/*
 * Return a tuple of the base class of a type that has no explicit super-type.
 */
static PyObject *getDefaultBase(void)
{
    static PyObject *default_base = NULL;

    /* Only do this once. */
    if (default_base == NULL)
    {
        if ((default_base = PyTuple_Pack(1, (PyObject *)&sipWrapper_Type)) == NULL)
            return NULL;
    }

    Py_INCREF(default_base);

    return default_base;
}


/*
 * Return a tuple of the base class of a simple type that has no explicit
 * super-type.
 */
static PyObject *getDefaultSimpleBase(void)
{
    static PyObject *default_simple_base = NULL;

    /* Only do this once. */
    if (default_simple_base == NULL)
    {
        if ((default_simple_base = PyTuple_Pack(1, (PyObject *)&sipSimpleWrapper_Type)) == NULL)
            return NULL;
    }

    Py_INCREF(default_simple_base);

    return default_simple_base;
}


/*
 * Return the dictionary of a type.
 */
static PyObject *getScopeDict(sipTypeDef *td, PyObject *mod_dict,
        sipExportedModuleDef *client)
{
    /*
     * Initialise the scoping type if necessary.  It will always be in the
     * same module if it needs doing.
     */
    if (sipTypeIsMapped(td))
    {
        if (createMappedType(client, (sipMappedTypeDef *)td, mod_dict) < 0)
            return NULL;

        /* Check that the mapped type can act as a container. */
        assert(sipTypeAsPyTypeObject(td) != NULL);
    }
    else
    {
        if (createClassType(client, (sipClassTypeDef *)td, mod_dict) < 0)
            return NULL;
    }

    return (sipTypeAsPyTypeObject(td))->tp_dict;
}


/*
 * Create a container type and return a borrowed reference to it.
 */
static PyObject *createContainerType(sipContainerDef *cod, sipTypeDef *td,
        PyObject *bases, PyObject *metatype, PyObject *mod_dict,
        PyObject *type_dict, sipExportedModuleDef *client)
{
    PyObject *py_type, *scope_dict, *name, *args;
    sipTypeDef *scope_td;

    /* Get the dictionary to place the type in. */
    if (cod->cod_scope.sc_flag)
    {
        scope_td = NULL;
        scope_dict = mod_dict;
    }
    else
    {
        scope_td = getGeneratedType(&cod->cod_scope, client);
        scope_dict = getScopeDict(scope_td, mod_dict, client);

        if (scope_dict == NULL)
            goto reterr;
    }

    /* Create an object corresponding to the type name. */
    if ((name = PyUnicode_FromString(sipPyNameOfContainer(cod, td))) == NULL)
        goto reterr;

    /* Create the type by calling the metatype. */
    if ((args = PyTuple_Pack(3, name, bases, type_dict)) == NULL)
        goto relname;

    /* Pass the type via the back door. */
    assert(currentType == NULL);
    currentType = td;
    py_type = PyObject_Call(metatype, args, NULL);
    currentType = NULL;

    if (py_type == NULL)
        goto relargs;

    /* Fix __qualname__ if there is a scope. */
    if (scope_td != NULL)
    {
        PyHeapTypeObject *ht;
        PyObject *qualname = get_qualname(scope_td, name);

        if (qualname == NULL)
            goto reltype;

        ht = (PyHeapTypeObject *)py_type;

        Py_CLEAR(ht->ht_qualname);
        ht->ht_qualname = qualname;
    }

    /* Add the type to the "parent" dictionary. */
    if (PyDict_SetItem(scope_dict, name, py_type) < 0)
        goto reltype;

    Py_DECREF(args);
    Py_DECREF(name);

    return py_type;

    /* Unwind on error. */

reltype:
    Py_DECREF(py_type);

relargs:
    Py_DECREF(args);

relname:
    Py_DECREF(name);

reterr:
    return NULL;
}


/*
 * Create a single class type object.
 */
static int createClassType(sipExportedModuleDef *client, sipClassTypeDef *ctd,
        PyObject *mod_dict)
{
    PyObject *bases, *metatype, *py_type, *type_dict;
    sipEncodedTypeDef *sup;
    int i;

    /* Handle the trivial case where we have already been initialised. */
    if (ctd->ctd_base.td_module != NULL)
        return 0;

    /* Set this up now to gain access to the string pool. */
    ctd->ctd_base.td_module = client;

    /* Create the tuple of super-types. */
    if ((sup = ctd->ctd_supers) == NULL)
    {
        if (ctd->ctd_supertype < 0)
        {
            bases = (sipTypeIsNamespace(&ctd->ctd_base) ? getDefaultSimpleBase() : getDefaultBase());
        }
        else
        {
            PyObject *supertype;
            const char *supertype_name = sipNameFromPool(client,
                    ctd->ctd_supertype);

            if ((supertype = findPyType(supertype_name)) == NULL)
                goto reterr;

            bases = PyTuple_Pack(1, supertype);
        }

        if (bases == NULL)
            goto reterr;
    }
    else
    {
        int nrsupers = 0;

        do
            ++nrsupers;
        while (!sup++->sc_flag);

        if ((bases = PyTuple_New(nrsupers)) == NULL)
            goto reterr;

        for (sup = ctd->ctd_supers, i = 0; i < nrsupers; ++i, ++sup)
        {
            PyObject *st;
            sipTypeDef *sup_td = getGeneratedType(sup, client);

            /*
             * Initialise the super-class if necessary.  It will always be in
             * the same module if it needs doing.
             */
            if (createClassType(client, (sipClassTypeDef *)sup_td, mod_dict) < 0)
                goto relbases;

            st = (PyObject *)sipTypeAsPyTypeObject(sup_td);

            Py_INCREF(st);
            PyTuple_SET_ITEM(bases, i, st);

            /*
             * Inherit any garbage collector code rather than look for it each
             * time it is needed.
             */
            if (ctd->ctd_traverse == NULL)
                ctd->ctd_traverse = ((sipClassTypeDef *)sup_td)->ctd_traverse;

            if (ctd->ctd_clear == NULL)
                ctd->ctd_clear = ((sipClassTypeDef *)sup_td)->ctd_clear;
        }
    }

    /*
     * Use the explicit meta-type if there is one, otherwise use the meta-type
     * of the first super-type.
     */
    if (ctd->ctd_metatype >= 0)
    {
        const char *metatype_name = sipNameFromPool(client, ctd->ctd_metatype);

        if ((metatype = findPyType(metatype_name)) == NULL)
            goto relbases;
    }
    else
        metatype = (PyObject *)Py_TYPE(PyTuple_GET_ITEM(bases, 0));

    /* Create the type dictionary and populate it with any non-lazy methods. */
    if ((type_dict = createTypeDict(client)) == NULL)
        goto relbases;

    if (sipTypeHasNonlazyMethod(&ctd->ctd_base))
    {
        PyMethodDef *pmd = ctd->ctd_container.cod_methods;

        for (i = 0; i < ctd->ctd_container.cod_nrmethods; ++i)
        {
            if (isNonlazyMethod(pmd) && addMethod(type_dict, pmd) < 0)
                goto reldict;

            ++pmd;
        }
    }

    if ((py_type = createContainerType(&ctd->ctd_container, (sipTypeDef *)ctd, bases, metatype, mod_dict, type_dict, client)) == NULL)
        goto reldict;

    if (ctd->ctd_pyslots != NULL)
        fix_slots((PyTypeObject *)py_type, ctd->ctd_pyslots);

    /* Handle the pickle function. */
    if (ctd->ctd_pickle != NULL)
    {
        static PyMethodDef md = {
            "_pickle_type", pickle_type, METH_NOARGS, NULL
        };

        if (setReduce((PyTypeObject *)py_type, &md) < 0)
            goto reltype;
    }

    /* We can now release our references. */
    Py_DECREF(bases);
    Py_DECREF(type_dict);

    return 0;

    /* Unwind after an error. */

reltype:
    Py_DECREF(py_type);

reldict:
    Py_DECREF(type_dict);

relbases:
    Py_DECREF(bases);

reterr:
    ctd->ctd_base.td_module = NULL;
    return -1;
}


/*
 * Create a single mapped type object.
 */
static int createMappedType(sipExportedModuleDef *client,
        sipMappedTypeDef *mtd, PyObject *mod_dict)
{
    PyObject *bases, *type_dict;

    /* Handle the trivial case where we have already been initialised. */
    if (mtd->mtd_base.td_module != NULL)
        return 0;

    /* Set this up now to gain access to the string pool. */
    mtd->mtd_base.td_module = client;

    /* Create the tuple of super-types. */
    if ((bases = getDefaultBase()) == NULL)
        goto reterr;

    /* Create the type dictionary. */
    if ((type_dict = createTypeDict(client)) == NULL)
        goto relbases;

    if (createContainerType(&mtd->mtd_container, (sipTypeDef *)mtd, bases, (PyObject *)&sipWrapperType_Type, mod_dict, type_dict, client) == NULL)
        goto reldict;

    /* We can now release our references. */
    Py_DECREF(bases);
    Py_DECREF(type_dict);

    return 0;

    /* Unwind after an error. */

reldict:
    Py_DECREF(type_dict);

relbases:
    Py_DECREF(bases);

reterr:
    mtd->mtd_base.td_module = NULL;
    return -1;
}


/*
 * Return the module definition for a named module.
 */
static sipExportedModuleDef *getModule(PyObject *mname_obj)
{
    PyObject *mod;
    sipExportedModuleDef *em;

    /* Make sure the module is imported. */
    if ((mod = PyImport_Import(mname_obj)) == NULL)
        return NULL;

    /* Find the module definition. */
    for (em = moduleList; em != NULL; em = em->em_next)
        if (PyUnicode_Compare(mname_obj, em->em_nameobj) == 0)
            break;

    Py_DECREF(mod);

    if (em == NULL)
        PyErr_Format(PyExc_SystemError, "unable to find to find module: %U",
                mname_obj);

    return em;
}


/*
 * The type unpickler.
 */
static PyObject *unpickle_type(PyObject *obj, PyObject *args)
{
    PyObject *mname_obj, *init_args;
    const char *tname;
    sipExportedModuleDef *em;
    int i;

    (void)obj;

    if (!PyArg_ParseTuple(args, "UsO!:_unpickle_type", &mname_obj, &tname, &PyTuple_Type, &init_args))
        return NULL;

    /* Get the module definition. */
    if ((em = getModule(mname_obj)) == NULL)
        return NULL;

    /* Find the class type object. */
    for (i = 0; i < em->em_nrtypes; ++i)
    {
        sipTypeDef *td = em->em_types[i];

        if (td != NULL && !sipTypeIsStub(td) && sipTypeIsClass(td))
        {
            const char *pyname = sipPyNameOfContainer(
                    &((sipClassTypeDef *)td)->ctd_container, td);

            if (strcmp(pyname, tname) == 0)
                return PyObject_CallObject((PyObject *)sipTypeAsPyTypeObject(td), init_args);
        }
    }

    PyErr_Format(PyExc_SystemError, "unable to find to find type: %s", tname);

    return NULL;
}


/*
 * The type pickler.
 */
static PyObject *pickle_type(PyObject *obj, PyObject *args)
{
    sipExportedModuleDef *em;

    (void)args;

    /* Find the type definition and defining module. */
    for (em = moduleList; em != NULL; em = em->em_next)
    {
        int i;

        for (i = 0; i < em->em_nrtypes; ++i)
        {
            sipTypeDef *td = em->em_types[i];

            if (td != NULL && !sipTypeIsStub(td) && sipTypeIsClass(td))
                if (sipTypeAsPyTypeObject(td) == Py_TYPE(obj))
                {
                    PyObject *init_args;
                    sipClassTypeDef *ctd = (sipClassTypeDef *)td;
                    const char *pyname = sipPyNameOfContainer(&ctd->ctd_container, td);

                    /*
                     * Ask the handwritten pickle code for the tuple of
                     * arguments that will recreate the object.
                     */
                    init_args = ctd->ctd_pickle(sip_api_get_cpp_ptr((sipSimpleWrapper *)obj, NULL));

                    if (init_args == NULL)
                        return NULL;

                    if (!PyTuple_Check(init_args))
                    {
                        PyErr_Format(PyExc_TypeError,
                                "%%PickleCode for type %s.%s did not return a tuple",
                                sipNameOfModule(em), pyname);

                        return NULL;
                    }

                    return Py_BuildValue("O(OsN)", type_unpickler,
                            em->em_nameobj, pyname, init_args);
                }
        }
    }

    /* We should never get here. */
    PyErr_Format(PyExc_SystemError, "attempt to pickle unknown type '%s'",
            Py_TYPE(obj)->tp_name);

    return NULL;
}


/*
 * The enum unpickler.
 */
static PyObject *unpickle_enum(PyObject *obj, PyObject *args)
{
    PyObject *mname_obj, *evalue_obj;
    const char *ename;
    sipExportedModuleDef *em;
    int i;

    (void)obj;

    if (!PyArg_ParseTuple(args, "UsO:_unpickle_enum", &mname_obj, &ename, &evalue_obj))
        return NULL;

    /* Get the module definition. */
    if ((em = getModule(mname_obj)) == NULL)
        return NULL;

    /* Find the enum type object. */
    for (i = 0; i < em->em_nrtypes; ++i)
    {
        sipTypeDef *td = em->em_types[i];

        if (td != NULL && !sipTypeIsStub(td) && sipTypeIsEnum(td))
            if (strcmp(sipPyNameOfEnum((sipEnumTypeDef *)td), ename) == 0)
                return PyObject_CallFunctionObjArgs((PyObject *)sipTypeAsPyTypeObject(td), evalue_obj, NULL);
    }

    PyErr_Format(PyExc_SystemError, "unable to find to find enum: %s", ename);

    return NULL;
}


/*
 * The enum pickler.
 */
static PyObject *pickle_enum(PyObject *obj, PyObject *args)
{
    sipTypeDef *td = ((sipEnumTypeObject *)Py_TYPE(obj))->type;

    (void)args;

    return Py_BuildValue("O(Osi)", enum_unpickler, td->td_module->em_nameobj,
            sipPyNameOfEnum((sipEnumTypeDef *)td), (int)PyLong_AS_LONG(obj));
}


/*
 * Set the __reduce__method for a type.
 */
static int setReduce(PyTypeObject *type, PyMethodDef *pickler)
{
    static PyObject *rstr = NULL;
    PyObject *descr;
    int rc;

    if (objectify("__reduce__", &rstr) < 0)
        return -1;

    /* Create the method descripter. */
    if ((descr = PyDescr_NewMethod(type, pickler)) == NULL)
        return -1;

    /*
     * Save the method.  Note that we don't use PyObject_SetAttr() as we want
     * to bypass any lazy attribute loading (which may not be safe yet).
     */
    rc = PyType_Type.tp_setattro((PyObject *)type, rstr, descr);

    Py_DECREF(descr);

    return rc;
}


/*
 * Create an enum object.
 */
static int createEnum(sipExportedModuleDef *client, sipEnumTypeDef *etd,
        int enum_nr, PyObject *mod_dict)
{
    int rc;
    PyObject *name, *dict, *enum_obj;

    etd->etd_base.td_module = client;

    /* Get the dictionary into which the type will be placed. */
    if (etd->etd_scope < 0)
        dict = mod_dict;
    else if ((dict = getScopeDict(client->em_types[etd->etd_scope], mod_dict, client)) == NULL)
        return -1;

    /* Create an object corresponding to the type name. */
    if ((name = PyUnicode_FromString(sipPyNameOfEnum(etd))) == NULL)
        return -1;

    /* Create the enum. */
    if (sipTypeIsEnum(&etd->etd_base))
        enum_obj = createUnscopedEnum(client, etd, name);
    else
        enum_obj = createScopedEnum(client, etd, enum_nr, name);

    if (enum_obj == NULL)
    {
        Py_DECREF(name);
        return -1;
    }

    /* Add the enum to the "parent" dictionary. */
    rc = PyDict_SetItem(dict, name, enum_obj);

    /* We can now release our remaining references. */
    Py_DECREF(name);
    Py_DECREF(enum_obj);

    return rc;
}


/*
 * Create an unscoped enum.
 */
static PyObject *createUnscopedEnum(sipExportedModuleDef *client,
        sipEnumTypeDef *etd, PyObject *name)
{
    static PyObject *bases = NULL;
    PyObject *type_dict, *args;
    sipEnumTypeObject *eto;

    /* Create the base type tuple if it hasn't already been done. */
    if (bases == NULL)
        if ((bases = PyTuple_Pack(1, (PyObject *)&PyLong_Type)) == NULL)
            return NULL;

    /* Create the type dictionary. */
    if ((type_dict = createTypeDict(client)) == NULL)
        return NULL;

    /* Create the type by calling the metatype. */
    args = PyTuple_Pack(3, name, bases, type_dict);

    Py_DECREF(type_dict);

    if (args == NULL)
        return NULL;

    /* Pass the type via the back door. */
    assert(currentType == NULL);
    currentType = &etd->etd_base;
    eto = (sipEnumTypeObject *)PyObject_Call((PyObject *)&sipEnumType_Type,
            args, NULL);
    currentType = NULL;

    Py_DECREF(args);

    if (eto == NULL)
        return NULL;

    if (etd->etd_pyslots != NULL)
        fix_slots((PyTypeObject *)eto, etd->etd_pyslots);

    /*
     * If the enum has a scope then the default __qualname__ will be incorrect.
     */
     if (etd->etd_scope >= 0)
     {
        /* Append the name of the enum to the scope's __qualname__. */
        Py_CLEAR(eto->super.ht_qualname);
        eto->super.ht_qualname = get_qualname(
                client->em_types[etd->etd_scope], name);

        if (eto->super.ht_qualname == NULL)
        {
            Py_DECREF((PyObject *)eto);
            return NULL;
        }
    }

    return (PyObject *)eto;
}


/*
 * Create a scoped enum.
 */
static PyObject *createScopedEnum(sipExportedModuleDef *client,
        sipEnumTypeDef *etd, int enum_nr, PyObject *name)
{
    static PyObject *enum_type = NULL, *module_arg = NULL;
    static PyObject *qualname_arg = NULL;
    int i, nr_members;
    sipEnumMemberDef *enm;
    PyObject *members, *enum_obj, *args, *kw_args;

    /* Get the enum type if we haven't done so already. */
    if (enum_type == NULL)
    {
        if ((enum_type = import_module_attr("enum", "IntEnum")) == NULL)
            goto ret_err;
    }

    /* Create a dict of the members. */
    if ((members = PyDict_New()) == NULL)
        goto ret_err;

    /*
     * Note that the current structures for defining scoped enums are not ideal
     * as we are re-using the ones used for unscoped enums (which are designed
     * to support lazy implementations).
     */
    if (etd->etd_scope < 0)
    {
        nr_members = client->em_nrenummembers;
        enm = client->em_enummembers;
    }
    else
    {
        const sipContainerDef *cod = get_container(client->em_types[etd->etd_scope]);

        nr_members = cod->cod_nrenummembers;
        enm = cod->cod_enummembers;
    }

    for (i = 0; i < nr_members; ++i)
    {
        if (enm->em_enum == enum_nr)
        {
            PyObject *val = PyLong_FromLong(enm->em_val);

            if (dict_set_and_discard(members, enm->em_name, val) < 0)
                goto rel_members;
        }

        ++enm;
    }

    if ((args = PyTuple_Pack(2, name, members)) == NULL)
        goto rel_members;

    if ((kw_args = PyDict_New()) == NULL)
        goto rel_args;

    if (objectify("module", &module_arg) < 0)
        goto rel_kw_args;

    if (PyDict_SetItem(kw_args, module_arg, client->em_nameobj) < 0)
        goto rel_kw_args;

    /*
     * If the enum has a scope then the default __qualname__ will be incorrect.
     */
     if (etd->etd_scope >= 0)
     {
        int rc;
        PyObject *qualname;

        if (objectify("qualname", &qualname_arg) < 0)
            goto rel_kw_args;

        if ((qualname = get_qualname(client->em_types[etd->etd_scope], name)) == NULL)
            goto rel_kw_args;

        rc = PyDict_SetItem(kw_args, qualname_arg, qualname);

        Py_DECREF(qualname);

        if (rc < 0)
            goto rel_kw_args;
    }

    if ((enum_obj = PyObject_Call(enum_type, args, kw_args)) == NULL)
        goto rel_kw_args;

    Py_DECREF(kw_args);
    Py_DECREF(args);
    Py_DECREF(members);

    /* Note that it isn't actually a PyTypeObject. */
    etd->etd_base.td_py_type = (PyTypeObject *)enum_obj;

    return enum_obj;

    /* Unwind on errors. */

rel_kw_args:
    Py_DECREF(kw_args);

rel_args:
    Py_DECREF(args);

rel_members:
    Py_DECREF(members);

ret_err:
    return NULL;
}


/*
 * Create a type dictionary for dynamic type being created in a module.
 */
static PyObject *createTypeDict(sipExportedModuleDef *em)
{
    static PyObject *mstr = NULL;
    PyObject *dict;

    if (objectify("__module__", &mstr) < 0)
        return NULL;

    /* Create the dictionary. */
    if ((dict = PyDict_New()) == NULL)
        return NULL;

    /* We need to set the module name as an attribute for dynamic types. */
    if (PyDict_SetItem(dict, mstr, em->em_nameobj) < 0)
    {
        Py_DECREF(dict);
        return NULL;
    }

    return dict;
}


/*
 * Convert an ASCII string to a Python object if it hasn't already been done.
 */
static int objectify(const char *s, PyObject **objp)
{
    if (*objp == NULL)
        if ((*objp = PyUnicode_FromString(s)) == NULL)
            return -1;

    return 0;
}


/*
 * Add a set of static instances to a dictionary.
 */
static int addInstances(PyObject *dict, sipInstancesDef *id)
{
    if (id->id_type != NULL && addTypeInstances(dict, id->id_type) < 0)
        return -1;

    if (id->id_voidp != NULL && addVoidPtrInstances(dict,id->id_voidp) < 0)
        return -1;

    if (id->id_char != NULL && addCharInstances(dict,id->id_char) < 0)
        return -1;

    if (id->id_string != NULL && addStringInstances(dict,id->id_string) < 0)
        return -1;

    if (id->id_int != NULL && addIntInstances(dict, id->id_int) < 0)
        return -1;

    if (id->id_long != NULL && addLongInstances(dict,id->id_long) < 0)
        return -1;

    if (id->id_ulong != NULL && addUnsignedLongInstances(dict, id->id_ulong) < 0)
        return -1;

    if (id->id_llong != NULL && addLongLongInstances(dict, id->id_llong) < 0)
        return -1;

    if (id->id_ullong != NULL && addUnsignedLongLongInstances(dict, id->id_ullong) < 0)
        return -1;

    if (id->id_double != NULL && addDoubleInstances(dict,id->id_double) < 0)
        return -1;

    return 0;
}


/*
 * Get "self" from the argument tuple for a method called as
 * Class.Method(self, ...) rather than self.Method(...).
 */
static int getSelfFromArgs(sipTypeDef *td, PyObject *args, int argnr,
        sipSimpleWrapper **selfp)
{
    PyObject *self;

    /* Get self from the argument tuple. */

    if (argnr >= PyTuple_GET_SIZE(args))
        return FALSE;

    self = PyTuple_GET_ITEM(args, argnr);

    if (!PyObject_TypeCheck(self, sipTypeAsPyTypeObject(td)))
        return FALSE;

    *selfp = (sipSimpleWrapper *)self;

    return TRUE;
}


/*
 * Return non-zero if a method is non-lazy, ie. it must be added to the type
 * when it is created.
 */
static int isNonlazyMethod(PyMethodDef *pmd)
{
    static const char *lazy[] = {
        "__getattribute__",
        "__getattr__",
        "__enter__",
        "__exit__",
        "__aenter__",
        "__aexit__",
        NULL
    };  

    const char **l;

    for (l = lazy; *l != NULL; ++l)
        if (strcmp(pmd->ml_name, *l) == 0)
            return TRUE;

    return FALSE;
}


/*
 * Add a method to a dictionary.
 */
static int addMethod(PyObject *dict, PyMethodDef *pmd)
{
    PyObject *descr = sipMethodDescr_New(pmd);

    return dict_set_and_discard(dict, pmd->ml_name, descr);
}


/*
 * Populate a container's type dictionary.
 */
static int add_lazy_container_attrs(sipTypeDef *td, sipContainerDef *cod,
        PyObject *dict)
{
    int i;
    PyMethodDef *pmd;
    sipEnumMemberDef *enm;
    sipVariableDef *vd;

    /* Do the methods. */
    for (pmd = cod->cod_methods, i = 0; i < cod->cod_nrmethods; ++i, ++pmd)
    {
        /* Non-lazy methods will already have been handled. */
        if (!sipTypeHasNonlazyMethod(td) || !isNonlazyMethod(pmd))
        {
            if (addMethod(dict, pmd) < 0)
                return -1;
        }
    }

    /* Do the unscoped enum members. */
    for (enm = cod->cod_enummembers, i = 0; i < cod->cod_nrenummembers; ++i, ++enm)
    {
        PyObject *val;

        if (enm->em_enum < 0)
        {
            /* It's an unnamed unscoped enum. */
            val = PyLong_FromLong(enm->em_val);
        }
        else
        {
            sipTypeDef *etd = td->td_module->em_types[enm->em_enum];

            if (sipTypeIsScopedEnum(etd))
                continue;

            val = sip_api_convert_from_enum(enm->em_val, etd);
        }

        if (dict_set_and_discard(dict, enm->em_name, val) < 0)
            return -1;
    }

    /* Do the variables. */
    for (vd = cod->cod_variables, i = 0; i < cod->cod_nrvariables; ++i, ++vd)
    {
        PyObject *descr;

        if (vd->vd_type == PropertyVariable)
            descr = create_property(vd);
        else
            descr = sipVariableDescr_New(vd, td, cod);

        if (dict_set_and_discard(dict, vd->vd_name, descr) < 0)
            return -1;
    }

    return 0;
}


/*
 * Create a Python property object from the SIP generated structure.
 */
static PyObject *create_property(sipVariableDef *vd)
{
    PyObject *descr, *fget, *fset, *fdel, *doc;

    descr = fget = fset = fdel = doc = NULL;

    if ((fget = create_function(vd->vd_getter)) == NULL)
        goto done;

    if ((fset = create_function(vd->vd_setter)) == NULL)
        goto done;

    if ((fdel = create_function(vd->vd_deleter)) == NULL)
        goto done;

    if (vd->vd_docstring == NULL)
    {
        doc = Py_None;
        Py_INCREF(doc);
    }
    else if ((doc = PyUnicode_FromString(vd->vd_docstring)) == NULL)
    {
        goto done;
    }

    descr = PyObject_CallFunctionObjArgs((PyObject *)&PyProperty_Type, fget,
            fset, fdel, doc, NULL);

done:
    Py_XDECREF(fget);
    Py_XDECREF(fset);
    Py_XDECREF(fdel);
    Py_XDECREF(doc);

    return descr;
}


/*
 * Return a PyCFunction as an object or Py_None if there isn't one.
 */
static PyObject *create_function(PyMethodDef *ml)
{
    if (ml != NULL)
        return PyCFunction_New(ml, NULL);

    Py_INCREF(Py_None);
    return Py_None;
}


/*
 * Populate a type dictionary with all lazy attributes if it hasn't already
 * been done.
 */
static int add_lazy_attrs(sipTypeDef *td)
{
    sipWrapperType *wt = (sipWrapperType *)sipTypeAsPyTypeObject(td);
    PyObject *dict;
    sipAttrGetter *ag;

    /* Handle the trivial case. */
    if (wt->wt_dict_complete)
        return 0;

    dict = ((PyTypeObject *)wt)->tp_dict;

    if (sipTypeIsMapped(td))
    {
        if (add_lazy_container_attrs(td, &((sipMappedTypeDef *)td)->mtd_container, dict) < 0)
            return -1;
    }
    else
    {
        sipClassTypeDef *nsx;

        /* Search the possible linked list of namespace extenders. */
        for (nsx = (sipClassTypeDef *)td; nsx != NULL; nsx = nsx->ctd_nsextender)
            if (add_lazy_container_attrs((sipTypeDef *)nsx, &nsx->ctd_container, dict) < 0)
                return -1;
    }

    /*
     * Get any lazy attributes from registered getters.  This must be done last
     * to allow any existing attributes to be replaced.
     */
    /* TODO: Deprecate this mechanism in favour of an event handler. */
    for (ag = sipAttrGetters; ag != NULL; ag = ag->next)
        if (ag->type == NULL || PyType_IsSubtype((PyTypeObject *)wt, ag->type))
            if (ag->getter(td, dict) < 0)
                return -1;

    wt->wt_dict_complete = TRUE;

    PyType_Modified((PyTypeObject *)wt);

    return 0;
}


/*
 * Populate the type dictionary and all its super-types.
 */
static int add_all_lazy_attrs(sipTypeDef *td)
{
    if (td == NULL)
        return 0;

    if (add_lazy_attrs(td) < 0)
        return -1;

    if (sipTypeIsClass(td))
    {
        sipClassTypeDef *ctd = (sipClassTypeDef *)td;
        sipEncodedTypeDef *sup;

        if ((sup = ctd->ctd_supers) != NULL)
            do
            {
                sipTypeDef *sup_td = getGeneratedType(sup, td->td_module);

                if (add_all_lazy_attrs(sup_td) < 0)
                    return -1;
            }
            while (!sup++->sc_flag);
    }

    return 0;
}


/*
 * Return the generated type structure corresponding to the given Python type
 * object.
 */
static const sipTypeDef *sip_api_type_from_py_type_object(PyTypeObject *py_type)
{
    if (PyObject_TypeCheck((PyObject *)py_type, &sipWrapperType_Type))
        return ((sipWrapperType *)py_type)->wt_td;

    if (PyObject_TypeCheck((PyObject *)py_type, &sipEnumType_Type))
        return ((sipEnumTypeObject *)py_type)->type;

    return NULL;
}


/*
 * Return the generated type structure corresponding to the scope of the given
 * type.
 */
static const sipTypeDef *sip_api_type_scope(const sipTypeDef *td)
{
    if (sipTypeIsEnum(td) || sipTypeIsScopedEnum(td))
    {
        const sipEnumTypeDef *etd = (const sipEnumTypeDef *)td;

        if (etd->etd_scope >= 0)
            return td->td_module->em_types[etd->etd_scope];
    }
    else
    {
        const sipContainerDef *cod;

        if (sipTypeIsMapped(td))
            cod = &((const sipMappedTypeDef *)td)->mtd_container;
        else
            cod = &((const sipClassTypeDef *)td)->ctd_container;

        if (!cod->cod_scope.sc_flag)
            return getGeneratedType(&cod->cod_scope, td->td_module);
    }

    return NULL;
}


/*
 * Return TRUE if an object can be converted to a named enum.
 */
static int sip_api_can_convert_to_enum(PyObject *obj, const sipTypeDef *td)
{
    assert(sipTypeIsEnum(td));

    /* If the object is an enum then it must be the right enum. */
    if (PyObject_TypeCheck((PyObject *)Py_TYPE(obj), &sipEnumType_Type))
        return (PyObject_TypeCheck(obj, sipTypeAsPyTypeObject(td)));

    return PyLong_Check(obj);
}


/*
 * Convert a Python object implementing a named enum to an integer value.
 */
static int sip_api_convert_to_enum(PyObject *obj, const sipTypeDef *td)
{
    return convert_to_enum(obj, td, TRUE);
}


/*
 * Convert a Python object implementing a named enum (or, optionally, an int)
 * to an integer value.
 */
static int convert_to_enum(PyObject *obj, const sipTypeDef *td, int allow_int)
{
    int val;

    assert(sipTypeIsEnum(td) || sipTypeIsScopedEnum(td));

    if (sipTypeIsScopedEnum(td))
    {
        static PyObject *value = NULL;
        PyObject *val_obj;

        if (PyObject_IsInstance(obj, (PyObject *)sipTypeAsPyTypeObject(td)) <= 0)
        {
            enum_expected(obj, td);
            return -1;
        }

        if (objectify("value", &value) < 0)
            return -1;

        if ((val_obj = PyObject_GetAttr(obj, value)) == NULL)
            return -1;

        /* This will never overflow. */
        val = long_as_nonoverflow_int(val_obj);

        Py_DECREF(val_obj);
    }
    else
    {
        if (PyObject_TypeCheck((PyObject *)Py_TYPE(obj), &sipEnumType_Type))
        {
            if (!PyObject_TypeCheck(obj, sipTypeAsPyTypeObject(td)))
            {
                enum_expected(obj, td);
                return -1;
            }

            /* This will never overflow. */
            val = long_as_nonoverflow_int(obj);
        }
        else if (allow_int && PyLong_Check(obj))
        {
            val = long_as_nonoverflow_int(obj);
        }
        else
        {
            enum_expected(obj, td);
            return -1;
        }
    }

    return val;
}


/*
 * Raise an exception when failing to convert an enum because of its type.
 */
static void enum_expected(PyObject *obj, const sipTypeDef *td)
{
    PyErr_Format(PyExc_TypeError, "a member of enum '%s' is expected not '%s'",
            sipPyNameOfEnum((sipEnumTypeDef *)td), Py_TYPE(obj)->tp_name);
}


/* Convert to a C/C++ int while checking for overflow. */
static int long_as_nonoverflow_int(PyObject *val_obj)
{
    int old_overflow, val;

    old_overflow = sip_api_enable_overflow_checking(TRUE);
    val = sip_api_long_as_int(val_obj);
    sip_api_enable_overflow_checking(old_overflow);

    return val;
}


/*
 * Create a Python object for a member of a named enum.
 */
static PyObject *sip_api_convert_from_enum(int eval, const sipTypeDef *td)
{
    assert(sipTypeIsEnum(td) || sipTypeIsScopedEnum(td));

    return PyObject_CallFunction((PyObject *)sipTypeAsPyTypeObject(td), "(i)",
            eval);
}


/*
 * Register a getter for unknown attributes.
 */
static int sip_api_register_attribute_getter(const sipTypeDef *td,
        sipAttrGetterFunc getter)
{
    sipAttrGetter *ag = sip_api_malloc(sizeof (sipAttrGetter));

    if (ag == NULL)
        return -1;

    ag->type = sipTypeAsPyTypeObject(td);
    ag->getter = getter;
    ag->next = sipAttrGetters;

    sipAttrGetters = ag;

    return 0;
}


/*
 * Register a proxy resolver.
 */
static int sip_api_register_proxy_resolver(const sipTypeDef *td,
        sipProxyResolverFunc resolver)
{
    sipProxyResolver *pr = sip_api_malloc(sizeof (sipProxyResolver));

    if (pr == NULL)
        return -1;

    pr->td = td;
    pr->resolver = resolver;
    pr->next = proxyResolvers;

    proxyResolvers = pr;

    return 0;
}


/*
 * Report a function with invalid argument types.
 */
static void sip_api_no_function(PyObject *parseErr, const char *func,
        const char *doc)
{
    sip_api_no_method(parseErr, NULL, func, doc);
}


/*
 * Report a method/function/signal with invalid argument types.
 */
static void sip_api_no_method(PyObject *parseErr, const char *scope,
        const char *method, const char *doc)
{
    const char *sep = ".";

    if (scope == NULL)
        scope = ++sep;

    if (parseErr == NULL)
    {
        /*
         * If we have got this far without trying a parse then there must be no
         * overloads.
         */
        PyErr_Format(PyExc_TypeError, "%s%s%s() is a private method", scope,
                sep, method);
    }
    else if (PyList_Check(parseErr))
    {
        PyObject *exc;

        /* There is an entry for each overload that was tried. */
        if (PyList_GET_SIZE(parseErr) == 1)
        {
            PyObject *detail = detail_FromFailure(
                    PyList_GET_ITEM(parseErr, 0));

            if (detail != NULL)
            {
                if (doc != NULL)
                {
                    PyObject *doc_obj = signature_FromDocstring(doc, 0);

                    if (doc_obj != NULL)
                    {
                        exc = PyUnicode_FromFormat("%U: %U", doc_obj, detail);
                        Py_DECREF(doc_obj);
                    }
                    else
                    {
                        exc = NULL;
                    }
                }
                else
                {
                    exc = PyUnicode_FromFormat("%s%s%s(): %U", scope, sep,
                            method, detail);
                }

                Py_DECREF(detail);
            }
            else
            {
                exc = NULL;
            }
        }
        else
        {
            static const char *summary = "arguments did not match any overloaded call:";

            Py_ssize_t i;

            if (doc != NULL)
                exc = PyUnicode_FromString(summary);
            else
                exc = PyUnicode_FromFormat("%s%s%s(): %s", scope, sep, method,
                        summary);

            for (i = 0; i < PyList_GET_SIZE(parseErr); ++i)
            {
                PyObject *failure;
                PyObject *detail = detail_FromFailure(
                        PyList_GET_ITEM(parseErr, i));

                if (detail != NULL)
                {
                    if (doc != NULL)
                    {
                        PyObject *doc_obj = signature_FromDocstring(doc, i);

                        if (doc_obj != NULL)
                        {
                            failure = PyUnicode_FromFormat("\n  %U: %U",
                                    doc_obj, detail);

                            Py_DECREF(doc_obj);
                        }
                        else
                        {
                            Py_XDECREF(exc);
                            exc = NULL;
                            break;
                        }
                    }
                    else
                    {
                        failure = PyUnicode_FromFormat("\n  overload %zd: %U",
                                i + 1, detail);
                    }

                    Py_DECREF(detail);

                    PyUnicode_AppendAndDel(&exc, failure);
                }
                else
                {
                    Py_XDECREF(exc);
                    exc = NULL;
                    break;
                }
            }
        }

        if (exc != NULL)
        {
            PyErr_SetObject(PyExc_TypeError, exc);
            Py_DECREF(exc);
        }
    }
    else
    {
        /*
         * None is used as a marker to say that an exception has already been
         * raised.
         */
        assert(parseErr == Py_None);
    }

    Py_XDECREF(parseErr);
}


/*
 * Return a string/unicode object extracted from a particular line of a
 * docstring.
 */
static PyObject *signature_FromDocstring(const char *doc, Py_ssize_t line)
{
    const char *eol;
    Py_ssize_t size = 0;

    /*
     * Find the start of the line.  If there is a non-default versioned
     * overload that has been enabled then it won't have an entry in the
     * docstring.  This means that the returned signature may be incorrect.
     */
    while (line-- > 0)
    {
        const char *next = strchr(doc, '\n');

        if (next == NULL)
            break;

        doc = next + 1;
    }

    /* Find the last closing parenthesis. */
    for (eol = doc; *eol != '\n' && *eol != '\0'; ++eol)
        if (*eol == ')')
            size = eol - doc + 1;

    return PyUnicode_FromStringAndSize(doc, size);
}


/*
 * Return a string/unicode object that describes the given failure.
 */
static PyObject *detail_FromFailure(PyObject *failure_obj)
{
    sipParseFailure *failure;
    PyObject *detail;

    failure = (sipParseFailure *)PyCapsule_GetPointer(failure_obj, NULL);

    switch (failure->reason)
    {
    case Unbound:
        detail = PyUnicode_FromFormat(
                "first argument of unbound method must have type '%s'",
                failure->detail_str);
        break;

    case TooFew:
        detail = PyUnicode_FromString("not enough arguments");
        break;

    case TooMany:
        detail = PyUnicode_FromString("too many arguments");
        break;

    case KeywordNotString:
        detail = PyUnicode_FromFormat(
                "%S keyword argument name is not a string",
                failure->detail_obj);
        break;

    case UnknownKeyword:
        detail = PyUnicode_FromFormat("'%U' is not a valid keyword argument",
                failure->detail_obj);
        break;

    case Duplicate:
        detail = PyUnicode_FromFormat(
                "'%U' has already been given as a positional argument",
                failure->detail_obj);
        break;

    case WrongType:
        if (failure->arg_nr >= 0)
            detail = bad_type_str(failure->arg_nr, failure->detail_obj);
        else
            detail = PyUnicode_FromFormat(
                    "argument '%s' has unexpected type '%s'",
                    failure->arg_name, Py_TYPE(failure->detail_obj)->tp_name);

        break;

    case Exception:
        detail = failure->detail_obj;

        if (detail)
        {
            Py_INCREF(detail);
            break;
        }

        /* Drop through. */

    default:
        detail = PyUnicode_FromString("unknown reason");
    }

    return detail;
}


/*
 * Report an abstract method called with an unbound self.
 */
static void sip_api_abstract_method(const char *classname, const char *method)
{
    PyErr_Format(PyExc_TypeError,
            "%s.%s() is abstract and cannot be called as an unbound method",
            classname, method);
}


/*
 * Report a deprecated class or method.
 */
int sip_api_deprecated(const char *classname, const char *method)
{
    char buf[100];

    if (classname == NULL)
        PyOS_snprintf(buf, sizeof (buf), "%s() is deprecated", method);
    else if (method == NULL)
        PyOS_snprintf(buf, sizeof (buf), "%s constructor is deprecated",
                classname);
    else
        PyOS_snprintf(buf, sizeof (buf), "%s.%s() is deprecated", classname,
                method);

    return PyErr_WarnEx(PyExc_DeprecationWarning, buf, 1);
}


/*
 * Report a bad operator argument.  Only a small subset of operators need to
 * be handled (those that don't return Py_NotImplemented).
 */
static void sip_api_bad_operator_arg(PyObject *self, PyObject *arg,
        sipPySlotType st)
{
    const char *sn = NULL;

    /* Try and get the text to match a Python exception. */

    switch (st)
    {
    case concat_slot:
    case iconcat_slot:
        PyErr_Format(PyExc_TypeError,
                "cannot concatenate '%s' and '%s' objects",
                Py_TYPE(self)->tp_name, Py_TYPE(arg)->tp_name);
        break;

    case repeat_slot:
        sn = "*";
        break;

    case irepeat_slot:
        sn = "*=";
        break;

    default:
        sn = "unknown";
    }

    if (sn != NULL)
        PyErr_Format(PyExc_TypeError,
                "unsupported operand type(s) for %s: '%s' and '%s'", sn,
                Py_TYPE(self)->tp_name, Py_TYPE(arg)->tp_name);
}


/*
 * Report a sequence length that does not match the length of a slice.
 */
static void sip_api_bad_length_for_slice(Py_ssize_t seqlen,
        Py_ssize_t slicelen)
{
    PyErr_Format(PyExc_ValueError,
            "attempt to assign sequence of size %zd to slice of size %zd",
            seqlen, slicelen);
}


/*
 * Report a Python object that cannot be converted to a particular class.
 */
static void sip_api_bad_class(const char *classname)
{
    PyErr_Format(PyExc_TypeError,
            "cannot convert Python object to an instance of %s", classname);
}


/*
 * Report a Python member function with an unexpected result.
 */
static void sip_api_bad_catcher_result(PyObject *method)
{
    PyObject *mname, *etype, *evalue, *etraceback;

    /*
     * Get the current exception object if there is one.  Its string
     * representation will be used as the detail of a new exception.
     */
    PyErr_Fetch(&etype, &evalue, &etraceback);
    PyErr_NormalizeException(&etype, &evalue, &etraceback);
    Py_XDECREF(etraceback);

    /*
     * This is part of the public API so we make no assumptions about the
     * method object.
     */
    if (!PyMethod_Check(method) ||
        PyMethod_GET_FUNCTION(method) == NULL ||
        !PyFunction_Check(PyMethod_GET_FUNCTION(method)) ||
        PyMethod_GET_SELF(method) == NULL)
    {
        PyErr_Format(PyExc_TypeError,
                "invalid argument to sipBadCatcherResult()");
        return;
    }

    mname = ((PyFunctionObject *)PyMethod_GET_FUNCTION(method))->func_name;

    if (evalue != NULL)
    {
        PyErr_Format(etype, "invalid result from %s.%U(), %S",
                Py_TYPE(PyMethod_GET_SELF(method))->tp_name, mname, evalue);
        Py_DECREF(evalue);
    }
    else
    {
        PyErr_Format(PyExc_TypeError, "invalid result from %s.%U()",
                Py_TYPE(PyMethod_GET_SELF(method))->tp_name, mname);
    }

    Py_XDECREF(etype);
}


/*
 * Transfer ownership of a class instance to Python from C/C++.
 */
static void sip_api_transfer_back(PyObject *self)
{
    if (self != NULL && PyObject_TypeCheck(self, (PyTypeObject *)&sipWrapper_Type))
    {
        sipSimpleWrapper *sw = (sipSimpleWrapper *)self;

        if (sipCppHasRef(sw))
        {
            sipResetCppHasRef(sw);
            Py_DECREF(sw);
        }
        else
        {
            removeFromParent((sipWrapper *)sw);
        }

        sipSetPyOwned(sw);
    }
}


/*
 * Break the association of a C++ owned Python object with any parent.  This is
 * deprecated because it is the equivalent of sip_api_transfer_to(self, NULL).
 */
static void sip_api_transfer_break(PyObject *self)
{
    if (self != NULL && PyObject_TypeCheck(self, (PyTypeObject *)&sipWrapper_Type))
    {
        sipSimpleWrapper *sw = (sipSimpleWrapper *)self;

        if (sipCppHasRef(sw))
        {
            sipResetCppHasRef(sw);
            Py_DECREF(sw);
        }
        else
        {
            removeFromParent((sipWrapper *)sw);
        }
    }
}


/*
 * Transfer ownership of a class instance to C/C++ from Python.
 */
static void sip_api_transfer_to(PyObject *self, PyObject *owner)
{
    /*
     * There is a legitimate case where we try to transfer a PyObject that
     * may not be a SIP generated class.  The virtual handler code calls
     * this function to keep the C/C++ instance alive when it gets rid of
     * the Python object returned by the Python method.  A class may have
     * handwritten code that converts a regular Python type - so we can't
     * assume that we can simply cast to sipWrapper.
     */
    if (self != NULL && PyObject_TypeCheck(self, (PyTypeObject *)&sipWrapper_Type))
    {
        sipSimpleWrapper *sw = (sipSimpleWrapper *)self;

        if (owner == NULL)
        {
            /* There is no owner. */

            if (sipCppHasRef(sw))
            {
                sipResetCppHasRef(sw);
            }
            else
            {
                Py_INCREF(sw);
                removeFromParent((sipWrapper *)sw);
                sipResetPyOwned(sw);
            }

            Py_DECREF(sw);
        }
        else if (owner == Py_None)
        {
            /*
             * The owner is a C++ instance and not a Python object (ie. there
             * is no parent) so there is an explicit extra reference to keep
             * this Python object alive.  Note that there is no way to
             * specify this from a .sip file - it is useful when embedding in
             * C/C++ applications.
             */

            if (!sipCppHasRef(sw))
            {
                Py_INCREF(sw);
                removeFromParent((sipWrapper *)sw);
                sipResetPyOwned(sw);

                sipSetCppHasRef(sw);
            }
        }
        else if (PyObject_TypeCheck(owner, (PyTypeObject *)&sipWrapper_Type))
        {
            /*
             * The owner is a Python object (ie. the C++ instance that the
             * Python object wraps).
             */

            if (sipCppHasRef(sw))
            {
                sipResetCppHasRef(sw);
            }
            else
            {
                Py_INCREF(sw);
                removeFromParent((sipWrapper *)sw);
                sipResetPyOwned(sw);
            }

            addToParent((sipWrapper *)sw, (sipWrapper *)owner);

            Py_DECREF(sw);
        }
    }
}


/*
 * Add a license to a dictionary.
 */
static int addLicense(PyObject *dict,sipLicenseDef *lc)
{
    int rc;
    PyObject *ldict, *proxy, *o;

    /* Convert the strings we use to objects if not already done. */

    if (objectify("__license__", &licenseName) < 0)
        return -1;

    if (objectify("Licensee", &licenseeName) < 0)
        return -1;

    if (objectify("Type", &typeName) < 0)
        return -1;

    if (objectify("Timestamp", &timestampName) < 0)
        return -1;

    if (objectify("Signature", &signatureName) < 0)
        return -1;

    /* We use a dictionary to hold the license information. */
    if ((ldict = PyDict_New()) == NULL)
        return -1;

    /* The license type is compulsory, the rest are optional. */
    if (lc->lc_type == NULL)
        goto deldict;

    if ((o = PyUnicode_FromString(lc->lc_type)) == NULL)
        goto deldict;

    rc = PyDict_SetItem(ldict,typeName,o);
    Py_DECREF(o);

    if (rc < 0)
        goto deldict;

    if (lc->lc_licensee != NULL)
    {
        if ((o = PyUnicode_FromString(lc->lc_licensee)) == NULL)
            goto deldict;

        rc = PyDict_SetItem(ldict,licenseeName,o);
        Py_DECREF(o);

        if (rc < 0)
            goto deldict;
    }

    if (lc->lc_timestamp != NULL)
    {
        if ((o = PyUnicode_FromString(lc->lc_timestamp)) == NULL)
            goto deldict;

        rc = PyDict_SetItem(ldict,timestampName,o);
        Py_DECREF(o);

        if (rc < 0)
            goto deldict;
    }

    if (lc->lc_signature != NULL)
    {
        if ((o = PyUnicode_FromString(lc->lc_signature)) == NULL)
            goto deldict;

        rc = PyDict_SetItem(ldict,signatureName,o);
        Py_DECREF(o);

        if (rc < 0)
            goto deldict;
    }

    /* Create a read-only proxy. */
    if ((proxy = PyDictProxy_New(ldict)) == NULL)
        goto deldict;

    Py_DECREF(ldict);

    rc = PyDict_SetItem(dict, licenseName, proxy);
    Py_DECREF(proxy);

    return rc;

deldict:
    Py_DECREF(ldict);

    return -1;
}


/*
 * Add the void pointer instances to a dictionary.
 */
static int addVoidPtrInstances(PyObject *dict,sipVoidPtrInstanceDef *vi)
{
    while (vi->vi_name != NULL)
    {
        PyObject *w = sip_api_convert_from_void_ptr(vi->vi_val);

        if (dict_set_and_discard(dict, vi->vi_name, w) < 0)
            return -1;

        ++vi;
    }

    return 0;
}


/*
 * Add the char instances to a dictionary.
 */
static int addCharInstances(PyObject *dict, sipCharInstanceDef *ci)
{
    while (ci->ci_name != NULL)
    {
        PyObject *w;

        switch (ci->ci_encoding)
        {
        case 'A':
            w = PyUnicode_DecodeASCII(&ci->ci_val, 1, NULL);
            break;

        case 'L':
            w = PyUnicode_DecodeLatin1(&ci->ci_val, 1, NULL);
            break;

        case '8':
            w = PyUnicode_FromStringAndSize(&ci->ci_val, 1);
            break;

        default:
            w = PyBytes_FromStringAndSize(&ci->ci_val, 1);
        }

        if (dict_set_and_discard(dict, ci->ci_name, w) < 0)
            return -1;

        ++ci;
    }

    return 0;
}


/*
 * Add the string instances to a dictionary.
 */
static int addStringInstances(PyObject *dict, sipStringInstanceDef *si)
{
    while (si->si_name != NULL)
    {
        PyObject *w;

        switch (si->si_encoding)
        {
        case 'A':
            w = PyUnicode_DecodeASCII(si->si_val, strlen(si->si_val), NULL);
            break;

        case 'L':
            w = PyUnicode_DecodeLatin1(si->si_val, strlen(si->si_val), NULL);
            break;

        case '8':
            w = PyUnicode_FromString(si->si_val);
            break;

        case 'w':
            /* The hack for wchar_t. */
#if defined(HAVE_WCHAR_H)
            w = PyUnicode_FromWideChar((const wchar_t *)si->si_val, 1);
            break;
#else
            raiseNoWChar();
            return -1;
#endif

        case 'W':
            /* The hack for wchar_t*. */
#if defined(HAVE_WCHAR_H)
            w = PyUnicode_FromWideChar((const wchar_t *)si->si_val,
                    wcslen((const wchar_t *)si->si_val));
            break;
#else
            raiseNoWChar();
            return -1;
#endif

        default:
            w = PyBytes_FromString(si->si_val);
        }

        if (dict_set_and_discard(dict, si->si_name, w) < 0)
            return -1;

        ++si;
    }

    return 0;
}


/*
 * Add the int instances to a dictionary.
 */
static int addIntInstances(PyObject *dict, sipIntInstanceDef *ii)
{
    while (ii->ii_name != NULL)
    {
        PyObject *w = PyLong_FromLong(ii->ii_val);

        if (dict_set_and_discard(dict, ii->ii_name, w) < 0)
            return -1;

        ++ii;
    }

    return 0;
}


/*
 * Add the long instances to a dictionary.
 */
static int addLongInstances(PyObject *dict,sipLongInstanceDef *li)
{
    while (li->li_name != NULL)
    {
        PyObject *w = PyLong_FromLong(li->li_val);

        if (dict_set_and_discard(dict, li->li_name, w) < 0)
            return -1;

        ++li;
    }

    return 0;
}


/*
 * Add the unsigned long instances to a dictionary.
 */
static int addUnsignedLongInstances(PyObject *dict, sipUnsignedLongInstanceDef *uli)
{
    while (uli->uli_name != NULL)
    {
        PyObject *w = PyLong_FromUnsignedLong(uli->uli_val);

        if (dict_set_and_discard(dict, uli->uli_name, w) < 0)
            return -1;

        ++uli;
    }

    return 0;
}


/*
 * Add the long long instances to a dictionary.
 */
static int addLongLongInstances(PyObject *dict, sipLongLongInstanceDef *lli)
{
    while (lli->lli_name != NULL)
    {
        PyObject *w;

#if defined(HAVE_LONG_LONG)
        w = PyLong_FromLongLong(lli->lli_val);
#else
        w = PyLong_FromLong(lli->lli_val);
#endif

        if (dict_set_and_discard(dict, lli->lli_name, w) < 0)
            return -1;

        ++lli;
    }

    return 0;
}


/*
 * Add the unsigned long long instances to a dictionary.
 */
static int addUnsignedLongLongInstances(PyObject *dict, sipUnsignedLongLongInstanceDef *ulli)
{
    while (ulli->ulli_name != NULL)
    {
        PyObject *w;

#if defined(HAVE_LONG_LONG)
        w = PyLong_FromUnsignedLongLong(ulli->ulli_val);
#else
        w = PyLong_FromUnsignedLong(ulli->ulli_val);
#endif

        if (dict_set_and_discard(dict, ulli->ulli_name, w) < 0)
            return -1;

        ++ulli;
    }

    return 0;
}


/*
 * Add the double instances to a dictionary.
 */
static int addDoubleInstances(PyObject *dict,sipDoubleInstanceDef *di)
{
    while (di->di_name != NULL)
    {
        PyObject *w = PyFloat_FromDouble(di->di_val);

        if (dict_set_and_discard(dict, di->di_name, w) < 0)
            return -1;

        ++di;
    }

    return 0;
}


/*
 * Wrap a set of type instances and add them to a dictionary.
 */
static int addTypeInstances(PyObject *dict, sipTypeInstanceDef *ti)
{
    while (ti->ti_name != NULL)
    {
        if (addSingleTypeInstance(dict, ti->ti_name, ti->ti_ptr, *ti->ti_type, ti->ti_flags) < 0)
            return -1;

        ++ti;
    }

    return 0;
}


/*
 * Wrap a single type instance and add it to a dictionary.
 */
static int addSingleTypeInstance(PyObject *dict, const char *name,
        void *cppPtr, const sipTypeDef *td, int initflags)
{
    PyObject *obj;

    if (sipTypeIsEnum(td) || sipTypeIsScopedEnum(td))
    {
        obj = sip_api_convert_from_enum(*(int *)cppPtr, td);
    }
    else
    {
        sipConvertFromFunc cfrom;

        cppPtr = resolve_proxy(td, cppPtr);

        cfrom = get_from_convertor(td);

        if (cfrom != NULL)
            obj = cfrom(cppPtr, NULL);
        else
            obj = wrap_simple_instance(cppPtr, td, NULL, initflags);
    }

    return dict_set_and_discard(dict, name, obj);
}


/*
 * Convert a type instance and add it to a dictionary.
 */
static int sip_api_add_type_instance(PyObject *dict, const char *name,
        void *cppPtr, const sipTypeDef *td)
{
    return addSingleTypeInstance(getDictFromObject(dict), name, cppPtr, td, 0);
}


/*
 * Return the instance dictionary for an object if it is a wrapped type.
 * Otherwise assume that it is a module dictionary.
 */
static PyObject *getDictFromObject(PyObject *obj)
{
    if (PyObject_TypeCheck(obj, (PyTypeObject *)&sipWrapperType_Type))
        obj = ((PyTypeObject *)obj)->tp_dict;

    return obj;
}


/*
 * Return a Python reimplementation corresponding to a C/C++ virtual function,
 * if any.  If one was found then the GIL is acquired.  This is deprecated, use
 * sip_api_is_python_method_12_8() instead.
 */
static PyObject *sip_api_is_py_method(sip_gilstate_t *gil, char *pymc,
        sipSimpleWrapper *sipSelf, const char *cname, const char *mname)
{
    return sip_api_is_py_method_12_8(gil, pymc, &sipSelf, cname, mname);
}


/*
 * Return a Python reimplementation corresponding to a C/C++ virtual function,
 * if any.  If one was found then the GIL is acquired.
 */
static PyObject *sip_api_is_py_method_12_8(sip_gilstate_t *gil, char *pymc,
        sipSimpleWrapper **sipSelfp, const char *cname, const char *mname)
{
    sipSimpleWrapper *sipSelf;
    PyObject *mname_obj, *reimp, *mro, *cls;
    Py_ssize_t i;

    /*
     * This is the most common case (where there is no Python reimplementation)
     * so we take a fast shortcut without acquiring the GIL.
     */
    if (*pymc != 0)
        return NULL;

    /* We might still have C++ going after the interpreter has gone. */
    if (sipInterpreter == NULL)
        return NULL;

#ifdef WITH_THREAD
    *gil = PyGILState_Ensure();
#endif

    /* Only read this when we have the GIL. */
    sipSelf = *sipSelfp;

    /*
     * It's possible that the Python object has been deleted but the underlying
     * C++ instance is still working and trying to handle virtual functions.
     * Alternatively, an instance has started handling virtual functions before
     * its ctor has returned.  In either case say there is no Python
     * reimplementation.
     */
    if (sipSelf != NULL)
        sipSelf = deref_mixin(sipSelf);

    if (sipSelf == NULL)
        goto release_gil;

    /*
     * It's possible that the object's type's tp_mro is NULL.  A possible
     * circumstance is when a type has been created dynamically and the only
     * reference to it is the single instance of the type which is in the
     * process of being garbage collected.
     */
    cls = (PyObject *)Py_TYPE(sipSelf);
    mro = ((PyTypeObject *)cls)->tp_mro;

    if (mro == NULL)
        goto release_gil;

    /* Get any reimplementation. */

    if ((mname_obj = PyUnicode_FromString(mname)) == NULL)
        goto release_gil;

    /*
     * We don't use PyObject_GetAttr() because that might find the generated
     * C function before a reimplementation defined in a mixin (ie. later in
     * the MRO).  However that means we must explicitly check that the class
     * hierarchy is fully initialised.
     */
    if (add_all_lazy_attrs(((sipWrapperType *)Py_TYPE(sipSelf))->wt_td) < 0)
    {
        Py_DECREF(mname_obj);
        goto release_gil;
    }

    if (sipSelf->dict != NULL)
    {
        /* Check the instance dictionary in case it has been monkey patched. */
        if ((reimp = PyDict_GetItem(sipSelf->dict, mname_obj)) != NULL && PyCallable_Check(reimp))
        {
            Py_DECREF(mname_obj);

            Py_INCREF(reimp);
            return reimp;
        }
    }

    assert(PyTuple_Check(mro));

    reimp = NULL;

    for (i = 0; i < PyTuple_GET_SIZE(mro); ++i)
    {
        PyObject *cls_dict, *cls_attr;

        cls = PyTuple_GET_ITEM(mro, i);

        cls_dict = ((PyTypeObject *)cls)->tp_dict;

        /*
         * Check any possible reimplementation is not the wrapped C++ method or
         * a default special method implementation.
         */
        if (cls_dict != NULL && (cls_attr = PyDict_GetItem(cls_dict, mname_obj)) != NULL && Py_TYPE(cls_attr) != &sipMethodDescr_Type && Py_TYPE(cls_attr) != &PyWrapperDescr_Type)
        {
            reimp = cls_attr;
            break;
        }
    }

    Py_DECREF(mname_obj);

    if (reimp != NULL)
    {
        /*
         * Emulate the behaviour of a descriptor to make sure we return a bound
         * method.
         */
        if (PyMethod_Check(reimp))
        {
            /* It's already a method but make sure it is bound. */
            if (PyMethod_GET_SELF(reimp) != NULL)
                Py_INCREF(reimp);
            else
                reimp = PyMethod_New(PyMethod_GET_FUNCTION(reimp),
                        (PyObject *)sipSelf);
        }
        else if (PyFunction_Check(reimp))
        {
            reimp = PyMethod_New(reimp, (PyObject *)sipSelf);
        }
        else if (Py_TYPE(reimp)->tp_descr_get)
        {
            /* It is a descriptor, so assume it will do the right thing. */
            reimp = Py_TYPE(reimp)->tp_descr_get(reimp, (PyObject *)sipSelf,
                    cls);
        }
        else
        {
            /*
             * We don't know what it is so just return and assume that an
             * appropriate exception will be raised later on.
             */
            Py_INCREF(reimp);
        }
    }
    else
    {
        /* Use the fast track in future. */
        *pymc = 1;

        if (cname != NULL)
        {
            /* Note that this will only be raised once per method. */
            PyErr_Format(PyExc_NotImplementedError,
                    "%s.%s() is abstract and must be overridden", cname,
                    mname);
            PyErr_Print();
        }

#ifdef WITH_THREAD
        PyGILState_Release(*gil);
#endif
    }

    return reimp;

release_gil:
#ifdef WITH_THREAD
    PyGILState_Release(*gil);
#endif
    return NULL;
}


/*
 * Convert a C/C++ pointer to the object that wraps it.
 */
static PyObject *sip_api_get_pyobject(void *cppPtr, const sipTypeDef *td)
{
    return (PyObject *)sipOMFindObject(&cppPyMap, cppPtr, td);
}


/*
 * The default access function.
 */
void *sip_api_get_address(sipSimpleWrapper *w)
{
    return (w->access_func != NULL) ? w->access_func(w, GuardedPointer) : w->data;
}


/*
 * The access function for handwritten access functions.
 */
static void *explicit_access_func(sipSimpleWrapper *sw, AccessFuncOp op)
{
    typedef void *(*explicitAccessFunc)(void);

    if (op == ReleaseGuard)
        return NULL;

    return ((explicitAccessFunc)(sw->data))();
}


/*
 * The access function for indirect access.
 */
static void *indirect_access_func(sipSimpleWrapper *sw, AccessFuncOp op)
{
    void *addr;

    switch (op)
    {
    case UnguardedPointer:
        addr = sw->data;
        break;

    case GuardedPointer:
        addr = *((void **)sw->data);
        break;

    default:
        addr = NULL;
    }

    return addr;
}


/*
 * Get the C/C++ pointer for a complex object.  Note that not casting the C++
 * pointer is a bug.  However this would only ever be called by PyQt3 signal
 * emitter code and PyQt doesn't contain anything that multiply inherits from
 * QObject.
 */
static void *sip_api_get_complex_cpp_ptr(sipSimpleWrapper *sw)
{
    return getComplexCppPtr(sw, NULL);
}


/*
 * Get the C/C++ pointer for a complex object and optionally cast it to the
 * required type.
 */
static void *getComplexCppPtr(sipSimpleWrapper *sw, const sipTypeDef *td)
{
    if (!sipIsDerived(sw))
    {
        PyErr_SetString(PyExc_RuntimeError,
                "no access to protected functions or signals for objects not created from Python");

        return NULL;
    }

    return sip_api_get_cpp_ptr(sw, td);
}


/*
 * Get the C/C++ pointer from a wrapper and optionally cast it to the required
 * type.
 */
void *sip_api_get_cpp_ptr(sipSimpleWrapper *sw, const sipTypeDef *td)
{
    void *ptr = sip_api_get_address(sw);

    if (checkPointer(ptr, sw) < 0)
        return NULL;

    if (td != NULL)
    {
        if (PyObject_TypeCheck((PyObject *)sw, sipTypeAsPyTypeObject(td)))
            ptr = cast_cpp_ptr(ptr, Py_TYPE(sw), td);
        else
            ptr = NULL;

        if (ptr == NULL)
            PyErr_Format(PyExc_TypeError, "could not convert '%s' to '%s'",
                    Py_TYPE(sw)->tp_name,
                    sipPyNameOfContainer(&((const sipClassTypeDef *)td)->ctd_container, td));
    }

    return ptr;
}


/*
 * Cast a C/C++ pointer from a source type to a destination type.
 */
static void *cast_cpp_ptr(void *ptr, PyTypeObject *src_type,
        const sipTypeDef *dst_type)
{
    sipCastFunc cast = ((const sipClassTypeDef *)((sipWrapperType *)src_type)->wt_td)->ctd_cast;

    /* C structures and base classes don't have cast functions. */
    if (cast != NULL)
        ptr = (*cast)(ptr, dst_type);

    return ptr;
}


/*
 * Check that a pointer is non-NULL.
 */
static int checkPointer(void *ptr, sipSimpleWrapper *sw)
{
    if (ptr == NULL)
    {
        PyErr_Format(PyExc_RuntimeError, (sipWasCreated(sw) ?
                        "wrapped C/C++ object of type %s has been deleted" :
                        "super-class __init__() of type %s was never called"),
                Py_TYPE(sw)->tp_name);
        return -1;
    }

    return 0;
}


/*
 * Keep an extra reference to an object.
 */
static void sip_api_keep_reference(PyObject *self, int key, PyObject *obj)
{
    PyObject *dict, *key_obj;

    /*
     * If there isn't a "self" to keep the extra reference for later garbage
     * collection then just take a reference and let it leak.
     */
    if (self == NULL)
    {
        Py_XINCREF(obj);
        return;
    }

    /* Create the extra references dictionary if needed. */
    if ((dict = ((sipSimpleWrapper *)self)->extra_refs) == NULL)
    {
        if ((dict = PyDict_New()) == NULL)
            return;

        ((sipSimpleWrapper *)self)->extra_refs = dict;
    }

    if ((key_obj = PyLong_FromLong(key)) != NULL)
    {
        /* This can happen if the argument was optional. */
        if (obj == NULL)
            obj = Py_None;

        PyDict_SetItem(dict, key_obj, obj);
        Py_DECREF(key_obj);
    }
}


/*
 * Get an object that has an extra reference.
 */
static PyObject *sip_api_get_reference(PyObject *self, int key)
{
    PyObject *dict, *key_obj, *obj;

    /* Get the extra references dictionary if there is one. */
    if ((dict = ((sipSimpleWrapper *)self)->extra_refs) == NULL)
        return NULL;

    if ((key_obj = PyLong_FromLong(key)) == NULL)
        return NULL;

    obj = PyDict_GetItem(dict, key_obj);
    Py_DECREF(key_obj);
    Py_XINCREF(obj);

    return obj;
}


/*
 * Return TRUE if an object is owned by Python.  Note that this isn't
 * implemented as a macro in sip.h because the position of the sw_flags field
 * is dependent on the version of Python.
 */
static int sip_api_is_owned_by_python(sipSimpleWrapper *sw)
{
    return sipIsPyOwned(sw);
}


/*
 * Return TRUE if the type of a C++ instance is a derived class.  Note that
 * this isn't implemented as a macro in sip.h because the position of the
 * sw_flags field is dependent on the version of Python.
 */
static int sip_api_is_derived_class(sipSimpleWrapper *sw)
{
    return sipIsDerived(sw);
}


/*
 * Get the user defined object from a wrapped object.  Note that this isn't
 * implemented as a macro in sip.h because the position of the user field is
 * dependent on the version of Python.
 */
static PyObject *sip_api_get_user_object(const sipSimpleWrapper *sw)
{
    return sw->user;
}


/*
 * Set the user defined object in a wrapped object.  Note that this isn't
 * implemented as a macro in sip.h because the position of the user field is
 * dependent on the version of Python.
 */
static void sip_api_set_user_object(sipSimpleWrapper *sw, PyObject *user)
{
    sw->user = user;
}


/*
 * Check to see if a Python object can be converted to a type.
 */
static int sip_api_can_convert_to_type(PyObject *pyObj, const sipTypeDef *td,
        int flags)
{
    int ok;

    assert(td == NULL || sipTypeIsClass(td) || sipTypeIsMapped(td));

    if (td == NULL)
    {
        /*
         * The type must be /External/ and the module that contains the
         * implementation hasn't been imported.
         */
        ok = FALSE;
    }
    else if (pyObj == Py_None)
    {
        /* If the type explicitly handles None then ignore the flags. */
        if (sipTypeAllowNone(td))
            ok = TRUE;
        else
            ok = ((flags & SIP_NOT_NONE) == 0);
    }
    else
    {
        sipConvertToFunc cto;

        if (sipTypeIsClass(td))
        {
            cto = ((const sipClassTypeDef *)td)->ctd_cto;

            if (cto == NULL || (flags & SIP_NO_CONVERTORS) != 0)
                ok = PyObject_TypeCheck(pyObj, sipTypeAsPyTypeObject(td));
            else
                ok = cto(pyObj, NULL, NULL, NULL);
        }
        else
        {
            cto = ((const sipMappedTypeDef *)td)->mtd_cto;
            ok = cto(pyObj, NULL, NULL, NULL);
        }
    }

    return ok;
}


/*
 * Convert a Python object to a C/C++ pointer, assuming a previous call to
 * sip_api_can_convert_to_type() has been successful.  Allow ownership to be
 * transferred and any type convertors to be disabled.
 */
static void *sip_api_convert_to_type(PyObject *pyObj, const sipTypeDef *td,
        PyObject *transferObj, int flags, int *statep, int *iserrp)
{
    void *cpp = NULL;
    int state = 0;

    assert(sipTypeIsClass(td) || sipTypeIsMapped(td));

    /* Don't convert if there has already been an error. */
    if (!*iserrp)
    {
        /* Do the conversion. */
        if (pyObj == Py_None && !sipTypeAllowNone(td))
            cpp = NULL;
        else
        {
            sipConvertToFunc cto;

            if (sipTypeIsClass(td))
            {
                cto = ((const sipClassTypeDef *)td)->ctd_cto;

                if (cto == NULL || (flags & SIP_NO_CONVERTORS) != 0)
                {
                    if ((cpp = sip_api_get_cpp_ptr((sipSimpleWrapper *)pyObj, td)) == NULL)
                        *iserrp = TRUE;
                    else if (transferObj != NULL)
                    {
                        if (transferObj == Py_None)
                            sip_api_transfer_back(pyObj);
                        else
                            sip_api_transfer_to(pyObj, transferObj);
                    }
                }
                else
                {
                    state = cto(pyObj, &cpp, iserrp, transferObj);
                }
            }
            else
            {
                cto = ((const sipMappedTypeDef *)td)->mtd_cto;
                state = cto(pyObj, &cpp, iserrp, transferObj);
            }
        }
    }

    if (statep != NULL)
        *statep = state;

    return cpp;
}


/*
 * Convert a Python object to a C/C++ pointer and raise an exception if it
 * can't be done.
 */
void *sip_api_force_convert_to_type(PyObject *pyObj, const sipTypeDef *td,
        PyObject *transferObj, int flags, int *statep, int *iserrp)
{
    /* Don't even try if there has already been an error. */
    if (*iserrp)
        return NULL;

    /* See if the object's type can be converted. */
    if (!sip_api_can_convert_to_type(pyObj, td, flags))
    {
        if (sipTypeIsMapped(td))
            PyErr_Format(PyExc_TypeError,
                    "%s cannot be converted to a C/C++ %s in this context",
                    Py_TYPE(pyObj)->tp_name, sipTypeName(td));
        else
            PyErr_Format(PyExc_TypeError,
                    "%s cannot be converted to %s.%s in this context",
                    Py_TYPE(pyObj)->tp_name, sipNameOfModule(td->td_module),
                    sipPyNameOfContainer(&((const sipClassTypeDef *)td)->ctd_container, td));

        if (statep != NULL)
            *statep = 0;

        *iserrp = TRUE;
        return NULL;
    }

    /* Do the conversion. */
    return sip_api_convert_to_type(pyObj, td, transferObj, flags, statep,
            iserrp);
}


/*
 * Release a possibly temporary C/C++ instance created by a type convertor.
 */
static void sip_api_release_type(void *cpp, const sipTypeDef *td, int state)
{
    /* See if there is something to release. */
    if (state & SIP_TEMPORARY)
        release(cpp, td, state);
}


/*
 * Release an instance.
 */
static void release(void *addr, const sipTypeDef *td, int state)
{
    sipReleaseFunc rel;

    if (sipTypeIsClass(td))
    {
        rel = ((const sipClassTypeDef *)td)->ctd_release;

        /*
         * If there is no release function then it must be a C structure and we
         * can just free it.
         */
        if (rel == NULL)
            sip_api_free(addr);
    }
    else if (sipTypeIsMapped(td))
        rel = ((const sipMappedTypeDef *)td)->mtd_release;
    else
        rel = NULL;

    if (rel != NULL)
        rel(addr, state);
}


/*
 * Convert a C/C++ instance to a Python instance.
 */
PyObject *sip_api_convert_from_type(void *cpp, const sipTypeDef *td,
        PyObject *transferObj)
{
    PyObject *py;
    sipConvertFromFunc cfrom;

    assert(sipTypeIsClass(td) || sipTypeIsMapped(td));

    /* Handle None. */
    if (cpp == NULL)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }

    cpp = resolve_proxy(td, cpp);

    cfrom = get_from_convertor(td);

    if (cfrom != NULL)
        return cfrom(cpp, transferObj);

    /*
     * See if we have already wrapped it.  Invoking sub-class code can be
     * expensive so we check the cache first, even though the sub-class code
     * might perform a down-cast.
     */
    if ((py = sip_api_get_pyobject(cpp, td)) == NULL && sipTypeHasSCC(td))
    {
        void *orig_cpp = cpp;
        const sipTypeDef *orig_td = td;

        /* Apply the sub-class convertor. */
        td = convertSubClass(td, &cpp);

        /*
         * If the sub-class convertor has done something then check the cache
         * again using the modified values.
         */
        if (cpp != orig_cpp || td != orig_td)
            py = sip_api_get_pyobject(cpp, td);
    }

    if (py != NULL)
        Py_INCREF(py);
    else if ((py = wrap_simple_instance(cpp, td, NULL, SIP_SHARE_MAP)) == NULL)
        return NULL;

    /* Handle any ownership transfer. */
    if (transferObj != NULL)
    {
        if (transferObj == Py_None)
            sip_api_transfer_back(py);
        else
            sip_api_transfer_to(py, transferObj);
    }

    return py;
}


/*
 * Convert a new C/C++ instance to a Python instance.
 */
static PyObject *sip_api_convert_from_new_type(void *cpp, const sipTypeDef *td,
        PyObject *transferObj)
{
    sipWrapper *owner;
    sipConvertFromFunc cfrom;

    /* Handle None. */
    if (cpp == NULL)
    {
        Py_INCREF(Py_None);
        return Py_None;
    }

    cpp = resolve_proxy(td, cpp);

    cfrom = get_from_convertor(td);

    if (cfrom != NULL)
    {
        PyObject *res = cfrom(cpp, transferObj);

        if (res != NULL)
        {
            /*
             * We no longer need the C/C++ instance so we release it (unless
             * its ownership is transferred).  This means this call is
             * semantically equivalent to the case where we are wrapping a
             * class.
             */
            if (transferObj == NULL || transferObj == Py_None)
                release(cpp, td, 0);
        }

        return res;
    }

    /* Apply any sub-class convertor. */
    if (sipTypeHasSCC(td))
        td = convertSubClass(td, &cpp);

    /* Handle any ownership transfer. */
    if (transferObj == NULL || transferObj == Py_None)
        owner = NULL;
    else
        owner = (sipWrapper *)transferObj;

    return wrap_simple_instance(cpp, td, owner,
            (owner == NULL ? SIP_PY_OWNED : 0));
}


/*
 * Implement the normal transfer policy for the result of %ConvertToTypeCode,
 * ie. it is temporary unless it is being transferred from Python.
 */
int sip_api_get_state(PyObject *transferObj)
{
    return (transferObj == NULL || transferObj == Py_None) ? SIP_TEMPORARY : 0;
}


/*
 * This is set by sip_api_find_type() before calling bsearch() on the types
 * table for the module.  This is a hack that works around the problem of
 * unresolved externally defined types.
 */
static sipExportedModuleDef *module_searched;


/*
 * The bsearch() helper function for searching the types table.
 */
static int compareTypeDef(const void *key, const void *el)
{
    const char *s1 = (const char *)key;
    const char *s2 = NULL;
    const sipTypeDef *td;
    char ch1, ch2;

    /* Allow for unresolved externally defined types. */
    td = *(const sipTypeDef **)el;

    if (td != NULL)
    {
        s2 = sipTypeName(td);
    }
    else
    {
        sipExternalTypeDef *etd = module_searched->em_external;

        assert(etd != NULL);

        /* Find which external type it is. */
        while (etd->et_nr >= 0)
        {
            const void *tdp = &module_searched->em_types[etd->et_nr];

            if (tdp == el)
            {
                s2 = etd->et_name;
                break;
            }

            ++etd;
        }

        assert(s2 != NULL);
    }

    /*
     * Compare while ignoring spaces so that we don't impose a rigorous naming
     * standard.  This only really affects template-based mapped types.
     */
    do
    {
        while ((ch1 = *s1++) == ' ')
            ;

        while ((ch2 = *s2++) == ' ')
            ;

        /* We might be looking for a pointer or a reference. */
        if ((ch1 == '*' || ch1 == '&' || ch1 == '\0') && ch2 == '\0')
            return 0;
    }
    while (ch1 == ch2);

    return (ch1 < ch2 ? -1 : 1);
}


/*
 * Return the type structure for a particular type.
 */
static const sipTypeDef *sip_api_find_type(const char *type)
{
    sipExportedModuleDef *em;

    for (em = moduleList; em != NULL; em = em->em_next)
    {
        sipTypeDef **tdp;

        /* The backdoor to the comparison helper. */
        module_searched = em;

        tdp = (sipTypeDef **)bsearch((const void *)type,
                (const void *)em->em_types, em->em_nrtypes,
                sizeof (sipTypeDef *), compareTypeDef);

        if (tdp != NULL)
        {
            /*
             * Note that this will be NULL for unresolved externally defined
             * types.
             */
            return *tdp;
        }
    }

    return NULL;
}


/*
 * Return the mapped type structure for a particular mapped type.  This is
 * deprecated.
 */
static const sipMappedType *sip_api_find_mapped_type(const char *type)
{
    const sipTypeDef *td = sip_api_find_type(type);

    if (td != NULL && sipTypeIsMapped(td))
        return (const sipMappedType *)td;

    return NULL;
}


/*
 * Return the type structure for a particular class.  This is deprecated.
 */
static sipWrapperType *sip_api_find_class(const char *type)
{
    const sipTypeDef *td = sip_api_find_type(type);

    if (td != NULL && sipTypeIsClass(td))
        return (sipWrapperType *)sipTypeAsPyTypeObject(td);

    return NULL;
}


/*
 * Return the type structure for a particular named unscoped enum.  This is
 * deprecated.
 */
static PyTypeObject *sip_api_find_named_enum(const char *type)
{
    const sipTypeDef *td = sip_api_find_type(type);

    if (td != NULL && sipTypeIsEnum(td))
        return sipTypeAsPyTypeObject(td);

    return NULL;
}


/*
 * Save the components of a Python method.
 */
void sipSaveMethod(sipPyMethod *pm, PyObject *meth)
{
    pm->mfunc = PyMethod_GET_FUNCTION(meth);
    pm->mself = PyMethod_GET_SELF(meth);
}


/*
 * Call a hook.
 */
static void sip_api_call_hook(const char *hookname)
{
    PyObject *dictofmods, *mod, *dict, *hook, *res;
 
    /* Get the dictionary of modules. */
    if ((dictofmods = PyImport_GetModuleDict()) == NULL)
        return;
 
    /* Get the builtins module. */
    if ((mod = PyDict_GetItemString(dictofmods, "builtins")) == NULL)
        return;
 
    /* Get it's dictionary. */
    if ((dict = PyModule_GetDict(mod)) == NULL)
        return;
 
    /* Get the function hook. */
    if ((hook = PyDict_GetItemString(dict, hookname)) == NULL)
        return;
 
    /* Call the hook and discard any result. */
    res = PyObject_Call(hook, empty_tuple, NULL);
 
    Py_XDECREF(res);
}


/*
 * Call any sub-class convertors for a given type returning a pointer to the
 * sub-type object, and possibly modifying the C++ address (in the case of
 * multiple inheritence).
 */
static const sipTypeDef *convertSubClass(const sipTypeDef *td, void **cppPtr)
{
    /* Handle the trivial case. */
    if (*cppPtr == NULL)
        return NULL;

    /* Try the conversions until told to stop. */
    while (convertPass(&td, cppPtr))
        ;

    return td;
}


/*
 * Do a single pass through the available convertors.
 */
static int convertPass(const sipTypeDef **tdp, void **cppPtr)
{
    PyTypeObject *py_type = sipTypeAsPyTypeObject(*tdp);
    sipExportedModuleDef *em;

    /*
     * Note that this code depends on the fact that a module appears in the
     * list of modules before any module it imports, ie. sub-class convertors
     * will be invoked for more specific types first.
     */
    for (em = moduleList; em != NULL; em = em->em_next)
    {
        sipSubClassConvertorDef *scc;

        if ((scc = em->em_convertors) == NULL)
            continue;

        while (scc->scc_convertor != NULL)
        {
            PyTypeObject *base_type = sipTypeAsPyTypeObject(scc->scc_basetype);

            /*
             * The base type is the "root" class that may have a number of
             * convertors each handling a "branch" of the derived tree of
             * classes.  The "root" normally implements the base function that
             * provides the RTTI used by the convertors and is re-implemented
             * by derived classes.  We therefore see if the target type is a
             * sub-class of the root, ie. see if the convertor might be able to
             * convert the target type to something more specific.
             */
            if (PyType_IsSubtype(py_type, base_type))
            {
                void *ptr;
                const sipTypeDef *sub_td;

                ptr = cast_cpp_ptr(*cppPtr, py_type, scc->scc_basetype);

                if ((sub_td = (*scc->scc_convertor)(&ptr)) != NULL)
                {
                    PyTypeObject *sub_type = sipTypeAsPyTypeObject(sub_td);

                    /*
                     * We are only interested in types that are not
                     * super-classes of the target.  This happens either
                     * because it is in an earlier convertor than the one that
                     * handles the type or it is in a later convertor that
                     * handles a different branch of the hierarchy.  Either
                     * way, the ordering of the modules ensures that there will
                     * be no more than one and that it will be the right one.
                     */
                    if (!PyType_IsSubtype(py_type, sub_type))
                    {
                        *tdp = sub_td;
                        *cppPtr = ptr;

                        /*
                         * Finally we allow the convertor to return a type that
                         * is apparently unrelated to the current convertor.
                         * This causes the whole process to be restarted with
                         * the new values.  The use case is PyQt's QLayoutItem.
                         */
                        return !PyType_IsSubtype(sub_type, base_type);
                    }
                }
            }

            ++scc;
        }
    }

    /*
     * We haven't found the exact type, so return the most specific type that
     * it must be.  This can happen legitimately if the wrapped library is
     * returning an internal class that is down-cast to a more generic class.
     * Also we want this function to be safe when a class doesn't have any
     * convertors.
     */
    return FALSE;
}


/*
 * The bsearch() helper function for searching a sorted string map table.
 */
static int compareStringMapEntry(const void *key,const void *el)
{
    return strcmp((const char *)key,((const sipStringTypeClassMap *)el)->typeString);
}


/*
 * A convenience function for %ConvertToSubClassCode for types represented as a
 * string.  Returns the Python class object or NULL if the type wasn't
 * recognised.  This is deprecated.
 */
static sipWrapperType *sip_api_map_string_to_class(const char *typeString,
        const sipStringTypeClassMap *map, int maplen)
{
    sipStringTypeClassMap *me;

    me = (sipStringTypeClassMap *)bsearch((const void *)typeString,
                          (const void *)map,maplen,
                          sizeof (sipStringTypeClassMap),
                          compareStringMapEntry);

        return ((me != NULL) ? *me->pyType : NULL);
}


/*
 * The bsearch() helper function for searching a sorted integer map table.
 */
static int compareIntMapEntry(const void *keyp,const void *el)
{
    int key = *(int *)keyp;

    if (key > ((const sipIntTypeClassMap *)el)->typeInt)
        return 1;

    if (key < ((const sipIntTypeClassMap *)el)->typeInt)
        return -1;

    return 0;
}


/*
 * A convenience function for %ConvertToSubClassCode for types represented as
 * an integer.  Returns the Python class object or NULL if the type wasn't
 * recognised.  This is deprecated.
 */
static sipWrapperType *sip_api_map_int_to_class(int typeInt,
        const sipIntTypeClassMap *map, int maplen)
{
    sipIntTypeClassMap *me;

    me = (sipIntTypeClassMap *)bsearch((const void *)&typeInt,
                       (const void *)map,maplen,
                       sizeof (sipIntTypeClassMap),
                       compareIntMapEntry);

        return ((me != NULL) ? *me->pyType : NULL);
}


/*
 * Raise an unknown exception.  Make no assumptions about the GIL.
 */
static void sip_api_raise_unknown_exception(void)
{
    static PyObject *mobj = NULL;

    SIP_BLOCK_THREADS

    objectify("unknown", &mobj);

    PyErr_SetObject(PyExc_Exception, mobj);

    SIP_UNBLOCK_THREADS
}


/*
 * Raise an exception implemented as a type.  Make no assumptions about the
 * GIL.
 */
static void sip_api_raise_type_exception(const sipTypeDef *td, void *ptr)
{
    PyObject *self;

    assert(sipTypeIsClass(td));

    SIP_BLOCK_THREADS

    self = wrap_simple_instance(ptr, td, NULL, SIP_PY_OWNED);

    PyErr_SetObject((PyObject *)sipTypeAsPyTypeObject(td), self);

    Py_XDECREF(self);

    SIP_UNBLOCK_THREADS
}


/*
 * Return the generated type structure of an encoded type.
 */
static sipTypeDef *getGeneratedType(const sipEncodedTypeDef *enc,
        sipExportedModuleDef *em)
{
    if (enc->sc_module == 255)
        return em->em_types[enc->sc_type];

    return em->em_imports[enc->sc_module].im_imported_types[enc->sc_type].it_td;
}


/*
 * Return the generated class type structure of a class's super-class.
 */
sipClassTypeDef *sipGetGeneratedClassType(const sipEncodedTypeDef *enc,
        const sipClassTypeDef *ctd)
{
    return (sipClassTypeDef *)getGeneratedType(enc, ctd->ctd_base.td_module);
}


/*
 * Find a particular slot function for a type.
 */
static void *findSlot(PyObject *self, sipPySlotType st)
{
    void *slot;
    PyTypeObject *py_type = Py_TYPE(self);

    /* See if it is a wrapper. */
    if (PyObject_TypeCheck((PyObject *)py_type, &sipWrapperType_Type))
    {
        const sipClassTypeDef *ctd;

        ctd = (sipClassTypeDef *)((sipWrapperType *)(py_type))->wt_td;

        slot = findSlotInClass(ctd, st);
    }
    else
    {
        sipEnumTypeDef *etd;

        /* If it is not a wrapper then it must be an enum. */
        assert(PyObject_TypeCheck((PyObject *)py_type, &sipEnumType_Type));

        etd = (sipEnumTypeDef *)((sipEnumTypeObject *)(py_type))->type;

        assert(etd->etd_pyslots != NULL);

        slot = findSlotInSlotList(etd->etd_pyslots, st);
    }

    return slot;
}


/*
 * Find a particular slot function in a class hierarchy.
 */
static void *findSlotInClass(const sipClassTypeDef *ctd, sipPySlotType st)
{
    void *slot;

    if (ctd->ctd_pyslots != NULL)
        slot = findSlotInSlotList(ctd->ctd_pyslots, st);
    else
        slot = NULL;

    if (slot == NULL)
    {
        sipEncodedTypeDef *sup;

        /* Search any super-types. */
        if ((sup = ctd->ctd_supers) != NULL)
        {
            do
            {
                const sipClassTypeDef *sup_ctd = sipGetGeneratedClassType(
                        sup, ctd);

                slot = findSlotInClass(sup_ctd, st);
            }
            while (slot == NULL && !sup++->sc_flag);
        }
    }

    return slot;
}


/*
 * Find a particular slot function in a particular type.
 */
static void *findSlotInSlotList(sipPySlotDef *psd, sipPySlotType st)
{
    while (psd->psd_func != NULL)
    {
        if (psd->psd_type == st)
            return psd->psd_func;

        ++psd;
    }

    return NULL;
}


/*
 * Return the C/C++ address and the generated class structure for a wrapper.
 */
static void *getPtrTypeDef(sipSimpleWrapper *self, const sipClassTypeDef **ctd)
{
    *ctd = (const sipClassTypeDef *)((sipWrapperType *)Py_TYPE(self))->wt_td;

    return (sipNotInMap(self) ? NULL : sip_api_get_address(self));
}


/*
 * Handle an objobjargproc slot.
 */
static int objobjargprocSlot(PyObject *self, PyObject *arg1, PyObject *arg2,
        sipPySlotType st)
{
    int (*f)(PyObject *, PyObject *);
    int res;

    f = (int (*)(PyObject *, PyObject *))findSlot(self, st);

    if (f != NULL)
    {
        PyObject *args;

        /*
         * Slot handlers require a single PyObject *.  The second argument is
         * optional.
         */
        if (arg2 == NULL)
        {
            args = arg1;
            Py_INCREF(args);
        }
        else if ((args = PyTuple_Pack(2, arg1, arg2)) == NULL)
        {
            return -1;
        }

        res = f(self, args);
        Py_DECREF(args);
    }
    else
    {
        PyErr_SetNone(PyExc_NotImplementedError);
        res = -1;
    }

    return res;
}


/*
 * Handle an ssizeobjargproc slot.
 */
static int ssizeobjargprocSlot(PyObject *self, Py_ssize_t arg1,
        PyObject *arg2, sipPySlotType st)
{
    int (*f)(PyObject *, PyObject *);
    int res;

    f = (int (*)(PyObject *, PyObject *))findSlot(self, st);

    if (f != NULL)
    {
        PyObject *args;

        /*
         * Slot handlers require a single PyObject *.  The second argument is
         * optional.
         */
        if (arg2 == NULL)
            args = PyLong_FromSsize_t(arg1);
        else
            args = Py_BuildValue("(nO)", arg1, arg2);

        if (args == NULL)
            return -1;

        res = f(self, args);
        Py_DECREF(args);
    }
    else
    {
        PyErr_SetNone(PyExc_NotImplementedError);
        res = -1;
    }

    return res;
}


/*
 * The metatype alloc slot.
 */
static PyObject *sipWrapperType_alloc(PyTypeObject *self, Py_ssize_t nitems)
{
    PyObject *o;

    /* Call the standard super-metatype alloc. */
    if ((o = PyType_Type.tp_alloc(self, nitems)) == NULL)
        return NULL;

    /*
     * Consume any extra type specific information and use it to initialise the
     * slots.  This only happens for directly wrapped classes (and not
     * programmer written sub-classes).  This must be done in the alloc
     * function because it is the only place we can break out of the default
     * new() function before PyType_Ready() is called.
     */
    if (currentType != NULL)
    {
        assert(!sipTypeIsEnum(currentType));

        ((sipWrapperType *)o)->wt_td = currentType;

        if (sipTypeIsClass(currentType))
        {
            const sipClassTypeDef *ctd = (const sipClassTypeDef *)currentType;
            const char *docstring = ctd->ctd_docstring;

            /*
             * Skip the marker that identifies the docstring as being
             * automatically generated.
             */
            if (docstring != NULL && *docstring == AUTO_DOCSTRING)
                ++docstring;

            ((PyTypeObject *)o)->tp_doc = docstring;

            addClassSlots((sipWrapperType *)o, ctd);

            /* Patch any mixin initialiser. */
            if (ctd->ctd_init_mixin != NULL)
                ((PyTypeObject *)o)->tp_init = ctd->ctd_init_mixin;
        }
    }

    return o;
}


/*
 * The metatype init slot.
 */
static int sipWrapperType_init(sipWrapperType *self, PyObject *args,
        PyObject *kwds)
{
    /* Call the standard super-metatype init. */
    if (PyType_Type.tp_init((PyObject *)self, args, kwds) < 0)
        return -1;

    /*
     * If we don't yet have any extra type specific information (because we are
     * a programmer defined sub-class) then get it from the (first) super-type.
     */
    if (self->wt_td == NULL)
    {
        PyTypeObject *base = ((PyTypeObject *)self)->tp_base;

        self->wt_user_type = TRUE;

        /*
         * We allow the class to use this as a meta-type without being derived
         * from a class that uses it.  This allows mixin classes that need
         * their own meta-type to work so long as their meta-type is derived
         * from this meta-type.  This condition is indicated by the pointer to
         * the generated type structure being NULL.
         */
        if (base != NULL && PyObject_TypeCheck((PyObject *)base, (PyTypeObject *)&sipWrapperType_Type))
        {
            /* TODO: Deprecate this mechanism in favour of an event handler. */
            sipNewUserTypeFunc new_user_type_handler;

            self->wt_td = ((sipWrapperType *)base)->wt_td;

            if (self->wt_td != NULL)
            {
                /* Call any new type handler. */
                new_user_type_handler = find_new_user_type_handler(
                        (sipWrapperType *)sipTypeAsPyTypeObject(self->wt_td));

                if (new_user_type_handler != NULL)
                    if (new_user_type_handler(self) < 0)
                        return -1;
            }
        }
    }
    else
    {
        /*
         * We must be a generated type so remember the type object in the
         * generated type structure.
         */
        assert(self->wt_td->td_py_type == NULL);

        self->wt_td->td_py_type = (PyTypeObject *)self;
    }

    return 0;
}


/*
 * The metatype getattro slot.
 */
static PyObject *sipWrapperType_getattro(PyObject *self, PyObject *name)
{
    if (add_all_lazy_attrs(((sipWrapperType *)self)->wt_td) < 0)
        return NULL;

    return PyType_Type.tp_getattro(self, name);
}


/*
 * The metatype setattro slot.
 */
static int sipWrapperType_setattro(PyObject *self, PyObject *name,
        PyObject *value)
{
    if (add_all_lazy_attrs(((sipWrapperType *)self)->wt_td) < 0)
        return -1;

    return PyType_Type.tp_setattro(self, name, value);
}


/*
 * The instance new slot.
 */
static PyObject *sipSimpleWrapper_new(sipWrapperType *wt, PyObject *args,
        PyObject *kwds)
{
    sipTypeDef *td = wt->wt_td;

    (void)args;
    (void)kwds;

    /* Check the base types are not being used directly. */
    if (wt == &sipSimpleWrapper_Type || wt == &sipWrapper_Type)
    {
        PyErr_Format(PyExc_TypeError,
                "the %s type cannot be instantiated or sub-classed",
                ((PyTypeObject *)wt)->tp_name);

        return NULL;
    }

    if (add_all_lazy_attrs(td) < 0)
        return NULL;

    /* See if it is a mapped type. */
    if (sipTypeIsMapped(td))
    {
        PyErr_Format(PyExc_TypeError,
                "%s.%s represents a mapped type and cannot be instantiated",
                sipNameOfModule(td->td_module),
                sipPyNameOfContainer(get_container(td), td));

        return NULL;
    }

    /* See if it is a namespace. */
    if (sipTypeIsNamespace(td))
    {
        PyErr_Format(PyExc_TypeError,
                "%s.%s represents a C++ namespace and cannot be instantiated",
                sipNameOfModule(td->td_module),
                sipPyNameOfContainer(get_container(td), td));

        return NULL;
    }

    /*
     * See if the object is being created explicitly rather than being wrapped.
     */
    if (!sipIsPending())
    {
        /*
         * See if it cannot be instantiated or sub-classed from Python, eg.
         * it's an opaque class.  Some restrictions might be overcome with
         * better SIP support.
         */
        if (((sipClassTypeDef *)td)->ctd_init == NULL)
        {
            PyErr_Format(PyExc_TypeError,
                    "%s.%s cannot be instantiated or sub-classed",
                    sipNameOfModule(td->td_module),
                    sipPyNameOfContainer(get_container(td), td));

            return NULL;
        }

        /* See if it is an abstract type. */
        if (sipTypeIsAbstract(td) && !wt->wt_user_type && ((sipClassTypeDef *)td)->ctd_init_mixin == NULL)
        {
            PyErr_Format(PyExc_TypeError,
                    "%s.%s represents a C++ abstract class and cannot be instantiated",
                    sipNameOfModule(td->td_module),
                    sipPyNameOfContainer(get_container(td), td));

            return NULL;
        }
    }

    /* Call the standard super-type new. */
    return PyBaseObject_Type.tp_new((PyTypeObject *)wt, empty_tuple, NULL);
}


/*
 * The instance init slot.
 */
static int sipSimpleWrapper_init(sipSimpleWrapper *self, PyObject *args,
        PyObject *kwds)
{
    void *sipNew;
    int sipFlags, from_cpp = TRUE;
    sipWrapper *owner;
    sipWrapperType *wt = (sipWrapperType *)Py_TYPE(self);
    sipTypeDef *td = wt->wt_td;
    sipClassTypeDef *ctd = (sipClassTypeDef *)td;
    PyObject *unused = NULL;
    sipFinalFunc final_func = find_finalisation(ctd);

    /* Check for an existing C++ instance waiting to be wrapped. */
    if (sipGetPending(&sipNew, &owner, &sipFlags) < 0)
        return -1;

    if (sipNew == NULL)
    {
        PyObject *parseErr = NULL, **unused_p = NULL;

        /* See if we are interested in any unused keyword arguments. */
        if (sipTypeCallSuperInit(&ctd->ctd_base) || final_func != NULL || kw_handler != NULL)
            unused_p = &unused;

        /* Call the C++ ctor. */
        owner = NULL;

        sipNew = ctd->ctd_init(self, args, kwds, unused_p, (PyObject **)&owner,
                &parseErr);

        if (sipNew != NULL)
        {
            sipFlags = SIP_DERIVED_CLASS;
        }
        else if (parseErr == NULL)
        {
            /*
             * The C++ ctor must have raised an exception which has been
             * translated to a Python exception.
             */
            return -1;
        }
        else
        {
            sipInitExtenderDef *ie = wt->wt_iextend;

            /*
             * If we have not found an appropriate overload then try any
             * extenders.
             */
            while (PyList_Check(parseErr) && ie != NULL)
            {
                sipNew = ie->ie_extender(self, args, kwds, &unused,
                        (PyObject **)&owner, &parseErr);

                if (sipNew != NULL)
                    break;

                ie = ie->ie_next;
            }

            if (sipNew == NULL)
            {
                const char *docstring = ctd->ctd_docstring;

                /*
                 * Use the docstring for errors if it was automatically
                 * generated.
                 */
                if (docstring != NULL)
                {
                    if (*docstring == AUTO_DOCSTRING)
                        ++docstring;
                    else
                        docstring = NULL;
                }

                sip_api_no_function(parseErr,
                        sipPyNameOfContainer(&ctd->ctd_container, td),
                        docstring);

                return -1;
            }

            sipFlags = 0;
        }

        if (owner == NULL)
            sipFlags |= SIP_PY_OWNED;
        else if ((PyObject *)owner == Py_None)
        {
            /* This is the hack that means that C++ owns the new instance. */
            sipFlags |= SIP_CPP_HAS_REF;
            Py_INCREF(self);
            owner = NULL;
        }

        /* The instance was created from Python. */
        from_cpp = FALSE;
    }

    /* Handler any owner if the type supports the concept. */
    if (PyObject_TypeCheck((PyObject *)self, (PyTypeObject *)&sipWrapper_Type))
    {
        /*
         * The application may be doing something very unadvisable (like
         * calling __init__() for a second time), so make sure we don't already
         * have a parent.
         */
        removeFromParent((sipWrapper *)self);

        if (owner != NULL)
        {
            assert(PyObject_TypeCheck((PyObject *)owner, (PyTypeObject *)&sipWrapper_Type));

            addToParent((sipWrapper *)self, (sipWrapper *)owner);
        }
    }

    self->data = sipNew;
    self->sw_flags = sipFlags | SIP_CREATED;

    /* Set the access function. */
    if (sipIsAccessFunc(self))
        self->access_func = explicit_access_func;
    else if (sipIsIndirect(self))
        self->access_func = indirect_access_func;
    else
        self->access_func = NULL;

    if (!sipNotInMap(self))
        sipOMAddObject(&cppPyMap, self);

    /* If we are wrapping an instance returned from C/C++ then we are done. */
    if (from_cpp)
    {
        /*
         * Invoke any event handlers for instances that are accessed directly.
         */
        if (self->access_func == NULL)
        {
            sipEventHandler *eh;

            for (eh = event_handlers[sipEventWrappedInstance]; eh != NULL; eh = eh->next)
            {
                if (is_subtype(ctd, eh->ctd))
                {
                    sipWrappedInstanceEventHandler handler = (sipWrappedInstanceEventHandler)eh->handler;

                    handler(sipNew);
                }
            }
        }

        return 0;
    }

    /* Call any finalisation code. */
    if (final_func != NULL)
    {
        PyObject *new_unused = NULL, **new_unused_p;

        if (unused == NULL || unused != kwds)
        {
            /*
             * There are no unused arguments or we have already created a dict
             * containing the unused sub-set, so there is no need to create
             * another.
             */
            new_unused_p = NULL;
        }
        else
        {
            /*
             * All of the keyword arguments are unused, so if some of them are
             * now going to be used then a new dict will be needed.
             */
            new_unused_p = &new_unused;
        }
            
        if (final_func((PyObject *)self, sipNew, unused, new_unused_p) < 0)
        {
            Py_XDECREF(unused);
            return -1;
        }

        if (new_unused != NULL)
        {
            Py_DECREF(unused);
            unused = new_unused;
        }
    }

    /* Call the handler if we have one. */
    if (kw_handler != NULL && unused != NULL && isQObject((PyObject *)self))
    {
        int rc = kw_handler((PyObject *)self, sipNew, unused);

        /*
         * A handler will always consume all unused keyword arguments (or raise
         * an exception) so discard the dict now.
         */
        Py_DECREF(unused);

        if (rc < 0)
            return -1;

        unused = NULL;
    }

    /* See if we should call the equivalent of super().__init__(). */
    if (sipTypeCallSuperInit(&ctd->ctd_base))
    {
        PyObject *next;

        /* Find the next type in the MRO. */
        next = next_in_mro((PyObject *)self,
                (PyObject *)&sipSimpleWrapper_Type);

        /*
         * If the next type in the MRO is object then take a shortcut by not
         * calling super().__init__() but emulating object.__init__() instead.
         * This will be the most common case and also allows us to generate a
         * better exception message if there are unused keyword arguments.  The
         * disadvantage is that the exception message will be different if
         * there is a mixin.
         */
        if (next != (PyObject *)&PyBaseObject_Type)
        {
            int rc = super_init((PyObject *)self, empty_tuple, unused, next);

            Py_XDECREF(unused);

            return rc;
        }
    }

    if (unused_backdoor != NULL)
    {
        /*
         * We are being called by a mixin's __init__ so save any unused
         * arguments for it to pass on to the main class's __init__.
         */
        *unused_backdoor = unused;
    }
    else if (unused != NULL)
    {
        /* We shouldn't have any unused keyword arguments. */
        if (PyDict_Size(unused) != 0)
        {
            PyObject *key, *value;
            Py_ssize_t pos = 0;

            /* Just report one of the unused arguments. */
            PyDict_Next(unused, &pos, &key, &value);

            PyErr_Format(PyExc_TypeError,
                    "'%S' is an unknown keyword argument", key);

            Py_DECREF(unused);

            return -1;
        }

        Py_DECREF(unused);
    }

    return 0;
}


/*
 * Get the C++ address of a mixin.
 */
static void *sip_api_get_mixin_address(sipSimpleWrapper *w,
        const sipTypeDef *td)
{
    PyObject *mixin;
    void *cpp;

    if ((mixin = PyObject_GetAttrString((PyObject *)w, sipTypeName(td))) == NULL)
    {
        PyErr_Clear();
        return NULL;
    }

    cpp = sip_api_get_address((sipSimpleWrapper *)mixin);

    Py_DECREF(mixin);

    return cpp;
}


/*
 * Initialise a mixin.
 */
static int sip_api_init_mixin(PyObject *self, PyObject *args, PyObject *kwds,
        const sipClassTypeDef *ctd)
{
    int rc;
    Py_ssize_t pos;
    PyObject *unused, *mixin, *mixin_name, *key, *value;
    PyTypeObject *self_wt = sipTypeAsPyTypeObject(((sipWrapperType *)Py_TYPE(self))->wt_td);
    PyTypeObject *wt = sipTypeAsPyTypeObject(&ctd->ctd_base);

    static PyObject *double_us = NULL;

    if (objectify("__", &double_us) < 0)
        return -1;

    /* If we are not a mixin to another wrapped class then behave as normal. */
    if (PyType_IsSubtype(self_wt, wt))
        return super_init(self, args, kwds, next_in_mro(self, (PyObject *)wt));

    /*
     * Create the mixin instance.  Retain the positional arguments for the
     * super-class.  Remember that, even though the mixin appears after the
     * main class in the MRO, it appears before sipWrapperType where the main
     * class's arguments are actually parsed.
     */
    unused = NULL;
    unused_backdoor = &unused;
    mixin = PyObject_Call((PyObject *)wt, empty_tuple, kwds);
    unused_backdoor = NULL;

    if (mixin == NULL)
        goto gc_unused;

    /* Make sure the mixin can find the main instance. */
    ((sipSimpleWrapper *)mixin)->mixin_main = self;
    Py_INCREF(self);

    if ((mixin_name = PyUnicode_FromString(sipTypeName(&ctd->ctd_base))) == NULL)
    {
        Py_DECREF(mixin);
        goto gc_unused;
    }

    rc = PyObject_SetAttr(self, mixin_name, mixin);
    Py_DECREF(mixin);

    if (rc < 0)
        goto gc_mixin_name;

    /* Add the mixin's useful attributes to the main class. */
    pos = 0;

    while (PyDict_Next(wt->tp_dict, &pos, &key, &value))
    {
        /* Don't replace existing values. */
        if (PyDict_Contains(Py_TYPE(self)->tp_dict, key) != 0)
            continue;

        /* Skip values with names that start with double underscore. */
        if (!PyUnicode_Check(key))
            continue;

        /*
         * Despite what the docs say this returns a Py_ssize_t - although the
         * docs are probably right.
         */
        rc = (int)PyUnicode_Tailmatch(key, double_us, 0, 2, -1);

        if (rc < 0)
            goto gc_mixin_name;

        if (rc > 0)
            continue;

        if (PyObject_IsInstance(value, (PyObject *)&sipMethodDescr_Type))
        {
            if ((value = sipMethodDescr_Copy(value, mixin_name)) == NULL)
                goto gc_mixin_name;
        }
        else if (PyObject_IsInstance(value, (PyObject *)&sipVariableDescr_Type))
        {
            if ((value = sipVariableDescr_Copy(value, mixin_name)) == NULL)
                goto gc_mixin_name;
        }
        else
        {
            Py_INCREF(value);
        }

        rc = PyDict_SetItem(Py_TYPE(self)->tp_dict, key, value);

        Py_DECREF(value);

        if (rc < 0)
            goto gc_mixin_name;
    }

    Py_DECREF(mixin_name);

    /* Call the super-class's __init__ with any remaining arguments. */
    rc = super_init(self, args, unused, next_in_mro(self, (PyObject *)wt));
    Py_XDECREF(unused);

    return rc;

gc_mixin_name:
    Py_DECREF(mixin_name);

gc_unused:
    Py_XDECREF(unused);

    return -1;
}


/*
 * Return the next in the MRO of an instance after a given type.
 */
static PyObject *next_in_mro(PyObject *self, PyObject *after)
{
    Py_ssize_t i;
    PyObject *mro;

    mro = Py_TYPE(self)->tp_mro;
    assert(PyTuple_Check(mro));

    for (i = 0; i < PyTuple_GET_SIZE(mro); ++i)
        if (PyTuple_GET_ITEM(mro, i) == after)
            break;

    /* Assert that we have found ourself and that we are not the last. */
    assert(i + 1 < PyTuple_GET_SIZE(mro));

    return PyTuple_GET_ITEM(mro, i + 1);
}


/*
 * Call the equivalent of super()__init__() of an instance.
 */
static int super_init(PyObject *self, PyObject *args, PyObject *kwds,
        PyObject *type)
{
    int i;
    PyObject *init, *init_args, *init_res;

    if ((init = PyObject_GetAttr(type, init_name)) == NULL)
        return -1;

    if ((init_args = PyTuple_New(1 + PyTuple_GET_SIZE(args))) == NULL)
    {
        Py_DECREF(init);
        return -1;
    }

    PyTuple_SET_ITEM(init_args, 0, self);
    Py_INCREF(self);

    for (i = 0; i < PyTuple_GET_SIZE(args); ++i)
    {
        PyObject *arg = PyTuple_GET_ITEM(args, i);

        PyTuple_SET_ITEM(init_args, 1 + i, arg);
        Py_INCREF(arg);
    }

    init_res = PyObject_Call(init, init_args, kwds);
    Py_DECREF(init_args);
    Py_DECREF(init);
    Py_XDECREF(init_res);

    return (init_res != NULL) ? 0 : -1;
}


/*
 * Find any finalisation function for a class, searching its super-classes if
 * necessary.
 */
static sipFinalFunc find_finalisation(sipClassTypeDef *ctd)
{
    sipEncodedTypeDef *sup;

    if (ctd->ctd_final != NULL)
        return ctd->ctd_final;

    if ((sup = ctd->ctd_supers) != NULL)
        do
        {
            sipClassTypeDef *sup_ctd = sipGetGeneratedClassType(sup, ctd);
            sipFinalFunc func;

            if ((func = find_finalisation(sup_ctd)) != NULL)
                return func;
        }
        while (!sup++->sc_flag);

    return NULL;
}


/*
 * Find any new user type handler function for a class, searching its
 * super-classes if necessary.
 */
static sipNewUserTypeFunc find_new_user_type_handler(sipWrapperType *wt)
{
    sipEncodedTypeDef *sup;
    sipClassTypeDef *ctd;

    if (wt->wt_new_user_type_handler != NULL)
        return wt->wt_new_user_type_handler;

    ctd = (sipClassTypeDef *)wt->wt_td;

    if ((sup = ctd->ctd_supers) != NULL)
    {
        do
        {
            sipTypeDef *sup_td = getGeneratedType(sup, ctd->ctd_base.td_module);
            sipNewUserTypeFunc func;

            wt = (sipWrapperType *)sipTypeAsPyTypeObject(sup_td);

            if ((func = find_new_user_type_handler(wt)) != NULL)
                return func;
        }
        while (!sup++->sc_flag);
    }

    return NULL;
}


/*
 * The instance traverse slot.
 */
static int sipSimpleWrapper_traverse(sipSimpleWrapper *self, visitproc visit,
        void *arg)
{
    int vret;
    void *ptr;
    const sipClassTypeDef *ctd;

    /* Call any handwritten traverse code. */
    if ((ptr = getPtrTypeDef(self, &ctd)) != NULL)
        if (ctd->ctd_traverse != NULL)
            if ((vret = ctd->ctd_traverse(ptr, visit, arg)) != 0)
                return vret;

    if (self->dict != NULL)
        if ((vret = visit(self->dict, arg)) != 0)
            return vret;

    if (self->extra_refs != NULL)
        if ((vret = visit(self->extra_refs, arg)) != 0)
            return vret;

    if (self->user != NULL)
        if ((vret = visit(self->user, arg)) != 0)
            return vret;

    if (self->mixin_main != NULL)
        if ((vret = visit(self->mixin_main, arg)) != 0)
            return vret;

    return 0;
}


/*
 * The instance clear slot.
 */
static int sipSimpleWrapper_clear(sipSimpleWrapper *self)
{
    int vret = 0;
    void *ptr;
    const sipClassTypeDef *ctd;
    PyObject *tmp;

    /* Call any handwritten clear code. */
    if ((ptr = getPtrTypeDef(self, &ctd)) != NULL)
        if (ctd->ctd_clear != NULL)
            vret = ctd->ctd_clear(ptr);

    /* Remove the instance dictionary. */
    tmp = self->dict;
    self->dict = NULL;
    Py_XDECREF(tmp);

    /* Remove any extra references dictionary. */
    tmp = self->extra_refs;
    self->extra_refs = NULL;
    Py_XDECREF(tmp);

    /* Remove any user object. */
    tmp = self->user;
    self->user = NULL;
    Py_XDECREF(tmp);

    /* Remove any mixin main. */
    tmp = self->mixin_main;
    self->mixin_main = NULL;
    Py_XDECREF(tmp);

    return vret;
}


/*
 * The instance get buffer slot.
 */
static int sipSimpleWrapper_getbuffer(sipSimpleWrapper *self, Py_buffer *buf,
        int flags)
{
    void *ptr;
    const sipClassTypeDef *ctd;

    if ((ptr = getPtrTypeDef(self, &ctd)) == NULL)
        return -1;

    if (sipTypeUseLimitedAPI(&ctd->ctd_base))
    {
        sipGetBufferFuncLimited getbuffer = (sipGetBufferFuncLimited)ctd->ctd_getbuffer;
        sipBufferDef bd;

        /*
         * Ensure all fields have a default value.  This means that extra
         * fields can be appended in the future that older handwritten code
         * doesn't know about.
         */
        memset(&bd, 0, sizeof(sipBufferDef));

        if (getbuffer((PyObject *)self, ptr, &bd) < 0)
            return -1;

        return PyBuffer_FillInfo(buf, (PyObject *)self, bd.bd_buffer,
                bd.bd_length, bd.bd_readonly, flags);
    }

    return ctd->ctd_getbuffer((PyObject *)self, ptr, buf, flags);
}


/*
 * The instance release buffer slot.
 */
static void sipSimpleWrapper_releasebuffer(sipSimpleWrapper *self,
        Py_buffer *buf)
{
    void *ptr;
    const sipClassTypeDef *ctd;

    if ((ptr = getPtrTypeDef(self, &ctd)) == NULL)
        return;

    if (sipTypeUseLimitedAPI(&ctd->ctd_base))
    {
        sipReleaseBufferFuncLimited releasebuffer = (sipReleaseBufferFuncLimited)ctd->ctd_releasebuffer;

        releasebuffer((PyObject *)self, ptr);

        return;
    }

    ctd->ctd_releasebuffer((PyObject *)self, ptr, buf);
}


/*
 * The instance dealloc slot.
 */
static void sipSimpleWrapper_dealloc(sipSimpleWrapper *self)
{
    PyObject *error_type, *error_value, *error_traceback;

    /* Save the current exception, if any. */
    PyErr_Fetch(&error_type, &error_value, &error_traceback);

    forgetObject(self);

    /*
     * Now that the C++ object no longer exists we can tidy up the Python
     * object.  We used to do this first but that meant lambda slots were
     * removed too soon (if they were connected to QObject.destroyed()).
     */
    sipSimpleWrapper_clear(self);

    /* Call the standard super-type dealloc. */
    PyBaseObject_Type.tp_dealloc((PyObject *)self);

    /* Restore the saved exception. */
    PyErr_Restore(error_type, error_value, error_traceback);
}


/*
 * The type call slot.
 */
static PyObject *slot_call(PyObject *self, PyObject *args, PyObject *kw)
{
    PyObject *(*f)(PyObject *, PyObject *, PyObject *);

    f = (PyObject *(*)(PyObject *, PyObject *, PyObject *))findSlot(self, call_slot);

    assert(f != NULL);

    return f(self, args, kw);
}


/*
 * The sequence type item slot.
 */
static PyObject *slot_sq_item(PyObject *self, Py_ssize_t n)
{
    PyObject *(*f)(PyObject *,PyObject *);
    PyObject *arg, *res;

    if ((arg = PyLong_FromSsize_t(n)) == NULL)
        return NULL;

    f = (PyObject *(*)(PyObject *,PyObject *))findSlot(self, getitem_slot);

    assert(f != NULL);

    res = f(self,arg);

    Py_DECREF(arg);

    return res;
}


/*
 * The mapping type assign subscript slot.
 */
static int slot_mp_ass_subscript(PyObject *self, PyObject *key,
        PyObject *value)
{
    return objobjargprocSlot(self, key, value,
            (value != NULL ? setitem_slot : delitem_slot));
}


/*
 * The sequence type assign item slot.
 */
static int slot_sq_ass_item(PyObject *self, Py_ssize_t i, PyObject *o)
{
    return ssizeobjargprocSlot(self, i, o,
            (o != NULL ? setitem_slot : delitem_slot));
}


/*
 * The type rich compare slot.
 */
static PyObject *slot_richcompare(PyObject *self, PyObject *arg, int op)
{
    PyObject *(*f)(PyObject *,PyObject *);
    sipPySlotType st;

    /* Convert the operation to a slot type. */
    switch (op)
    {
    case Py_LT:
        st = lt_slot;
        break;

    case Py_LE:
        st = le_slot;
        break;

    case Py_EQ:
        st = eq_slot;
        break;

    case Py_NE:
        st = ne_slot;
        break;

    case Py_GT:
        st = gt_slot;
        break;

    case Py_GE:
        st = ge_slot;
        break;
    }

    /* It might not exist if not all the above have been implemented. */
    if ((f = (PyObject *(*)(PyObject *,PyObject *))findSlot(self, st)) == NULL)
    {
        Py_INCREF(Py_NotImplemented);
        return Py_NotImplemented;
    }

    return f(self, arg);
}


/*
 * The __dict__ getter.
 */
static PyObject *sipSimpleWrapper_get_dict(sipSimpleWrapper *sw, void *closure)
{
    (void)closure;

    /* Create the dictionary if needed. */
    if (sw->dict == NULL)
    {
        sw->dict = PyDict_New();

        if (sw->dict == NULL)
            return NULL;
    }

    Py_INCREF(sw->dict);
    return sw->dict;
}


/*
 * The __dict__ setter.
 */
static int sipSimpleWrapper_set_dict(sipSimpleWrapper *sw, PyObject *value,
        void *closure)
{
    (void)closure;

    /* Check that any new value really is a dictionary. */
    if (value != NULL && !PyDict_Check(value))
    {
        PyErr_Format(PyExc_TypeError,
                "__dict__ must be set to a dictionary, not a '%s'",
                Py_TYPE(value)->tp_name);
        return -1;
    }

    Py_XDECREF(sw->dict);
    
    Py_XINCREF(value);
    sw->dict = value;

    return 0;
}


/*
 * The table of getters and setters.
 */
static PyGetSetDef sipSimpleWrapper_getset[] = {
    {(char *)"__dict__", (getter)sipSimpleWrapper_get_dict,
            (setter)sipSimpleWrapper_set_dict, NULL, NULL},
    {NULL, NULL, NULL, NULL, NULL}
};


/*
 * The type data structure.  Note that we pretend to be a mapping object and a
 * sequence object at the same time.  Python will choose one over another,
 * depending on the context, but we implement as much as we can and don't make
 * assumptions about which Python will choose.
 */
sipWrapperType sipSimpleWrapper_Type = {
#if !defined(STACKLESS)
    {
#endif
        {
            PyVarObject_HEAD_INIT(&sipWrapperType_Type, 0)
            "sip.simplewrapper",    /* tp_name */
            sizeof (sipSimpleWrapper),  /* tp_basicsize */
            0,              /* tp_itemsize */
            (destructor)sipSimpleWrapper_dealloc,   /* tp_dealloc */
            0,              /* tp_print */
            0,              /* tp_getattr */
            0,              /* tp_setattr */
            0,              /* tp_as_async (Python v3.5), tp_compare (Python v2) */
            0,              /* tp_repr */
            0,              /* tp_as_number */
            0,              /* tp_as_sequence */
            0,              /* tp_as_mapping */
            0,              /* tp_hash */
            0,              /* tp_call */
            0,              /* tp_str */
            0,              /* tp_getattro */
            0,              /* tp_setattro */
            0,              /* tp_as_buffer */
            Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC,  /* tp_flags */
            0,              /* tp_doc */
            (traverseproc)sipSimpleWrapper_traverse,    /* tp_traverse */
            (inquiry)sipSimpleWrapper_clear,    /* tp_clear */
            0,              /* tp_richcompare */
            0,              /* tp_weaklistoffset */
            0,              /* tp_iter */
            0,              /* tp_iternext */
            0,              /* tp_methods */
            0,              /* tp_members */
            sipSimpleWrapper_getset,    /* tp_getset */
            0,              /* tp_base */
            0,              /* tp_dict */
            0,              /* tp_descr_get */
            0,              /* tp_descr_set */
            offsetof(sipSimpleWrapper, dict),   /* tp_dictoffset */
            (initproc)sipSimpleWrapper_init,    /* tp_init */
            0,              /* tp_alloc */
            (newfunc)sipSimpleWrapper_new,  /* tp_new */
            0,              /* tp_free */
            0,              /* tp_is_gc */
            0,              /* tp_bases */
            0,              /* tp_mro */
            0,              /* tp_cache */
            0,              /* tp_subclasses */
            0,              /* tp_weaklist */
            0,              /* tp_del */
            0,              /* tp_version_tag */
            0,              /* tp_finalize */
#if PY_VERSION_HEX >= 0x03080000
            0,              /* tp_vectorcall */
#endif
        },
        {
            0,              /* am_await */
            0,              /* am_aiter */
            0,              /* am_anext */
        },
        {
            0,              /* nb_add */
            0,              /* nb_subtract */
            0,              /* nb_multiply */
            0,              /* nb_remainder */
            0,              /* nb_divmod */
            0,              /* nb_power */
            0,              /* nb_negative */
            0,              /* nb_positive */
            0,              /* nb_absolute */
            0,              /* nb_bool */
            0,              /* nb_invert */
            0,              /* nb_lshift */
            0,              /* nb_rshift */
            0,              /* nb_and */
            0,              /* nb_xor */
            0,              /* nb_or */
            0,              /* nb_int */
            0,              /* nb_reserved */
            0,              /* nb_float */
            0,              /* nb_inplace_add */
            0,              /* nb_inplace_subtract */
            0,              /* nb_inplace_multiply */
            0,              /* nb_inplace_remainder */
            0,              /* nb_inplace_power */
            0,              /* nb_inplace_lshift */
            0,              /* nb_inplace_rshift */
            0,              /* nb_inplace_and */
            0,              /* nb_inplace_xor */
            0,              /* nb_inplace_or */
            0,              /* nb_floor_divide */
            0,              /* nb_true_divide */
            0,              /* nb_inplace_floor_divide */
            0,              /* nb_inplace_true_divide */
            0,              /* nb_index */
            0,              /* nb_matrix_multiply */
            0,              /* nb_inplace_matrix_multiply */
        },
        {
            0,              /* mp_length */
            0,              /* mp_subscript */
            0,              /* mp_ass_subscript */
        },
        {
            0,              /* sq_length */
            0,              /* sq_concat */
            0,              /* sq_repeat */
            0,              /* sq_item */
            0,              /* was_sq_slice */
            0,              /* sq_ass_item */
            0,              /* was_sq_ass_slice */
            0,              /* sq_contains */
            0,              /* sq_inplace_concat */
            0,              /* sq_inplace_repeat */
        },
        {
            0,              /* bf_getbuffer */
            0,              /* bf_releasebuffer */
        },
        0,                  /* ht_name */
        0,                  /* ht_slots */
        0,                  /* ht_qualname */
        0,                  /* ht_cached_keys */
#if PY_VERSION_HEX >= 0x03090000
        0,                  /* ht_module */
#endif
#if !defined(STACKLESS)
    },
#endif
    0,                      /* wt_user_type */
    0,                      /* wt_dict_complete */
    0,                      /* wt_unused */
    0,                      /* wt_td */
    0,                      /* wt_iextend */
    0,                      /* wt_new_user_type_handler */
    0,                      /* wt_user_data */
};


/*
 * The wrapper clear slot.
 */
static int sipWrapper_clear(sipWrapper *self)
{
    int vret;
    sipSimpleWrapper *sw = (sipSimpleWrapper *)self;

    vret = sipSimpleWrapper_clear(sw);

    /* Remove any slots connected via a proxy. */
    if (sipQtSupport != NULL && sipPossibleProxy(sw) && !sipNotInMap(sw))
    {
        void *tx = sip_api_get_address(sw);

        if (tx != NULL)
        {
            sipSlot *slot;
            void *context = NULL;

            assert (sipQtSupport->qt_find_sipslot);

            while ((slot = sipQtSupport->qt_find_sipslot(tx, &context)) != NULL)
            {
                sip_api_clear_any_slot_reference(slot);

                if (context == NULL)
                    break;
            }
        }
    }

    /* Detach any children (which will be owned by C/C++). */
    detachChildren(self);

    return vret;
}


/*
 * The wrapper dealloc slot.
 */
static void sipWrapper_dealloc(sipWrapper *self)
{
    PyObject *error_type, *error_value, *error_traceback;

    /* Save the current exception, if any. */
    PyErr_Fetch(&error_type, &error_value, &error_traceback);

    /*
     * We can't simply call the super-type because things have to be done in a
     * certain order.  The first thing is to get rid of the wrapped instance.
     */
    forgetObject((sipSimpleWrapper *)self);

    sipWrapper_clear(self);

    /* Skip the super-type's dealloc. */
    PyBaseObject_Type.tp_dealloc((PyObject *)self);

    /* Restore the saved exception. */
    PyErr_Restore(error_type, error_value, error_traceback);
}


/*
 * The wrapper traverse slot.
 */
static int sipWrapper_traverse(sipWrapper *self, visitproc visit, void *arg)
{
    int vret;
    sipSimpleWrapper *sw = (sipSimpleWrapper *)self;
    sipWrapper *w;

    if ((vret = sipSimpleWrapper_traverse(sw, visit, arg)) != 0)
        return vret;

    /*
     * This should be handwritten code in PyQt.  The map check is a bit of a
     * hack to work around PyQt4 problems with qApp and a user created
     * instance.  qt_find_sipslot() will return the same slot information for
     * both causing the gc module to trigger assert() failures.
     */
    if (sipQtSupport != NULL && sipQtSupport->qt_find_sipslot && !sipNotInMap(sw))
    {
        void *tx = sip_api_get_address(sw);

        if (tx != NULL)
        {
            sipSlot *slot;
            void *context = NULL;

            while ((slot = sipQtSupport->qt_find_sipslot(tx, &context)) != NULL)
            {
                if ((vret = sip_api_visit_slot(slot, visit, arg)) != 0)
                    return vret;

                if (context == NULL)
                    break;
            }
        }
    }

    for (w = self->first_child; w != NULL; w = w->sibling_next)
    {
        /*
         * We don't traverse if the wrapper is a child of itself.  We do this
         * so that wrapped objects returned by virtual methods with the
         * /Factory/ don't have those objects collected.  This then means that
         * plugins implemented in Python have a chance of working.
         */
        if (w != self)
            if ((vret = visit((PyObject *)w, arg)) != 0)
                return vret;
    }

    return 0;
}


/*
 * Add the slots for a class type and all its super-types.
 */
static void addClassSlots(sipWrapperType *wt, const sipClassTypeDef *ctd)
{
    PyHeapTypeObject *heap_to = &wt->super;
    PyBufferProcs *bp = &heap_to->as_buffer;

    /* Add the buffer interface. */
    if (ctd->ctd_getbuffer != NULL)
        bp->bf_getbuffer = (getbufferproc)sipSimpleWrapper_getbuffer;

    if (ctd->ctd_releasebuffer != NULL)
        bp->bf_releasebuffer = (releasebufferproc)sipSimpleWrapper_releasebuffer;

    /* Add the slots for this type. */
    if (ctd->ctd_pyslots != NULL)
        addTypeSlots(heap_to, ctd->ctd_pyslots);
}


/*
 * Add the slot handler for each slot present in the type.
 */
static void addTypeSlots(PyHeapTypeObject *heap_to, sipPySlotDef *slots)
{
    PyTypeObject *to;
    PyNumberMethods *nb;
    PySequenceMethods *sq;
    PyMappingMethods *mp;
    PyAsyncMethods *am;
    void *f;

    to = &heap_to->ht_type;
    nb = &heap_to->as_number;
    sq = &heap_to->as_sequence;
    mp = &heap_to->as_mapping;
    am = &heap_to->as_async;

    while ((f = slots->psd_func) != NULL)
        switch (slots++->psd_type)
        {
        case str_slot:
            to->tp_str = (reprfunc)f;
            break;

        case int_slot:
            nb->nb_int = (unaryfunc)f;
            break;

        case float_slot:
            nb->nb_float = (unaryfunc)f;
            break;

        case len_slot:
            mp->mp_length = (lenfunc)f;
            sq->sq_length = (lenfunc)f;
            break;

        case contains_slot:
            sq->sq_contains = (objobjproc)f;
            break;

        case add_slot:
            nb->nb_add = (binaryfunc)f;
            break;

        case concat_slot:
            sq->sq_concat = (binaryfunc)f;
            break;

        case sub_slot:
            nb->nb_subtract = (binaryfunc)f;
            break;

        case mul_slot:
            nb->nb_multiply = (binaryfunc)f;
            break;

        case repeat_slot:
            sq->sq_repeat = (ssizeargfunc)f;
            break;

        case div_slot:
            nb->nb_true_divide = (binaryfunc)f;
            break;

        case mod_slot:
            nb->nb_remainder = (binaryfunc)f;
            break;

        case floordiv_slot:
            nb->nb_floor_divide = (binaryfunc)f;
            break;

        case truediv_slot:
            nb->nb_true_divide = (binaryfunc)f;
            break;

        case and_slot:
            nb->nb_and = (binaryfunc)f;
            break;

        case or_slot:
            nb->nb_or = (binaryfunc)f;
            break;

        case xor_slot:
            nb->nb_xor = (binaryfunc)f;
            break;

        case lshift_slot:
            nb->nb_lshift = (binaryfunc)f;
            break;

        case rshift_slot:
            nb->nb_rshift = (binaryfunc)f;
            break;

        case iadd_slot:
            nb->nb_inplace_add = (binaryfunc)f;
            break;

        case iconcat_slot:
            sq->sq_inplace_concat = (binaryfunc)f;
            break;

        case isub_slot:
            nb->nb_inplace_subtract = (binaryfunc)f;
            break;

        case imul_slot:
            nb->nb_inplace_multiply = (binaryfunc)f;
            break;

        case irepeat_slot:
            sq->sq_inplace_repeat = (ssizeargfunc)f;
            break;

        case idiv_slot:
            nb->nb_inplace_true_divide = (binaryfunc)f;
            break;

        case imod_slot:
            nb->nb_inplace_remainder = (binaryfunc)f;
            break;

        case ifloordiv_slot:
            nb->nb_inplace_floor_divide = (binaryfunc)f;
            break;

        case itruediv_slot:
            nb->nb_inplace_true_divide = (binaryfunc)f;
            break;

        case iand_slot:
            nb->nb_inplace_and = (binaryfunc)f;
            break;

        case ior_slot:
            nb->nb_inplace_or = (binaryfunc)f;
            break;

        case ixor_slot:
            nb->nb_inplace_xor = (binaryfunc)f;
            break;

        case ilshift_slot:
            nb->nb_inplace_lshift = (binaryfunc)f;
            break;

        case irshift_slot:
            nb->nb_inplace_rshift = (binaryfunc)f;
            break;

        case invert_slot:
            nb->nb_invert = (unaryfunc)f;
            break;

        case call_slot:
            to->tp_call = slot_call;
            break;

        case getitem_slot:
            mp->mp_subscript = (binaryfunc)f;
            sq->sq_item = slot_sq_item;
            break;

        case setitem_slot:
        case delitem_slot:
            mp->mp_ass_subscript = slot_mp_ass_subscript;
            sq->sq_ass_item = slot_sq_ass_item;
            break;

        case lt_slot:
        case le_slot:
        case eq_slot:
        case ne_slot:
        case gt_slot:
        case ge_slot:
            to->tp_richcompare = slot_richcompare;
            break;

        case bool_slot:
            nb->nb_bool = (inquiry)f;
            break;

        case neg_slot:
            nb->nb_negative = (unaryfunc)f;
            break;

        case repr_slot:
            to->tp_repr = (reprfunc)f;
            break;

        case hash_slot:
            to->tp_hash = (hashfunc)f;
            break;

        case pos_slot:
            nb->nb_positive = (unaryfunc)f;
            break;

        case abs_slot:
            nb->nb_absolute = (unaryfunc)f;
            break;

        case index_slot:
            nb->nb_index = (unaryfunc)f;
            break;

        case iter_slot:
            to->tp_iter = (getiterfunc)f;
            break;

        case next_slot:
            to->tp_iternext = (iternextfunc)f;
            break;

        case setattr_slot:
            to->tp_setattro = (setattrofunc)f;
            break;

        case matmul_slot:
            nb->nb_matrix_multiply = (binaryfunc)f;
            break;

        case imatmul_slot:
            nb->nb_inplace_matrix_multiply = (binaryfunc)f;
            break;

        case await_slot:
            am->am_await = (unaryfunc)f;
            break;

        case aiter_slot:
            am->am_aiter = (unaryfunc)f;
            break;

        case anext_slot:
            am->am_anext = (unaryfunc)f;
            break;

        /* Suppress a compiler warning. */
        default:
            ;
        }
}


/*
 * Remove the object from the map and call the C/C++ dtor if we own the
 * instance.
 */
static void forgetObject(sipSimpleWrapper *sw)
{
    sipEventHandler *eh;
    const sipClassTypeDef *ctd = (const sipClassTypeDef *)((sipWrapperType *)Py_TYPE(sw))->wt_td;

    /* Invoke any event handlers. */
    for (eh = event_handlers[sipEventCollectingWrapper]; eh != NULL; eh = eh->next)
    {
        if (is_subtype(ctd, eh->ctd))
        {
            sipCollectingWrapperEventHandler handler = (sipCollectingWrapperEventHandler)eh->handler;

            handler(sw);
        }
    }

    /*
     * This is needed because we might release the GIL when calling a C++ dtor.
     * Without it the cyclic garbage collector can be invoked from another
     * thread resulting in a crash.
     */
    PyObject_GC_UnTrack((PyObject *)sw);

    /*
     * Remove the object from the map before calling the class specific dealloc
     * code.  This code calls the C++ dtor and may result in further calls that
     * pass the instance as an argument.  If this is still in the map then it's
     * reference count would be increased (to one) and bad things happen when
     * it drops back to zero again.  (An example is PyQt events generated
     * during the dtor call being passed to an event filter implemented in
     * Python.)  By removing it from the map first we ensure that a new Python
     * object is created.
     */
    sipOMRemoveObject(&cppPyMap, sw);

    if (sipInterpreter != NULL || destroy_on_exit)
    {
        const sipClassTypeDef *ctd;

        if (getPtrTypeDef(sw, &ctd) != NULL && ctd->ctd_dealloc != NULL)
            ctd->ctd_dealloc(sw);
    }

    clear_access_func(sw);
}


/*
 * If the given name is that of a typedef then the corresponding type is
 * returned.
 */
static const char *sip_api_resolve_typedef(const char *name)
{
    const sipExportedModuleDef *em;

    /*
     * Note that if the same name is defined as more than one type (which is
     * possible if more than one completely independent modules are being
     * used) then we might pick the wrong one.
     */
    for (em = moduleList; em != NULL; em = em->em_next)
    {
        if (em->em_nrtypedefs > 0)
        {
            sipTypedefDef *tdd;

            tdd = (sipTypedefDef *)bsearch(name, em->em_typedefs,
                    em->em_nrtypedefs, sizeof (sipTypedefDef),
                    compareTypedefName);

            if (tdd != NULL)
                return tdd->tdd_type_name;
        }
    }

    return NULL;
}


/*
 * The bsearch() helper function for searching a sorted typedef table.
 */
static int compareTypedefName(const void *key, const void *el)
{
    return strcmp((const char *)key, ((const sipTypedefDef *)el)->tdd_name);
}


/*
 * Add the given Python object to the given list.  Return 0 if there was no
 * error.
 */
static int addPyObjectToList(sipPyObject **head, PyObject *object)
{
    sipPyObject *po;

    if ((po = sip_api_malloc(sizeof (sipPyObject))) == NULL)
        return -1;

    po->object = object;
    po->next = *head;

    *head = po;

    return 0;
}


/*
 * Register a symbol with a name.  A negative value is returned if the name was
 * already registered.
 */
static int sip_api_export_symbol(const char *name, void *sym)
{
    sipSymbol *ss;

    if (sip_api_import_symbol(name) != NULL)
        return -1;

    if ((ss = sip_api_malloc(sizeof (sipSymbol))) == NULL)
        return -1;

    ss->name = name;
    ss->symbol = sym;
    ss->next = sipSymbolList;

    sipSymbolList = ss;

    return 0;
}


/*
 * Return the symbol registered with the given name.  NULL is returned if the
 * name was not registered.
 */
static void *sip_api_import_symbol(const char *name)
{
    sipSymbol *ss;

    for (ss = sipSymbolList; ss != NULL; ss = ss->next)
        if (strcmp(ss->name, name) == 0)
            return ss->symbol;

    return NULL;
}


/*
 * Visit a slot connected to an object for the cyclic garbage collector.  This
 * would only be called externally by PyQt3.
 */
static int sip_api_visit_slot(sipSlot *slot, visitproc visit, void *arg)
{
    /* See if the slot has an extra reference. */
    if (slot->weakSlot == Py_True && slot->pyobj != Py_None)
        return visit(slot->pyobj, arg);

    return 0;
}


/*
 * Clear a slot if it has an extra reference to keep it alive.  This would only
 * be called externally by PyQt3.
 */
static void sip_api_clear_any_slot_reference(sipSlot *slot)
{
    if (slot->weakSlot == Py_True)
    {
        PyObject *xref = slot->pyobj;

        /*
         * Replace the slot with None.  We don't use NULL as this has another
         * meaning.
         */
        Py_INCREF(Py_None);
        slot->pyobj = Py_None;

        Py_DECREF(xref);
    }
}


/*
 * Convert a Python object to a character and raise an exception if there was
 * an error.
 */
static char sip_api_bytes_as_char(PyObject *obj)
{
    char ch;

    if (parseBytes_AsChar(obj, &ch) < 0)
    {
        PyErr_Format(PyExc_TypeError, "bytes of length 1 expected not '%s'",
                Py_TYPE(obj)->tp_name);

        return '\0';
    }

    return ch;
}


/*
 * Convert a Python object to a string and raise an exception if there was
 * an error.
 */
static const char *sip_api_bytes_as_string(PyObject *obj)
{
    const char *a;

    if (parseBytes_AsString(obj, &a) < 0)
    {
        PyErr_Format(PyExc_TypeError, "bytes expected not '%s'",
                Py_TYPE(obj)->tp_name);

        return NULL;
    }

    return a;
}


/*
 * Convert a Python ASCII string object to a character and raise an exception
 * if there was an error.
 */
static char sip_api_string_as_ascii_char(PyObject *obj)
{
    char ch;

    if (parseString_AsASCIIChar(obj, &ch) < 0)
        ch = '\0';

    return ch;
}


/*
 * Parse an ASCII character and return it.
 */
static int parseString_AsASCIIChar(PyObject *obj, char *ap)
{
    if (parseString_AsEncodedChar(PyUnicode_AsASCIIString(obj), obj, ap) < 0)
    {
        /* Use the exception set if it was an encoding error. */
        if (!PyUnicode_Check(obj) || PyUnicode_GET_LENGTH(obj) != 1)
            PyErr_SetString(PyExc_TypeError,
                    "bytes or ASCII string of length 1 expected");

        return -1;
    }

    return 0;
}


/*
 * Convert a Python Latin-1 string object to a character and raise an exception
 * if there was an error.
 */
static char sip_api_string_as_latin1_char(PyObject *obj)
{
    char ch;

    if (parseString_AsLatin1Char(obj, &ch) < 0)
        ch = '\0';

    return ch;
}


/*
 * Parse a Latin-1 character and return it via a pointer.
 */
static int parseString_AsLatin1Char(PyObject *obj, char *ap)
{
    if (parseString_AsEncodedChar(PyUnicode_AsLatin1String(obj), obj, ap) < 0)
    {
        /* Use the exception set if it was an encoding error. */
        if (!PyUnicode_Check(obj) || PyUnicode_GET_LENGTH(obj) != 1)
            PyErr_SetString(PyExc_TypeError,
                    "bytes or Latin-1 string of length 1 expected");

        return -1;
    }

    return 0;
}


/*
 * Convert a Python UTF-8 string object to a character and raise an exception
 * if there was an error.
 */
static char sip_api_string_as_utf8_char(PyObject *obj)
{
    char ch;

    if (parseString_AsUTF8Char(obj, &ch) < 0)
        ch = '\0';

    return ch;
}


/*
 * Parse a UTF-8 character and return it.
 */
static int parseString_AsUTF8Char(PyObject *obj, char *ap)
{
    if (parseString_AsEncodedChar(PyUnicode_AsUTF8String(obj), obj, ap) < 0)
    {
        /* Use the exception set if it was an encoding error. */
        if (!PyUnicode_Check(obj) || PyUnicode_GET_LENGTH(obj) != 1)
            PyErr_SetString(PyExc_TypeError,
                    "bytes or UTF-8 string of length 1 expected");

        return -1;
    }

    return 0;
}


/*
 * Parse an encoded character and return it.
 */
static int parseString_AsEncodedChar(PyObject *bytes, PyObject *obj, char *ap)
{
    Py_ssize_t size;

    if (bytes == NULL)
    {
        PyErr_Clear();

        return parseBytes_AsChar(obj, ap);
    }

    size = PyBytes_GET_SIZE(bytes);

    if (size != 1)
    {
        Py_DECREF(bytes);
        return -1;
    }

    if (ap != NULL)
        *ap = *PyBytes_AS_STRING(bytes);

    Py_DECREF(bytes);

    return 0;
}


/*
 * Convert a Python ASCII string object to a string and raise an exception if
 * there was an error.  The object is updated with the one that owns the
 * string.  Note that None is considered an error.
 */
static const char *sip_api_string_as_ascii_string(PyObject **obj)
{
    PyObject *s = *obj;
    const char *a;

    if (s == Py_None || (*obj = parseString_AsASCIIString(s, &a)) == NULL)
    {
        /* Use the exception set if it was an encoding error. */
        if (!PyUnicode_Check(s))
            PyErr_Format(PyExc_TypeError,
                    "bytes or ASCII string expected not '%s'",
                    Py_TYPE(s)->tp_name);

        return NULL;
    }

    return a;
}


/*
 * Parse an ASCII string and return it and a new reference to the object that
 * owns the string.
 */
static PyObject *parseString_AsASCIIString(PyObject *obj, const char **ap)
{
    return parseString_AsEncodedString(PyUnicode_AsASCIIString(obj), obj, ap);
}


/*
 * Convert a Python Latin-1 string object to a string and raise an exception if
 * there was an error.  The object is updated with the one that owns the
 * string.  Note that None is considered an error.
 */
static const char *sip_api_string_as_latin1_string(PyObject **obj)
{
    PyObject *s = *obj;
    const char *a;

    if (s == Py_None || (*obj = parseString_AsLatin1String(s, &a)) == NULL)
    {
        /* Use the exception set if it was an encoding error. */
        if (!PyUnicode_Check(s))
            PyErr_Format(PyExc_TypeError,
                    "bytes or Latin-1 string expected not '%s'",
                    Py_TYPE(s)->tp_name);

        return NULL;
    }

    return a;
}


/*
 * Parse a Latin-1 string and return it and a new reference to the object that
 * owns the string.
 */
static PyObject *parseString_AsLatin1String(PyObject *obj, const char **ap)
{
    return parseString_AsEncodedString(PyUnicode_AsLatin1String(obj), obj, ap);
}


/*
 * Convert a Python UTF-8 string object to a string and raise an exception if
 * there was an error.  The object is updated with the one that owns the
 * string.  Note that None is considered an error.
 */
static const char *sip_api_string_as_utf8_string(PyObject **obj)
{
    PyObject *s = *obj;
    const char *a;

    if (s == Py_None || (*obj = parseString_AsUTF8String(s, &a)) == NULL)
    {
        /* Use the exception set if it was an encoding error. */
        if (!PyUnicode_Check(s))
            PyErr_Format(PyExc_TypeError,
                    "bytes or UTF-8 string expected not '%s'",
                    Py_TYPE(s)->tp_name);

        return NULL;
    }

    return a;
}


/*
 * Parse a UTF-8 string and return it and a new reference to the object that
 * owns the string.
 */
static PyObject *parseString_AsUTF8String(PyObject *obj, const char **ap)
{
    return parseString_AsEncodedString(PyUnicode_AsUTF8String(obj), obj, ap);
}


/*
 * Parse an encoded string and return it and a new reference to the object that
 * owns the string.
 */
static PyObject *parseString_AsEncodedString(PyObject *bytes, PyObject *obj,
        const char **ap)
{
    if (bytes != NULL)
    {
        *ap = PyBytes_AS_STRING(bytes);

        return bytes;
    }

    /* Don't try anything else if there was an encoding error. */
    if (PyUnicode_Check(obj))
        return NULL;

    PyErr_Clear();

    if (parseBytes_AsString(obj, ap) < 0)
        return NULL;

    Py_INCREF(obj);

    return obj;
}


/*
 * Parse a character array and return it's address and length.
 */
static int parseBytes_AsCharArray(PyObject *obj, const char **ap,
        Py_ssize_t *aszp)
{
    const char *a;
    Py_ssize_t asz;

    if (obj == Py_None)
    {
        a = NULL;
        asz = 0;
    }
    else if (PyBytes_Check(obj))
    {
        a = PyBytes_AS_STRING(obj);
        asz = PyBytes_GET_SIZE(obj);
    }
    else
    {
        Py_buffer view;

        if (PyObject_GetBuffer(obj, &view, PyBUF_SIMPLE) < 0)
            return -1;

        a = view.buf;
        asz = view.len;

        PyBuffer_Release(&view);
    }

    if (ap != NULL)
        *ap = a;

    if (aszp != NULL)
        *aszp = asz;

    return 0;
}


/*
 * Parse a character and return it.
 */
static int parseBytes_AsChar(PyObject *obj, char *ap)
{
    const char *chp;
    Py_ssize_t sz;

    if (PyBytes_Check(obj))
    {
        chp = PyBytes_AS_STRING(obj);
        sz = PyBytes_GET_SIZE(obj);
    }
    else
    {
        Py_buffer view;

        if (PyObject_GetBuffer(obj, &view, PyBUF_SIMPLE) < 0)
            return -1;

        chp = view.buf;
        sz = view.len;

        PyBuffer_Release(&view);
    }

    if (sz != 1)
        return -1;

    if (ap != NULL)
        *ap = *chp;

    return 0;
}


/*
 * Parse a character string and return it.
 */
static int parseBytes_AsString(PyObject *obj, const char **ap)
{
    const char *a;
    Py_ssize_t sz;

    if (parseBytes_AsCharArray(obj, &a, &sz) < 0)
        return -1;

    if (ap != NULL)
        *ap = a;

    return 0;
}


#if defined(HAVE_WCHAR_H)
/*
 * Convert a Python object to a wide character.
 */
static wchar_t sip_api_unicode_as_wchar(PyObject *obj)
{
    wchar_t ch;

    if (parseWChar(obj, &ch) < 0)
    {
        PyErr_Format(PyExc_ValueError,
                "string of length 1 expected, not %s", Py_TYPE(obj)->tp_name);

        return L'\0';
    }

    return ch;
}


/*
 * Convert a Python object to a wide character string on the heap.
 */
static wchar_t *sip_api_unicode_as_wstring(PyObject *obj)
{
    wchar_t *p;

    if (parseWCharString(obj, &p) < 0)
    {
        PyErr_Format(PyExc_ValueError,
                "string expected, not %s", Py_TYPE(obj)->tp_name);

        return NULL;
    }

    return p;
}


/*
 * Parse a wide character array and return it's address and length.
 */
static int parseWCharArray(PyObject *obj, wchar_t **ap, Py_ssize_t *aszp)
{
    wchar_t *a;
    Py_ssize_t asz;

    if (obj == Py_None)
    {
        a = NULL;
        asz = 0;
    }
    else if (PyUnicode_Check(obj))
    {
        if (convertToWCharArray(obj, &a, &asz) < 0)
            return -1;
    }
    else
    {
        return -1;
    }

    if (ap != NULL)
        *ap = a;

    if (aszp != NULL)
        *aszp = asz;

    return 0;
}


/*
 * Convert a Unicode object to a wide character array and return it's address
 * and length.
 */
static int convertToWCharArray(PyObject *obj, wchar_t **ap, Py_ssize_t *aszp)
{
    Py_ssize_t ulen;
    wchar_t *wc;

    ulen = PyUnicode_GET_LENGTH(obj);

    if ((wc = sip_api_malloc(ulen * sizeof (wchar_t))) == NULL)
        return -1;

    if ((ulen = PyUnicode_AsWideChar(obj, wc, ulen)) < 0)
    {
        sip_api_free(wc);
        return -1;
    }

    *ap = wc;
    *aszp = ulen;

    return 0;
}


/*
 * Parse a wide character and return it.
 */
static int parseWChar(PyObject *obj, wchar_t *ap)
{
    wchar_t a;

    if (PyUnicode_Check(obj))
    {
        if (convertToWChar(obj, &a) < 0)
            return -1;
    }
    else
    {
        return -1;
    }

    if (ap != NULL)
        *ap = a;

    return 0;
}


/*
 * Convert a Unicode object to a wide character and return it.
 */
static int convertToWChar(PyObject *obj, wchar_t *ap)
{
    if (PyUnicode_GET_LENGTH(obj) != 1)
        return -1;

    if (PyUnicode_AsWideChar(obj, ap, 1) != 1)
        return -1;

    return 0;
}


/*
 * Parse a wide character string and return a copy on the heap.
 */
static int parseWCharString(PyObject *obj, wchar_t **ap)
{
    wchar_t *a;

    if (obj == Py_None)
    {
        a = NULL;
    }
    else if (PyUnicode_Check(obj))
    {
        if (convertToWCharString(obj, &a) < 0)
            return -1;
    }
    else
    {
        return -1;
    }

    if (ap != NULL)
        *ap = a;

    return 0;
}


/*
 * Convert a Unicode object to a wide character string and return a copy on
 * the heap.
 */
static int convertToWCharString(PyObject *obj, wchar_t **ap)
{
    Py_ssize_t ulen;
    wchar_t *wc;

    ulen = PyUnicode_GET_LENGTH(obj);

    if ((wc = sip_api_malloc((ulen + 1) * sizeof (wchar_t))) == NULL)
        return -1;

    if ((ulen = PyUnicode_AsWideChar(obj, wc, ulen)) < 0)
    {
        sip_api_free(wc);
        return -1;
    }

    wc[ulen] = L'\0';

    *ap = wc;

    return 0;
}

#else

/*
 * Convert a Python object to a wide character.
 */
static int sip_api_unicode_as_wchar(PyObject *obj)
{
    raiseNoWChar();

    return 0;
}


/*
 * Convert a Python object to a wide character.
 */
static int *sip_api_unicode_as_wstring(PyObject *obj)
{
    raiseNoWChar();

    return NULL;
}


/*
 * Report the need for absent wide character support.
 */
static void raiseNoWChar()
{
    PyErr_SetString(PyExc_SystemError, "sip built without wchar_t support");
}

#endif


/*
 * The enum type alloc slot.
 */
static PyObject *sipEnumType_alloc(PyTypeObject *self, Py_ssize_t nitems)
{
    sipEnumTypeObject *py_type;
    sipPySlotDef *psd;

    if (currentType == NULL)
    {
        PyErr_SetString(PyExc_TypeError, "enums cannot be sub-classed");
        return NULL;
    }

    assert(sipTypeIsEnum(currentType));

    /* Call the standard super-metatype alloc. */
    if ((py_type = (sipEnumTypeObject *)PyType_Type.tp_alloc(self, nitems)) == NULL)
        return NULL;

    /*
     * Set the links between the Python type object and the generated type
     * structure.  Strictly speaking this doesn't need to be done here.
     */
    py_type->type = currentType;
    currentType->td_py_type = (PyTypeObject *)py_type;

    /*
     * Initialise any slots.  This must be done here, after the type is
     * allocated but before PyType_Ready() is called.
     */
    if ((psd = ((sipEnumTypeDef *)currentType)->etd_pyslots) != NULL)
        addTypeSlots(&py_type->super, psd);

    return (PyObject *)py_type;
}


/*
 * The enum type getattro slot.
 */
static PyObject *sipEnumType_getattro(PyObject *self, PyObject *name)
{
    PyObject *res;
    sipEnumTypeDef *etd;
    sipExportedModuleDef *client;
    const sipEnumMemberDef *enm, *emd;
    int enum_nr, nr_members, m;
    const char *name_str;

    /*
     * Try a generic lookup first.  This has the side effect of checking the
     * type of the name object.
     */
    if ((res = PyObject_GenericGetAttr(self, name)) != NULL)
        return res;

    if (!PyErr_ExceptionMatches(PyExc_AttributeError))
        return NULL;

    PyErr_Clear();

    /* Get the member name. */
    if ((name_str = PyUnicode_AsUTF8(name)) == NULL)
        return NULL;

    etd = (sipEnumTypeDef *)((sipEnumTypeObject *)self)->type;
    client = ((sipTypeDef *)etd)->td_module;

    /* Find the number of this enum. */
    for (enum_nr = 0; enum_nr < client->em_nrtypes; ++enum_nr)
        if (client->em_types[enum_nr] == (sipTypeDef *)etd)
            break;

    /* Get the enum members in the same scope. */
    if (etd->etd_scope < 0)
    {
        nr_members = client->em_nrenummembers;
        enm = client->em_enummembers;
    }
    else
    {
        const sipContainerDef *cod = get_container(client->em_types[etd->etd_scope]);

        nr_members = cod->cod_nrenummembers;
        enm = cod->cod_enummembers;
    }

    /* Find the enum member. */
    for (emd = enm, m = 0; m < nr_members; ++m, ++emd)
        if (emd->em_enum == enum_nr && strcmp(emd->em_name, name_str) == 0)
            return sip_api_convert_from_enum(emd->em_val, (sipTypeDef *)etd);

    PyErr_Format(PyExc_AttributeError,
            "sip.enumtype object '%s' has no member '%s'",
            sipPyNameOfEnum(etd), name_str);

    return NULL;
}


/*
 * Check if an object is of the right type to convert to an encoded string.
 */
static int check_encoded_string(PyObject *obj)
{
    Py_buffer view;

    if (obj == Py_None)
        return 0;

    if (PyUnicode_Check(obj))
        return 0;

    if (PyBytes_Check(obj))
        return 0;

    if (PyObject_GetBuffer(obj, &view, PyBUF_SIMPLE) < 0)
    {
        PyErr_Clear();
    }
    else
    {
        PyBuffer_Release(&view);
        return 0;
    }

    return -1;
}


/*
 * This is called by the atexit module.
 */
static PyObject *sip_exit(PyObject *self, PyObject *args)
{
    (void)self;
    (void)args;

    /* Disable all Python reimplementations of virtuals. */
    sipInterpreter = NULL;

    Py_INCREF(Py_None);
    return Py_None;
}


/*
 * Register an exit notifier with the atexit module.
 */
static int sip_api_register_exit_notifier(PyMethodDef *md)
{
    static PyObject *register_func = NULL;
    PyObject *notifier, *res;

    if (register_func == NULL && (register_func = import_module_attr("atexit", "register")) == NULL)
        return -1;

    if ((notifier = PyCFunction_New(md, NULL)) == NULL)
        return -1;

    res = PyObject_CallFunctionObjArgs(register_func, notifier, NULL);

    Py_DECREF(notifier);

    if (res == NULL)
        return -1;

    Py_DECREF(res);

    return 0;
}


/*
 * Return the function that converts a C++ instance to a Python object.
 */
static sipConvertFromFunc get_from_convertor(const sipTypeDef *td)
{
    if (sipTypeIsMapped(td))
        return ((const sipMappedTypeDef *)td)->mtd_cfrom;

    assert(sipTypeIsClass(td));

    if (autoconversion_disabled(td) != NULL)
        return NULL;

    return ((const sipClassTypeDef *)td)->ctd_cfrom;
}


/*
 * Enable or disable the auto-conversion.  Returns the previous enabled state
 * or -1 on error.
 */
static int sip_api_enable_autoconversion(const sipTypeDef *td, int enable)
{
    sipPyObject **pop;

    assert(sipTypeIsClass(td));

    pop = autoconversion_disabled(td);

    /* See if there is anything to do. */
    if (pop == NULL && enable)
        return TRUE;

    if (pop != NULL && !enable)
        return FALSE;

    if (pop != NULL)
    {
        /* Remove it from the list. */
        sipPyObject *po = *pop;

        *pop = po->next;
        sip_api_free(po);
    }
    else
    {
        /* Add it to the list. */
        if (addPyObjectToList(&sipDisabledAutoconversions, (PyObject *)sipTypeAsPyTypeObject(td)) < 0)
            return -1;
    }

    return !enable;
}


/*
 * Return a pointer to the entry in the list of disabled auto-conversions for a
 * type.
 */
static sipPyObject **autoconversion_disabled(const sipTypeDef *td)
{
    PyObject *type = (PyObject *)sipTypeAsPyTypeObject(td);
    sipPyObject **pop;

    for (pop = &sipDisabledAutoconversions; *pop != NULL; pop = &(*pop)->next)
        if ((*pop)->object == type)
            return pop;

    return NULL;
}


/*
 * Enable or disable auto-conversion of a class that supports it.
 */
static PyObject *enableAutoconversion(PyObject *self, PyObject *args)
{
    sipWrapperType *wt;
    int enable;

    (void)self;

    if (PyArg_ParseTuple(args, "O!i:enableautoconversion", &sipWrapperType_Type, &wt, &enable))
    {
        sipTypeDef *td = wt->wt_td;
        int was_enabled;
        PyObject *res;

        if (!sipTypeIsClass(td) || ((sipClassTypeDef *)td)->ctd_cfrom == NULL)
        {
            PyErr_Format(PyExc_TypeError,
                    "%s is not a wrapped class that supports optional auto-conversion", ((PyTypeObject *)wt)->tp_name);

            return NULL;
        }

        if ((was_enabled = sip_api_enable_autoconversion(td, enable)) < 0)
            return NULL;

        res = (was_enabled ? Py_True : Py_False);

        Py_INCREF(res);
        return res;
    }

    return NULL;
}


/*
 * Python copies the nb_inplace_add slot to the sq_inplace_concat slot and vice
 * versa if either are missing.  This is a bug because they don't have the same
 * API.  We therefore reverse this.
 */
static void fix_slots(PyTypeObject *py_type, sipPySlotDef *psd)
{
    while (psd->psd_func != NULL)
    {
        if (psd->psd_type == iadd_slot && py_type->tp_as_sequence != NULL)
            py_type->tp_as_sequence->sq_inplace_concat = NULL;

        if (psd->psd_type == iconcat_slot && py_type->tp_as_number != NULL)
            py_type->tp_as_number->nb_inplace_add = NULL;

        ++psd;
    }
}


/*
 * Return the main instance for an object if it is a mixin.
 */
static sipSimpleWrapper *deref_mixin(sipSimpleWrapper *w)
{
    return w->mixin_main != NULL ? (sipSimpleWrapper *)w->mixin_main : w;
}


/*
 * Convert a new C/C++ pointer to a Python instance.
 */
static PyObject *wrap_simple_instance(void *cpp, const sipTypeDef *td,
        sipWrapper *owner, int flags)
{
    return sipWrapInstance(cpp, sipTypeAsPyTypeObject(td), empty_tuple, owner,
            flags);
}


/*
 * Resolve a proxy, if applicable.
 */
static void *resolve_proxy(const sipTypeDef *td, void *proxy)
{
    sipProxyResolver *pr;

    /* TODO: Deprecate this mechanism in favour of an event handler. */
    for (pr = proxyResolvers; pr != NULL; pr = pr->next)
        if (pr->td == td)
            proxy = pr->resolver(proxy);

    return proxy;
}


/*
 * Clear a simple wrapper.
 */
static void clear_wrapper(sipSimpleWrapper *sw)
{
    if (PyObject_TypeCheck((PyObject *)sw, (PyTypeObject *)&sipWrapper_Type))
        removeFromParent((sipWrapper *)sw);

    /*
     * Transfer ownership to C++ so we don't try to release it when the
     * Python object is garbage collected.
     */
    sipResetPyOwned(sw);

    sipOMRemoveObject(&cppPyMap, sw);

    clear_access_func(sw);
}


/*
 * Set the handler to invoke when a new user Python sub-class is defined and
 * return the old handler.
 */
static sipNewUserTypeFunc sip_api_set_new_user_type_handler(
        const sipTypeDef *td, sipNewUserTypeFunc handler)
{
    sipWrapperType *wt = (sipWrapperType *)sipTypeAsPyTypeObject(td);
    sipNewUserTypeFunc old_handler = wt->wt_new_user_type_handler;;

    wt->wt_new_user_type_handler = handler;

    return old_handler;
}


/*
 * Set the user-specific type data.
 */
static void sip_api_set_type_user_data(sipWrapperType *wt, void *data)
{
    wt->wt_user_data = data;
}


/*
 * Get the user-specific type data.
 */
static void *sip_api_get_type_user_data(const sipWrapperType *wt)
{
    return wt->wt_user_data;
}


/*
 * Get the dict of a Python type (on behalf of the limited API).
 */
static PyObject *sip_api_py_type_dict(const PyTypeObject *py_type)
{
#if PY_VERSION_HEX >= 0x030c0000
    return PyType_GetDict(py_type);
#else
    return py_type->tp_dict;
#endif
}


/*
 * Get the name of a Python type (on behalf of the limited API).
 */
static const char *sip_api_py_type_name(const PyTypeObject *py_type)
{
    return py_type->tp_name;
}


/*
 * Check an object is a method and return TRUE and its component parts if it
 * is.
 */
static int sip_api_get_method(PyObject *obj, sipMethodDef *method)
{
    if (!PyMethod_Check(obj))
        return FALSE;

    if (method != NULL)
    {
        method->pm_self = PyMethod_GET_SELF(obj);
        method->pm_function = PyMethod_GET_FUNCTION(obj);
    }

    return TRUE;
}


/*
 * Create a method from its component parts.
 */
static PyObject *sip_api_from_method(const sipMethodDef *method)
{
    return PyMethod_New(method->pm_function, method->pm_self);
}


/*
 * Check an object is a C function and return TRUE and its component parts if
 * it is.
 */
static int sip_api_get_c_function(PyObject *obj, sipCFunctionDef *c_function)
{
    if (!PyCFunction_Check(obj))
        return FALSE;

    if (c_function != NULL)
    {
        c_function->cf_function = ((PyCFunctionObject *)obj)->m_ml;
        c_function->cf_self = PyCFunction_GET_SELF(obj);
    }

    return TRUE;
}


/*
 * Check an object is a date and return TRUE and its component parts if it is.
 */
static int sip_api_get_date(PyObject *obj, sipDateDef *date)
{
    if (!PyDateTimeAPI)
        PyDateTime_IMPORT;

    if (!PyDate_Check(obj))
        return FALSE;

    if (date != NULL)
    {
        date->pd_year = PyDateTime_GET_YEAR(obj);
        date->pd_month = PyDateTime_GET_MONTH(obj);
        date->pd_day = PyDateTime_GET_DAY(obj);
    }

    return TRUE;
}


/*
 * Create a date from its component parts.
 */
static PyObject *sip_api_from_date(const sipDateDef *date)
{
    if (!PyDateTimeAPI)
        PyDateTime_IMPORT;

    return PyDate_FromDate(date->pd_year, date->pd_month, date->pd_day);
}


/*
 * Check an object is a datetime and return TRUE and its component parts if it
 * is.
 */
static int sip_api_get_datetime(PyObject *obj, sipDateDef *date,
        sipTimeDef *time)
{
    if (!PyDateTimeAPI)
        PyDateTime_IMPORT;

    if (!PyDateTime_Check(obj))
        return FALSE;

    if (date != NULL)
    {
        date->pd_year = PyDateTime_GET_YEAR(obj);
        date->pd_month = PyDateTime_GET_MONTH(obj);
        date->pd_day = PyDateTime_GET_DAY(obj);
    }

    if (time != NULL)
    {
        time->pt_hour = PyDateTime_DATE_GET_HOUR(obj);
        time->pt_minute = PyDateTime_DATE_GET_MINUTE(obj);
        time->pt_second = PyDateTime_DATE_GET_SECOND(obj);
        time->pt_microsecond = PyDateTime_DATE_GET_MICROSECOND(obj);
    }

    return TRUE;
}


/*
 * Create a datetime from its component parts.
 */
static PyObject *sip_api_from_datetime(const sipDateDef *date,
        const sipTimeDef *time)
{
    if (!PyDateTimeAPI)
        PyDateTime_IMPORT;

    return PyDateTime_FromDateAndTime(date->pd_year, date->pd_month,
            date->pd_day, time->pt_hour, time->pt_minute, time->pt_second,
            time->pt_microsecond);
}


/*
 * Check an object is a time and return TRUE and its component parts if it is.
 */
static int sip_api_get_time(PyObject *obj, sipTimeDef *time)
{
    if (!PyDateTimeAPI)
        PyDateTime_IMPORT;

    if (!PyTime_Check(obj))
        return FALSE;

    if (time != NULL)
    {
        time->pt_hour = PyDateTime_TIME_GET_HOUR(obj);
        time->pt_minute = PyDateTime_TIME_GET_MINUTE(obj);
        time->pt_second = PyDateTime_TIME_GET_SECOND(obj);
        time->pt_microsecond = PyDateTime_TIME_GET_MICROSECOND(obj);
    }

    return TRUE;
}


/*
 * Create a time from its component parts.
 */
static PyObject *sip_api_from_time(const sipTimeDef *time)
{
    if (!PyDateTimeAPI)
        PyDateTime_IMPORT;

    return PyTime_FromTime(time->pt_hour, time->pt_minute, time->pt_second,
            time->pt_microsecond);
}


/*
 * See if a type is user defined.
 */
static int sip_api_is_user_type(const sipWrapperType *wt)
{
    return wt->wt_user_type;
}


/*
 * Return a frame from the execution stack.  Note that we use 'struct _frame'
 * rather than PyFrameObject because the latter wasn't exposed to the limited
 * API until Python v3.9.
 */
static struct _frame *sip_api_get_frame(int depth)
{
#if defined(PYPY_VERSION)
    /* PyPy only supports a depth of 0. */
    return NULL;
#else
    struct _frame *frame = PyEval_GetFrame();

    while (frame != NULL && depth > 0)
    {
#if PY_VERSION_HEX < 0x03090000
        frame = frame->f_back;
#else
        frame = PyFrame_GetBack(frame);

        /* Historically we return a borrowed reference. */
        Py_XDECREF(frame);
#endif
        --depth;
    }

    return frame;
#endif
}


/*
 * Check if a type was generated using the given plugin.  Note that, although
 * this is part of the public API it is undocumented on purpose.
 */
static int sip_api_check_plugin_for_type(const sipTypeDef *td,
        const char *name)
{
    /*
     * The current thinking on plugins is that SIP v5 will look for a plugin
     * with a name derived from the name as the current module in the same
     * directory as the .sip defining the module (ie. no %Plugin directive).  A
     * module hierachy may have multiple plugins but they must co-operate.  If
     * a plugin generates user data then it should include a void* (and a
     * run-time API) so that other plugins can extend it further.  This
     * approach means that a plugin's user data structure can be opaque.
     */

    sipExportedModuleDef *em = td->td_module;
    sipImportedModuleDef *im;

    if (strcmp(sipNameOfModule(em), name) == 0)
        return TRUE;

    if ((im = em->em_imports) == NULL)
        return FALSE;

    while (im->im_name != NULL)
    {
        if (strcmp(im->im_name, name) == 0)
            return TRUE;

        ++im;
    }

    return FALSE;
}


/*
 * Create a new Unicode object and return the character size and buffer.
 */
static PyObject *sip_api_unicode_new(Py_ssize_t len, unsigned maxchar,
        int *kind, void **data)
{
    PyObject *obj;

    if ((obj = PyUnicode_New(len, maxchar)) != NULL)
    {
        *kind = PyUnicode_KIND(obj);
        *data = PyUnicode_DATA(obj);
    }

    return obj;
}


/*
 * Update a new Unicode object with a new character.
 */
static void sip_api_unicode_write(int kind, void *data, int index,
        unsigned value)
{
    PyUnicode_WRITE(kind, data, index, value);
}


/*
 * Get the address of the contents of a Unicode object, the character size and
 * the length.
 */
static void *sip_api_unicode_data(PyObject *obj, int *char_size,
        Py_ssize_t *len)
{
    void *data;

    /* Assume there will be an error. */
    *char_size = -1;

    if (PyUnicode_READY(obj) < 0)
        return NULL;

    *len = PyUnicode_GET_LENGTH(obj);

    switch (PyUnicode_KIND(obj))
    {
    case PyUnicode_1BYTE_KIND:
        *char_size = 1;
        data = PyUnicode_1BYTE_DATA(obj);
        break;

    case PyUnicode_2BYTE_KIND:
        *char_size = 2;
        data = PyUnicode_2BYTE_DATA(obj);
        break;

    case PyUnicode_4BYTE_KIND:
        *char_size = 4;
        data = PyUnicode_4BYTE_DATA(obj);
        break;

    default:
        data = NULL;
    }

    return data;
}


/*
 * Get the buffer information supplied by an object that supports the buffer
 * protocol.
 */
static int sip_api_get_buffer_info(PyObject *obj, sipBufferInfoDef *bi)
{
    int rc;
    Py_buffer *buffer;

    if (!PyObject_CheckBuffer(obj))
        return 0;

    if (bi == NULL)
        return 1;

    if ((bi->bi_internal = sip_api_malloc(sizeof (Py_buffer))) == NULL)
        return -1;

    buffer = (Py_buffer *)bi->bi_internal;

    if (PyObject_GetBuffer(obj, buffer, PyBUF_FORMAT) < 0)
        return -1;

    if (buffer->ndim == 1)
    {
        bi->bi_buf = buffer->buf;
        bi->bi_obj = buffer->obj;
        bi->bi_len = buffer->len;
        bi->bi_format = buffer->format;

        rc = 1;
    }
    else
    {
        PyErr_SetString(PyExc_TypeError, "a 1-dimensional buffer is required");
        PyBuffer_Release(buffer);
        rc = -1;
    }

    return rc;
}


/*
 * Release the buffer information obtained from a previous call to
 * sipGetBufferInfo().
 */
static void sip_api_release_buffer_info(sipBufferInfoDef *bi)
{
    if (bi->bi_internal != NULL)
    {
        PyBuffer_Release((Py_buffer *)bi->bi_internal);
        sip_api_free(bi->bi_internal);
        bi->bi_internal = NULL;
    }
}


/*
 * Import all the required types from an imported module.
 */
static int importTypes(sipExportedModuleDef *client, sipImportedModuleDef *im,
        sipExportedModuleDef *em)
{
    const char *name;
    int i, e;

    /*
     * Look for each required type in turn.  Both tables are sorted so a single
     * pass will find them all.
     */
    for (i = e = 0; (name = im->im_imported_types[i].it_name) != NULL; ++i)
    {
        sipTypeDef *td = NULL;

        do
        {
            sipTypeDef *e_td;

            if (e >= em->em_nrtypes)
            {
                PyErr_Format(PyExc_RuntimeError,
                        "%s cannot import type '%s' from %s",
                        sipNameOfModule(client), name, sipNameOfModule(em));

                return -1;
            }

            e_td = em->em_types[e++];

            /* Ignore unresolved external types. */
            if (e_td != NULL && strcmp(name, sipTypeName(e_td)) == 0)
                td = e_td;
        }
        while (td == NULL);

        im->im_imported_types[i].it_td = td;
    }

    return 0;
}


/*
 * Import all the required virtual error handlers from an imported module.
 */
static int importErrorHandlers(sipExportedModuleDef *client,
        sipImportedModuleDef *im, sipExportedModuleDef *em)
{
    const char *name;
    int i;

    for (i = 0; (name = im->im_imported_veh[i].iveh_name) != NULL; ++i)
    {
        sipVirtErrorHandlerDef *veh = em->em_virterrorhandlers;
        sipVirtErrorHandlerFunc handler = NULL;

        if (veh != NULL)
        {
            while (veh->veh_name != NULL)
            {
                if (strcmp(veh->veh_name, name) == 0)
                {
                    handler = veh->veh_handler;
                    break;
                }

                ++veh;
            }
        }

        if (handler == NULL)
        {
            PyErr_Format(PyExc_RuntimeError,
                    "%s cannot import virtual error handler '%s' from %s",
                    sipNameOfModule(client), name, sipNameOfModule(em));

            return -1;
        }

        im->im_imported_veh[i].iveh_handler = handler;
    }

    return 0;
}


/*
 * Import all the required exceptions from an imported module.
 */
static int importExceptions(sipExportedModuleDef *client,
        sipImportedModuleDef *im, sipExportedModuleDef *em)
{
    const char *name;
    int i;

    for (i = 0; (name = im->im_imported_exceptions[i].iexc_name) != NULL; ++i)
    {
        PyObject **exc = em->em_exceptions;
        PyObject *exception = NULL;

        if (exc != NULL)
        {
            while (*exc != NULL)
            {
                if (strcmp(((PyTypeObject *)(*exc))->tp_name, name) == 0)
                {
                    exception = *exc;
                    break;
                }

                ++exc;
            }
        }

        if (exception == NULL)
        {
            PyErr_Format(PyExc_RuntimeError,
                    "%s cannot import exception '%s' from %s",
                    sipNameOfModule(client), name, sipNameOfModule(em));

            return -1;
        }

        im->im_imported_exceptions[i].iexc_object = exception;
    }

    return 0;
}


/*
 * Enable or disable the garbage collector.  Return the previous state or -1 if
 * there was an error.
 */
static int sip_api_enable_gc(int enable)
{
    static PyObject *enable_func = NULL, *disable_func, *isenabled_func;
    PyObject *result;
    int was_enabled;

    /*
     * This may be -ve in the highly unusual event that a previous call failed.
     */
    if (enable < 0)
        return -1;

    /* Get the functions if we haven't already got them. */
    if (enable_func == NULL)
    {
        PyObject *gc_module;

        if ((gc_module = PyImport_ImportModule("gc")) == NULL)
            return -1;

        if ((enable_func = PyObject_GetAttrString(gc_module, "enable")) == NULL)
        {
            Py_DECREF(gc_module);
            return -1;
        }

        if ((disable_func = PyObject_GetAttrString(gc_module, "disable")) == NULL)
        {
            Py_DECREF(enable_func);
            Py_DECREF(gc_module);
            return -1;
        }

        if ((isenabled_func = PyObject_GetAttrString(gc_module, "isenabled")) == NULL)
        {
            Py_DECREF(disable_func);
            Py_DECREF(enable_func);
            Py_DECREF(gc_module);
            return -1;
        }

        Py_DECREF(gc_module);
    }

    /* Get the current state. */
    if ((result = PyObject_Call(isenabled_func, empty_tuple, NULL)) == NULL)
        return -1;

    was_enabled = PyObject_IsTrue(result);
    Py_DECREF(result);

    if (was_enabled < 0)
        return -1;

    /* See if the state needs changing. */
    if (!was_enabled != !enable)
    {
        /* Enable or disable as required. */
        result = PyObject_Call((enable ? enable_func : disable_func),
                empty_tuple, NULL);

        Py_XDECREF(result);

        if (result != Py_None)
            return -1;
    }

    return was_enabled;
}


/*
 * A thin wrapper around PyObject_Print() usually used when debugging with the
 * limited API.
 */
static void sip_api_print_object(PyObject *o)
{
    PyObject_Print(o, stdout, 0);
}


/*
 * Register a handler for a particular event.
 */
static int sip_api_register_event_handler(sipEventType type,
        const sipTypeDef *td, void *handler)
{
    sipEventHandler *eh;

    assert(sipTypeIsClass(td));

    if ((eh = sip_api_malloc(sizeof (sipEventHandler))) == NULL)
        return -1;

    eh->ctd = (const sipClassTypeDef *)td;
    eh->handler = handler;

    eh->next = event_handlers[(int)type];
    event_handlers[(int)type] = eh;

    return 0;
}


/*
 * Returns TRUE if a generated class type is a sub-class of a base generated
 * class type.
 */
static int is_subtype(const sipClassTypeDef *ctd,
        const sipClassTypeDef *base_ctd)
{
    const sipEncodedTypeDef *sup;

    /* Handle the trivial cases. */
    if (ctd == base_ctd)
        return TRUE;

    if ((sup = ctd->ctd_supers) == NULL)
        return FALSE;

    /* Search the super-types. */
    do
    {
        const sipClassTypeDef *sup_ctd = sipGetGeneratedClassType(sup, ctd);

        if (is_subtype(sup_ctd, base_ctd))
            return TRUE;
    }
    while (!sup++->sc_flag);

    return FALSE;
}


/*
 * Return an attribute of an imported module.
 */
static PyObject *import_module_attr(const char *module, const char *attr)
{
    PyObject *mod_obj, *attr_obj;

    if ((mod_obj = PyImport_ImportModule(module)) == NULL)
        return NULL;

    attr_obj = PyObject_GetAttrString(mod_obj, attr);

    Py_DECREF(mod_obj);

    return attr_obj;
}


/*
 * Get the container for a generated type.
 */
static const sipContainerDef *get_container(const sipTypeDef *td)
{
    if (sipTypeIsMapped(td))
        return &((const sipMappedTypeDef *)td)->mtd_container;

    return &((const sipClassTypeDef *)td)->ctd_container;
}


/*
 * Get the __qualname__ of an object based on its enclosing scope.
 */
static PyObject *get_qualname(const sipTypeDef *td, PyObject *name)
{
    PyTypeObject *scope_type;

    /* Get the type that is the scope. */
    scope_type = sipTypeAsPyTypeObject(td);

    return PyUnicode_FromFormat("%U.%U",
            ((PyHeapTypeObject *)scope_type)->ht_qualname, name);
}


/*
 * Implement PySlice_GetIndicesEx() (or its subsequent replacement).
 */
int sip_api_convert_from_slice_object(PyObject *slice, Py_ssize_t length,
        Py_ssize_t *start, Py_ssize_t *stop, Py_ssize_t *step,
        Py_ssize_t *slicelength)
{
    if (PySlice_Unpack(slice, start, stop, step) < 0)
        return -1;

    *slicelength = PySlice_AdjustIndices(length, start, stop, *step);

    return 0;
}


/*
 * Call a visitor function for every wrapped object.
 */
static void sip_api_visit_wrappers(sipWrapperVisitorFunc visitor,
        void *closure)
{
    const sipHashEntry *he;
    unsigned long i;

    for (he = cppPyMap.hash_array, i = 0; i < cppPyMap.size; ++i, ++he)
    {
        if (he->key != NULL)
        {
            sipSimpleWrapper *sw;

            for (sw = he->first; sw != NULL; sw = sw->next)
                visitor(sw, closure);
        }
    }
}




/*
 * Return the next exception handler.  The order is undefined.
 */
sipExceptionHandler sip_api_next_exception_handler(void **statep)
{
    sipExportedModuleDef *em = *(sipExportedModuleDef **)statep;

    if (em != NULL)
        em = em->em_next;
    else
        em = moduleList;

    while (em->em_exception_handler == NULL)
        if ((em = em->em_next) == NULL)
            return NULL;

    *statep = em;

    return em->em_exception_handler;
}
