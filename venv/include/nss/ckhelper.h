/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * ckhelper.h
 *
 * This file contains some helper utilities for interaction with cryptoki.
 */

#ifndef CKHELPER_H
#define CKHELPER_H

PR_BEGIN_EXTERN_C

/* Some globals to keep from constantly redeclaring common cryptoki
 * attribute types on the stack.
 */

/* Boolean values */
NSS_EXTERN_DATA const NSSItem g_ck_true;
NSS_EXTERN_DATA const NSSItem g_ck_false;

/* Object classes */
NSS_EXTERN_DATA const NSSItem g_ck_class_cert;
NSS_EXTERN_DATA const NSSItem g_ck_class_pubkey;
NSS_EXTERN_DATA const NSSItem g_ck_class_privkey;

#define NSS_CK_TEMPLATE_START(_template, attr, size) \
    attr = _template;                                \
    size = 0;

#define NSS_CK_SET_ATTRIBUTE_ITEM(pattr, kind, item) \
    (pattr)->type = kind;                            \
    (pattr)->pValue = (CK_VOID_PTR)(item)->data;     \
    (pattr)->ulValueLen = (CK_ULONG)(item)->size;    \
    (pattr)++;

#define NSS_CK_SET_ATTRIBUTE_UTF8(pattr, kind, utf8)          \
    (pattr)->type = kind;                                     \
    (pattr)->pValue = (CK_VOID_PTR)utf8;                      \
    (pattr)->ulValueLen = (CK_ULONG)nssUTF8_Size(utf8, NULL); \
    if ((pattr)->ulValueLen)                                  \
        ((pattr)->ulValueLen)--;                              \
    (pattr)++;

#define NSS_CK_SET_ATTRIBUTE_VAR(pattr, kind, var) \
    (pattr)->type = kind;                          \
    (pattr)->pValue = (CK_VOID_PTR)&var;           \
    (pattr)->ulValueLen = (CK_ULONG)sizeof(var);   \
    (pattr)++;

#define NSS_CK_SET_ATTRIBUTE_NULL(pattr, kind) \
    (pattr)->type = kind;                      \
    (pattr)->pValue = (CK_VOID_PTR)NULL;       \
    (pattr)->ulValueLen = 0;                   \
    (pattr)++;

#define NSS_CK_TEMPLATE_FINISH(_template, attr, size) \
    size = (attr) - (_template);                      \
    PR_ASSERT(size <= sizeof(_template) / sizeof(_template[0]));

/* NSS_CK_ATTRIBUTE_TO_ITEM(attrib, item)
 *
 * Convert a CK_ATTRIBUTE to an NSSItem.
 */
#define NSS_CK_ATTRIBUTE_TO_ITEM(attrib, item)         \
    if ((CK_LONG)(attrib)->ulValueLen > 0) {           \
        (item)->data = (void *)(attrib)->pValue;       \
        (item)->size = (PRUint32)(attrib)->ulValueLen; \
    } else {                                           \
        (item)->data = 0;                              \
        (item)->size = 0;                              \
    }

#define NSS_CK_ATTRIBUTE_TO_BOOL(attrib, boolvar)         \
    if ((attrib)->ulValueLen > 0) {                       \
        if (*((CK_BBOOL *)(attrib)->pValue) == CK_TRUE) { \
            boolvar = PR_TRUE;                            \
        } else {                                          \
            boolvar = PR_FALSE;                           \
        }                                                 \
    }

#define NSS_CK_ATTRIBUTE_TO_ULONG(attrib, ulongvar) \
    if ((attrib)->ulValueLen > 0) {                 \
        ulongvar = *((CK_ULONG *)(attrib)->pValue); \
    }

/* NSS_CK_ATTRIBUTE_TO_UTF8(attrib, str)
 *
 * Convert a CK_ATTRIBUTE to a string.
 */
#define NSS_CK_ATTRIBUTE_TO_UTF8(attrib, str) \
    str = (NSSUTF8 *)((attrib)->pValue);

/* NSS_CK_ITEM_TO_ATTRIBUTE(item, attrib)
 *
 * Convert an NSSItem to a  CK_ATTRIBUTE.
 */
#define NSS_CK_ITEM_TO_ATTRIBUTE(item, attrib)    \
    (attrib)->pValue = (CK_VOID_PTR)(item)->data; \
    (attrib)->ulValueLen = (CK_ULONG)(item)->size;

/* Get an array of attributes from an object. */
NSS_EXTERN PRStatus
nssCKObject_GetAttributes(
    CK_OBJECT_HANDLE object,
    CK_ATTRIBUTE_PTR obj_template,
    CK_ULONG count,
    NSSArena *arenaOpt,
    nssSession *session,
    NSSSlot *slot);

/* Get a single attribute as an item. */
NSS_EXTERN PRStatus
nssCKObject_GetAttributeItem(
    CK_OBJECT_HANDLE object,
    CK_ATTRIBUTE_TYPE attribute,
    NSSArena *arenaOpt,
    nssSession *session,
    NSSSlot *slot,
    NSSItem *rvItem);

NSS_EXTERN PRBool
nssCKObject_IsAttributeTrue(
    CK_OBJECT_HANDLE object,
    CK_ATTRIBUTE_TYPE attribute,
    nssSession *session,
    NSSSlot *slot,
    PRStatus *rvStatus);

NSS_EXTERN PRStatus
nssCKObject_SetAttributes(
    CK_OBJECT_HANDLE object,
    CK_ATTRIBUTE_PTR obj_template,
    CK_ULONG count,
    nssSession *session,
    NSSSlot *slot);

NSS_EXTERN PRBool
nssCKObject_IsTokenObjectTemplate(
    CK_ATTRIBUTE_PTR objectTemplate,
    CK_ULONG otsize);

PR_END_EXTERN_C

#endif /* CKHELPER_H */
