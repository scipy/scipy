/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * These functions provide support for a number of other functions
 * by creating and manipulating data structures used by those functions.
 *
 */

#ifndef _PKIX_UTIL_H
#define _PKIX_UTIL_H

#include "pkixt.h"

#ifdef __cplusplus
extern "C" {
#endif

/* General
 *
 * Please refer to the libpkix Programmer's Guide for detailed information
 * about how to use the libpkix library. Certain key warnings and notices from
 * that document are repeated here for emphasis.
 *
 * All identifiers in this file (and all public identifiers defined in
 * libpkix) begin with "PKIX_". Private identifiers only intended for use
 * within the library begin with "pkix_".
 *
 * A function returns NULL upon success, and a PKIX_Error pointer upon failure.
 *
 * Unless otherwise noted, for all accessor (gettor) functions that return a
 * PKIX_PL_Object pointer, callers should assume that this pointer refers to a
 * shared object. Therefore, the caller should treat this shared object as
 * read-only and should not modify this shared object. When done using the
 * shared object, the caller should release the reference to the object by
 * using the PKIX_PL_Object_DecRef function.
 *
 * While a function is executing, if its arguments (or anything referred to by
 * its arguments) are modified, free'd, or destroyed, the function's behavior
 * is undefined.
 *
 */

/* PKIX_Logger
 *
 * PKIX_Loggers provide a standard way for the caller to insert custom logging
 * facilities. These are used by libpkix to log errors, debug information,
 * status, etc. The LogCallback allows custom logging to take place.
 * Additionally, a Logger can be initialized with a loggerContext, which is
 * where the caller can specify configuration data such as the name of a
 * logfile or database. Note that this loggerContext must be a PKIX_PL_Object,
 * allowing it to be reference-counted and allowing it to provide the standard
 * PKIX_PL_Object functions (Equals, Hashcode, ToString, Compare, Duplicate).
 *
 * Once the caller has created the Logger object(s) (and set the loggerContext
 * (if any) and the Log callback), the caller then registers these Loggers
 * with the system by calling PKIX_SetLoggers or PKIX_AddLogger. All log
 * entries will then be logged using the specified Loggers. If multiple
 * Loggers are specified, every log entry will be logged with each of them.
 *
 * XXX Maybe give some guidance somewhere on how much detail each logging
 * level should have and where component boundaries should be. Maybe in
 * Implementor's Guide or Programmer's Guide.
 */

#define PKIX_LOGGER_LEVEL_TRACE                5
#define PKIX_LOGGER_LEVEL_DEBUG                4
#define PKIX_LOGGER_LEVEL_WARNING              3
#define PKIX_LOGGER_LEVEL_ERROR                2
#define PKIX_LOGGER_LEVEL_FATALERROR           1

#define PKIX_LOGGER_LEVEL_MAX                  5

/*
 * FUNCTION: PKIX_Logger_LogCallback
 * DESCRIPTION:
 *
 *  This callback function logs a log entry containing the String pointed to
 *  by "message", the integer value of logLevel, and the String pointed to by
 *  "logComponent". A log entry can be associated with a particular log
 *  level (i.e. level 3) and a particular log component (i.e. "CertStore").
 *  For example, someone reading the log may only be interested in very general
 *  log entries so they look only for log level 1. Similarly, they may only be
 *  interested in log entries pertaining to the CertStore component so they
 *  look only for that log component. This function can be used before calling
 *  PKIX_Initialize.
 *
 * PARAMETERS:
 *  "logger"
 *      Address of logger whose LogCallback is to be used. Must be non-NULL.
 *  "message"
 *      Address of String that is to be logged used "logger". Must be non-NULL.
 *  "logLevel"
 *      Integer value representing the log level for this entry. The higher the
 *      level, the more detail. Must be non-NULL.
 *  "logComponent"
 *      PKIXERRORNUM value (defined in pkixt.h) designating the log component
 *      for this entry.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe
 *
 *  Multiple threads must be able to safely call this function without
 *  worrying about conflicts, even if they're operating on the same objects.
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
typedef PKIX_Error *
(*PKIX_Logger_LogCallback)(
        PKIX_Logger *logger,
        PKIX_PL_String *message,
        PKIX_UInt32 logLevel,
        PKIX_ERRORCLASS logComponent,
        void *plContext);

/*
 * FUNCTION: PKIX_Logger_Create
 * DESCRIPTION:
 *
 *  Creates a new Logger using the Object pointed to by "loggerContext"
 *  (if any) and stores it at "pLogger". The new Logger uses the LogCallback
 *  pointed to by "callback". The Logger's maximum logging level is initially
 *  set to a very high level and its logging component is set to NULL (all
 *  components).
 *
 * PARAMETERS:
 *  "callback"
 *      The LogCallback function to be used. Must be non-NULL.
 *  "loggerContext"
 *      Address of Object representing the Logger's context (if any).
 *  "pLogger"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Logger_Create(
        PKIX_Logger_LogCallback callback,
        PKIX_PL_Object *loggerContext,
        PKIX_Logger **pLogger,
        void *plContext);

/*
 * FUNCTION: PKIX_Logger_GetLogCallback
 * DESCRIPTION:
 *
 *  Retrieves a pointer to "logger's" Log callback function and puts it in
 *  "pCallback".
 *
 * PARAMETERS:
 *  "logger"
 *      Address of Logger whose Log callback is desired. Must be non-NULL.
 *  "pCallback"
 *      Address where Log callback function pointer will be stored.
 *      Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Logger_GetLogCallback(
        PKIX_Logger *logger,
        PKIX_Logger_LogCallback *pCallback,
        void *plContext);

/*
 * FUNCTION: PKIX_Logger_GetLoggerContext
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a PKIX_PL_Object representing the context (if any)
 *  of the Logger pointed to by "logger" and stores it at "pLoggerContext".
 *
 * PARAMETERS:
 *  "logger"
 *      Address of Logger whose context is to be stored. Must be non-NULL.
 *  "pLoggerContext"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Logger_GetLoggerContext(
        PKIX_Logger *logger,
        PKIX_PL_Object **pLoggerContext,
        void *plContext);

/*
 * FUNCTION: PKIX_Logger_GetMaxLoggingLevel
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a PKIX_UInt32 representing the maximum logging
 *  level of the Logger pointed to by "logger" and stores it at "pLevel". Only
 *  log entries whose log level is less than or equal to this maximum logging
 *  level will be logged.
 *
 * PARAMETERS:
 *  "logger"
 *      Address of Logger whose maximum logging level is to be stored.
 *      Must be non-NULL.
 *  "pLevel"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Logger_GetMaxLoggingLevel(
        PKIX_Logger *logger,
        PKIX_UInt32 *pLevel,
        void *plContext);

/*
 * FUNCTION: PKIX_Logger_SetMaxLoggingLevel
 * DESCRIPTION:
 *
 *  Sets the maximum logging level of the Logger pointed to by "logger" with
 *  the integer value of "level".
 *
 * PARAMETERS:
 *  "logger"
 *      Address of Logger whose maximum logging level is to be set.
 *      Must be non-NULL.
 *  "level"
 *      Maximum logging level to be set
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "logger"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Logger_SetMaxLoggingLevel(
        PKIX_Logger *logger,
        PKIX_UInt32 level,
        void *plContext);

/*
 * FUNCTION: PKIX_Logger_GetLoggingComponent
 * DESCRIPTION:
 *
 *  Retrieves a pointer to a String representing the logging component of the
 *  Logger pointed to by "logger" and stores it at "pComponent". Only log
 *  entries whose log component matches the specified logging component will
 *  be logged.
 *
 * PARAMETERS:
 *  "logger"
 *      Address of Logger whose logging component is to be stored.
 *      Must be non-NULL.
 *  "pComponent"
 *      Address where PKIXERRORNUM will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Logger_GetLoggingComponent(
        PKIX_Logger *logger,
        PKIX_ERRORCLASS *pComponent,
        void *plContext);

/*
 * FUNCTION: PKIX_Logger_SetLoggingComponent
 * DESCRIPTION:
 *
 *  Sets the logging component of the Logger pointed to by "logger" with the
 *  PKIXERRORNUM pointed to by "component". To match a small set of components,
 *  create a Logger for each.
 *
 * PARAMETERS:
 *  "logger"
 *      Address of Logger whose logging component is to be set.
 *      Must be non-NULL.
 *  "component"
 *      PKIXERRORNUM value representing logging component to be set.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "logger"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Logger_SetLoggingComponent(
        PKIX_Logger *logger,
        PKIX_ERRORCLASS component,
        void *plContext);

/*
 * FUNCTION: PKIX_GetLoggers
 * DESCRIPTION:
 *
 *  Retrieves a pointer to the List of Loggers (if any) being used for logging
 *  by libpkix and stores it at "pLoggers". If no loggers are being used, this
 *  function stores an empty List at "pLoggers".
 *
 *  Note that the List returned by this function is immutable.
 *
 * PARAMETERS:
 *  "pLoggers"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_GetLoggers(
        PKIX_List **pLoggers,  /* list of PKIX_Logger */
        void *plContext);

/*
 * FUNCTION: PKIX_SetLoggers
 * DESCRIPTION:
 *
 *  Sets the Loggers to be used by libpkix to the List of Loggers pointed to
 *  by "loggers". If "loggers" is NULL, no Loggers will be used.
 *
 * PARAMETERS:
 *  "loggers"
 *      Address of List of Loggers to be set. NULL for no Loggers.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_SetLoggers(
        PKIX_List *loggers,  /* list of PKIX_Logger */
        void *plContext);

/*
 * FUNCTION: PKIX_AddLogger
 * DESCRIPTION:
 *
 *  Adds the Logger pointed to by "logger" to the List of Loggers used by
 *  libpkix.
 *
 * PARAMETERS:
 *  "logger"
 *      Address of Logger to be added. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Logger Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_AddLogger(
        PKIX_Logger *logger,
        void *plContext);

/* Functions pertaining to the PKIX_Error type */

/* Error
 *
 * An Error object is returned by a function upon encountering some error
 * condition. Each Error is associated with an errorCode specified in pkixt.h.
 * The remaining components of an Error are optional. An Error's description
 * specifies a text message describing the Error. An Error's supplementary info
 * specifies additional information that might be useful. Finally, an Error's
 * cause specifies the underlying Error (if any) that resulted in this Error
 * being returned, thereby allowing Errors to be chained so that an entire
 * "error stack trace" can be represented. Once created, an Error is immutable.
 *
 * Note that the Error's supplementary info must be an Object (although any
 * object type), allowing it to be reference-counted and allowing it to
 * provide the standard Object functions (Equals, Hashcode, ToString, Compare,
 * Duplicate).
 *
 * Errors are classified as either being fatal or non-fatal. If a function
 * fails in an unrecoverable way, it returns an Error whose errorCode is
 * PKIX_FATAL_ERROR. If such an error is encountered, the caller should
 * not attempt to recover since something seriously wrong has happened
 * (e.g. corrupted memory, memory finished, etc.). All other errorCodes
 * are considered non-fatal errors and can be handled by the caller as they
 * see fit.
 */

/*
 * FUNCTION: PKIX_Error_Create
 * DESCRIPTION:
 *
 *  Creates a new Error using the value of "errorCode", the Error pointed to by
 *  "cause" (if any), the Object pointed to by "info" (if any), and the String
 *  pointed to by "desc" and stores it at "pError". If any error occurs during
 *  error allocation, it will be returned without chaining, since new errors
 *  cannot be created. Once created, an Error is immutable.
 *
 * PARAMETERS:
 *  "errorCode"
 *      Value of error code.
 *  "cause"
 *      Address of Error representing error's cause.
 *      NULL if none or unspecified.
 *  "info"
 *      Address of Object representing error's supplementary information.
 *      NULL if none.
 *  "desc"
 *      Address of String representing error's description. NULL if none.
 *  "pError"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Error Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Error_Create(
        PKIX_ERRORCLASS errClass,
        PKIX_Error *cause,
        PKIX_PL_Object *info,
        PKIX_ERRORCODE errCode,
        PKIX_Error **pError,
        void *plContext);

/*
 * FUNCTION: PKIX_Error_GetErrorClass
 * DESCRIPTION:
 *
 *  Retrieves the error class of the Error pointed to by "error" and 
 *  stores it at "pClass". Supported error codes are defined in pkixt.h.
 *
 * PARAMETERS:
 *  "error"
 *      Address of Error whose error code is desired. Must be non-NULL.
 *  "pClass"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Error Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Error_GetErrorClass(
        PKIX_Error *error,
        PKIX_ERRORCLASS *pClass,
        void *plContext);

/*
 * FUNCTION: PKIX_Error_GetErrorCode
 * DESCRIPTION:
 *
 *  Retrieves the error code of the Error pointed to by "error" and 
 *  stores it at "pCode". Supported error codes are defined in pkixt.h.
 *
 * PARAMETERS:
 *  "error"
 *      Address of Error whose error code is desired. Must be non-NULL.
 *  "pCode"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Error Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Error_GetErrorCode(
        PKIX_Error *error,
        PKIX_ERRORCODE *pCode,
        void *plContext);

/*
 * FUNCTION: PKIX_Error_GetCause
 * DESCRIPTION:
 *
 *  Retrieves the cause of the Error pointed to by "error" and stores it at
 *  "pCause". If no cause was specified, NULL will be stored at "pCause".
 *
 * PARAMETERS:
 *  "error"
 *      Address of Error whose cause is desired. Must be non-NULL.
 *  "pCause"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Error Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Error_GetCause(
        PKIX_Error *error,
        PKIX_Error **pCause,
        void *plContext);

/*
 * FUNCTION: PKIX_Error_GetSupplementaryInfo
 * DESCRIPTION:
 *
 *  Retrieves the supplementary info of the Error pointed to by "error" and
 *  stores it at "pInfo".
 *
 * PARAMETERS:
 *  "error"
 *      Address of Error whose info is desired. Must be non-NULL.
 *  "pInfo"
 *      Address where info pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Error Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Error_GetSupplementaryInfo(
        PKIX_Error *error,
        PKIX_PL_Object **pInfo,
        void *plContext);

/*
 * FUNCTION: PKIX_Error_GetDescription
 * DESCRIPTION:
 *
 *  Retrieves the description of the Error pointed to by "error" and stores it
 *  at "pDesc". If no description was specified, NULL will be stored at
 *  "pDesc".
 *
 * PARAMETERS:
 *  "error"
 *      Address of Error whose description is desired. Must be non-NULL.
 *  "pDesc"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns an Error Error if the function fails in a non-fatal way.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_Error_GetDescription(
        PKIX_Error *error,
        PKIX_PL_String **pDesc,
        void *plContext);

/* PKIX_List
 *
 * Represents a collection of items. NULL is considered a valid item.
 */

/*
 * FUNCTION: PKIX_List_Create
 * DESCRIPTION:
 *
 *  Creates a new List and stores it at "pList". The List is initially empty
 *  and holds no items. To initially add items to the List, use
 *  PKIX_List_AppendItem
 *
 * PARAMETERS:
 *  "pList"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Thread Safe (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_Create(
        PKIX_List **pList,
        void *plContext);

/*
 * FUNCTION: PKIX_List_SetImmutable
 * DESCRIPTION:
 *
 *  Sets the List pointed to by "list" to be immutable. If a caller tries to
 *  change a List after it has been marked immutable (i.e. by calling
 *  PKIX_List_AppendItem, PKIX_List_InsertItem, PKIX_List_SetItem, or
 *  PKIX_List_DeleteItem), an Error is returned.
 *
 * PARAMETERS:
 *  "list"
 *      Address of List to be marked immutable. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "list"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_SetImmutable(
        PKIX_List *list,
        void *plContext);

/*
 * FUNCTION: PKIX_List_IsImmutable
 * DESCRIPTION:
 *
 *  Checks whether the List pointed to by "list" is immutable and stores
 *  the Boolean result at "pImmutable". If a caller tries to change a List
 *  after it has been marked immutable (i.e. by calling PKIX_List_AppendItem,
 *  PKIX_List_InsertItem, PKIX_List_SetItem, or PKIX_List_DeleteItem), an
 *  Error is returned.
 *
 * PARAMETERS:
 *  "list"
 *      Address of List whose immutability is to be determined.
 *      Must be non-NULL.
 *  "pImmutable"
 *      Address where PKIX_Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_IsImmutable(
        PKIX_List *list,
        PKIX_Boolean *pImmutable,
        void *plContext);

/*
 * FUNCTION: PKIX_List_GetLength
 * DESCRIPTION:
 *
 *  Retrieves the length of the List pointed to by "list" and stores it at
 *  "pLength".
 *
 * PARAMETERS:
 *  "list"
 *      Address of List whose length is desired. Must be non-NULL.
 *  "pLength"
 *      Address where PKIX_UInt32 will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_GetLength(
        PKIX_List *list,
        PKIX_UInt32 *pLength,
        void *plContext);

/*
 * FUNCTION: PKIX_List_IsEmpty
 * DESCRIPTION:
 *
 *  Checks whether the List pointed to by "list" is empty and stores
 *  the Boolean result at "pEmpty".
 *
 * PARAMETERS:
 *  "list"
 *      Address of List whose emptiness is to be determined. Must be non-NULL.
 *  "pEmpty"
 *      Address where PKIX_Boolean will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_IsEmpty(
        PKIX_List *list,
        PKIX_Boolean *pEmpty,
        void *plContext);

/*
 * FUNCTION: PKIX_List_AppendItem
 * DESCRIPTION:
 *
 *  Appends the Object pointed to by "item" after the last non-NULL item in
 *  List pointed to by "list", if any. Note that a List may validly contain
 *  NULL items. Appending "c" into the List ("a", NULL, "b", NULL) will result
 *  in ("a", NULL, "b", "c").
 *
 * PARAMETERS:
 *  "list"
 *      Address of List to append to. Must be non-NULL.
 *  "item"
 *      Address of new item to append.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "list"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_AppendItem(
        PKIX_List *list,
        PKIX_PL_Object *item,
        void *plContext);

/*
 * FUNCTION: PKIX_List_InsertItem
 * DESCRIPTION:
 *
 *  Inserts the Object pointed to by "item" into the List pointed to by "list"
 *  at the given "index". The index counts from zero and must be less than the
 *  List's length. Existing list entries at or after this index will be moved
 *  to the next highest index.
 *
 *  XXX why not allow equal to length which would be equivalent to AppendItem?
 *
 * PARAMETERS:
 *  "list"
 *      Address of List to insert into. Must be non-NULL.
 *  "index"
 *      Position to insert into. Must be less than List's length.
 *  "item"
 *      Address of new item to append.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "list"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_InsertItem(
        PKIX_List *list,
        PKIX_UInt32 index,
        PKIX_PL_Object *item,
        void *plContext);

/*
 * FUNCTION: PKIX_List_GetItem
 * DESCRIPTION:
 *
 *  Copies the "list"'s item at "index" into "pItem". The index counts from
 *  zero and must be less than the list's length. Increments the reference
 *  count on the returned object, if non-NULL.
 *
 * PARAMETERS:
 *  "list"
 *      Address of List to get item from. Must be non-NULL.
 *  "index"
 *      Index of list to get item from. Must be less than List's length.
 *  "pItem"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_GetItem(
        PKIX_List *list,
        PKIX_UInt32 index,
        PKIX_PL_Object **pItem,
        void *plContext);

/*
 * FUNCTION: PKIX_List_SetItem
 * DESCRIPTION:
 *
 *  Sets the item at "index" of the List pointed to by "list" with the Object
 *  pointed to by "item". The index counts from zero and must be less than the
 *  List's length. The previous entry at this index will have its reference
 *  count decremented and the new entry will have its reference count
 *  incremented.
 *
 * PARAMETERS:
 *  "list"
 *      Address of List to modify. Must be non-NULL.
 *  "index"
 *      Position in List to set. Must be less than List's length.
 *  "item"
 *      Address of Object to set at "index".
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "list"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_SetItem(
        PKIX_List *list,
        PKIX_UInt32 index,
        PKIX_PL_Object *item,
        void *plContext);

/*
 * FUNCTION: PKIX_List_DeleteItem
 *
 *  Deletes the item at "index" from the List pointed to by "list". The index
 *  counts from zero and must be less than the List's length. Note that this
 *  function does not destroy the List. It simply decrements the reference
 *  count of the item at "index" in the List, deletes that item from the list
 *  and moves all subsequent entries to a lower index in the list. If there is
 *  only a single element in the List and that element is deleted, then the
 *  List will be empty.
 *
 * PARAMETERS:
 *  "list"
 *      Address of List to delete from. Must be non-NULL.
 *  "index"
 *      Position in List to delete. Must be less than List's length.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Not Thread Safe - assumes exclusive access to "list"
 *  (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_DeleteItem(
        PKIX_List *list,
        PKIX_UInt32 index,
        void *plContext);

/*
 * FUNCTION: PKIX_List_ReverseList
 * DESCRIPTION:
 *
 *  Creates a new List whose elements are in the reverse order as the elements
 *  of the Object pointed to by "list" and stores the copy at "pReversedList".
 *  If "list" is empty, the new reversed List will be a copy of "list".
 *  Changes to the new object will not affect the original and vice versa.
 *
 * PARAMETERS:
 *  "list"
 *      Address of List whose elements are to be reversed. Must be non-NULL.
 *  "pReversedList"
 *      Address where object pointer will be stored. Must be non-NULL.
 *  "plContext"
 *      Platform-specific context pointer.
 * THREAD SAFETY:
 *  Conditionally Thread Safe
 *      (see Thread Safety Definitions in Programmer's Guide)
 * RETURNS:
 *  Returns NULL if the function succeeds.
 *  Returns a Fatal Error if the function fails in an unrecoverable way.
 */
PKIX_Error *
PKIX_List_ReverseList(
        PKIX_List *list,
        PKIX_List **pReversedList,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_UTIL_H */
