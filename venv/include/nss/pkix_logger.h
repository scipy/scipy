/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_logger.h
 *
 * Logger Object Type Definition
 *
 */

#ifndef _PKIX_LOGGER_H
#define _PKIX_LOGGER_H

#include "pkix_tools.h"

#ifdef __cplusplus
extern "C" {
#endif

extern PKIX_List *pkixLoggers;
extern PKIX_List *pkixLoggersErrors;
extern PKIX_List *pkixLoggersDebugTrace;

struct PKIX_LoggerStruct {
        PKIX_Logger_LogCallback callback;
        PKIX_PL_Object *context;
        PKIX_UInt32 maxLevel;
        PKIX_ERRORCLASS logComponent;
};

PKIX_Error *
pkix_Logger_Check(
        PKIX_List *pkixLoggersList,
        const char *message,
        const char *message2,
        PKIX_ERRORCLASS logComponent,
        PKIX_UInt32 maxLevel,
        void *plContext);

PKIX_Error *
pkix_Logger_CheckWithCode(
        PKIX_List *pkixLoggersList,
        PKIX_UInt32 errorCode,
        const char *message2,
        PKIX_ERRORCLASS logComponent,
        PKIX_UInt32 maxLevel,
        void *plContext);

/* see source file for function documentation */

PKIX_Error *pkix_Logger_RegisterSelf(void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_LOGGER_H */
