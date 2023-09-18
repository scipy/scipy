/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
/*
 * pkix_pl_socket.h
 *
 * Socket Object Type Definition
 *
 */

#ifndef _PKIX_PL_SOCKET_H
#define _PKIX_PL_SOCKET_H

#include <errno.h>
#include "pkix_pl_common.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
        SOCKET_BOUND,
        SOCKET_LISTENING,
        SOCKET_ACCEPTPENDING,
        SOCKET_UNCONNECTED,
        SOCKET_CONNECTPENDING,
        SOCKET_CONNECTED,
        SOCKET_SENDPENDING,
        SOCKET_RCVPENDING,
        SOCKET_SENDRCVPENDING,
        SOCKET_SHUTDOWN
} SockStatus;

/* This is the default port number, if none is supplied to CreateByName. */
#define LDAP_PORT 389

/*
 * These callbacks allow a user to substitute a counterfeit socket in places
 * where a PKIX_PL_Socket is expected. A conforming usage will use the
 * ListenCallback function instead of Listen, AcceptCallback instead of Accept,
 * etc. The counterfeit socket may have special capabilites such as the
 * ability to do proxy authentication, etc.
 */

typedef PKIX_Error *
(*pkix_pl_Socket_ListenCallback)(
        PKIX_PL_Socket *socket,
        PKIX_UInt32 backlog,
        void *plContext);

typedef PKIX_Error *
(*pkix_pl_Socket_AcceptCallback)(
        PKIX_PL_Socket *socket,
        PKIX_PL_Socket **pRendezvousSock,
        void *plContext);

typedef PKIX_Error *
(*pkix_pl_Socket_ConnectContinueCallback)(
        PKIX_PL_Socket *socket,
        PRErrorCode *pStatus,
        void *plContext);

typedef PKIX_Error *
(*pkix_pl_Socket_SendCallback)(
        PKIX_PL_Socket *sendSock,
        void *buf,
        PKIX_UInt32 bytesToWrite,
        PKIX_Int32 *pBytesWritten,
        void *plContext);

typedef PKIX_Error *
(*pkix_pl_Socket_RecvCallback)(
        PKIX_PL_Socket *rcvSock,
        void *buf,
        PKIX_UInt32 capacity,
        PKIX_Int32 *pBytesRead,
        void *plContext);

typedef PKIX_Error *
(*pkix_pl_Socket_PollCallback)(
        PKIX_PL_Socket *sock,
        PKIX_Int32 *pBytesWritten,
        PKIX_Int32 *pBytesRead,
        void *plContext);

typedef PKIX_Error *
(*pkix_pl_Socket_ShutdownCallback)(
        PKIX_PL_Socket *socket, void *plContext);

typedef struct PKIX_PL_Socket_CallbackStruct {
        pkix_pl_Socket_ListenCallback listenCallback;
        pkix_pl_Socket_AcceptCallback acceptCallback;
        pkix_pl_Socket_ConnectContinueCallback connectcontinueCallback;
        pkix_pl_Socket_SendCallback sendCallback;
        pkix_pl_Socket_RecvCallback recvCallback;
        pkix_pl_Socket_PollCallback pollCallback;
        pkix_pl_Socket_ShutdownCallback shutdownCallback;
} PKIX_PL_Socket_Callback;

struct PKIX_PL_SocketStruct {
        PKIX_Boolean isServer;
        PRIntervalTime timeout; /* zero for non-blocking I/O */
        SockStatus status;
        PRFileDesc *clientSock;
        PRFileDesc *serverSock;
        void *readBuf;
        void *writeBuf;
        PKIX_UInt32 readBufSize;
        PKIX_UInt32 writeBufSize;
        PRNetAddr *netAddr;
        PKIX_PL_Socket_Callback callbackList;
};

/* see source file for function documentation */

PKIX_Error *pkix_pl_Socket_RegisterSelf(void *plContext);

PKIX_Error *
pkix_pl_Socket_Create(
        PKIX_Boolean isServer,
        PRIntervalTime timeout, /* zero for non-blocking I/O */
        PRNetAddr *netAddr,
        PRErrorCode *status,
        PKIX_PL_Socket **pSocket,
        void *plContext);

PKIX_Error *
pkix_pl_Socket_CreateByName(
        PKIX_Boolean isServer,
        PRIntervalTime timeout,
        char *serverName,
        PRErrorCode *pStatus,
        PKIX_PL_Socket **pSocket,
        void *plContext);

PKIX_Error *
pkix_pl_Socket_CreateByHostAndPort(
        PKIX_Boolean isServer,
        PRIntervalTime timeout,
        char *hostname,
        PRUint16 portnum,
        PRErrorCode *pStatus,
        PKIX_PL_Socket **pSocket,
        void *plContext);

/* Do not use these functions directly; use their callback variants instead
 *      static PKIX_Error *
 *      pkix_pl_Socket_Listen(
 *              PKIX_PL_Socket *socket,
 *              PKIX_UInt32 backlog,
 *              void *plContext);
 *      
 *      static PKIX_Error *
 *      pkix_pl_Socket_Accept(
 *              PKIX_PL_Socket *socket,
 *              PKIX_PL_Socket **pRendezvousSock,
 *              void *plContext);
 *      
 *      static PKIX_Error *
 *      pkix_pl_Socket_ConnectContinue(
 *              PKIX_PL_Socket *socket,
 *              PRErrorCode *pStatus,
 *              void *plContext);
 *      
 *      static PKIX_Error *
 *      pkix_pl_Socket_Send(
 *              PKIX_PL_Socket *sendSock,
 *              void *buf,
 *              PKIX_UInt32 bytesToWrite,
 *              PKIX_Int32 *pBytesWritten,
 *              void *plContext);
 *      
 *      static PKIX_Error *
 *      pkix_pl_Socket_Recv(
 *              PKIX_PL_Socket *rcvSock,
 *              void *buf,
 *              PKIX_UInt32 capacity,
 *              PKIX_Int32 *pBytesRead,
 *              void *plContext);
 *      
 *      static PKIX_Error *
 *      pkix_pl_Socket_Poll(
 *              PKIX_PL_Socket *sock,
 *              PKIX_Int32 *pBytesWritten,
 *              PKIX_Int32 *pBytesRead,
 *              void *plContext);
 *      
 *      static PKIX_Error *
 *      pkix_pl_Socket_Shutdown(
 *              PKIX_PL_Socket *socket, void *plContext);
 */

PKIX_Error *
pkix_pl_Socket_GetCallbackList(
        PKIX_PL_Socket *socket,
        PKIX_PL_Socket_Callback **pCallbackList,
        void *plContext);

PKIX_Error *
pkix_pl_Socket_GetPRFileDesc(
        PKIX_PL_Socket *socket,
        PRFileDesc **pDesc,
        void *plContext);

#ifdef __cplusplus
}
#endif

#endif /* _PKIX_PL_SOCKET_H */
