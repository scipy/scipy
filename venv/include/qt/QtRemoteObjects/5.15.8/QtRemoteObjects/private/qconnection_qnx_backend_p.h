/****************************************************************************
**
** Copyright (C) 2017-2016 Ford Motor Company
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtRemoteObjects module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QCONNECTIONQNXBACKEND_P_H
#define QCONNECTIONQNXBACKEND_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include "qconnectionfactories_p.h"
#include "qconnection_qnx_qiodevices.h"
#include "qconnection_qnx_server.h"

QT_BEGIN_NAMESPACE

/*!
    QtRO provides ClientIoDevice, ServerIoDevice and QConnectionAbstractServer
    as abstract interfaces to allow different backends to be used by QtRO. The
    concept behind these classes is that there needs to be a Host node, which
    has an address that can be connected to. Then there is a client object,
    which can be publicly constructed, and can connect to the server. When the
    server gets a connection request, it creates the server side of the
    connection, which communicates directly with the client. There are thus
    three abstractions, one for the server, one for the client-side of the
    connection, and the third for the server-side of the connection. The later
    two need to inherit from QIODevice.

    Creating a backend for something that is already implemented in Qt is a
    matter of creating the three needed abstractions. In the case of creating a
    QNX backend using QNX's Native Messaging, the backend needs to create the
    Server (which has an address for accepting connections), the client
    QIODevice, and the server side QIODevice. Since Native Messaging is one
    way, and recommends using pulses to support two-way communication, the
    logic for the client-side and server-side QIODevice are very different.
    Thus, three additional backend classes are needed as well.

    QnxClientIo implements the QtRO ClientIoDevice wrapper around the QNX
    specific QQnxNativeIo QIODevice (the client-side QIODevice).

    QnxServerIo implements the QtRO ServerIoDevice wrapper around the QNX
    specific QIOQnxSource QIODevice (the server-side QIODevice).

    QnxServerImpl implements the QtRO QConnectionAbstractServer wrapper around
    the QNX specific QQnxNativeServer, which is the server object listening for
    connections.

    Not sure if it is of interest to the Qt community, but it seems like
    QQnxNativeIo, QIOQnxSource and QQnxNativeServer could used as an optimized
    QLocalServer/QLocalSocket QPA for QNX.
*/

class QnxClientIo final : public ClientIoDevice
{
    Q_OBJECT

public:
    explicit QnxClientIo(QObject *parent = nullptr);
    ~QnxClientIo() override;

    QIODevice *connection() const override;
    void connectToServer() override;
    bool isOpen() const override;

public Q_SLOTS:
    void onError(QAbstractSocket::SocketError error);
    void onStateChanged(QAbstractSocket::SocketState state);

protected:
    void doClose() override;
    void doDisconnectFromServer() override;
private:
    QQnxNativeIo *m_socket;
};

class QnxServerIo final : public ServerIoDevice
{
public:
    explicit QnxServerIo(QSharedPointer<QIOQnxSource> conn, QObject *parent = nullptr);

    QIODevice *connection() const override;
protected:
    void doClose() override;

private:
    //TODO Source or Replica
    QSharedPointer<QIOQnxSource> m_connection;
};

class QnxServerImpl final : public QConnectionAbstractServer
{
    Q_OBJECT

public:
    explicit QnxServerImpl(QObject *parent);
    ~QnxServerImpl() override;

    bool hasPendingConnections() const override;
    ServerIoDevice *configureNewConnection() override;
    QUrl address() const override;
    bool listen(const QUrl &address) override;
    QAbstractSocket::SocketError serverError() const override;
    void close() override;

private:
    QQnxNativeServer m_server;
};

QT_END_NAMESPACE

#endif // QCONNECTIONQNXBACKEND_P_H

