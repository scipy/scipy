/****************************************************************************
**
** Copyright (C) 2017-2015 Ford Motor Company
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

#ifndef QCONNECTIONTCPIPBACKEND_P_H
#define QCONNECTIONTCPIPBACKEND_P_H

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

#include <QtNetwork/qtcpserver.h>
#include <QtNetwork/qtcpsocket.h>

QT_BEGIN_NAMESPACE

class TcpClientIo final : public ClientIoDevice
{
    Q_OBJECT

public:
    explicit TcpClientIo(QObject *parent = nullptr);
    ~TcpClientIo() override;

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
    QTcpSocket *m_socket;
};

class TcpServerIo final : public ServerIoDevice
{
    Q_OBJECT
public:
    explicit TcpServerIo(QTcpSocket *conn, QObject *parent = nullptr);

    QIODevice *connection() const override;
protected:
    void doClose() override;

private:
    QTcpSocket *m_connection;
};

class TcpServerImpl final : public QConnectionAbstractServer
{
    Q_OBJECT
    Q_DISABLE_COPY(TcpServerImpl)

public:
    explicit TcpServerImpl(QObject *parent);
    ~TcpServerImpl() override;

    bool hasPendingConnections() const override;
    ServerIoDevice *configureNewConnection() override;
    QUrl address() const override;
    bool listen(const QUrl &address) override;
    QAbstractSocket::SocketError serverError() const override;
    void close() override;

private:
    QTcpServer m_server;
    QUrl m_originalUrl; // necessary because of a QHostAddress bug
};

QT_END_NAMESPACE
#endif // QCONNECTIONTCPIPBACKEND_P_H
