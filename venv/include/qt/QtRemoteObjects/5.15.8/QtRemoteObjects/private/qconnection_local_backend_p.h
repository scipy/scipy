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

#ifndef QCONNECTIONCLIENTFACTORY_P_H
#define QCONNECTIONCLIENTFACTORY_P_H

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

#include <QtNetwork/qlocalserver.h>
#include <QtNetwork/qlocalsocket.h>

QT_BEGIN_NAMESPACE

class LocalClientIo final : public ClientIoDevice
{
    Q_OBJECT

public:
    explicit LocalClientIo(QObject *parent = nullptr);
    ~LocalClientIo() override;

    QIODevice *connection() const override;
    void connectToServer() override;
    bool isOpen() const override;

public Q_SLOTS:
    void onError(QLocalSocket::LocalSocketError error);
    void onStateChanged(QLocalSocket::LocalSocketState state);

protected:
    void doClose() override;
    void doDisconnectFromServer() override;
private:
    QLocalSocket* m_socket;
};

class LocalServerIo final : public ServerIoDevice
{
    Q_OBJECT
public:
    explicit LocalServerIo(QLocalSocket *conn, QObject *parent = nullptr);

    QIODevice *connection() const override;
protected:
    void doClose() override;

private:
    QLocalSocket *m_connection;
};

class LocalServerImpl final : public QConnectionAbstractServer
{
    Q_OBJECT
    Q_DISABLE_COPY(LocalServerImpl)

public:
    explicit LocalServerImpl(QObject *parent);
    ~LocalServerImpl() override;

    bool hasPendingConnections() const override;
    ServerIoDevice *configureNewConnection() override;
    QUrl address() const override;
    bool listen(const QUrl &address) override;
    QAbstractSocket::SocketError serverError() const override;
    void close() override;

private:
    QLocalServer m_server;
};

QT_END_NAMESPACE

#endif
