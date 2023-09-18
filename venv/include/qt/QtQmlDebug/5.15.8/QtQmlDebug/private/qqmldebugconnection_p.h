/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQml module of the Qt Toolkit.
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

#ifndef QQMLDEBUGCONNECTION_P_H
#define QQMLDEBUGCONNECTION_P_H

#include <QtCore/qobject.h>
#include <QtNetwork/qabstractsocket.h>

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

QT_BEGIN_NAMESPACE

class QQmlDebugClient;
class QQmlDebugConnectionPrivate;
class QQmlDebugConnection : public QObject
{
    Q_OBJECT
    Q_DISABLE_COPY(QQmlDebugConnection)
    Q_DECLARE_PRIVATE(QQmlDebugConnection)
public:
    QQmlDebugConnection(QObject *parent = nullptr);
    ~QQmlDebugConnection();

    void connectToHost(const QString &hostName, quint16 port);
    void startLocalServer(const QString &fileName);

    int currentDataStreamVersion() const;
    void setMaximumDataStreamVersion(int maximumVersion);

    bool isConnected() const;
    bool isConnecting() const;

    void close();
    bool waitForConnected(int msecs = 30000);

    QQmlDebugClient *client(const QString &name) const;
    bool addClient(const QString &name, QQmlDebugClient *client);
    bool removeClient(const QString &name);

    float serviceVersion(const QString &serviceName) const;
    bool sendMessage(const QString &name, const QByteArray &message);

signals:
    void connected();
    void disconnected();
    void socketError(QAbstractSocket::SocketError socketError);
    void socketStateChanged(QAbstractSocket::SocketState socketState);

private:
    void newConnection();
    void socketConnected();
    void socketDisconnected();
    void protocolReadyRead();
    void handshakeTimeout();
};

QT_END_NAMESPACE

#endif // QQMLDEBUGCONNECTION_P_H
