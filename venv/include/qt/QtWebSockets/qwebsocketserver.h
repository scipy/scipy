/****************************************************************************
**
** Copyright (C) 2016 Kurt Pattyn <pattyn.kurt@gmail.com>.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtWebSockets module of the Qt Toolkit.
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

#ifndef QWEBSOCKETSERVER_H
#define QWEBSOCKETSERVER_H

#include "QtWebSockets/qwebsockets_global.h"
#include "QtWebSockets/qwebsocketprotocol.h"

#include <QtCore/QObject>
#include <QtCore/QString>
#include <QtNetwork/QHostAddress>

#ifndef QT_NO_SSL
#include <QtNetwork/QSslConfiguration>
#include <QtNetwork/QSslError>
#endif

#if QT_HAS_INCLUDE(<chrono>)
#include <chrono>
#endif

QT_BEGIN_NAMESPACE

class QTcpSocket;
class QWebSocketServerPrivate;
class QWebSocket;
class QWebSocketCorsAuthenticator;

class Q_WEBSOCKETS_EXPORT QWebSocketServer : public QObject
{
    Q_OBJECT
    Q_DISABLE_COPY(QWebSocketServer)
    Q_DECLARE_PRIVATE(QWebSocketServer)

public:
    enum SslMode {
#ifndef QT_NO_SSL
        SecureMode = 0,
#endif
        NonSecureMode = 1
    };
    Q_ENUM(SslMode)

    explicit QWebSocketServer(const QString &serverName, SslMode secureMode,
                              QObject *parent = nullptr);
    ~QWebSocketServer() override;

    bool listen(const QHostAddress &address = QHostAddress::Any, quint16 port = 0);
    void close();

    bool isListening() const;

    void setMaxPendingConnections(int numConnections);
    int maxPendingConnections() const;

#if QT_HAS_INCLUDE(<chrono>) || defined(Q_CLANG_QDOC)
    void setHandshakeTimeout(std::chrono::milliseconds msec)
    {
        setHandshakeTimeout(int(msec.count()));
    }
    std::chrono::milliseconds handshakeTimeout() const
    {
        return std::chrono::milliseconds(handshakeTimeoutMS());
    }
#endif
    void setHandshakeTimeout(int msec);
    int handshakeTimeoutMS() const;

    quint16 serverPort() const;
    QHostAddress serverAddress() const;
    QUrl serverUrl() const;

    SslMode secureMode() const;

#if (QT_VERSION >= QT_VERSION_CHECK(6, 0, 0))
    bool setSocketDescriptor(qintptr socketDescriptor);
    qintptr socketDescriptor() const;
    bool setNativeDescriptor(qintptr descriptor) { return setSocketDescriptor(descriptor); }
    qintptr nativeDescriptor() const { return socketDescriptor(); }
#else // ### Qt 6: Remove leftovers
    Q_DECL_DEPRECATED_X("Use setNativeDescriptor") bool setSocketDescriptor(int socketDescriptor);
    Q_DECL_DEPRECATED_X("Use nativeDescriptor") int socketDescriptor() const;
    bool setNativeDescriptor(qintptr descriptor);
    qintptr nativeDescriptor() const;
#endif // (QT_VERSION >= QT_VERSION_CHECK(6, 0, 0))


    bool hasPendingConnections() const;
    virtual QWebSocket *nextPendingConnection();

    QWebSocketProtocol::CloseCode error() const;
    QString errorString() const;

    void pauseAccepting();
    void resumeAccepting();

    void setServerName(const QString &serverName);
    QString serverName() const;

#ifndef QT_NO_NETWORKPROXY
    void setProxy(const QNetworkProxy &networkProxy);
    QNetworkProxy proxy() const;
#endif
#ifndef QT_NO_SSL
    void setSslConfiguration(const QSslConfiguration &sslConfiguration);
    QSslConfiguration sslConfiguration() const;
#endif

    QList<QWebSocketProtocol::Version> supportedVersions() const;

    void handleConnection(QTcpSocket *socket) const;

Q_SIGNALS:
    void acceptError(QAbstractSocket::SocketError socketError);
    void serverError(QWebSocketProtocol::CloseCode closeCode);
    //TODO: should use a delegate iso of a synchronous signal
    //see also QTBUG-16251
    void originAuthenticationRequired(QWebSocketCorsAuthenticator *pAuthenticator);
    void newConnection();
#ifndef QT_NO_SSL
    void peerVerifyError(const QSslError &error);
    void sslErrors(const QList<QSslError> &errors);
    void preSharedKeyAuthenticationRequired(QSslPreSharedKeyAuthenticator *authenticator);
#endif
    void closed();
};

QT_END_NAMESPACE

#endif // QWEBSOCKETSERVER_H
