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

#ifndef QWEBSOCKET_H
#define QWEBSOCKET_H

#include <QtCore/QUrl>
#include <QtNetwork/QAbstractSocket>
#include <QtNetwork/QNetworkRequest>
#ifndef QT_NO_NETWORKPROXY
#include <QtNetwork/QNetworkProxy>
#endif
#ifndef QT_NO_SSL
#include <QtNetwork/QSslError>
#include <QtNetwork/QSslConfiguration>
#endif
#include "QtWebSockets/qwebsockets_global.h"
#include "QtWebSockets/qwebsocketprotocol.h"

QT_BEGIN_NAMESPACE

class QTcpSocket;
class QWebSocketPrivate;
class QMaskGenerator;

class Q_WEBSOCKETS_EXPORT QWebSocket : public QObject
{
    Q_OBJECT
    Q_DISABLE_COPY(QWebSocket)
    Q_DECLARE_PRIVATE(QWebSocket)

public:
    explicit QWebSocket(const QString &origin = QString(),
                        QWebSocketProtocol::Version version = QWebSocketProtocol::VersionLatest,
                        QObject *parent = nullptr);
    ~QWebSocket() override;

    void abort();
    QAbstractSocket::SocketError error() const;
    QString errorString() const;
    bool flush();
    bool isValid() const;
    QHostAddress localAddress() const;
    quint16 localPort() const;
    QAbstractSocket::PauseModes pauseMode() const;
    QHostAddress peerAddress() const;
    QString peerName() const;
    quint16 peerPort() const;
#ifndef QT_NO_NETWORKPROXY
    QNetworkProxy proxy() const;
    void setProxy(const QNetworkProxy &networkProxy);
#endif
    void setMaskGenerator(const QMaskGenerator *maskGenerator);
    const QMaskGenerator *maskGenerator() const;
    qint64 readBufferSize() const;
    void setReadBufferSize(qint64 size);

    void resume();
    void setPauseMode(QAbstractSocket::PauseModes pauseMode);

    QAbstractSocket::SocketState state() const;

    QWebSocketProtocol::Version version() const;
    QString resourceName() const;
    QUrl requestUrl() const;
    QNetworkRequest request() const;
    QString origin() const;
    QWebSocketProtocol::CloseCode closeCode() const;
    QString closeReason() const;

    qint64 sendTextMessage(const QString &message);
    qint64 sendBinaryMessage(const QByteArray &data);

#ifndef QT_NO_SSL
    void ignoreSslErrors(const QList<QSslError> &errors);
    void setSslConfiguration(const QSslConfiguration &sslConfiguration);
    QSslConfiguration sslConfiguration() const;
#endif

    qint64 bytesToWrite() const;

    void setMaxAllowedIncomingFrameSize(quint64 maxAllowedIncomingFrameSize);
    quint64 maxAllowedIncomingFrameSize() const;
    void setMaxAllowedIncomingMessageSize(quint64 maxAllowedIncomingMessageSize);
    quint64 maxAllowedIncomingMessageSize() const;
    static quint64 maxIncomingMessageSize();
    static quint64 maxIncomingFrameSize();

    void setOutgoingFrameSize(quint64 outgoingFrameSize);
    quint64 outgoingFrameSize() const;
    static quint64 maxOutgoingFrameSize();

public Q_SLOTS:
    void close(QWebSocketProtocol::CloseCode closeCode = QWebSocketProtocol::CloseCodeNormal,
               const QString &reason = QString());
    void open(const QUrl &url);
    void open(const QNetworkRequest &request);
    void ping(const QByteArray &payload = QByteArray());
#ifndef QT_NO_SSL
    void ignoreSslErrors();
#endif

Q_SIGNALS:
    void aboutToClose();
    void connected();
    void disconnected();
    void stateChanged(QAbstractSocket::SocketState state);
#ifndef QT_NO_NETWORKPROXY
    void proxyAuthenticationRequired(const QNetworkProxy &proxy, QAuthenticator *pAuthenticator);
#endif
    void readChannelFinished();
    void textFrameReceived(const QString &frame, bool isLastFrame);
    void binaryFrameReceived(const QByteArray &frame, bool isLastFrame);
    void textMessageReceived(const QString &message);
    void binaryMessageReceived(const QByteArray &message);
    void error(QAbstractSocket::SocketError error);
    void pong(quint64 elapsedTime, const QByteArray &payload);
    void bytesWritten(qint64 bytes);

#ifndef QT_NO_SSL
    void sslErrors(const QList<QSslError> &errors);
    void preSharedKeyAuthenticationRequired(QSslPreSharedKeyAuthenticator *authenticator);
#endif

private:
    QWebSocket(QTcpSocket *pTcpSocket, QWebSocketProtocol::Version version,
               QObject *parent = nullptr);
};

QT_END_NAMESPACE

#endif // QWEBSOCKET_H
