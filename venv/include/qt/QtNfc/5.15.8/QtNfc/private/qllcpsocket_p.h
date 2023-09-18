/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtNfc module of the Qt Toolkit.
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

#ifndef QLLCPSOCKET_H
#define QLLCPSOCKET_H

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

#include <QtCore/QIODevice>
#include <QtNetwork/QAbstractSocket>
#include <QtNfc/qtnfcglobal.h>

QT_BEGIN_NAMESPACE

class QNearFieldTarget;
class QLlcpSocketPrivate;

class Q_NFC_EXPORT QLlcpSocket : public QIODevice
{
    Q_OBJECT

    Q_DECLARE_PRIVATE(QLlcpSocket)

    friend class QLlcpServerPrivate;

public:
    enum SocketState {
        UnconnectedState = QAbstractSocket::UnconnectedState,
        ConnectingState = QAbstractSocket::ConnectingState,
        ConnectedState = QAbstractSocket::ConnectedState,
        ClosingState = QAbstractSocket::ClosingState,
        BoundState = QAbstractSocket::BoundState,
        ListeningState = QAbstractSocket::ListeningState
    };
    Q_ENUM(SocketState)

    enum SocketError {
        UnknownSocketError = QAbstractSocket::UnknownSocketError,
        RemoteHostClosedError = QAbstractSocket::RemoteHostClosedError,
        SocketAccessError = QAbstractSocket::SocketAccessError,
        SocketResourceError = QAbstractSocket::SocketResourceError
    };
    Q_ENUM(SocketError)

    explicit QLlcpSocket(QObject *parent = 0);
    ~QLlcpSocket();

    void connectToService(QNearFieldTarget *target, const QString &serviceUri);
    void disconnectFromService();

    void close();

    bool bind(quint8 port);

    bool hasPendingDatagrams() const;
    qint64 pendingDatagramSize() const;

    qint64 writeDatagram(const char *data, qint64 size);
    qint64 writeDatagram(const QByteArray &datagram);

    qint64 readDatagram(char *data, qint64 maxSize,
                        QNearFieldTarget **target = 0, quint8 *port = 0);
    qint64 writeDatagram(const char *data, qint64 size,
                         QNearFieldTarget *target, quint8 port);
    qint64 writeDatagram(const QByteArray &datagram, QNearFieldTarget *target, quint8 port);

    SocketError error() const;
    SocketState state() const;

    qint64 bytesAvailable() const;
    bool canReadLine() const;

    bool waitForReadyRead(int msecs = 30000);
    bool waitForBytesWritten(int msecs = 30000);
    virtual bool waitForConnected(int msecs = 30000);
    virtual bool waitForDisconnected(int msecs = 30000);
    bool isSequential() const;

signals:
    void connected();
    void disconnected();
    void error(QLlcpSocket::SocketError socketError);
    void stateChanged(QLlcpSocket::SocketState socketState);

protected:
    qint64 readData(char *data, qint64 maxlen);
    qint64 writeData(const char *data, qint64 len);

private:
    QLlcpSocket(QLlcpSocketPrivate *d, QObject *parent);

    QLlcpSocketPrivate *d_ptr;
};

QT_END_NAMESPACE

#endif // QLLCPSOCKET_H
