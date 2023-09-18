/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtBluetooth module of the Qt Toolkit.
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

#ifndef QBLUETOOTHSOCKET_WINRT_P_H
#define QBLUETOOTHSOCKET_WINRT_P_H

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

#include "qbluetoothsocket.h"
#include "qbluetoothsocketbase_p.h"
#include <QtGlobal>

QT_FORWARD_DECLARE_CLASS(SocketWorker)

namespace ABI {
    namespace Windows {
        namespace  Networking {
            struct IHostName;
        }
    }
}

namespace Microsoft {
    namespace WRL {
        template <typename T> class ComPtr;
    }
}
QT_BEGIN_NAMESPACE

class QBluetoothSocketPrivateWinRT final: public QBluetoothSocketBasePrivate
{
    Q_OBJECT
    friend class QBluetoothServerPrivate;

public:
    QBluetoothSocketPrivateWinRT();
    ~QBluetoothSocketPrivateWinRT() override;

    void connectToServiceHelper(const QBluetoothAddress &address,
                          quint16 port,
                          QIODevice::OpenMode openMode) override;

    void connectToService(const QBluetoothServiceInfo &service,
                          QIODevice::OpenMode openMode) override;
    void connectToService(const QBluetoothAddress &address, const QBluetoothUuid &uuid,
                          QIODevice::OpenMode openMode) override;
    void connectToService(const QBluetoothAddress &address, quint16 port,
                          QIODevice::OpenMode openMode) override;

    bool ensureNativeSocket(QBluetoothServiceInfo::Protocol type) override;

    QString localName() const override;
    QBluetoothAddress localAddress() const override;
    quint16 localPort() const override;

    QString peerName() const override;
    QBluetoothAddress peerAddress() const override;
    quint16 peerPort() const override;

    void abort() override;
    void close() override;

    qint64 writeData(const char *data, qint64 maxSize) override;
    qint64 readData(char *data, qint64 maxSize) override;

    bool setSocketDescriptor(Microsoft::WRL::ComPtr<ABI::Windows::Networking::Sockets::IStreamSocket> socket,
                             QBluetoothServiceInfo::Protocol socketType,
                             QBluetoothSocket::SocketState socketState = QBluetoothSocket::ConnectedState,
                             QBluetoothSocket::OpenMode openMode = QBluetoothSocket::ReadWrite) override;

    bool setSocketDescriptor(int socketDescriptor, QBluetoothServiceInfo::Protocol socketType,
                             QBluetoothSocket::SocketState socketState = QBluetoothSocket::ConnectedState,
                             QBluetoothSocket::OpenMode openMode = QBluetoothSocket::ReadWrite) override;

    qint64 bytesAvailable() const override;
    bool canReadLine() const override;
    qint64 bytesToWrite() const override;

    SocketWorker *m_worker;

    Microsoft::WRL::ComPtr<ABI::Windows::Networking::Sockets::IStreamSocket> m_socketObject;
    Microsoft::WRL::ComPtr<ABI::Windows::Foundation::IAsyncAction> m_connectOp;

    QMutex m_readMutex;

    // Protected by m_readMutex. Written in addToPendingData (native callback)
    QVector<QByteArray> m_pendingData;

    Q_INVOKABLE void addToPendingData(const QVector<QByteArray> &data);

private slots:
    void handleNewData(const QVector<QByteArray> &data);
    void handleError(QBluetoothSocket::SocketError error);

private:
    void connectToService(Microsoft::WRL::ComPtr<ABI::Windows::Networking::IHostName> hostName,
                          const QString &serviceName, QIODevice::OpenMode openMode);
    HRESULT handleConnectOpFinished(ABI::Windows::Foundation::IAsyncAction *action,
                                    ABI::Windows::Foundation::AsyncStatus status);

    QIODevice::OpenMode requestedOpenMode = QIODevice::NotOpen;
};

QT_END_NAMESPACE

#endif
