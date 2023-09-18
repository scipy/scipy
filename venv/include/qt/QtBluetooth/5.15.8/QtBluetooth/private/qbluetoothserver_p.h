/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
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

#ifndef QBLUETOOTHSERVER_P_H
#define QBLUETOOTHSERVER_P_H

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

#include <QtGlobal>
#include <QList>
#include <QtBluetooth/QBluetoothSocket>
#include "qbluetoothserver.h"
#include "qbluetooth.h"

#if QT_CONFIG(bluez) || defined(QT_WIN_BLUETOOTH)
QT_FORWARD_DECLARE_CLASS(QSocketNotifier)
#endif

#ifdef QT_ANDROID_BLUETOOTH
#include <QtAndroidExtras/QAndroidJniEnvironment>
#include <QtAndroidExtras/QAndroidJniObject>
#include <QtBluetooth/QBluetoothUuid>

class ServerAcceptanceThread;
#endif

#ifdef QT_WINRT_BLUETOOTH
#include <QtCore/QMutex>

#include <wrl.h>
// No forward declares because QBluetoothServerPrivate::listener does not work with them
#include <windows.networking.sockets.h>
#endif

#ifdef QT_OSX_BLUETOOTH

#include "osx/btdelegates_p.h"
#include "osx/btraii_p.h"

#include <QtCore/qvector.h>

#endif // QT_OSX_BLUETOOTH

QT_BEGIN_NAMESPACE

class QBluetoothAddress;
class QBluetoothSocket;
class QBluetoothServer;

class QBluetoothServerPrivate
#ifdef QT_OSX_BLUETOOTH
        : public DarwinBluetooth::SocketListener
#endif
{
    Q_DECLARE_PUBLIC(QBluetoothServer)

public:
    QBluetoothServerPrivate(QBluetoothServiceInfo::Protocol serverType, QBluetoothServer *parent);
    ~QBluetoothServerPrivate();

#if QT_CONFIG(bluez)
    void _q_newConnection();
    void setSocketSecurityLevel(QBluetooth::SecurityFlags requestedSecLevel, int *errnoCode);
    QBluetooth::SecurityFlags socketSecurityLevel() const;
    static QBluetoothSocket *createSocketForServer(
                QBluetoothServiceInfo::Protocol socketType = QBluetoothServiceInfo::RfcommProtocol);
#endif
#if defined(QT_WIN_BLUETOOTH)
    void _q_newConnection();
#endif

public:
    QBluetoothSocket *socket;

    int maxPendingConnections;
    QBluetooth::SecurityFlags securityFlags;
    QBluetoothServiceInfo::Protocol serverType;

protected:
    QBluetoothServer *q_ptr;

private:
    QBluetoothServer::Error m_lastError;
#if QT_CONFIG(bluez) || defined(QT_WIN_BLUETOOTH)
    QSocketNotifier *socketNotifier = nullptr;
#elif defined(QT_ANDROID_BLUETOOTH)
    ServerAcceptanceThread *thread;
    QString m_serviceName;
    QBluetoothUuid m_uuid;
public:
    bool isListening() const;
    bool initiateActiveListening(const QBluetoothUuid& uuid, const QString &serviceName);
    bool deactivateActiveListening();
#elif defined(QT_WINRT_BLUETOOTH)
    EventRegistrationToken connectionToken {-1};

    mutable QMutex pendingConnectionsMutex;
    QVector<Microsoft::WRL::ComPtr<ABI::Windows::Networking::Sockets::IStreamSocket>> pendingConnections;

    Microsoft::WRL::ComPtr<ABI::Windows::Networking::Sockets::IStreamSocketListener> socketListener;
    HRESULT handleClientConnection(ABI::Windows::Networking::Sockets::IStreamSocketListener *listener,
                                   ABI::Windows::Networking::Sockets::IStreamSocketListenerConnectionReceivedEventArgs *args);

public:
    bool isListening() const;
    Microsoft::WRL::ComPtr<ABI::Windows::Networking::Sockets::IStreamSocketListener> listener() { return socketListener; }
    bool initiateActiveListening(const QString &serviceName);
    bool deactivateActiveListening();
#endif

#ifdef QT_OSX_BLUETOOTH

public:

    friend class QBluetoothServer;
    friend class QBluetoothServiceInfoPrivate;

private:
    bool startListener(quint16 realPort);
    void stopListener();
    bool isListening() const;

    // SocketListener (delegate):
    void openNotifyRFCOMM(void *channel) override;
    void openNotifyL2CAP(void *channel) override;

    // Either a "temporary" channelID/PSM assigned by QBluetoothServer::listen,
    // or a real channelID/PSM returned by IOBluetooth after we've registered
    // a service.
    quint16 port;

    DarwinBluetooth::StrongReference listener;

    // These static functions below
    // deal with differences between bluetooth sockets
    // (bluez and QtBluetooth's API) and IOBluetooth, where it's not possible
    // to have a real PSM/channelID _before_ a service is registered,
    // the solution - "fake" ports.
    // These functions require external locking - using channelMapMutex.
    static QMutex &channelMapMutex();

    static bool channelIsBusy(quint16 channelID);
    static quint16 findFreeChannel();

    static bool psmIsBusy(quint16 psm);
    static quint16 findFreePSM();

    static void registerServer(QBluetoothServerPrivate *server, quint16 port);
    static QBluetoothServerPrivate *registeredServer(quint16 port, QBluetoothServiceInfo::Protocol protocol);
    static void unregisterServer(QBluetoothServerPrivate *server);

    using PendingConnection = DarwinBluetooth::StrongReference;
    QVector<PendingConnection> pendingConnections;

#endif // QT_OSX_BLUETOOTH
};

QT_END_NAMESPACE

#endif
