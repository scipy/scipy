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

#ifndef HCIMANAGER_P_H
#define HCIMANAGER_P_H

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

#include <QObject>
#include <QtCore/QSet>
#include <QtCore/QSocketNotifier>
#include <QtBluetooth/QBluetoothAddress>
#include <QVector>
#include "bluez/bluez_data_p.h"

QT_BEGIN_NAMESPACE

class QLowEnergyConnectionParameters;

class HciManager : public QObject
{
    Q_OBJECT
public:
    enum HciEvent {
        EncryptChangeEvent = EVT_ENCRYPT_CHANGE,
        CommandCompleteEvent = EVT_CMD_COMPLETE,
        LeMetaEvent = 0x3e,
    };

    explicit HciManager(const QBluetoothAddress &deviceAdapter, QObject *parent = nullptr);
    ~HciManager();

    bool isValid() const;
    bool monitorEvent(HciManager::HciEvent event);
    bool monitorAclPackets();
    bool sendCommand(OpCodeGroupField ogf, OpCodeCommandField ocf, const QByteArray &parameters);

    void stopEvents();
    QBluetoothAddress addressForConnectionHandle(quint16 handle) const;

    // active connections
    QVector<quint16> activeLowEnergyConnections() const;

    bool sendConnectionUpdateCommand(quint16 handle, const QLowEnergyConnectionParameters &params);
    bool sendConnectionParameterUpdateRequest(quint16 handle,
                                              const QLowEnergyConnectionParameters &params);

signals:
    void encryptionChangedEvent(const QBluetoothAddress &address, bool wasSuccess);
    void commandCompleted(quint16 opCode, quint8 status, const QByteArray &data);
    void connectionComplete(quint16 handle);
    void connectionUpdate(quint16 handle, const QLowEnergyConnectionParameters &parameters);
    void signatureResolvingKeyReceived(quint16 connHandle, bool remoteKey, const quint128 &csrk);

private slots:
    void _q_readNotify();

private:
    int hciForAddress(const QBluetoothAddress &deviceAdapter);
    void handleHciEventPacket(const quint8 *data, int size);
    void handleHciAclPacket(const quint8 *data, int size);
    void handleLeMetaEvent(const quint8 *data);

    int hciSocket;
    int hciDev;
    quint8 sigPacketIdentifier = 0;
    QSocketNotifier *notifier = nullptr;
    QSet<HciManager::HciEvent> runningEvents;
};

QT_END_NAMESPACE

#endif // HCIMANAGER_P_H
