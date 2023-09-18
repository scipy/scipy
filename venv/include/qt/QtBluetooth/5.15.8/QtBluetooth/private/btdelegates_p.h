/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
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

#ifndef BTDELEGATES_P_H
#define BTDELEGATES_P_H

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

#include "qbluetoothdevicediscoveryagent.h"
#include "qlowenergycontroller.h"
#include "qbluetooth.h"

#include <QtCore/qsharedpointer.h>
#include <QtCore/qglobal.h>

#if defined(Q_OS_MACOS)

#include <IOKit/IOReturn.h>

#include <cstdint>

QT_BEGIN_NAMESPACE

class QLowEnergyServicePrivate;
class QBluetoothAddress;
class QByteArray;

namespace DarwinBluetooth {

class DeviceInquiryDelegate
{
public:
    virtual ~DeviceInquiryDelegate();

    virtual void inquiryFinished() = 0;
    virtual void error(IOReturn error) = 0;
    virtual void classicDeviceFound(void *ioBluetoothDevice) = 0;
};

class PairingDelegate
{
public:
    using BluetoothNumericValue = uint32_t;
    using BluetoothPasskey = BluetoothNumericValue;

    virtual ~PairingDelegate();

    virtual void connecting(void *pair) = 0;
    virtual void requestPIN(void *pair) = 0;
    virtual void requestUserConfirmation(void *pair,
                                         BluetoothNumericValue) = 0;
    virtual void passkeyNotification(void *pair,
                                     BluetoothPasskey passkey) = 0;
    virtual void error(void *pair, IOReturn errorCode) = 0;
    virtual void pairingFinished(void *pair) = 0;
};

class SDPInquiryDelegate {
public:
    virtual ~SDPInquiryDelegate();

    virtual void SDPInquiryFinished(void *ioBluetoothDevice) = 0;
    virtual void SDPInquiryError(void *ioBluetoothDevice, IOReturn errorCode) = 0;
};

// L2CAP and RFCOMM.
class ChannelDelegate
{
public:
    virtual ~ChannelDelegate();

    virtual void setChannelError(IOReturn errorCode) = 0;
    virtual void channelOpenComplete() = 0;
    virtual void channelClosed() = 0;

    virtual void readChannelData(void *data, std::size_t size) = 0;
    virtual void writeComplete() = 0;
};

class ConnectionMonitor {
public:
    virtual ~ConnectionMonitor();

    virtual void deviceConnected(const QBluetoothAddress &address) = 0;
    virtual void deviceDisconnected(const QBluetoothAddress &address) = 0;
};

class SocketListener
{
public:
    virtual ~SocketListener();

    virtual void openNotifyRFCOMM(void *rfcommChannel) = 0;
    virtual void openNotifyL2CAP(void *l2capChannel) = 0;
};


} // namespace DarwinBluetooth

QT_END_NAMESPACE

#endif // Q_OS_MACOS

#endif // DARWINBTDELEGATES_P_H
