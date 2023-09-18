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

#ifndef OSXBTCENTRALMANAGER_P_H
#define OSXBTCENTRALMANAGER_P_H

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

#include "qlowenergycontroller.h"
#include "qlowenergyservice.h"
#include "osxbtgcdtimer_p.h"
#include "qbluetoothuuid.h"
#include "osxbtutility_p.h"
#include "osxbluetooth_p.h"

#include <QtCore/qbytearray.h>
#include <QtCore/qglobal.h>
#include <QtCore/qqueue.h>
#include <QtCore/qhash.h>

#include <Foundation/Foundation.h>

@class QT_MANGLE_NAMESPACE(OSXBTCentralManager);

QT_BEGIN_NAMESPACE

class QLowEnergyServicePrivate;

namespace OSXBluetooth {

class LECBManagerNotifier;

enum CentralManagerState
{
    // QLowEnergyController already has some of these states,
    // but it's not enough and we need more special states here.
    CentralManagerIdle,
    // Required by CBCentralManager to avoid problems with dangled pointers.
    CentralManagerUpdating,
    CentralManagerConnecting,
    CentralManagerDiscovering,
    CentralManagerDisconnecting
};

// In Qt we work with handles and UUIDs. Core Bluetooth
// has NSArrays (and nested NSArrays inside servces/characteristics).
// To simplify a navigation, I need a simple way to map from a handle
// to a Core Bluetooth object. These are weak pointers,
// will probably require '__weak' with ARC.
typedef QHash<QLowEnergyHandle, CBService *> ServiceHash;
typedef QHash<QLowEnergyHandle, CBCharacteristic *> CharHash;
typedef QHash<QLowEnergyHandle, CBDescriptor *> DescHash;

// Descriptor/charactesirstic read/write requests
// - we have to serialize 'concurrent' requests.
struct LERequest {
    enum RequestType {
        CharRead,
        CharWrite,
        DescRead,
        DescWrite,
        ClientConfiguration,
        Invalid
    };

    LERequest() : type(Invalid),
                  withResponse(false),
                  handle(0)
    {}

    RequestType type;
    bool withResponse;
    QLowEnergyHandle handle;
    QByteArray value;
};

typedef QQueue<LERequest> RequestQueue;

// Core Bluetooth's write confirmation does not provide
// the updated value (characteristic or descriptor).
// But the Qt's Bluetooth API ('write with response')
// expects this updated value as a response, so we have
// to cache this write value and report it back.
// 'NSObject *' will require '__weak' with ARC.
typedef QHash<NSObject *, QByteArray> ValueHash;

}

QT_END_NAMESPACE

@interface QT_MANGLE_NAMESPACE(OSXBTCentralManager) : NSObject<CBCentralManagerDelegate,
                                                               CBPeripheralDelegate,
                                                               QT_MANGLE_NAMESPACE(GCDTimerDelegate)>
- (id)initWith:(QT_PREPEND_NAMESPACE(OSXBluetooth)::LECBManagerNotifier *)notifier;
- (void)dealloc;

- (CBPeripheral *)peripheral;

// IMPORTANT: _all_ these methods are to be executed on qt_LE_queue,
// when passing parameters - C++ objects _must_ be copied (see the controller's code).
- (void)connectToDevice:(const QT_PREPEND_NAMESPACE(QBluetoothUuid) &)aDeviceUuid;

- (void)disconnectFromDevice;

- (void)discoverServices;
- (void)discoverServiceDetails:(const QT_PREPEND_NAMESPACE(QBluetoothUuid) &)serviceUuid;

- (void)setNotifyValue:(const QT_PREPEND_NAMESPACE(QByteArray) &)value
        forCharacteristic:(QT_PREPEND_NAMESPACE(QLowEnergyHandle))charHandle
        onService:(const QT_PREPEND_NAMESPACE(QBluetoothUuid) &)serviceUuid;

- (void)readCharacteristic:(QT_PREPEND_NAMESPACE(QLowEnergyHandle))charHandle
        onService:(const QT_PREPEND_NAMESPACE(QBluetoothUuid) &)serviceUuid;

- (void)write:(const QT_PREPEND_NAMESPACE(QByteArray) &)value
        charHandle:(QT_PREPEND_NAMESPACE(QLowEnergyHandle))charHandle
        onService:(const QT_PREPEND_NAMESPACE(QBluetoothUuid) &)serviceUuid
        withResponse:(bool)writeWithResponse;

- (void)readDescriptor:(QT_PREPEND_NAMESPACE(QLowEnergyHandle))descHandle
        onService:(const QT_PREPEND_NAMESPACE(QBluetoothUuid) &)serviceUuid;

- (void)write:(const QT_PREPEND_NAMESPACE(QByteArray) &)value
        descHandle:(QT_PREPEND_NAMESPACE(QLowEnergyHandle))descHandle
        onService:(const QT_PREPEND_NAMESPACE(QBluetoothUuid) &)serviceUuid;

- (void)detach;

@end

#endif
