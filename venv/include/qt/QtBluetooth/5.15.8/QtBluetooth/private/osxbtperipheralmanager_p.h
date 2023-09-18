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

#ifndef OSXBTPERIPHERALMANAGER_P_H
#define OSXBTPERIPHERALMANAGER_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API. It exists for the convenience
// of internal files. This header file may change from version to version
// without notice, or even be removed.
//
// We mean it.
//

#include "osxbtutility_p.h"

#include "qlowenergyadvertisingparameters.h"
#include "qlowenergyserviceprivate_p.h"
#include "qbluetoothuuid.h"
#include "qbluetooth.h"

#include <QtCore/qsharedpointer.h>
#include <QtCore/qbytearray.h>
#include <QtCore/qsysinfo.h>
#include <QtCore/qglobal.h>
#include <QtCore/qpair.h>
#include <QtCore/qmap.h>

#include <vector>
#include <deque>
#include <map>

#include <Foundation/Foundation.h>

#include "osxbluetooth_p.h"

QT_BEGIN_NAMESPACE

class QLowEnergyServiceData;

namespace OSXBluetooth
{

class LECBManagerNotifier;

}

QT_END_NAMESPACE


// Exposing names in a header is ugly, but constant QT_PREPEND_NAMESPACE is even worse ...
// After all, this header is to be included only in its own and controller's *.mm files.

QT_USE_NAMESPACE

using namespace OSXBluetooth;


template<class Type>
using GenericLEMap = QMap<QLowEnergyHandle, Type>;

enum class PeripheralState
{
    idle,
    waitingForPowerOn,
    advertising,
    connected
};

struct UpdateRequest
{
    UpdateRequest() = default;
    UpdateRequest(QLowEnergyHandle handle, const ObjCStrongReference<NSData> &val)
        : charHandle(handle),
          value(val)
    {
    }

    QLowEnergyHandle charHandle = {};
    ObjCStrongReference<NSData> value;
};

using ValueRange = QPair<NSUInteger, NSUInteger>;

@interface QT_MANGLE_NAMESPACE(OSXBTPeripheralManager) : NSObject<CBPeripheralManagerDelegate>

- (id)initWith:(LECBManagerNotifier *)notifier;
- (void)dealloc;

- (QSharedPointer<QLowEnergyServicePrivate>)addService:(const QLowEnergyServiceData &)data;
- (void) setParameters:(const QLowEnergyAdvertisingParameters &)parameters
         data:(const QLowEnergyAdvertisingData &)data
         scanResponse:(const QLowEnergyAdvertisingData &)scanResponse;

// To be executed on the Qt's special BTLE dispatch queue.
- (void)startAdvertising;
- (void)stopAdvertising;
- (void)detach;

- (void)write:(const QByteArray &)value
        charHandle:(QLowEnergyHandle)charHandle;


// CBPeripheralManagerDelegate's callbacks (BTLE queue).
- (void)peripheralManagerDidUpdateState:(CBPeripheralManager *)peripheral;
- (void)peripheralManager:(CBPeripheralManager *)peripheral
        willRestoreState:(NSDictionary *)dict;
- (void)peripheralManagerDidStartAdvertising:(CBPeripheralManager *)peripheral
        error:(NSError *)error;
- (void)peripheralManager:(CBPeripheralManager *)peripheral
        didAddService:(CBService *)service error:(NSError *)error;
- (void)peripheralManager:(CBPeripheralManager *)peripheral central:(CBCentral *)central
        didSubscribeToCharacteristic:(CBCharacteristic *)characteristic;
- (void)peripheralManager:(CBPeripheralManager *)peripheral central:(CBCentral *)central
        didUnsubscribeFromCharacteristic:(CBCharacteristic *)characteristic;
- (void)peripheralManager:(CBPeripheralManager *)peripheral
        didReceiveReadRequest:(CBATTRequest *)request;
- (void)peripheralManager:(CBPeripheralManager *)peripheral
        didReceiveWriteRequests:(NSArray *)requests;
- (void)peripheralManagerIsReadyToUpdateSubscribers:(CBPeripheralManager *)peripheral;

@end

#endif
