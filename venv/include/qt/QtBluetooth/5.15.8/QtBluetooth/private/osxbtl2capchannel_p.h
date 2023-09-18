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

#ifndef OSXBTL2CAPCHANNEL_P_H
#define OSXBTL2CAPCHANNEL_P_H

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

#include "osxbluetooth_p.h"

#include <QtCore/qglobal.h>

#include <Foundation/Foundation.h>

#include <cstddef>

QT_BEGIN_NAMESPACE

class QBluetoothAddress;

namespace DarwinBluetooth {

class ChannelDelegate;

}

QT_END_NAMESPACE

@interface QT_MANGLE_NAMESPACE(OSXBTL2CAPChannel) : NSObject<IOBluetoothL2CAPChannelDelegate>

- (id)initWithDelegate:(QT_PREPEND_NAMESPACE(DarwinBluetooth)::ChannelDelegate *)aDelegate;
- (id)initWithDelegate:(QT_PREPEND_NAMESPACE(DarwinBluetooth)::ChannelDelegate *)aDelegate
      channel:(IOBluetoothL2CAPChannel *)aChannel;

- (void)dealloc;

// Async. connection (connect can be called only once).
- (IOReturn)connectAsyncToDevice:(const QT_PREPEND_NAMESPACE(QBluetoothAddress) &)address
            withPSM:(BluetoothL2CAPChannelID)psm;

// IOBluetoothL2CAPChannelDelegate:
- (void)l2capChannelData:(IOBluetoothL2CAPChannel*)l2capChannel
        data:(void *)dataPointer length:(size_t)dataLength;
- (void)l2capChannelOpenComplete:(IOBluetoothL2CAPChannel*)
        l2capChannel status:(IOReturn)error;
- (void)l2capChannelClosed:(IOBluetoothL2CAPChannel*)l2capChannel;
- (void)l2capChannelReconfigured:(IOBluetoothL2CAPChannel*)l2capChannel;
- (void)l2capChannelWriteComplete:(IOBluetoothL2CAPChannel*)l2capChannel
        refcon:(void*)refcon status:(IOReturn)error;
- (void)l2capChannelQueueSpaceAvailable:(IOBluetoothL2CAPChannel*)l2capChannel;

//
- (BluetoothL2CAPPSM)getPSM;
- (BluetoothDeviceAddress)peerAddress;
- (NSString *)peerName;
- (bool)isConnected;

// Writes the given data synchronously over the target L2CAP channel to the remote
// device.
// The length of the data may not exceed the L2CAP channel's outgoing MTU.
// This method will block until the data has been successfully sent to the
// hardware for transmission (or an error occurs).
- (IOReturn) writeSync:(void*)data length:(UInt16)length;

// The length of the data may not exceed the L2CAP channel's outgoing MTU.
// When the data has been successfully passed to the hardware to be transmitted,
// the delegate method -l2capChannelWriteComplete:refcon:status: will be called.
// Returns kIOReturnSuccess if the data was buffered successfully.
- (IOReturn) writeAsync:(void*)data length:(UInt16)length;

@end

#endif
