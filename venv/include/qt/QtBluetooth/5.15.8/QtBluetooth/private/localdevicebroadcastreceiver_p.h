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

#include "android/androidbroadcastreceiver_p.h"
#include <QtBluetooth/QBluetoothAddress>
#include <QtBluetooth/QBluetoothLocalDevice>

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

#ifndef LOCALDEVICEBROADCASTRECEIVER_H
#define LOCALDEVICEBROADCASTRECEIVER_H

QT_BEGIN_NAMESPACE

class LocalDeviceBroadcastReceiver : public AndroidBroadcastReceiver
{
    Q_OBJECT
public:
    explicit LocalDeviceBroadcastReceiver(QObject *parent = nullptr);
    virtual ~LocalDeviceBroadcastReceiver() {}
    virtual void onReceive(JNIEnv *env, jobject context, jobject intent);
    virtual void onReceiveLeScan(JNIEnv *, jobject, jint, jbyteArray) {}
    bool pairingConfirmation(bool accept);

signals:
    void hostModeStateChanged(QBluetoothLocalDevice::HostMode state);
    void pairingStateChanged(const QBluetoothAddress &address, QBluetoothLocalDevice::Pairing pairing);
    void connectDeviceChanges(const QBluetoothAddress &address, bool isConnectEvent);
    void pairingDisplayConfirmation(const QBluetoothAddress &address, const QString& pin);
    void pairingDisplayPinCode(const QBluetoothAddress &address, const QString& pin);
private:
    int previousScanMode;
    QAndroidJniObject pairingDevice;

    int bondingModePreset[3];
    int hostModePreset[3];
};

QT_END_NAMESPACE

#endif // LOCALDEVICEBROADCASTRECEIVER_H
