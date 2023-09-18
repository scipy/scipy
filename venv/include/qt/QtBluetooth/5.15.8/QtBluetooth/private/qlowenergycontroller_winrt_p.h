/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
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

#ifndef QLOWENERGYCONTROLLERPRIVATEWINRT_P_H
#define QLOWENERGYCONTROLLERPRIVATEWINRT_P_H

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

#include <qglobal.h>
#include <QtCore/QQueue>
#include <QtCore/QVector>
#include <QtBluetooth/qbluetooth.h>
#include <QtBluetooth/qlowenergycharacteristic.h>
#include <QtBluetooth/qlowenergyservicedata.h>
#include "qlowenergycontroller.h"
#include "qlowenergycontrollerbase_p.h"

#include <wrl.h>
#include <windows.devices.bluetooth.h>

#include <functional>

QT_BEGIN_NAMESPACE

class QLowEnergyServiceData;
class QTimer;
class QWinRTLowEnergyServiceHandler;

extern void registerQLowEnergyControllerMetaType();

class QLowEnergyControllerPrivateWinRT final : public QLowEnergyControllerPrivate
{
    Q_OBJECT
public:
    QLowEnergyControllerPrivateWinRT();
    ~QLowEnergyControllerPrivateWinRT() override;

    void init() override;

    void connectToDevice() override;
    void disconnectFromDevice() override;

    void discoverServices() override;
    void discoverServiceDetails(const QBluetoothUuid &service) override;

    void startAdvertising(const QLowEnergyAdvertisingParameters &params,
                          const QLowEnergyAdvertisingData &advertisingData,
                          const QLowEnergyAdvertisingData &scanResponseData) override;
    void stopAdvertising() override;

    void requestConnectionUpdate(const QLowEnergyConnectionParameters &params) override;

    // read data
    void readCharacteristic(const QSharedPointer<QLowEnergyServicePrivate> service,
                            const QLowEnergyHandle charHandle) override;
    void readDescriptor(const QSharedPointer<QLowEnergyServicePrivate> service,
                        const QLowEnergyHandle charHandle,
                        const QLowEnergyHandle descriptorHandle) override;

    // write data
    void writeCharacteristic(const QSharedPointer<QLowEnergyServicePrivate> service,
                             const QLowEnergyHandle charHandle,
                             const QByteArray &newValue, QLowEnergyService::WriteMode mode) override;
    void writeDescriptor(const QSharedPointer<QLowEnergyServicePrivate> service,
                         const QLowEnergyHandle charHandle,
                         const QLowEnergyHandle descriptorHandle,
                         const QByteArray &newValue) override;

    void addToGenericAttributeList(const QLowEnergyServiceData &service,
                                   QLowEnergyHandle startHandle) override;

signals:
    void characteristicChanged(quint16 charHandle, const QByteArray &data);

private slots:
    void handleCharacteristicChanged(quint16 charHandle, const QByteArray &data);
    void doConnectToDevice();

private:
    Microsoft::WRL::ComPtr<ABI::Windows::Devices::Bluetooth::IBluetoothLEDevice> mDevice;
    EventRegistrationToken mStatusChangedToken;
    struct ValueChangedEntry {
        ValueChangedEntry() {}
        ValueChangedEntry(Microsoft::WRL::ComPtr<ABI::Windows::Devices::Bluetooth::GenericAttributeProfile::IGattCharacteristic> c,
                          EventRegistrationToken t)
            : characteristic(c)
            , token(t)
        {
        }

        Microsoft::WRL::ComPtr<ABI::Windows::Devices::Bluetooth::GenericAttributeProfile::IGattCharacteristic> characteristic;
        EventRegistrationToken token;
    };
    QVector<ValueChangedEntry> mValueChangedTokens;

    Microsoft::WRL::ComPtr<ABI::Windows::Devices::Bluetooth::GenericAttributeProfile::IGattDeviceService> getNativeService(const QBluetoothUuid &serviceUuid);
    Microsoft::WRL::ComPtr<ABI::Windows::Devices::Bluetooth::GenericAttributeProfile::IGattCharacteristic> getNativeCharacteristic(const QBluetoothUuid &serviceUuid, const QBluetoothUuid &charUuid);

    void registerForValueChanges(const QBluetoothUuid &serviceUuid, const QBluetoothUuid &charUuid);
    void unregisterFromValueChanges();

    void obtainIncludedServices(QSharedPointer<QLowEnergyServicePrivate> servicePointer,
        Microsoft::WRL::ComPtr<ABI::Windows::Devices::Bluetooth::GenericAttributeProfile::IGattDeviceService> nativeService);

};

QT_END_NAMESPACE

#endif // QLOWENERGYCONTROLLERPRIVATEWINRT_P_H
