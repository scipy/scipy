/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Copyright (C) 2014 Denis Shienkov <denis.shienkov@gmail.com>
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

#ifndef QBLUETOOTHLOCALDEVICE_P_H
#define QBLUETOOTHLOCALDEVICE_P_H

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

#include <QtBluetooth/qtbluetoothglobal.h>

#include "qbluetoothlocaldevice.h"

#if QT_CONFIG(bluez)
#include <QObject>
#include <QDBusContext>
#include <QDBusObjectPath>
#include <QDBusMessage>
#include <QSet>
#include "bluez/bluez5_helper_p.h"

class OrgBluezAdapterInterface;
class OrgBluezAdapter1Interface;
class OrgFreedesktopDBusPropertiesInterface;
class OrgFreedesktopDBusObjectManagerInterface;
class OrgBluezAgentAdaptor;
class OrgBluezDeviceInterface;
class OrgBluezDevice1Interface;
class OrgBluezManagerInterface;

QT_BEGIN_NAMESPACE
class QDBusPendingCallWatcher;
QT_END_NAMESPACE
#endif

#ifdef QT_ANDROID_BLUETOOTH
#include <jni.h>
#include <QtAndroidExtras/QAndroidJniEnvironment>
#include <QtAndroidExtras/QAndroidJniObject>
#include <QtCore/QPair>
#endif

#ifdef QT_WINRT_BLUETOOTH
#include <wrl.h>

namespace ABI {
    namespace Windows {
        namespace Devices {
            namespace Bluetooth {
                struct IBluetoothDeviceStatics;
                struct IBluetoothLEDeviceStatics;
            }
        }
    }
}
#endif

QT_BEGIN_NAMESPACE

extern void registerQBluetoothLocalDeviceMetaType();

class QBluetoothAddress;

#ifdef QT_ANDROID_BLUETOOTH
class LocalDeviceBroadcastReceiver;
class QBluetoothLocalDevicePrivate : public QObject
{
    Q_OBJECT
public:
    QBluetoothLocalDevicePrivate(
        QBluetoothLocalDevice *q, const QBluetoothAddress &address = QBluetoothAddress());
    ~QBluetoothLocalDevicePrivate();

    QAndroidJniObject *adapter();
    void initialize(const QBluetoothAddress &address);
    static bool startDiscovery();
    static bool cancelDiscovery();
    static bool isDiscovering();
    bool isValid() const;

private slots:
    void processHostModeChange(QBluetoothLocalDevice::HostMode newMode);
    void processPairingStateChanged(const QBluetoothAddress &address,
                                    QBluetoothLocalDevice::Pairing pairing);
    void processConnectDeviceChanges(const QBluetoothAddress &address, bool isConnectEvent);
    void processDisplayConfirmation(const QBluetoothAddress &address, const QString &pin);

private:
    QBluetoothLocalDevice *q_ptr;
    QAndroidJniObject *obj = nullptr;

    int pendingPairing(const QBluetoothAddress &address);

public:
    LocalDeviceBroadcastReceiver *receiver;
    bool pendingHostModeTransition = false;
    QList<QPair<QBluetoothAddress, bool> > pendingPairings;

    QList<QBluetoothAddress> connectedDevices;
};

#elif QT_CONFIG(bluez)
class QBluetoothLocalDevicePrivate : public QObject, protected QDBusContext
{
    Q_OBJECT
    Q_DECLARE_PUBLIC(QBluetoothLocalDevice)
public:
    QBluetoothLocalDevicePrivate(QBluetoothLocalDevice *q,
                                 QBluetoothAddress localAddress = QBluetoothAddress());
    ~QBluetoothLocalDevicePrivate();

    QSet<OrgBluezDeviceInterface *> devices;
    QSet<QBluetoothAddress> connectedDevicesSet;
    OrgBluezAdapterInterface *adapter = nullptr; //Bluez 4
    OrgBluezAdapter1Interface *adapterBluez5 = nullptr; //Bluez 5
    OrgFreedesktopDBusPropertiesInterface *adapterProperties = nullptr; //Bluez 5
    OrgFreedesktopDBusObjectManagerInterface *managerBluez5 = nullptr; //Bluez 5
    QMap<QString, OrgFreedesktopDBusPropertiesInterface *> deviceChangeMonitors; //Bluez 5
    OrgBluezAgentAdaptor *agent = nullptr;
    OrgBluezManagerInterface *manager = nullptr;

    QList<QBluetoothAddress> connectedDevices() const;

    QString agent_path;
    QBluetoothAddress localAddress;
    QBluetoothAddress address;
    QBluetoothLocalDevice::Pairing pairing;
    OrgBluezDevice1Interface *pairingTarget = nullptr;
    QTimer *pairingDiscoveryTimer = nullptr;
    QBluetoothLocalDevice::HostMode currentMode;
    int pendingHostModeChange;

public slots:
    void Authorize(const QDBusObjectPath &in0, const QString &in1);
    void Cancel();
    void ConfirmModeChange(const QString &in0);
    void DisplayPasskey(const QDBusObjectPath &in0, uint in1, uchar in2);
    void Release();
    uint RequestPasskey(const QDBusObjectPath &in0);

    void RequestConfirmation(const QDBusObjectPath &in0, uint in1);
    QString RequestPinCode(const QDBusObjectPath &in0);

    void pairingCompleted(QDBusPendingCallWatcher *);

    void PropertyChanged(QString, QDBusVariant);
    void _q_deviceCreated(const QDBusObjectPath &device);
    void _q_deviceRemoved(const QDBusObjectPath &device);
    void _q_devicePropertyChanged(const QString &property, const QDBusVariant &value);
    bool isValid() const;
    void adapterRemoved(const QDBusObjectPath &device);

    void requestPairingBluez5(const QBluetoothAddress &address,
                              QBluetoothLocalDevice::Pairing targetPairing);

private Q_SLOTS:
    void PropertiesChanged(const QString &interface,
                           const QVariantMap &changed_properties,
                           const QStringList &invalidated_properties,
                           const QDBusMessage &signal);
    void InterfacesAdded(const QDBusObjectPath &object_path,
                         InterfaceList interfaces_and_properties);
    void InterfacesRemoved(const QDBusObjectPath &object_path,
                           const QStringList &interfaces);
    void processPairingBluez5(const QString &objectPath,
                              QBluetoothLocalDevice::Pairing target);
    void pairingDiscoveryTimedOut();

private:
    void createCache();
    void connectDeviceChanges();

    QDBusMessage msgConfirmation;
    QDBusConnection *msgConnection = nullptr;
    QString deviceAdapterPath;

    QBluetoothLocalDevice *q_ptr;

    void initializeAdapter();
    void initializeAdapterBluez5();
};

#elif defined(QT_WIN_BLUETOOTH)

class QBluetoothLocalDevicePrivate : public QObject
{
    Q_OBJECT
    Q_DECLARE_PUBLIC(QBluetoothLocalDevice)
public:
    QBluetoothLocalDevicePrivate(QBluetoothLocalDevice *q,
                                 const QBluetoothAddress &address = QBluetoothAddress());

    ~QBluetoothLocalDevicePrivate();
    bool isValid() const;
    void initialize(const QBluetoothAddress &address);

    static QList<QBluetoothHostInfo> localAdapters();

    QBluetoothAddress deviceAddress;
    QString deviceName;
    bool deviceValid;
private:
    QBluetoothLocalDevice *q_ptr;
};
#elif defined(QT_WINRT_BLUETOOTH)
class QBluetoothLocalDevicePrivate : public QObject
{
    Q_DECLARE_PUBLIC(QBluetoothLocalDevice)
public:
    QBluetoothLocalDevicePrivate(QBluetoothLocalDevice *q,
                                 QBluetoothAddress = QBluetoothAddress());
    ~QBluetoothLocalDevicePrivate();

    bool isValid() const;

private:
    QBluetoothLocalDevice *q_ptr;
    Microsoft::WRL::ComPtr<ABI::Windows::Devices::Bluetooth::IBluetoothDeviceStatics> mStatics;
    Microsoft::WRL::ComPtr<ABI::Windows::Devices::Bluetooth::IBluetoothLEDeviceStatics> mLEStatics;
};
#elif !defined(QT_OSX_BLUETOOTH) // dummy backend
class QBluetoothLocalDevicePrivate : public QObject
{
public:
    QBluetoothLocalDevicePrivate(QBluetoothLocalDevice * = nullptr,
                                 QBluetoothAddress = QBluetoothAddress())
    {
    }

    bool isValid() const
    {
        return false;
    }
};
#endif

QT_END_NAMESPACE

#endif // QBLUETOOTHLOCALDEVICE_P_H
