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

#ifndef QBLUETOOTHDEVICEDISCOVERYAGENT_P_H
#define QBLUETOOTHDEVICEDISCOVERYAGENT_P_H

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
#ifdef QT_ANDROID_BLUETOOTH
#include <QtAndroidExtras/QAndroidJniObject>
#include "android/devicediscoverybroadcastreceiver_p.h"
#include <QtCore/QTimer>
#endif

#ifdef Q_OS_DARWIN
#include "osx/btdelegates_p.h"
#include "osx/btraii_p.h"
#endif // Q_OS_DARWIN

#include <QtCore/QVariantMap>

#include <QtBluetooth/QBluetoothAddress>
#include <QtBluetooth/QBluetoothLocalDevice>

#if QT_CONFIG(bluez)
#include "bluez/bluez5_helper_p.h"

class OrgBluezManagerInterface;
class OrgBluezAdapterInterface;
class OrgFreedesktopDBusObjectManagerInterface;
class OrgFreedesktopDBusPropertiesInterface;
class OrgBluezAdapter1Interface;
class OrgBluezDevice1Interface;

QT_BEGIN_NAMESPACE
class QDBusVariant;
QT_END_NAMESPACE
#endif

#ifdef QT_WIN_BLUETOOTH
QT_BEGIN_NAMESPACE
class QThread;

class ThreadWorkerDeviceDiscovery : public QObject
{
    Q_OBJECT
signals:
    void discoveryCompleted(const QVariant res);
};

QT_END_NAMESPACE

#elif defined(QT_WINRT_BLUETOOTH)
#include <QtCore/QPointer>
#include <QtCore/QTimer>

using ManufacturerData = QHash<quint16, QByteArray>;
Q_DECLARE_METATYPE(ManufacturerData)
#endif

QT_BEGIN_NAMESPACE

#ifdef QT_WINRT_BLUETOOTH
class QWinRTBluetoothDeviceDiscoveryWorker;
#endif

class QBluetoothDeviceDiscoveryAgentPrivate
#if defined(QT_ANDROID_BLUETOOTH) || defined(QT_WINRT_BLUETOOTH) || defined(QT_WIN_BLUETOOTH) \
            || defined(Q_OS_DARWIN)
    : public QObject
#if defined(Q_OS_MACOS)
    , public DarwinBluetooth::DeviceInquiryDelegate
#endif // Q_OS_MACOS
{
    Q_OBJECT
#else // BlueZ
{
#endif
    Q_DECLARE_PUBLIC(QBluetoothDeviceDiscoveryAgent)
public:
    QBluetoothDeviceDiscoveryAgentPrivate(
            const QBluetoothAddress &deviceAdapter,
            QBluetoothDeviceDiscoveryAgent *parent);
    ~QBluetoothDeviceDiscoveryAgentPrivate();

    void start(QBluetoothDeviceDiscoveryAgent::DiscoveryMethods methods);
    void stop();
    bool isActive() const;

#if QT_CONFIG(bluez)
    void _q_deviceFound(const QString &address, const QVariantMap &dict);
    void _q_propertyChanged(const QString &name, const QDBusVariant &value);
    void _q_InterfacesAdded(const QDBusObjectPath &object_path,
                            InterfaceList interfaces_and_properties);
    void _q_discoveryFinished();
    void _q_discoveryInterrupted(const QString &path);
    void _q_PropertiesChanged(const QString &interface,
                              const QString &path,
                              const QVariantMap &changed_properties,
                              const QStringList &invalidated_properties);
    void _q_extendedDeviceDiscoveryTimeout();
#endif

private:
    QList<QBluetoothDeviceInfo> discoveredDevices;
    QBluetoothDeviceDiscoveryAgent::InquiryType inquiryType;

    QBluetoothDeviceDiscoveryAgent::Error lastError;
    QString errorString;

#ifdef QT_ANDROID_BLUETOOTH
private slots:
    void processSdpDiscoveryFinished();
    void processDiscoveredDevices(const QBluetoothDeviceInfo &info, bool isLeResult);
    friend void QtBluetoothLE_leScanResult(JNIEnv *, jobject, jlong, jobject);
    void stopLowEnergyScan();

private:
    void startLowEnergyScan();

    DeviceDiscoveryBroadcastReceiver *receiver;
    QBluetoothAddress m_adapterAddress;
    short m_active;
    QAndroidJniObject adapter;
    QAndroidJniObject leScanner;
    QTimer *leScanTimeout;

    bool pendingCancel, pendingStart;
#elif QT_CONFIG(bluez)
    QBluetoothAddress m_adapterAddress;
    bool pendingCancel;
    bool pendingStart;
    OrgBluezManagerInterface *manager = nullptr;
    OrgBluezAdapterInterface *adapter = nullptr;
    OrgFreedesktopDBusObjectManagerInterface *managerBluez5 = nullptr;
    OrgBluezAdapter1Interface *adapterBluez5 = nullptr;
    QTimer *discoveryTimer = nullptr;
    QList<OrgFreedesktopDBusPropertiesInterface *> propertyMonitors;

    void deviceFoundBluez5(const QString &devicePath, const QVariantMap &properties);
    void startBluez5(QBluetoothDeviceDiscoveryAgent::DiscoveryMethods methods);

    bool useExtendedDiscovery;
    QTimer extendedDiscoveryTimer;
    QMap<QString, QVariantMap> devicesProperties;
#endif

#ifdef QT_WIN_BLUETOOTH
public:
    static QString discoveredLeDeviceSystemPath(const QBluetoothAddress &deviceAddress);

private:
    void cancelDiscovery();
    void restartDiscovery();
    void finishDiscovery(QBluetoothDeviceDiscoveryAgent::Error errorCode, const QString &errorText);

    void startLeDevicesDiscovery();
    void completeLeDevicesDiscovery(const QVariant &res);
    void startClassicDevicesDiscovery(Qt::HANDLE hSearch = nullptr);
    void completeClassicDevicesDiscovery(const QVariant &res);

    void processDiscoveredDevice(const QBluetoothDeviceInfo &foundDevice);

    QBluetoothAddress adapterAddress;
    bool pendingCancel;
    bool pendingStart;
    bool active;

    QThread *threadLE = nullptr;
    QThread *threadClassic = nullptr;
    ThreadWorkerDeviceDiscovery *threadWorkerLE = nullptr;
    ThreadWorkerDeviceDiscovery *threadWorkerClassic = nullptr;
#endif

#ifdef QT_WINRT_BLUETOOTH
private slots:
    void registerDevice(const QBluetoothDeviceInfo &info);
    void updateDeviceData(const QBluetoothAddress &address, QBluetoothDeviceInfo::Fields fields,
                          qint16 rssi, ManufacturerData manufacturerData);
    void onErrorOccured(QBluetoothDeviceDiscoveryAgent::Error e);
    void onScanFinished();

private:
    void disconnectAndClearWorker();
    QPointer<QWinRTBluetoothDeviceDiscoveryWorker> worker;
    QTimer *leScanTimer;
#endif

#ifdef Q_OS_DARWIN

    void startLE();

#ifdef Q_OS_MACOS

    void startClassic();

    // Classic (IOBluetooth) inquiry delegate's methods:
    void inquiryFinished() override;
    void error(IOReturn error) override;
    void classicDeviceFound(void *device) override;
    // Classic (IOBluetooth) errors:
    void setError(IOReturn error, const QString &text = QString());

#endif // Q_OS_MACOS

    // LE scan delegates (CoreBluetooth, all Darwin OSes):
    void LEinquiryFinished();
    void LEinquiryError(QBluetoothDeviceDiscoveryAgent::Error error);
    void LEnotSupported();

    // LE errors:
    void setError(QBluetoothDeviceDiscoveryAgent::Error,
                  const QString &text = QString());

    // Both LE and Classic devices go there:
    void deviceFound(const QBluetoothDeviceInfo &newDeviceInfo);

    enum AgentState {
        NonActive,
        ClassicScan, // macOS (IOBluetooth) only
        LEScan
    } agentState;

    QBluetoothAddress adapterAddress;

    bool startPending;
    bool stopPending;

#ifdef Q_OS_MACOS

    DarwinBluetooth::ScopedPointer controller;
    DarwinBluetooth::ScopedPointer inquiry;

#endif // Q_OS_MACOS

    DarwinBluetooth::ScopedPointer inquiryLE;

#endif // Q_OS_DARWIN

    int lowEnergySearchTimeout;
    QBluetoothDeviceDiscoveryAgent::DiscoveryMethods requestedMethods;
    QBluetoothDeviceDiscoveryAgent *q_ptr;
};

QT_END_NAMESPACE

#endif
