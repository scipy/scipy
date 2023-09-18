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

#ifndef QBLUETOOTHDEVICEDISCOVERYAGENT_H
#define QBLUETOOTHDEVICEDISCOVERYAGENT_H

#include <QtBluetooth/qtbluetoothglobal.h>

#include <QtCore/QObject>
#include <QtBluetooth/QBluetoothDeviceInfo>
#include <QtBluetooth/QBluetoothAddress>

QT_BEGIN_NAMESPACE

class QBluetoothDeviceDiscoveryAgentPrivate;

class Q_BLUETOOTH_EXPORT QBluetoothDeviceDiscoveryAgent : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QBluetoothDeviceDiscoveryAgent::InquiryType inquiryType
               READ inquiryType WRITE setInquiryType)

public:
    // FIXME: add more errors
    // FIXME: add bluez error handling
    enum Error {
        NoError,
        InputOutputError,
        PoweredOffError,
        InvalidBluetoothAdapterError,
        UnsupportedPlatformError,
        UnsupportedDiscoveryMethod,
        UnknownError = 100 // New errors must be added before Unknown error
    };
    Q_ENUM(Error)

    enum InquiryType {
        GeneralUnlimitedInquiry,
        LimitedInquiry
    };
    Q_ENUM(InquiryType)

    enum DiscoveryMethod
    {
        NoMethod = 0x0,
        ClassicMethod = 0x01,
        LowEnergyMethod = 0x02,
    };
    Q_DECLARE_FLAGS(DiscoveryMethods, DiscoveryMethod)
    Q_FLAG(DiscoveryMethods)

    explicit QBluetoothDeviceDiscoveryAgent(QObject *parent = nullptr);
    explicit QBluetoothDeviceDiscoveryAgent(const QBluetoothAddress &deviceAdapter,
                                            QObject *parent = nullptr);
    ~QBluetoothDeviceDiscoveryAgent();

    // TODO Remove inquiry type in Qt 6 -> not really used anywhere
    QBluetoothDeviceDiscoveryAgent::InquiryType inquiryType() const;
    void setInquiryType(QBluetoothDeviceDiscoveryAgent::InquiryType type);

    bool isActive() const;

    Error error() const;
    QString errorString() const;

    QList<QBluetoothDeviceInfo> discoveredDevices() const;

    void setLowEnergyDiscoveryTimeout(int msTimeout);
    int lowEnergyDiscoveryTimeout() const;

    static DiscoveryMethods supportedDiscoveryMethods();
public Q_SLOTS:
    void start();
    void start(DiscoveryMethods method);
    void stop();

Q_SIGNALS:
    void deviceDiscovered(const QBluetoothDeviceInfo &info);
    void deviceUpdated(const QBluetoothDeviceInfo &info, QBluetoothDeviceInfo::Fields updatedFields);
    void finished();
    void error(QBluetoothDeviceDiscoveryAgent::Error error);
    void canceled();

private:
    Q_DECLARE_PRIVATE(QBluetoothDeviceDiscoveryAgent)
    QBluetoothDeviceDiscoveryAgentPrivate *d_ptr;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QBluetoothDeviceDiscoveryAgent::DiscoveryMethods)

QT_END_NAMESPACE

#endif
