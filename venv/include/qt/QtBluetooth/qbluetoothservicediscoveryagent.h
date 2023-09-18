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

#ifndef QBLUETOOTHSERVICEDISCOVERYAGENT_H
#define QBLUETOOTHSERVICEDISCOVERYAGENT_H

#include <QtBluetooth/qtbluetoothglobal.h>

#include <QtCore/QObject>
#include <QtCore/QVariant>

#include <QtBluetooth/QBluetoothServiceInfo>
#include <QtBluetooth/QBluetoothUuid>
#include <QtBluetooth/QBluetoothDeviceDiscoveryAgent>

#if QT_CONFIG(bluez)
#include <QtCore/qprocess.h>
#endif

QT_BEGIN_NAMESPACE

class QBluetoothAddress;
class QBluetoothServiceDiscoveryAgentPrivate;

class Q_BLUETOOTH_EXPORT QBluetoothServiceDiscoveryAgent : public QObject
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QBluetoothServiceDiscoveryAgent)

public:
    enum Error {
        NoError =  QBluetoothDeviceDiscoveryAgent::NoError,
        InputOutputError = QBluetoothDeviceDiscoveryAgent::InputOutputError,
        PoweredOffError = QBluetoothDeviceDiscoveryAgent::PoweredOffError,
        InvalidBluetoothAdapterError = QBluetoothDeviceDiscoveryAgent::InvalidBluetoothAdapterError,
        UnknownError = QBluetoothDeviceDiscoveryAgent::UnknownError //=100
        //New Errors must be added after Unknown Error the space before UnknownError is reserved
        //for future device discovery errors
    };
    Q_ENUM(Error)

    enum DiscoveryMode {
        MinimalDiscovery,
        FullDiscovery
    };
    Q_ENUM(DiscoveryMode)

    explicit QBluetoothServiceDiscoveryAgent(QObject *parent = nullptr);
    explicit QBluetoothServiceDiscoveryAgent(const QBluetoothAddress &deviceAdapter, QObject *parent = nullptr);
    ~QBluetoothServiceDiscoveryAgent();

    bool isActive() const;

    Error error() const;
    QString errorString() const;

    QList<QBluetoothServiceInfo> discoveredServices() const;

    void setUuidFilter(const QList<QBluetoothUuid> &uuids);
    void setUuidFilter(const QBluetoothUuid &uuid);
    QList<QBluetoothUuid> uuidFilter() const;
    bool setRemoteAddress(const QBluetoothAddress &address);
    QBluetoothAddress remoteAddress() const;

public Q_SLOTS:
    void start(DiscoveryMode mode = MinimalDiscovery);
    void stop();
    void clear();

Q_SIGNALS:
    void serviceDiscovered(const QBluetoothServiceInfo &info);
    void finished();
    void canceled();
    void error(QBluetoothServiceDiscoveryAgent::Error error);

private:
    QBluetoothServiceDiscoveryAgentPrivate *d_ptr;
};

QT_END_NAMESPACE

#endif
