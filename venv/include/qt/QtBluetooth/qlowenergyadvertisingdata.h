/***************************************************************************
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

#ifndef QLOWENERGYADVERTISINGDATA_H
#define QLOWENERGYADVERTISINGDATA_H

#include <QtBluetooth/qtbluetoothglobal.h>
#include <QtBluetooth/qbluetoothuuid.h>
#include <QtCore/qshareddata.h>

QT_BEGIN_NAMESPACE

class QLowEnergyAdvertisingDataPrivate;

class Q_BLUETOOTH_EXPORT QLowEnergyAdvertisingData
{
    friend Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyAdvertisingData &data1,
                                              const QLowEnergyAdvertisingData &data2);
public:
    QLowEnergyAdvertisingData();
    QLowEnergyAdvertisingData(const QLowEnergyAdvertisingData &other);
    ~QLowEnergyAdvertisingData();

    QLowEnergyAdvertisingData &operator=(const QLowEnergyAdvertisingData &other);

    void setLocalName(const QString &name);
    QString localName() const;

    static quint16 invalidManufacturerId() { return 0xffff; }
    void setManufacturerData(quint16 id, const QByteArray &data);
    quint16 manufacturerId() const;
    QByteArray manufacturerData() const;

    void setIncludePowerLevel(bool doInclude);
    bool includePowerLevel() const;

    enum Discoverability {
        DiscoverabilityNone, DiscoverabilityLimited, DiscoverabilityGeneral
    };
    void setDiscoverability(Discoverability mode);
    Discoverability discoverability() const;

    void setServices(const QList<QBluetoothUuid> &services);
    QList<QBluetoothUuid> services() const;

    // TODO: BR/EDR capability flag?

    void setRawData(const QByteArray &data);
    QByteArray rawData() const;

    void swap(QLowEnergyAdvertisingData &other) Q_DECL_NOTHROW { qSwap(d, other.d); }

private:
    QSharedDataPointer<QLowEnergyAdvertisingDataPrivate> d;
};

Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyAdvertisingData &data1,
                                   const QLowEnergyAdvertisingData &data2);
inline bool operator!=(const QLowEnergyAdvertisingData &data1,
                       const QLowEnergyAdvertisingData &data2)
{
    return !(data1 == data2);
}

Q_DECLARE_SHARED(QLowEnergyAdvertisingData)

QT_END_NAMESPACE

#endif // Include guard
