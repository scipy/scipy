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

#ifndef QLOWENERGYSERVICEDATA_H
#define QLOWENERGYSERVICEDATA_H

#include <QtBluetooth/qtbluetoothglobal.h>
#include <QtCore/qshareddata.h>

QT_BEGIN_NAMESPACE

class QBluetoothUuid;
class QLowEnergyCharacteristicData;
class QLowEnergyService;
struct QLowEnergyServiceDataPrivate;

class Q_BLUETOOTH_EXPORT QLowEnergyServiceData
{
    friend Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyServiceData &sd1,
                                              const QLowEnergyServiceData &sd2);
public:
    QLowEnergyServiceData();
    QLowEnergyServiceData(const QLowEnergyServiceData &other);
    ~QLowEnergyServiceData();

    QLowEnergyServiceData &operator=(const QLowEnergyServiceData &other);

    enum ServiceType { ServiceTypePrimary = 0x2800, ServiceTypeSecondary = 0x2801 };
    ServiceType type() const;
    void setType(ServiceType type);

    QBluetoothUuid uuid() const;
    void setUuid(const QBluetoothUuid &uuid);

    QList<QLowEnergyService *> includedServices() const;
    void setIncludedServices(const QList<QLowEnergyService *> &services);
    void addIncludedService(QLowEnergyService *service);

    QList<QLowEnergyCharacteristicData> characteristics() const;
    void setCharacteristics(const QList<QLowEnergyCharacteristicData> &characteristics);
    void addCharacteristic(const QLowEnergyCharacteristicData &characteristic);

    bool isValid() const;

    void swap(QLowEnergyServiceData &other) Q_DECL_NOTHROW { qSwap(d, other.d); }

private:
    QSharedDataPointer<QLowEnergyServiceDataPrivate> d;
};

Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyServiceData &sd1,
                                   const QLowEnergyServiceData &sd2);
inline bool operator!=(const QLowEnergyServiceData &sd1, const QLowEnergyServiceData &sd2)
{
    return !(sd1 == sd2);
}

Q_DECLARE_SHARED(QLowEnergyServiceData)

QT_END_NAMESPACE

#endif // Include guard.
