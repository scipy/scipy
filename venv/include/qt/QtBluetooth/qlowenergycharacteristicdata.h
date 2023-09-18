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
#ifndef QLOWENERGYCHARACTERISTICDATA_H
#define QLOWENERGYCHARACTERISTICDATA_H

#include <QtBluetooth/qlowenergycharacteristic.h>
#include <QtCore/qshareddata.h>

QT_BEGIN_NAMESPACE

class QLowEnergyDescriptorData;
struct QLowEnergyCharacteristicDataPrivate;
class Q_BLUETOOTH_EXPORT QLowEnergyCharacteristicData
{
    friend Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyCharacteristicData &cd1,
                                              const QLowEnergyCharacteristicData &cd2);
public:
    QLowEnergyCharacteristicData();
    QLowEnergyCharacteristicData(const QLowEnergyCharacteristicData &other);
    ~QLowEnergyCharacteristicData();

    QLowEnergyCharacteristicData &operator=(const QLowEnergyCharacteristicData &other);

    QBluetoothUuid uuid() const;
    void setUuid(const QBluetoothUuid &uuid);

    QByteArray value() const;
    void setValue(const QByteArray &value);

    QLowEnergyCharacteristic::PropertyTypes properties() const;
    void setProperties(QLowEnergyCharacteristic::PropertyTypes properties);

    QList<QLowEnergyDescriptorData> descriptors() const;
    void setDescriptors(const QList<QLowEnergyDescriptorData> &descriptors);
    void addDescriptor(const QLowEnergyDescriptorData &descriptor);

    void setReadConstraints(QBluetooth::AttAccessConstraints constraints);
    QBluetooth::AttAccessConstraints readConstraints() const;

    void setWriteConstraints(QBluetooth::AttAccessConstraints constraints);
    QBluetooth::AttAccessConstraints writeConstraints() const;

    void setValueLength(int minimum, int maximum);
    int minimumValueLength() const;
    int maximumValueLength() const;

    bool isValid() const;

    void swap(QLowEnergyCharacteristicData &other) Q_DECL_NOTHROW { qSwap(d, other.d); }

private:
    QSharedDataPointer<QLowEnergyCharacteristicDataPrivate> d;
};

Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyCharacteristicData &cd1,
                                   const QLowEnergyCharacteristicData &cd2);
inline bool operator!=(const QLowEnergyCharacteristicData &cd1,
                       const QLowEnergyCharacteristicData &cd2)
{
    return !(cd1 == cd2);
}

Q_DECLARE_SHARED(QLowEnergyCharacteristicData)

QT_END_NAMESPACE

#endif // Include guard.
