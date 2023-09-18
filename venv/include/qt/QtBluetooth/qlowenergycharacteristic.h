/***************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Copyright (C) 2016 BlackBerry Limited all rights reserved
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

#ifndef QLOWENERGYCHARACTERISTIC_H
#define QLOWENERGYCHARACTERISTIC_H
#include <QtCore/QSharedPointer>
#include <QtCore/QObject>
#include <QtBluetooth/qbluetooth.h>
#include <QtBluetooth/QBluetoothUuid>
#include <QtBluetooth/QLowEnergyDescriptor>

QT_BEGIN_NAMESPACE

class QBluetoothUuid;
class QLowEnergyServicePrivate;
struct QLowEnergyCharacteristicPrivate;
class Q_BLUETOOTH_EXPORT QLowEnergyCharacteristic
{
public:

    enum PropertyType {
        Unknown = 0x00,
        Broadcasting = 0x01,
        Read = 0x02,
        WriteNoResponse = 0x04,
        Write = 0x08,
        Notify = 0x10,
        Indicate = 0x20,
        WriteSigned = 0x40,
        ExtendedProperty = 0x80
    };
    Q_DECLARE_FLAGS(PropertyTypes, PropertyType)

    QLowEnergyCharacteristic();
    QLowEnergyCharacteristic(const QLowEnergyCharacteristic &other);
    ~QLowEnergyCharacteristic();

    QLowEnergyCharacteristic &operator=(const QLowEnergyCharacteristic &other);
    bool operator==(const QLowEnergyCharacteristic &other) const;
    bool operator!=(const QLowEnergyCharacteristic &other) const;

    QString name() const;

    QBluetoothUuid uuid() const;

    QByteArray value() const;

    QLowEnergyCharacteristic::PropertyTypes properties() const;
    QLowEnergyHandle handle() const;

    QLowEnergyDescriptor descriptor(const QBluetoothUuid &uuid) const;
    QList<QLowEnergyDescriptor> descriptors() const;

    bool isValid() const;

protected:
    QLowEnergyHandle attributeHandle() const;

    QSharedPointer<QLowEnergyServicePrivate> d_ptr;

    friend class QLowEnergyService;
    friend class QLowEnergyControllerPrivate;
    friend class QLowEnergyControllerPrivateAndroid;
    friend class QLowEnergyControllerPrivateBluez;
    friend class QLowEnergyControllerPrivateBluezDBus;
    friend class QLowEnergyControllerPrivateCommon;
    friend class QLowEnergyControllerPrivateWin32;
    friend class QLowEnergyControllerPrivateDarwin;
    friend class QLowEnergyControllerPrivateWinRT;
    friend class QLowEnergyControllerPrivateWinRTNew;
    QLowEnergyCharacteristicPrivate *data = nullptr;
    QLowEnergyCharacteristic(QSharedPointer<QLowEnergyServicePrivate> p,
                             QLowEnergyHandle handle);
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QLowEnergyCharacteristic::PropertyTypes)

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QLowEnergyCharacteristic)

#endif // QLOWENERGYCHARACTERISTIC_H
