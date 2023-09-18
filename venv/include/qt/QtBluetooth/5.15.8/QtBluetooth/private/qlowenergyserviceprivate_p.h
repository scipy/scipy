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

#ifndef QLOWENERGYSERVICEPRIVATE_P_H
#define QLOWENERGYSERVICEPRIVATE_P_H

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

#include <QtCore/QObject>
#include <QtCore/QPointer>
#include <QtBluetooth/qbluetooth.h>
#include <QtBluetooth/QLowEnergyService>
#include <QtBluetooth/QLowEnergyCharacteristic>

#if defined(QT_ANDROID_BLUETOOTH)
#include <QtAndroidExtras/QAndroidJniObject>
#endif
#if defined(QT_WIN_BLUETOOTH)
#include <qt_windows.h>
#endif

QT_BEGIN_NAMESPACE

class QLowEnergyControllerPrivate;

class QLowEnergyServicePrivate : public QObject
{
    Q_OBJECT
public:
    explicit QLowEnergyServicePrivate(QObject *parent = nullptr);
    ~QLowEnergyServicePrivate();

    struct DescData {
        QByteArray value;
        QBluetoothUuid uuid;
    };

    struct CharData {
        QLowEnergyHandle valueHandle;
        QBluetoothUuid uuid;
        QLowEnergyCharacteristic::PropertyTypes properties;
        QByteArray value;
        QHash<QLowEnergyHandle, DescData> descriptorList;
#ifdef QT_WIN_BLUETOOTH
        Qt::HANDLE hValueChangeEvent;
#endif
    };

    enum GattAttributeTypes {
        PrimaryService = 0x2800,
        SecondaryService = 0x2801,
        IncludeAttribute = 0x2802,
        Characteristic = 0x2803
    };

    void setController(QLowEnergyControllerPrivate* control);
    void setError(QLowEnergyService::ServiceError newError);
    void setState(QLowEnergyService::ServiceState newState);

signals:
    void stateChanged(QLowEnergyService::ServiceState newState);
    void error(QLowEnergyService::ServiceError error);
    void characteristicChanged(const QLowEnergyCharacteristic &characteristic,
                               const QByteArray &newValue);
    void characteristicRead(const QLowEnergyCharacteristic &info,
                            const QByteArray &value);
    void characteristicWritten(const QLowEnergyCharacteristic &characteristic,
                               const QByteArray &newValue);
    void descriptorRead(const QLowEnergyDescriptor &info,
                        const QByteArray &value);
    void descriptorWritten(const QLowEnergyDescriptor &descriptor,
                           const QByteArray &newValue);

public:
    QLowEnergyHandle startHandle;
    QLowEnergyHandle endHandle;

    QBluetoothUuid uuid;
    QList<QBluetoothUuid> includedServices;
    QLowEnergyService::ServiceTypes type;
    QLowEnergyService::ServiceState state;
    QLowEnergyService::ServiceError lastError;

    QHash<QLowEnergyHandle, CharData> characteristicList;

    QPointer<QLowEnergyControllerPrivate> controller;

#if defined(QT_ANDROID_BLUETOOTH)
    // reference to the BluetoothGattService object
    QAndroidJniObject androidService;
#endif
#if defined(QT_WIN_BLUETOOTH)
    Qt::HANDLE hService = nullptr;
#endif

};

typedef QHash<QLowEnergyHandle, QLowEnergyServicePrivate::CharData> CharacteristicDataMap;
typedef QHash<QLowEnergyHandle, QLowEnergyServicePrivate::DescData> DescriptorDataMap;

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QSharedPointer<QLowEnergyServicePrivate>)

#endif // QLOWENERGYSERVICEPRIVATE_P_H
