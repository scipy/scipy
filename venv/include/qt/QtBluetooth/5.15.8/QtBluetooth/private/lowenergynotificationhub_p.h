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

#ifndef LOWENERGYNOTIFICATIONHUB_H
#define LOWENERGYNOTIFICATIONHUB_H

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
#include <QtCore/QReadWriteLock>
#include <QtCore/private/qjnihelpers_p.h>
#include <QtAndroidExtras/QAndroidJniObject>
#include <QtBluetooth/QBluetoothAddress>
#include <QtBluetooth/QLowEnergyController>
#include <QtBluetooth/QLowEnergyService>
#include <jni.h>

#include <QtBluetooth/QLowEnergyCharacteristic>

QT_BEGIN_NAMESPACE

class LowEnergyNotificationHub : public QObject
{
    Q_OBJECT
public:
    explicit LowEnergyNotificationHub(const QBluetoothAddress &remote, bool isPeripheral,
                                      QObject *parent = nullptr);
    ~LowEnergyNotificationHub();

    static void lowEnergy_connectionChange(JNIEnv*, jobject, jlong qtObject,
                                           jint errorCode, jint newState);
    static void lowEnergy_servicesDiscovered(JNIEnv*, jobject, jlong qtObject,
                                             jint errorCode, jobject uuidList);
    static void lowEnergy_serviceDetailsDiscovered(JNIEnv *, jobject,
                                                   jlong qtObject, jobject uuid,
                                                   jint startHandle, jint endHandle);
    static void lowEnergy_characteristicRead(JNIEnv*env, jobject, jlong qtObject,
                                             jobject serviceUuid,
                                             jint handle, jobject charUuid,
                                             jint properties, jbyteArray data);
    static void lowEnergy_descriptorRead(JNIEnv *env, jobject, jlong qtObject,
                                         jobject sUuid, jobject cUuid,
                                         jint handle, jobject dUuid, jbyteArray data);
    static void lowEnergy_characteristicWritten(JNIEnv *, jobject, jlong qtObject,
                                                jint charHandle, jbyteArray data,
                                                jint errorCode);
    static void lowEnergy_descriptorWritten(JNIEnv *, jobject, jlong qtObject,
                                            jint descHandle, jbyteArray data,
                                            jint errorCode);
    static void lowEnergy_serverDescriptorWritten(JNIEnv *, jobject, jlong qtObject,
                                                  jobject descriptor, jbyteArray newValue);
    static void lowEnergy_characteristicChanged(JNIEnv *, jobject, jlong qtObject,
                                                jint charHandle, jbyteArray data);
    static void lowEnergy_serverCharacteristicChanged(JNIEnv *, jobject, jlong qtObject,
                                                jobject characteristic, jbyteArray newValue);
    static void lowEnergy_serviceError(JNIEnv *, jobject, jlong qtObject,
                                       jint attributeHandle, int errorCode);
    static void lowEnergy_advertisementError(JNIEnv *, jobject, jlong qtObject,
                                               jint status);

    QAndroidJniObject javaObject()
    {
        return jBluetoothLe;
    }

signals:
    void connectionUpdated(QLowEnergyController::ControllerState newState,
            QLowEnergyController::Error errorCode);
    void servicesDiscovered(QLowEnergyController::Error errorCode, const QString &uuids);
    void serviceDetailsDiscoveryFinished(const QString& serviceUuid,
            int startHandle, int endHandle);
    void characteristicRead(const QBluetoothUuid &serviceUuid,
            int handle, const QBluetoothUuid &charUuid,
            int properties, const QByteArray &data);
    void descriptorRead(const QBluetoothUuid &serviceUuid, const QBluetoothUuid &charUuid,
            int handle, const QBluetoothUuid &descUuid, const QByteArray &data);
    void characteristicWritten(int charHandle, const QByteArray &data,
                               QLowEnergyService::ServiceError errorCode);
    void descriptorWritten(int descHandle, const QByteArray &data,
                           QLowEnergyService::ServiceError errorCode);
    void serverDescriptorWritten(const QAndroidJniObject &descriptor, const QByteArray &newValue);
    void characteristicChanged(int charHandle, const QByteArray &data);
    void serverCharacteristicChanged(const QAndroidJniObject &characteristic, const QByteArray &newValue);
    void serviceError(int attributeHandle, QLowEnergyService::ServiceError errorCode);
    void advertisementError(int status);

public slots:
private:
    static QReadWriteLock lock;

    QAndroidJniObject jBluetoothLe;
    long javaToCtoken;

};

QT_END_NAMESPACE

#endif // LOWENERGYNOTIFICATIONHUB_H

