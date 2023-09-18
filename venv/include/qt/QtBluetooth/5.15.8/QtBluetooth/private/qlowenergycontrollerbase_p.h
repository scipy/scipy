/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
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

#ifndef QLOWENERGYCONTROLLERPRIVATEBASE_P_H
#define QLOWENERGYCONTROLLERPRIVATEBASE_P_H

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

#include <qglobal.h>
#include <QtCore/qobject.h>

#include <QtBluetooth/qlowenergycontroller.h>

#include "qlowenergyserviceprivate_p.h"

QT_BEGIN_NAMESPACE

typedef QMap<QBluetoothUuid, QSharedPointer<QLowEnergyServicePrivate> > ServiceDataMap;

class QLowEnergyControllerPrivate : public QObject
{
    Q_OBJECT
public:
    // This class is required to enable selection of multiple
    // alternative QLowEnergyControllerPrivate implementations on BlueZ.
    // Bluez has a low level ATT protocol stack implementation and a DBus
    // implementation.

    QLowEnergyControllerPrivate();
    virtual ~QLowEnergyControllerPrivate();

    // interface definition
    virtual void init() = 0;
    virtual void connectToDevice() = 0;
    virtual void disconnectFromDevice() = 0;

    virtual void discoverServices() = 0;
    virtual void discoverServiceDetails(const QBluetoothUuid &/*service*/) = 0;

    virtual void readCharacteristic(
                        const QSharedPointer<QLowEnergyServicePrivate> /*service*/,
                        const QLowEnergyHandle /*charHandle*/) = 0;
    virtual void readDescriptor(
                        const QSharedPointer<QLowEnergyServicePrivate> /*service*/,
                        const QLowEnergyHandle /*charHandle*/,
                        const QLowEnergyHandle /*descriptorHandle*/) = 0;

    virtual void writeCharacteristic(
                        const QSharedPointer<QLowEnergyServicePrivate> /*service*/,
                        const QLowEnergyHandle /*charHandle*/,
                        const QByteArray &/*newValue*/,
                        QLowEnergyService::WriteMode /*writeMode*/) = 0;
    virtual void writeDescriptor(
                        const QSharedPointer<QLowEnergyServicePrivate> /*service*/,
                        const QLowEnergyHandle /*charHandle*/,
                        const QLowEnergyHandle /*descriptorHandle*/,
                        const QByteArray &/*newValue*/) = 0;

    virtual void startAdvertising(
                        const QLowEnergyAdvertisingParameters &/* params */,
                        const QLowEnergyAdvertisingData &/* advertisingData */,
                        const QLowEnergyAdvertisingData &/* scanResponseData */) = 0;
    virtual void stopAdvertising() = 0;

    virtual void requestConnectionUpdate(
                        const QLowEnergyConnectionParameters & /* params */) = 0;
    virtual void addToGenericAttributeList(
                        const QLowEnergyServiceData &/* service */,
                        QLowEnergyHandle /* startHandle */) = 0;


    virtual QLowEnergyService *addServiceHelper(
                        const QLowEnergyServiceData &service);

    // common backend methods
    bool isValidLocalAdapter();
    void setError(QLowEnergyController::Error newError);
    void setState(QLowEnergyController::ControllerState newState);

    // public variables
    QLowEnergyController::Role role;
    QLowEnergyController::RemoteAddressType addressType;

    // list of all found service uuids on remote device
    ServiceDataMap serviceList;
    // list of all found service uuids on local peripheral device
    ServiceDataMap localServices;

    //common helper functions
    QSharedPointer<QLowEnergyServicePrivate> serviceForHandle(QLowEnergyHandle handle);
    QLowEnergyCharacteristic characteristicForHandle(QLowEnergyHandle handle);
    QLowEnergyDescriptor descriptorForHandle(QLowEnergyHandle handle);
    quint16 updateValueOfCharacteristic(QLowEnergyHandle charHandle,
                                 const QByteArray &value,
                                 bool appendValue);
    quint16 updateValueOfDescriptor(QLowEnergyHandle charHandle,
                                 QLowEnergyHandle descriptorHandle,
                                 const QByteArray &value,
                                 bool appendValue);
    void invalidateServices();

protected:
    QLowEnergyController::ControllerState state = QLowEnergyController::UnconnectedState;
    QLowEnergyController::Error error = QLowEnergyController::NoError;
    QString errorString;

    QBluetoothAddress remoteDevice;
    QBluetoothAddress localAdapter;

    QLowEnergyHandle lastLocalHandle{};

    QString remoteName; // device name of the remote
    QBluetoothUuid deviceUuid; // quite useless anywhere but Darwin (CoreBluetooth).

    Q_DECLARE_PUBLIC(QLowEnergyController)
    QLowEnergyController *q_ptr;
};

QT_END_NAMESPACE

#endif // QLOWENERGYCONTROLLERPRIVATEBASE_P_H
