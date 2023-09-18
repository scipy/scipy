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

#ifndef QLOWENERGYADVERTISINGPARAMETERS_H
#define QLOWENERGYADVERTISINGPARAMETERS_H

#include <QtBluetooth/qtbluetoothglobal.h>
#include <QtBluetooth/qbluetoothaddress.h>
#include <QtBluetooth/qlowenergycontroller.h>
#include <QtCore/qlist.h>
#include <QtCore/qshareddata.h>

QT_BEGIN_NAMESPACE

class QLowEnergyAdvertisingParametersPrivate;

class Q_BLUETOOTH_EXPORT QLowEnergyAdvertisingParameters
{
    friend Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyAdvertisingParameters &p1,
                                              const QLowEnergyAdvertisingParameters &p2);
public:
    QLowEnergyAdvertisingParameters();
    QLowEnergyAdvertisingParameters(const QLowEnergyAdvertisingParameters &other);
    ~QLowEnergyAdvertisingParameters();

    QLowEnergyAdvertisingParameters &operator=(const QLowEnergyAdvertisingParameters &other);

    enum Mode { AdvInd = 0x0, AdvScanInd = 0x2, AdvNonConnInd = 0x3 };
    void setMode(Mode mode);
    Mode mode() const;

    struct AddressInfo {
        AddressInfo(const QBluetoothAddress &addr, QLowEnergyController::RemoteAddressType t)
            : address(addr), type(t) {}
        AddressInfo() : type(QLowEnergyController::PublicAddress) {}

        QBluetoothAddress address;
        QLowEnergyController::RemoteAddressType type;
    };
    enum FilterPolicy {
        IgnoreWhiteList = 0x00,
        UseWhiteListForScanning = 0x01,
        UseWhiteListForConnecting = 0x02,
        UseWhiteListForScanningAndConnecting = 0x03,
    };
    void setWhiteList(const QList<AddressInfo> &whiteList, FilterPolicy policy);
    QList<AddressInfo> whiteList() const;
    FilterPolicy filterPolicy() const;

    void setInterval(quint16 minimum, quint16 maximum);
    int minimumInterval() const;
    int maximumInterval() const;

    // TODO: own address type
    // TODO: For ADV_DIRECT_IND: peer address + peer address type

    void swap(QLowEnergyAdvertisingParameters &other) Q_DECL_NOTHROW { qSwap(d, other.d); }

private:
    QSharedDataPointer<QLowEnergyAdvertisingParametersPrivate> d;
};

inline bool operator==(const QLowEnergyAdvertisingParameters::AddressInfo &ai1,
                       const QLowEnergyAdvertisingParameters::AddressInfo &ai2)
{
    return ai1.address == ai2.address && ai1.type == ai2.type;
}

Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyAdvertisingParameters &p1,
                                   const QLowEnergyAdvertisingParameters &p2);
inline bool operator!=(const QLowEnergyAdvertisingParameters &p1,
                       const QLowEnergyAdvertisingParameters &p2)
{
    return !(p1 == p2);
}

Q_DECLARE_SHARED(QLowEnergyAdvertisingParameters)

QT_END_NAMESPACE

#endif // Include guard
