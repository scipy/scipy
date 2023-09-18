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

#ifndef QBLUETOOTHSERVICEINFO_P_H
#define QBLUETOOTHSERVICEINFO_P_H

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

#include "qbluetoothuuid.h"
#include "qbluetoothaddress.h"
#include "qbluetoothdeviceinfo.h"
#include "qbluetoothserviceinfo.h"

#include <QMap>
#include <QVariant>

#ifdef Q_OS_MACOS
#include "osx/btraii_p.h"
#endif

class OrgBluezServiceInterface;
class OrgBluezProfileManager1Interface;

#ifdef QT_WINRT_BLUETOOTH
#include <wrl.h>

namespace ABI {
    namespace Windows {
        namespace Devices {
            namespace Bluetooth {
                namespace Rfcomm {
                    struct IRfcommServiceProvider;
                }
            }
        }
    }
}
#endif

#ifdef QT_WIN_BLUETOOTH
#include <winsock2.h>
#include <ws2bth.h>
#endif

QT_BEGIN_NAMESPACE

class QBluetoothServiceInfo;


class QBluetoothServiceInfoPrivate
    : public QObject
{
    Q_OBJECT
public:
    QBluetoothServiceInfoPrivate();
    ~QBluetoothServiceInfoPrivate();

    bool registerService(const QBluetoothAddress &localAdapter = QBluetoothAddress());

    bool isRegistered() const;

    bool unregisterService();

    QBluetoothDeviceInfo deviceInfo;
    QMap<quint16, QVariant> attributes;

    QBluetoothServiceInfo::Sequence protocolDescriptor(QBluetoothUuid::ProtocolUuid protocol) const;
    int serverChannel() const;
private:
#if QT_CONFIG(bluez)
    bool ensureSdpConnection(const QBluetoothAddress &localAdapter = QBluetoothAddress());

    OrgBluezServiceInterface *service = nullptr;
    OrgBluezProfileManager1Interface *serviceBluez5 = nullptr;
    quint32 serviceRecord;
    QBluetoothAddress currentLocalAdapter;
    QString profilePath;
#endif

#ifdef QT_WINRT_BLUETOOTH
    Microsoft::WRL::ComPtr<ABI::Windows::Devices::Bluetooth::Rfcomm::IRfcommServiceProvider> serviceProvider;

    bool writeSdpAttributes();
#endif

#ifdef QT_WIN_BLUETOOTH
    SOCKADDR_BTH sockaddr = {};
    CSADDR_INFO addrinfo = {};
    WSAQUERYSET regInfo = {};
    QVector<WCHAR> serviceName;
    QVector<WCHAR> serviceDescription;
#endif

#if QT_OSX_BLUETOOTH
public:
    bool registerService(const QBluetoothServiceInfo &info);

private:

    using SDPRecord = DarwinBluetooth::ScopedPointer;
    SDPRecord serviceRecord;
    quint32 serviceRecordHandle = 0;
#endif // QT_OSX_BLUETOOTH

    mutable bool registered = false;
};

QT_END_NAMESPACE

#endif // QBLUETOOTHSERVICEINFO_P_H
