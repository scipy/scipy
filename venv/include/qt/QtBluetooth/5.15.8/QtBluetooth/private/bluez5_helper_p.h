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

#ifndef BLUEZ5_HELPER_H
#define BLUEZ5_HELPER_H

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
#include <QtDBus/QtDBus>
#include <QtBluetooth/QBluetoothAddress>
#include <QtBluetooth/private/qtbluetoothglobal_p.h>

typedef QMap<QString, QVariantMap> InterfaceList;
typedef QMap<QDBusObjectPath, InterfaceList> ManagedObjectList;
typedef QMap<quint16, QDBusVariant> ManufacturerDataList;

Q_DECLARE_METATYPE(InterfaceList)
Q_DECLARE_METATYPE(ManufacturerDataList)
Q_DECLARE_METATYPE(ManagedObjectList)

QT_BEGIN_NAMESPACE

bool isBluez5();

// exported for unit test purposes
Q_BLUETOOTH_PRIVATE_EXPORT QVersionNumber bluetoothdVersion();

QString sanitizeNameForDBus(const QString& text);

QString findAdapterForAddress(const QBluetoothAddress &wantedAddress, bool *ok);

class QtBluezDiscoveryManagerPrivate;
class QtBluezDiscoveryManager : public QObject
{
    Q_OBJECT
public:
    QtBluezDiscoveryManager(QObject* parent = nullptr);
    ~QtBluezDiscoveryManager();
    static QtBluezDiscoveryManager *instance();

    bool registerDiscoveryInterest(const QString &adapterPath);
    void unregisterDiscoveryInterest(const QString &adapterPath);

    //void dumpState() const;

signals:
    void discoveryInterrupted(const QString &adapterPath);

private slots:
    void InterfacesRemoved(const QDBusObjectPath &object_path,
                           const QStringList &interfaces);
    void PropertiesChanged(const QString &interface,
                           const QVariantMap &changed_properties,
                           const QStringList &invalidated_properties,
                           const QDBusMessage &msg);

private:
    void removeAdapterFromMonitoring(const QString &dbusPath);

    QtBluezDiscoveryManagerPrivate *d;
};

QT_END_NAMESPACE

#endif
