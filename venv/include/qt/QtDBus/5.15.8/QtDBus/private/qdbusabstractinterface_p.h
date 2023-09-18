/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Copyright (C) 2016 Intel Corporation.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtDBus module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the public API.  This header file may
// change from version to version without notice, or even be
// removed.
//
// We mean it.
//
//

#ifndef QDBUSABSTRACTINTERFACE_P_H
#define QDBUSABSTRACTINTERFACE_P_H

#include <QtDBus/private/qtdbusglobal_p.h>
#include <qdbusabstractinterface.h>
#include <qdbusconnection.h>
#include <qdbuserror.h>
#include "qdbusconnection_p.h"
#include "private/qobject_p.h"

#define ANNOTATION_NO_WAIT      "org.freedesktop.DBus.Method.NoReply"

#ifndef QT_NO_DBUS

QT_BEGIN_NAMESPACE

class QDBusAbstractInterfacePrivate : public QObjectPrivate
{
public:
    Q_DECLARE_PUBLIC(QDBusAbstractInterface)

    mutable QDBusConnection connection; // mutable because we want to make calls from const functions
    QString service;
    QString currentOwner;
    QString path;
    QString interface;
    mutable QDBusError lastError;
    int timeout;

    // this is set during creation and never changed
    // it can't be const because QDBusInterfacePrivate has one more check
    bool isValid;

    QDBusAbstractInterfacePrivate(const QString &serv, const QString &p,
                                  const QString &iface, const QDBusConnection& con, bool dynamic);
    virtual ~QDBusAbstractInterfacePrivate() { }
    void initOwnerTracking();
    bool canMakeCalls() const;

    // these functions do not check if the property is valid
    bool property(const QMetaProperty &mp, void *returnValuePtr) const;
    bool setProperty(const QMetaProperty &mp, const QVariant &value);

    // return conn's d pointer
    inline QDBusConnectionPrivate *connectionPrivate() const
    { return QDBusConnectionPrivate::d(connection); }

    void _q_serviceOwnerChanged(const QString &name, const QString &oldOwner, const QString &newOwner);

    static void finishDisconnectNotify(QDBusAbstractInterface *iface, int signalId);
};

QT_END_NAMESPACE

#endif // QT_NO_DBUS
#endif
