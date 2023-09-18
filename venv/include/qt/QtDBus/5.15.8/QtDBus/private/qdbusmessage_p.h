/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
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

#ifndef QDBUSMESSAGE_P_H
#define QDBUSMESSAGE_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of the QLibrary class.  This header file may change from
// version to version without notice, or even be removed.
//
// We mean it.
//

#include <QtDBus/private/qtdbusglobal_p.h>
#include <qatomic.h>
#include <qstring.h>
#include <qdbusmessage.h>
#include <qdbusconnection.h>

struct DBusMessage;

#ifndef QT_NO_DBUS

QT_BEGIN_NAMESPACE

class QDBusConnectionPrivate;

class QDBusMessagePrivate
{
public:
    QDBusMessagePrivate();
    ~QDBusMessagePrivate();

    QList<QVariant> arguments;

    // the following parameters are "const": they are not changed after the constructors
    // the parametersValidated member below controls whether they've been validated already
    QString service, path, interface, name, message, signature;

    DBusMessage *msg;
    DBusMessage *reply;
    mutable QDBusMessage *localReply;
    QAtomicInt ref;
    QDBusMessage::MessageType type;

    mutable uint delayedReply : 1;
    uint localMessage : 1;
    mutable uint parametersValidated : 1;
    uint autoStartService : 1;
    uint interactiveAuthorizationAllowed : 1;

    static void setParametersValidated(QDBusMessage &msg, bool enable)
    { msg.d_ptr->parametersValidated = enable; }

    static DBusMessage *toDBusMessage(const QDBusMessage &message, QDBusConnection::ConnectionCapabilities capabilities,
                                      QDBusError *error);
    static QDBusMessage fromDBusMessage(DBusMessage *dmsg, QDBusConnection::ConnectionCapabilities capabilities);

    static bool isLocal(const QDBusMessage &msg);
    static QDBusMessage makeLocal(const QDBusConnectionPrivate &conn,
                                  const QDBusMessage &asSent);
    static QDBusMessage makeLocalReply(const QDBusConnectionPrivate &conn,
                                       const QDBusMessage &asSent);
};

QT_END_NAMESPACE

#endif // QT_NO_DBUS
#endif
