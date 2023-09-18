/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
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

#ifndef PROFILECONTEXT_P_H
#define PROFILECONTEXT_P_H

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

#include <QtCore/qobject.h>
#include <QtDBus/qdbuscontext.h>
#include <QtDBus/qdbusextratypes.h>
#include <QtDBus/qdbusunixfiledescriptor.h>

QT_BEGIN_NAMESPACE

class OrgBluezProfile1ContextInterface : public QObject, protected QDBusContext
{
    Q_OBJECT
    Q_CLASSINFO("D-Bus Interface", "org.bluez.Profile1")
    Q_CLASSINFO("D-Bus Introspection", ""
"  <interface name=\"org.bluez.Profile1\">\n"
"    <method name=\"NewConnection\">\n"
"      <arg direction=\"in\" type=\"o\"/>\n"
"      <arg direction=\"in\" type=\"h\"/>\n"
"      <arg direction=\"in\" type=\"a{sv}\"/>\n"
"      <annotation value=\"QVariantMap\" name=\"org.qtproject.QtDBus.QtTypeName.In2\"/>\n"
"    </method>\n"
"    <method name=\"RequestDisconnection\">\n"
"      <arg direction=\"in\" type=\"o\"/>\n"
"    </method>\n"
"    <method name=\"Release\">\n"
"      <annotation value=\"true\" name=\"org.freedesktop.DBus.Method.NoReply\"/>\n"
"    </method>\n"
"  </interface>\n"
        "")
public:
    explicit OrgBluezProfile1ContextInterface(QObject *parent = nullptr);

Q_SIGNALS:
    void newConnection(const QDBusUnixFileDescriptor &fd);

public Q_SLOTS:
    void NewConnection(const QDBusObjectPath &, const QDBusUnixFileDescriptor &,
                       const QVariantMap &);
    void RequestDisconnection(const QDBusObjectPath &);
    Q_NOREPLY void Release();
};

QT_END_NAMESPACE

#endif // PROFILECONTEXT_P_H
