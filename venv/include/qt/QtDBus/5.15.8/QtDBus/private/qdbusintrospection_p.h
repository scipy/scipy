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

#ifndef QDBUSINTROSPECTION_P_H
#define QDBUSINTROSPECTION_P_H

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
#include <QtCore/qstring.h>
#include <QtCore/qvector.h>
#include <QtCore/qstringlist.h>
#include <QtCore/qmap.h>
#include <QtCore/qpair.h>
#include <QtCore/qshareddata.h>

#ifndef QT_NO_DBUS

QT_BEGIN_NAMESPACE

class Q_DBUS_EXPORT QDBusIntrospection
{
public:
    // forward declarations
    struct Argument;
    struct Method;
    struct Signal;
    struct Property;
    struct Interface;
    struct Object;
    struct ObjectTree;

    // typedefs
    typedef QMap<QString, QString> Annotations;
    typedef QVector<Argument> Arguments;
    typedef QMultiMap<QString, Method> Methods;
    typedef QMultiMap<QString, Signal> Signals;
    typedef QMap<QString, Property> Properties;
    typedef QMap<QString, QSharedDataPointer<Interface> > Interfaces;
    typedef QMap<QString, QSharedDataPointer<ObjectTree> > Objects;

public:
    // the structs

    struct Argument
    {
        QString type;
        QString name;

        inline bool operator==(const Argument& other) const
        { return name == other.name && type == other.type; }
    };

    struct Method
    {
        QString name;
        Arguments inputArgs;
        Arguments outputArgs;
        Annotations annotations;

        inline bool operator==(const Method& other) const
        { return name == other.name && annotations == other.annotations &&
                inputArgs == other.inputArgs && outputArgs == other.outputArgs; }
    };

    struct Signal
    {
        QString name;
        Arguments outputArgs;
        Annotations annotations;

        inline bool operator==(const Signal& other) const
        { return name == other.name && annotations == other.annotations &&
                outputArgs == other.outputArgs; }
    };

    struct Property
    {
        enum Access { Read, Write, ReadWrite };
        QString name;
        QString type;
        Access access;
        Annotations annotations;

        inline bool operator==(const Property& other) const
        { return access == other.access && name == other.name &&
                annotations == other.annotations && type == other.type; }
    };

    struct Interface: public QSharedData
    {
        QString name;
        QString introspection;

        Annotations annotations;
        Methods methods;
        Signals signals_;
        Properties properties;

        inline bool operator==(const Interface &other) const
        { return !name.isEmpty() && name == other.name; }
    };

    struct Object: public QSharedData
    {
        QString service;
        QString path;

        QStringList interfaces;
        QStringList childObjects;
    };

public:
    static Interface parseInterface(const QString &xml);
    static Interfaces parseInterfaces(const QString &xml);
    static Object parseObject(const QString &xml, const QString &service = QString(),
                              const QString &path = QString());

private:
    QDBusIntrospection();
};

QT_END_NAMESPACE

#endif // QT_NO_DBUS
#endif
