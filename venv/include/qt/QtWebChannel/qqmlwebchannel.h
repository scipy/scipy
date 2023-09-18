/****************************************************************************
**
** Copyright (C) 2016 Klarälvdalens Datakonsult AB, a KDAB Group company, info@kdab.com, author Milian Wolff <milian.wolff@kdab.com>
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtWebChannel module of the Qt Toolkit.
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

#ifndef QQMLWEBCHANNEL_H
#define QQMLWEBCHANNEL_H

#include <QtWebChannel/QWebChannel>
#include <QtWebChannel/qwebchannelglobal.h>

#include <QtQml/qqml.h>
#include <QtQml/QQmlListProperty>

QT_BEGIN_NAMESPACE

class QQmlWebChannelPrivate;
class QQmlWebChannelAttached;
class Q_WEBCHANNEL_EXPORT QQmlWebChannel : public QWebChannel
{
    Q_OBJECT
    Q_DISABLE_COPY(QQmlWebChannel)

    Q_PROPERTY( QQmlListProperty<QObject> transports READ transports )
    Q_PROPERTY( QQmlListProperty<QObject> registeredObjects READ registeredObjects )

public:
    explicit QQmlWebChannel(QObject *parent = Q_NULLPTR);
    virtual ~QQmlWebChannel();

    Q_INVOKABLE void registerObjects(const QVariantMap &objects);
    QQmlListProperty<QObject> registeredObjects();

    QQmlListProperty<QObject> transports();

    static QQmlWebChannelAttached *qmlAttachedProperties(QObject *obj);

    Q_INVOKABLE void connectTo(QObject *transport);
    Q_INVOKABLE void disconnectFrom(QObject *transport);

private:
    Q_DECLARE_PRIVATE(QQmlWebChannel)
    Q_PRIVATE_SLOT(d_func(), void _q_objectIdChanged(const QString &newId))

    static void registeredObjects_append(QQmlListProperty<QObject> *prop, QObject *item);
    static int registeredObjects_count(QQmlListProperty<QObject> *prop);
    static QObject *registeredObjects_at(QQmlListProperty<QObject> *prop, int index);
    static void registeredObjects_clear(QQmlListProperty<QObject> *prop);

    static void transports_append(QQmlListProperty<QObject> *prop, QObject *item);
    static int transports_count(QQmlListProperty<QObject> *prop);
    static QObject *transports_at(QQmlListProperty<QObject> *prop, int index);
    static void transports_clear(QQmlListProperty<QObject> *prop);
};

QT_END_NAMESPACE

QML_DECLARE_TYPE( QQmlWebChannel )
QML_DECLARE_TYPEINFO( QQmlWebChannel, QML_HAS_ATTACHED_PROPERTIES )

#endif // QQMLWEBCHANNEL_H
