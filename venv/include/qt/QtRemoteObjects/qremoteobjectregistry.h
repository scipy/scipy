/****************************************************************************
**
** Copyright (C) 2017 Ford Motor Company
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtRemoteObjects module of the Qt Toolkit.
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

#ifndef QREMOTEOBJECTREGISTRY_P_H
#define QREMOTEOBJECTREGISTRY_P_H

#include <QtRemoteObjects/qremoteobjectreplica.h>

QT_BEGIN_NAMESPACE

class QRemoteObjectRegistryPrivate;

class Q_REMOTEOBJECTS_EXPORT QRemoteObjectRegistry : public QRemoteObjectReplica
{
    Q_OBJECT
    Q_CLASSINFO(QCLASSINFO_REMOTEOBJECT_TYPE, "Registry")

    Q_PROPERTY(QRemoteObjectSourceLocations sourceLocations READ sourceLocations)

public:
    ~QRemoteObjectRegistry() override;
    static void registerMetatypes();

    QRemoteObjectSourceLocations sourceLocations() const;

Q_SIGNALS:
    void remoteObjectAdded(const QRemoteObjectSourceLocation &entry);
    void remoteObjectRemoved(const QRemoteObjectSourceLocation &entry);

protected Q_SLOTS:
    void addSource(const QRemoteObjectSourceLocation &entry);
    void removeSource(const QRemoteObjectSourceLocation &entry);
    void pushToRegistryIfNeeded();

private:
    void initialize() override;

    explicit QRemoteObjectRegistry(QObject *parent = nullptr);
    explicit QRemoteObjectRegistry(QRemoteObjectNode *node, const QString &name, QObject *parent = nullptr);

    Q_DECLARE_PRIVATE(QRemoteObjectRegistry)
    friend class QT_PREPEND_NAMESPACE(QRemoteObjectNode);
};

QT_END_NAMESPACE

#endif
