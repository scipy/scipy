/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Templates 2 module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QQUICKACTIONGROUP_P_H
#define QQUICKACTIONGROUP_P_H

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
#include <QtQuickTemplates2/private/qtquicktemplates2global_p.h>
#include <QtQml/qqml.h>

QT_BEGIN_NAMESPACE

class QQuickAction;
class QQuickActionGroupPrivate;
class QQuickActionGroupAttached;
class QQuickActionGroupAttachedPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickActionGroup : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQuickAction *checkedAction READ checkedAction WRITE setCheckedAction NOTIFY checkedActionChanged FINAL)
    Q_PROPERTY(QQmlListProperty<QQuickAction> actions READ actions NOTIFY actionsChanged FINAL)
    Q_PROPERTY(bool exclusive READ isExclusive WRITE setExclusive NOTIFY exclusiveChanged FINAL)
    Q_PROPERTY(bool enabled READ isEnabled WRITE setEnabled NOTIFY enabledChanged FINAL)
    Q_CLASSINFO("DefaultProperty", "actions")

public:
    explicit QQuickActionGroup(QObject *parent = nullptr);
    ~QQuickActionGroup();

    static QQuickActionGroupAttached *qmlAttachedProperties(QObject *object);

    QQuickAction *checkedAction() const;
    void setCheckedAction(QQuickAction *checkedAction);

    QQmlListProperty<QQuickAction> actions();

    bool isExclusive() const;
    void setExclusive(bool exclusive);

    bool isEnabled() const;
    void setEnabled(bool enabled);

public Q_SLOTS:
    void addAction(QQuickAction *action);
    void removeAction(QQuickAction *action);

Q_SIGNALS:
    void checkedActionChanged();
    void actionsChanged();
    void exclusiveChanged();
    void enabledChanged();
    void triggered(QQuickAction *action);

private:
    Q_DISABLE_COPY(QQuickActionGroup)
    Q_DECLARE_PRIVATE(QQuickActionGroup)

    Q_PRIVATE_SLOT(d_func(), void _q_updateCurrent())
};

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickActionGroupAttached : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQuickActionGroup *group READ group WRITE setGroup NOTIFY groupChanged FINAL)

public:
    explicit QQuickActionGroupAttached(QObject *parent = nullptr);

    QQuickActionGroup *group() const;
    void setGroup(QQuickActionGroup *group);

Q_SIGNALS:
    void groupChanged();

private:
    Q_DISABLE_COPY(QQuickActionGroupAttached)
    Q_DECLARE_PRIVATE(QQuickActionGroupAttached)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickActionGroup)
QML_DECLARE_TYPEINFO(QQuickActionGroup, QML_HAS_ATTACHED_PROPERTIES)

#endif // QQUICKACTIONGROUP_P_H
