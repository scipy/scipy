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

#ifndef QQUICKACTION_P_P_H
#define QQUICKACTION_P_P_H

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

#include <QtCore/private/qobject_p.h>
#include <QtCore/qvariant.h>
#include <QtCore/qstring.h>
#if QT_CONFIG(shortcut)
#  include <QtGui/qkeysequence.h>
#endif
#include <QtQuick/private/qquickitemchangelistener_p.h>

QT_BEGIN_NAMESPACE

class QShortcutEvent;
class QQuickActionGroup;

class QQuickActionPrivate : public QObjectPrivate, public QQuickItemChangeListener
{
    Q_DECLARE_PUBLIC(QQuickAction)

public:
    static QQuickActionPrivate *get(QQuickAction *action)
    {
        return action->d_func();
    }

#if QT_CONFIG(shortcut)
    QVariant shortcut() const;
    void setShortcut(const QVariant &shortcut);
#endif

    void setEnabled(bool enable);

    bool watchItem(QQuickItem *item);
    bool unwatchItem(QQuickItem *item);

    void registerItem(QQuickItem *item);
    void unregisterItem(QQuickItem *item);

    void itemVisibilityChanged(QQuickItem *item) override;
    void itemDestroyed(QQuickItem *item) override;

    bool handleShortcutEvent(QObject *object, QShortcutEvent *event);

    void trigger(QObject*, bool doToggle);

#if QT_CONFIG(shortcut)
    class ShortcutEntry
    {
    public:
        explicit ShortcutEntry(QObject *target);
        ~ShortcutEntry();

        QObject *target() const;
        int shortcutId() const;

        void grab(const QKeySequence &vshortcut, bool enabled);
        void ungrab();

        void setEnabled(bool enabled);

    private:
        int m_shortcutId = 0;
        QObject *m_target = nullptr;
    };

    ShortcutEntry *findShortcutEntry(QObject *target) const;
    void updateDefaultShortcutEntry();
#endif // QT_CONFIG(shortcut)

    bool explicitEnabled = false;
    bool enabled = true;
    bool checked = false;
    bool checkable = false;
    QString text;
    QQuickIcon icon;
#if QT_CONFIG(shortcut)
    QKeySequence keySequence;
    QVariant vshortcut;
    ShortcutEntry *defaultShortcutEntry = nullptr;
    QVector<ShortcutEntry *> shortcutEntries;
#endif
    QQuickActionGroup *group = nullptr;
};

QT_END_NAMESPACE

#endif // QQUICKACTION_P_P_H
