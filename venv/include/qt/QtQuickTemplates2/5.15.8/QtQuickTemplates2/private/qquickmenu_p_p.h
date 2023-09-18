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

#ifndef QQUICKMENU_P_P_H
#define QQUICKMENU_P_P_H

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

#include <QtCore/qvector.h>
#include <QtCore/qpointer.h>

#include <QtQuickTemplates2/private/qquickmenu_p.h>
#include <QtQuickTemplates2/private/qquickpopup_p_p.h>

QT_BEGIN_NAMESPACE

class QQuickAction;
class QQmlComponent;
class QQmlObjectModel;
class QQuickMenuItem;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickMenuPrivate : public QQuickPopupPrivate
{
    Q_DECLARE_PUBLIC(QQuickMenu)

public:
    QQuickMenuPrivate();

    static QQuickMenuPrivate *get(QQuickMenu *menu)
    {
        return menu->d_func();
    }

    void init();

    QQuickItem *itemAt(int index) const;
    void insertItem(int index, QQuickItem *item);
    void moveItem(int from, int to);
    void removeItem(int index, QQuickItem *item);

    QQuickItem *beginCreateItem();
    void completeCreateItem();

    QQuickItem *createItem(QQuickMenu *menu);
    QQuickItem *createItem(QQuickAction *action);

    void resizeItem(QQuickItem *item);
    void resizeItems();

    void itemChildAdded(QQuickItem *item, QQuickItem *child) override;
    void itemSiblingOrderChanged(QQuickItem *item) override;
    void itemParentChanged(QQuickItem *item, QQuickItem *parent) override;
    void itemDestroyed(QQuickItem *item) override;
    void itemGeometryChanged(QQuickItem *, QQuickGeometryChange change, const QRectF &diff) override;

    QQuickPopupPositioner *getPositioner() override;
    bool prepareEnterTransition() override;
    bool prepareExitTransition() override;
    bool blockInput(QQuickItem *item, const QPointF &point) const override;

    void onItemHovered();
    void onItemTriggered();
    void onItemActiveFocusChanged();

    QQuickMenu *currentSubMenu() const;
    void setParentMenu(QQuickMenu *parent);
    void resolveParentItem();

    void propagateKeyEvent(QKeyEvent *event);

    void startHoverTimer();
    void stopHoverTimer();

    void setCurrentIndex(int index, Qt::FocusReason reason);
    bool activateNextItem();
    bool activatePreviousItem();

    QQuickMenuItem *firstEnabledMenuItem() const;

    static void contentData_append(QQmlListProperty<QObject> *prop, QObject *obj);
    static int contentData_count(QQmlListProperty<QObject> *prop);
    static QObject *contentData_at(QQmlListProperty<QObject> *prop, int index);
    static void contentData_clear(QQmlListProperty<QObject> *prop);

    bool cascade = false;
    int hoverTimer = 0;
    int currentIndex = -1;
    qreal overlap = 0;
    QPointer<QQuickMenu> parentMenu;
    QPointer<QQuickMenuItem> currentItem;
    QQuickItem *contentItem = nullptr; // TODO: cleanup
    QVector<QObject *> contentData;
    QQmlObjectModel *contentModel;
    QQmlComponent *delegate = nullptr;
    QString title;
};

QT_END_NAMESPACE

#endif // QQUICKMENU_P_P_H
