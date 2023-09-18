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

#ifndef QQUICKSWIPEVIEW_P_H
#define QQUICKSWIPEVIEW_P_H

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

#include <QtQuickTemplates2/private/qquickcontainer_p.h>

QT_BEGIN_NAMESPACE

class QQuickSwipeViewAttached;
class QQuickSwipeViewPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickSwipeView : public QQuickContainer
{
    Q_OBJECT
    // 2.1 (Qt 5.8)
    Q_PROPERTY(bool interactive READ isInteractive WRITE setInteractive NOTIFY interactiveChanged FINAL REVISION 1)
    // 2.2 (Qt 5.9)
    Q_PROPERTY(Qt::Orientation orientation READ orientation WRITE setOrientation NOTIFY orientationChanged FINAL REVISION 2)
    // 2.3 (Qt 5.10)
    Q_PROPERTY(bool horizontal READ isHorizontal NOTIFY orientationChanged FINAL REVISION 3)
    Q_PROPERTY(bool vertical READ isVertical NOTIFY orientationChanged FINAL REVISION 3)

public:
    explicit QQuickSwipeView(QQuickItem *parent = nullptr);

    static QQuickSwipeViewAttached *qmlAttachedProperties(QObject *object);

    // 2.1 (Qt 5.8)
    bool isInteractive() const;
    void setInteractive(bool interactive);

    // 2.2 (Qt 5.9)
    Qt::Orientation orientation() const;
    void setOrientation(Qt::Orientation orientation);

    // 2.3 (Qt 5.10)
    bool isHorizontal() const;
    bool isVertical() const;

Q_SIGNALS:
    // 2.1 (Qt 5.8)
    Q_REVISION(1) void interactiveChanged();
    // 2.2 (Qt 5.9)
    Q_REVISION(2) void orientationChanged();

protected:
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;
    void itemAdded(int index, QQuickItem *item) override;
    void itemMoved(int index, QQuickItem *item) override;
    void itemRemoved(int index, QQuickItem *item) override;

#if QT_CONFIG(accessibility)
    QAccessible::Role accessibleRole() const override;
#endif

private:
    Q_DISABLE_COPY(QQuickSwipeView)
    Q_DECLARE_PRIVATE(QQuickSwipeView)
};

class QQuickSwipeViewAttachedPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickSwipeViewAttached : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int index READ index NOTIFY indexChanged FINAL)
    Q_PROPERTY(bool isCurrentItem READ isCurrentItem NOTIFY isCurrentItemChanged FINAL)
    Q_PROPERTY(QQuickSwipeView *view READ view NOTIFY viewChanged FINAL)
    // 2.1 (Qt 5.8)
    Q_PROPERTY(bool isNextItem READ isNextItem NOTIFY isNextItemChanged FINAL REVISION 1)
    Q_PROPERTY(bool isPreviousItem READ isPreviousItem NOTIFY isPreviousItemChanged FINAL REVISION 1)

public:
    explicit QQuickSwipeViewAttached(QObject *parent = nullptr);

    int index() const;
    bool isCurrentItem() const;
    QQuickSwipeView *view() const;

    // 2.1 (Qt 5.8)
    bool isNextItem() const;
    bool isPreviousItem() const;

Q_SIGNALS:
    void indexChanged();
    void isCurrentItemChanged();
    void viewChanged();
    // 2.1 (Qt 5.8)
    /*Q_REVISION(1)*/ void isNextItemChanged();
    /*Q_REVISION(1)*/ void isPreviousItemChanged();

private:
    Q_DISABLE_COPY(QQuickSwipeViewAttached)
    Q_DECLARE_PRIVATE(QQuickSwipeViewAttached)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickSwipeView)
QML_DECLARE_TYPEINFO(QQuickSwipeView, QML_HAS_ATTACHED_PROPERTIES)

#endif // QQUICKSWIPEVIEW_P_H
