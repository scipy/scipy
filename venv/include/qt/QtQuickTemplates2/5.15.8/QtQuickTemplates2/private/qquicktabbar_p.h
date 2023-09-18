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

#ifndef QQUICKTABBAR_P_H
#define QQUICKTABBAR_P_H

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

class QQuickTabBarPrivate;
class QQuickTabBarAttached;
class QQuickTabBarAttachedPrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickTabBar : public QQuickContainer
{
    Q_OBJECT
    Q_PROPERTY(Position position READ position WRITE setPosition NOTIFY positionChanged FINAL)
    // 2.2 (Qt 5.9)
    Q_PROPERTY(qreal contentWidth READ contentWidth WRITE setContentWidth RESET resetContentWidth NOTIFY contentWidthChanged FINAL REVISION 2) // re-declare QQuickContainer::contentWidth (REV 5)
    Q_PROPERTY(qreal contentHeight READ contentHeight WRITE setContentHeight RESET resetContentHeight NOTIFY contentHeightChanged FINAL REVISION 2) // re-declare QQuickContainer::contentHeight (REV 5)

public:
    explicit QQuickTabBar(QQuickItem *parent = nullptr);

    enum Position {
        Header,
        Footer
    };
    Q_ENUM(Position)

    Position position() const;
    void setPosition(Position position);

    static QQuickTabBarAttached *qmlAttachedProperties(QObject *object);

Q_SIGNALS:
    void positionChanged();

protected:
    void updatePolish() override;
    void componentComplete() override;
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;
    bool isContent(QQuickItem *item) const override;
    void itemAdded(int index, QQuickItem *item) override;
    void itemMoved(int index, QQuickItem *item) override;
    void itemRemoved(int index, QQuickItem *item) override;

    QFont defaultFont() const override;
    QPalette defaultPalette() const override;

#if QT_CONFIG(accessibility)
    QAccessible::Role accessibleRole() const override;
#endif

private:
    Q_DISABLE_COPY(QQuickTabBar)
    Q_DECLARE_PRIVATE(QQuickTabBar)
};

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickTabBarAttached : public QObject
{
    Q_OBJECT
    Q_PROPERTY(int index READ index NOTIFY indexChanged FINAL)
    Q_PROPERTY(QQuickTabBar *tabBar READ tabBar NOTIFY tabBarChanged FINAL)
    Q_PROPERTY(QQuickTabBar::Position position READ position NOTIFY positionChanged FINAL)

public:
    explicit QQuickTabBarAttached(QObject *parent = nullptr);

    int index() const;
    QQuickTabBar *tabBar() const;
    QQuickTabBar::Position position() const;

Q_SIGNALS:
    void indexChanged();
    void tabBarChanged();
    void positionChanged();

private:
    Q_DISABLE_COPY(QQuickTabBarAttached)
    Q_DECLARE_PRIVATE(QQuickTabBarAttached)
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuickTabBar)
QML_DECLARE_TYPEINFO(QQuickTabBar, QML_HAS_ATTACHED_PROPERTIES)

#endif // QQUICKTABBAR_P_H
