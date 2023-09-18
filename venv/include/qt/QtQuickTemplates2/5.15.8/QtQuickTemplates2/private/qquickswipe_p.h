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

#ifndef QQUICKSWIPE_P_H
#define QQUICKSWIPE_P_H

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
#include <QtQuickTemplates2/private/qquickswipedelegate_p.h>

QT_BEGIN_NAMESPACE

class QQmlComponent;
class QQuickItem;
class QQuickTransition;
class QQuickSwipePrivate;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickSwipe : public QObject
{
    Q_OBJECT
    Q_PROPERTY(qreal position READ position WRITE setPosition NOTIFY positionChanged FINAL)
    Q_PROPERTY(bool complete READ isComplete NOTIFY completeChanged FINAL)
    Q_PROPERTY(QQmlComponent *left READ left WRITE setLeft NOTIFY leftChanged FINAL)
    Q_PROPERTY(QQmlComponent *behind READ behind WRITE setBehind NOTIFY behindChanged FINAL)
    Q_PROPERTY(QQmlComponent *right READ right WRITE setRight NOTIFY rightChanged FINAL)
    Q_PROPERTY(QQuickItem *leftItem READ leftItem NOTIFY leftItemChanged FINAL)
    Q_PROPERTY(QQuickItem *behindItem READ behindItem NOTIFY behindItemChanged FINAL)
    Q_PROPERTY(QQuickItem *rightItem READ rightItem NOTIFY rightItemChanged FINAL)
    // 2.2 (Qt 5.9)
    Q_PROPERTY(bool enabled READ isEnabled WRITE setEnabled NOTIFY enabledChanged FINAL) // REVISION 2
    Q_PROPERTY(QQuickTransition *transition READ transition WRITE setTransition NOTIFY transitionChanged FINAL) // REVISION 2

public:
    explicit QQuickSwipe(QQuickSwipeDelegate *control);

    qreal position() const;
    void setPosition(qreal position);

    bool isComplete() const;
    void setComplete(bool complete);

    QQmlComponent *left() const;
    void setLeft(QQmlComponent *left);

    QQmlComponent *behind() const;
    void setBehind(QQmlComponent *behind);

    QQmlComponent *right() const;
    void setRight(QQmlComponent *right);

    QQuickItem *leftItem() const;
    void setLeftItem(QQuickItem *item);

    QQuickItem *behindItem() const;
    void setBehindItem(QQuickItem *item);

    QQuickItem *rightItem() const;
    void setRightItem(QQuickItem *item);

    // 2.1 (Qt 5.8)
    Q_REVISION(1) Q_INVOKABLE void close();

    // 2.2 (Qt 5.9)
    bool isEnabled() const;
    void setEnabled(bool enabled);

    QQuickTransition *transition() const;
    void setTransition(QQuickTransition *transition);

    Q_REVISION(2) Q_INVOKABLE void open(QQuickSwipeDelegate::Side side);

Q_SIGNALS:
    void positionChanged();
    void completeChanged();
    void leftChanged();
    void behindChanged();
    void rightChanged();
    void leftItemChanged();
    void behindItemChanged();
    void rightItemChanged();
    // 2.1 (Qt 5.8)
    /*Q_REVISION(1)*/ void completed();
    // 2.2 (Qt 5.9)
    /*Q_REVISION(2)*/ void opened();
    /*Q_REVISION(2)*/ void closed();
    /*Q_REVISION(2)*/ void enabledChanged();
    /*Q_REVISION(2)*/ void transitionChanged();

private:
    Q_DISABLE_COPY(QQuickSwipe)
    Q_DECLARE_PRIVATE(QQuickSwipe)
};

QT_END_NAMESPACE

#endif // QQUICKSWIPE_P_H
