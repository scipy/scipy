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

#ifndef QQUICKTUMBLER_P_P_H
#define QQUICKTUMBLER_P_P_H

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

#include <QtQuickTemplates2/private/qquickcontrol_p_p.h>
#include <QtQuickTemplates2/private/qquicktumbler_p.h>

QT_BEGIN_NAMESPACE

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickTumblerPrivate : public QQuickControlPrivate
{
    Q_DECLARE_PUBLIC(QQuickTumbler)

public:
    enum ContentItemType {
        NoContentItem,
        UnsupportedContentItemType,
        PathViewContentItem,
        ListViewContentItem
    };

    QQuickItem *determineViewType(QQuickItem *contentItem);
    void resetViewData();
    QList<QQuickItem *> viewContentItemChildItems() const;

    static QQuickTumblerPrivate *get(QQuickTumbler *tumbler);

    QVariant model;
    QQmlComponent *delegate = nullptr;
    int visibleItemCount = 5;
    bool wrap = true;
    bool explicitWrap = false;
    bool modelBeingSet = false;
    bool currentIndexSetDuringModelChange = false;
    QQuickItem *view = nullptr;
    QQuickItem *viewContentItem = nullptr;
    ContentItemType viewContentItemType = UnsupportedContentItemType;
    union {
        qreal viewOffset; // PathView
        qreal viewContentY; // ListView
    };
    int currentIndex = -1;
    int pendingCurrentIndex = -1;
    bool ignoreCurrentIndexChanges = false;
    int count = 0;
    bool ignoreSignals = false;

    void _q_updateItemHeights();
    void _q_updateItemWidths();
    void _q_onViewCurrentIndexChanged();
    void _q_onViewCountChanged();
    void _q_onViewOffsetChanged();
    void _q_onViewContentYChanged();

    void calculateDisplacements();

    void disconnectFromView();
    void setupViewData(QQuickItem *newControlContentItem);
    void warnAboutIncorrectContentItem();
    void syncCurrentIndex();
    void setPendingCurrentIndex(int index);

    enum PropertyChangeReason {
        UserChange,
        InternalChange
    };

    static QString propertyChangeReasonToString(PropertyChangeReason changeReason);

    void setCurrentIndex(int newCurrentIndex, PropertyChangeReason changeReason = InternalChange);
    void setCount(int newCount);
    void setWrapBasedOnCount();
    void setWrap(bool shouldWrap, bool isExplicit);
    void beginSetModel();
    void endSetModel();

    void itemChildAdded(QQuickItem *, QQuickItem *) override;
    void itemChildRemoved(QQuickItem *, QQuickItem *) override;
    void itemGeometryChanged(QQuickItem *, QQuickGeometryChange , const QRectF &) override;
};

class QQuickTumblerAttachedPrivate : public QObjectPrivate
{
    Q_DECLARE_PUBLIC(QQuickTumblerAttached)

public:
    static QQuickTumblerAttachedPrivate *get(QQuickTumblerAttached *attached)
    {
        return attached->d_func();
    }

    void init(QQuickItem *delegateItem);

    void calculateDisplacement();
    void emitIfDisplacementChanged(qreal oldDisplacement, qreal newDisplacement);

    // The Tumbler that contains the delegate. Required to calculated the displacement.
    QPointer<QQuickTumbler> tumbler;
    // The index of the delegate. Used to calculate the displacement.
    int index = -1;
    // The displacement for our delegate.
    qreal displacement = 0;
};

QT_END_NAMESPACE

#endif // QQUICKTUMBLER_P_P_H
