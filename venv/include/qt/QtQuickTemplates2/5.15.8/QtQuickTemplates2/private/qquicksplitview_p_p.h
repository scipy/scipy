/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
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

#ifndef QQUICKSPLITVIEW_P_P_H
#define QQUICKSPLITVIEW_P_P_H

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

#include <QtQuickTemplates2/private/qquickcontainer_p_p.h>

QT_BEGIN_NAMESPACE

class QQuickSplitView;
class QQuickSplitViewAttached;
class QQuickSplitHandleAttached;

class Q_QUICKTEMPLATES2_PRIVATE_EXPORT QQuickSplitViewPrivate : public QQuickContainerPrivate
{
    Q_DECLARE_PUBLIC(QQuickSplitView)

public:
    void updateFillIndex();
    void layoutResizeSplitItems(qreal &usedWidth, qreal &usedHeight, int &indexBeingResizedDueToDrag);
    void layoutResizeFillItem(QQuickItem *fillItem, qreal &usedWidth, qreal &usedHeight, int indexBeingResizedDueToDrag);
    void layoutPositionItems(const QQuickItem *fillItem);
    void requestLayout();
    void layout();
    void createHandles();
    void createHandleItem(int index);
    void removeExcessHandles();
    void destroyHandles();
    void resizeHandle(QQuickItem *handleItem);
    void resizeHandles();
    void updateHandleVisibilities();
    void updateHoveredHandle(QQuickItem *hoveredItem);
    void setResizing(bool resizing);

    bool isHorizontal() const;
    qreal accumulatedSize(int firstIndex, int lastIndex) const;

    struct EffectiveSizeData {
        qreal effectiveMinimumWidth;
        qreal effectiveMinimumHeight;
        qreal effectivePreferredWidth;
        qreal effectivePreferredHeight;
        qreal effectiveMaximumWidth;
        qreal effectiveMaximumHeight;
    };

    EffectiveSizeData effectiveSizeData(const QQuickItemPrivate *itemPrivate,
        const QQuickSplitViewAttached *attached) const;

    int handleIndexForSplitIndex(int splitIndex) const;

    QQuickItem *getContentItem() override;
    void handlePress(const QPointF &point) override;
    void handleMove(const QPointF &point) override;
    void handleRelease(const QPointF &point) override;

    void itemVisibilityChanged(QQuickItem *item) override;
    void itemImplicitWidthChanged(QQuickItem *item) override;
    void itemImplicitHeightChanged(QQuickItem *item) override;

    void updatePolish() override;

    static QQuickSplitViewPrivate *get(QQuickSplitView *splitView);

    Qt::Orientation m_orientation = Qt::Horizontal;
    QQmlComponent *m_handle = nullptr;
    QVector<QQuickItem*> m_handleItems;
    int m_hoveredHandleIndex = -1;
    int m_pressedHandleIndex = -1;
    int m_nextVisibleIndexAfterPressedHandle = -1;
    QPointF m_pressPos;
    QPointF m_mousePos;
    QPointF m_handlePosBeforePress;
    qreal m_leftOrTopItemSizeBeforePress = 0.0;
    qreal m_rightOrBottomItemSizeBeforePress = 0.0;
    int m_fillIndex = -1;
    bool m_layingOut = false;
    bool m_ignoreNextLayoutRequest = false;
    bool m_resizing = false;
};

class QQuickSplitViewAttachedPrivate : public QObjectPrivate
{
    Q_DECLARE_PUBLIC(QQuickSplitViewAttached)

public:
    QQuickSplitViewAttachedPrivate();

    void setView(QQuickSplitView *newView);
    void requestLayoutView();

    static QQuickSplitViewAttachedPrivate *get(QQuickSplitViewAttached *attached);
    static const QQuickSplitViewAttachedPrivate *get(const QQuickSplitViewAttached *attached);

    QQuickItem *m_splitItem = nullptr;
    QQuickSplitView *m_splitView = nullptr;

    unsigned m_fillWidth : 1;
    unsigned m_fillHeight : 1;
    unsigned m_isFillWidthSet : 1;
    unsigned m_isFillHeightSet : 1;
    unsigned m_isMinimumWidthSet : 1;
    unsigned m_isMinimumHeightSet : 1;
    unsigned m_isPreferredWidthSet : 1;
    unsigned m_isPreferredHeightSet : 1;
    unsigned m_isMaximumWidthSet : 1;
    unsigned m_isMaximumHeightSet : 1;
    qreal m_minimumWidth;
    qreal m_minimumHeight;
    qreal m_preferredWidth;
    qreal m_preferredHeight;
    qreal m_maximumWidth;
    qreal m_maximumHeight;
};

class QQuickSplitHandleAttachedPrivate : public QObjectPrivate
{
    Q_DECLARE_PUBLIC(QQuickSplitHandleAttached)

public:
    QQuickSplitHandleAttachedPrivate();

    void setHovered(bool hovered);
    void setPressed(bool pressed);

    static QQuickSplitHandleAttachedPrivate *get(QQuickSplitHandleAttached *attached);
    static const QQuickSplitHandleAttachedPrivate *get(const QQuickSplitHandleAttached *attached);

    unsigned m_hovered : 1;
    unsigned m_pressed : 1;
};

QT_END_NAMESPACE

#endif // QQUICKSPLITVIEW_P_P_H
