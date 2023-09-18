/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt Charts module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

//  W A R N I N G
//  -------------
//
// This file is not part of the Qt Chart API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef SCATTERCHARTITEM_H
#define SCATTERCHARTITEM_H

#include <QtCharts/QChartGlobal>
#include <private/xychart_p.h>
#include <QtWidgets/QGraphicsEllipseItem>
#include <QtGui/QPen>
#include <QtWidgets/QGraphicsSceneMouseEvent>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class QScatterSeries;

class Q_CHARTS_PRIVATE_EXPORT ScatterChartItem : public XYChart
{
    Q_OBJECT
    Q_INTERFACES(QGraphicsItem)
public:
    explicit ScatterChartItem(QScatterSeries *series, QGraphicsItem *item = 0);

public:
    //from QGraphicsItem
    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);

    void setPen(const QPen &pen);
    void setBrush(const QBrush &brush);

    void markerSelected(QGraphicsItem *item);
    void markerHovered(QGraphicsItem *item, bool state);
    void markerPressed(QGraphicsItem *item);
    void markerReleased(QGraphicsItem *item);
    void markerDoubleClicked(QGraphicsItem *item);

    void setMousePressed(bool pressed = true) {m_mousePressed = pressed;}
    bool mousePressed() {return m_mousePressed;}


public Q_SLOTS:
    void handleUpdated();

private:
    void createPoints(int count);
    void deletePoints(int count);

protected:
    void updateGeometry();

private:
    QScatterSeries *m_series;
    QGraphicsItemGroup m_items;
    bool m_visible;
    int m_shape;
    int m_size;
    QRectF m_rect;
    QMap<QGraphicsItem *, QPointF> m_markerMap;

    bool m_pointLabelsVisible;
    QString m_pointLabelsFormat;
    QFont m_pointLabelsFont;
    QColor m_pointLabelsColor;
    bool m_pointLabelsClipping;

    bool m_mousePressed;
};

class Q_CHARTS_PRIVATE_EXPORT CircleMarker: public QGraphicsEllipseItem
{

public:
    CircleMarker(qreal x, qreal y, qreal w, qreal h, ScatterChartItem *parent)
        : QGraphicsEllipseItem(x, y, w, h, parent),
          m_parent(parent)
    {
        setAcceptHoverEvents(true);
        setFlag(QGraphicsItem::ItemIsSelectable);
    }

protected:
    void mousePressEvent(QGraphicsSceneMouseEvent *event)
    {
        QGraphicsEllipseItem::mousePressEvent(event);
        m_parent->markerPressed(this);
        m_parent->setMousePressed();
    }
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event)
    {
        QGraphicsEllipseItem::hoverEnterEvent(event);
        m_parent->markerHovered(this, true);
    }
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
    {
        QGraphicsEllipseItem::hoverLeaveEvent(event);
        m_parent->markerHovered(this, false);
    }
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
    {
        QGraphicsEllipseItem::mouseReleaseEvent(event);
        m_parent->markerReleased(this);
        if (m_parent->mousePressed())
            m_parent->markerSelected(this);
        m_parent->setMousePressed(false);
    }
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event)
    {
        QGraphicsEllipseItem::mouseDoubleClickEvent(event);
        m_parent->markerDoubleClicked(this);
    }

private:
    ScatterChartItem *m_parent;
};

class Q_CHARTS_PRIVATE_EXPORT RectangleMarker: public QGraphicsRectItem
{

public:
    RectangleMarker(qreal x, qreal y, qreal w, qreal h, ScatterChartItem *parent)
        : QGraphicsRectItem(x, y, w, h, parent),
          m_parent(parent)
    {
        setAcceptHoverEvents(true);
        setFlag(QGraphicsItem::ItemIsSelectable);
    }

protected:
    void mousePressEvent(QGraphicsSceneMouseEvent *event)
    {
        QGraphicsRectItem::mousePressEvent(event);
        m_parent->markerPressed(this);
        m_parent->setMousePressed();
    }
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event)
    {
        QGraphicsRectItem::hoverEnterEvent(event);
        m_parent->markerHovered(this, true);
    }
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
    {
        QGraphicsRectItem::hoverLeaveEvent(event);
        m_parent->markerHovered(this, false);
    }
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event)
    {
        QGraphicsRectItem::mouseReleaseEvent(event);
        m_parent->markerReleased(this);
        if (m_parent->mousePressed())
            m_parent->markerSelected(this);
        m_parent->setMousePressed(false);
    }
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event)
    {
        QGraphicsRectItem::mouseDoubleClickEvent(event);
        m_parent->markerDoubleClicked(this);
    }

private:
    ScatterChartItem *m_parent;
};

QT_CHARTS_END_NAMESPACE

#endif // SCATTERPRESENTER_H
