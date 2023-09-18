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

#ifndef BOXWHISKERS_H
#define BOXWHISKERS_H

#include <private/boxwhiskersdata_p.h>
#include <QtCharts/QChartGlobal>
#include <QtCharts/private/qchartglobal_p.h>
#include <private/abstractdomain_p.h>
#include <QtCharts/QBoxSet>
#include <QtWidgets/QGraphicsRectItem>
#include <QtWidgets/QGraphicsLineItem>
#include <QtWidgets/QGraphicsLayoutItem>
#include <QtGui/QPainterPath>

QT_CHARTS_BEGIN_NAMESPACE

class QBarSet;

class Q_CHARTS_PRIVATE_EXPORT BoxWhiskers : public QGraphicsObject
{
    Q_OBJECT

public:
    BoxWhiskers(QBoxSet *set, AbstractDomain *domain, QGraphicsObject *parent);
    ~BoxWhiskers();

    void setBrush(const QBrush &brush);
    void setPen(const QPen &pen);
    void setLayout(const BoxWhiskersData &data);
    void setBoxOutlined(const bool outlined) { m_boxOutlined = outlined; }
    void setBoxWidth(const qreal width);

    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event);
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event);

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = 0);

    void updateGeometry(AbstractDomain *domain);
protected:
    QSizeF sizeHint(Qt::SizeHint which, const QSizeF &constraint) const;
    void setGeometry(const QRectF &rect);

Q_SIGNALS:
    void clicked(QBoxSet *boxset);
    void hovered(bool status, QBoxSet *boxset);
    void pressed(QBoxSet *boxset);
    void released(QBoxSet *boxset);
    void doubleClicked(QBoxSet *boxset);

private:
    friend class BoxPlotChartItem;
    friend class BoxPlotAnimation;

    QBoxSet *m_boxSet;
    AbstractDomain *m_domain;
    QPainterPath    m_boxPath;
    QRectF m_boundingRect;
    bool m_hovering;
    bool m_validData;
    QBrush  m_brush;
    QPen m_pen;
    QPen m_medianPen;
    QPen m_outlinePen;
    bool m_boxOutlined;
    qreal m_boxWidth;
    BoxWhiskersData m_data;
    QSizeF m_domainSize;
    QRectF m_middleBox;
    qreal m_geometryMedian;
    qreal m_geometryLeft;
    qreal m_geometryRight;

    bool m_mousePressed;
};

QT_CHARTS_END_NAMESPACE

#endif // BOXWHISKERS_H
