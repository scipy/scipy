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

#ifndef CANDLESTICK_P_H
#define CANDLESTICK_P_H

#include <QtGui/QBrush>
#include <QtGui/QPainterPath>
#include <QtGui/QPen>
#include <QtWidgets/QGraphicsObject>
#include <private/candlestickdata_p.h>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class AbstractDomain;
class QCandlestickSet;

class Q_CHARTS_PRIVATE_EXPORT Candlestick : public QGraphicsObject
{
    Q_OBJECT

public:
    Candlestick(QCandlestickSet *set, AbstractDomain *domain, QGraphicsObject *parent);
    ~Candlestick();

    void setTimePeriod(qreal timePeriod);
    void setMaximumColumnWidth(qreal maximumColumnWidth);
    void setMinimumColumnWidth(qreal minimumColumnWidth);
    void setBodyWidth(qreal bodyWidth);
    void setBodyOutlineVisible(bool bodyOutlineVisible);
    void setCapsWidth(qreal capsWidth);
    void setCapsVisible(bool capsVisible);
    void setIncreasingColor(const QColor &color);
    void setDecreasingColor(const QColor &color);
    void setBrush(const QBrush &brush);
    void setPen(const QPen &pen);
    void setLayout(const CandlestickData &data);

    void mousePressEvent(QGraphicsSceneMouseEvent *event);
    void hoverEnterEvent(QGraphicsSceneHoverEvent *event);
    void hoverLeaveEvent(QGraphicsSceneHoverEvent *event);
    void mouseReleaseEvent(QGraphicsSceneMouseEvent *event);
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent *event);

    QRectF boundingRect() const;
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,QWidget *widget = nullptr);

Q_SIGNALS:
    void clicked(QCandlestickSet *set);
    void hovered(bool status, QCandlestickSet *set);
    void pressed(QCandlestickSet *set);
    void released(QCandlestickSet *set);
    void doubleClicked(QCandlestickSet *set);

private:
    void updateGeometry(AbstractDomain *domain);

private:
    QCandlestickSet *m_set;
    AbstractDomain *m_domain;
    qreal m_timePeriod;
    qreal m_maximumColumnWidth;
    qreal m_minimumColumnWidth;
    qreal m_bodyWidth;
    bool m_bodyOutlineVisible;
    qreal m_capsWidth;
    bool m_capsVisible;
    QColor m_increasingColor;
    QColor m_decreasingColor;
    QBrush m_brush;
    QPen m_pen;
    CandlestickData m_data;
    bool m_hovering;
    bool m_mousePressed;
    QRectF m_boundingRect;
    QRectF m_bodyRect;
    QPainterPath m_wicksPath;
    QPainterPath m_capsPath;

    friend class CandlestickAnimation;
    friend class CandlestickChartItem;
};

QT_CHARTS_END_NAMESPACE

#endif // CANDLESTICK_P_H
