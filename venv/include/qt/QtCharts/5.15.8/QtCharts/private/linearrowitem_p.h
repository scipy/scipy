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

#ifndef LINEARROWITEM_P_H
#define LINEARROWITEM_P_H

#include <private/chartaxiselement_p.h>
#include <private/qabstractaxis_p.h>
#include <QtWidgets/QGraphicsLineItem>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class Q_CHARTS_PRIVATE_EXPORT LineArrowItem: public QGraphicsLineItem
{
public:
    explicit LineArrowItem(ChartAxisElement *axis, QGraphicsItem *parent = 0)
        : QGraphicsLineItem(parent),
          m_axis(axis),
          m_axisOrientation(axis->axis()->orientation())
    {
    }

protected:
    void mousePressEvent(QGraphicsSceneMouseEvent *event)
    {
        m_axis->axisSelected();
        QGraphicsLineItem::mousePressEvent(event);
    }

    QRectF boundingRect() const
    {
        return shape().boundingRect();
    }

    QPainterPath shape() const
    {
        QPainterPath path = QGraphicsLineItem::shape();
        QRectF rect = path.boundingRect();
        path.addRect(rect.adjusted(0, 0, m_axisOrientation != Qt::Horizontal ? 8 : 0, m_axisOrientation != Qt::Vertical ? 8 : 0));
        return path;
    }

private:
    ChartAxisElement *m_axis;
    Qt::Orientation m_axisOrientation;
};

QT_CHARTS_END_NAMESPACE

#endif /* LINEARROWITEM_P_H */
