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

#ifndef POLARCHARTAXISANGULAR_P_H
#define POLARCHARTAXISANGULAR_P_H

#include <QtCharts/qvalueaxis.h>
#include <private/polarchartaxis_p.h>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class Q_CHARTS_PRIVATE_EXPORT PolarChartAxisAngular : public PolarChartAxis
{
    Q_OBJECT

public:
    PolarChartAxisAngular(QAbstractAxis *axis, QGraphicsItem *item, bool intervalAxis = false);
    ~PolarChartAxisAngular();

    Qt::Orientation orientation() const;
    QSizeF sizeHint(Qt::SizeHint which, const QSizeF &constraint = QSizeF()) const;

    virtual void updateGeometry();
    virtual void createItems(int count);
    virtual void updateMinorTickItems();

    qreal preferredAxisRadius(const QSizeF &maxSize);

public Q_SLOTS:
    virtual void handleArrowPenChanged(const QPen &pen);
    virtual void handleGridPenChanged(const QPen &pen);
    virtual void handleMinorArrowPenChanged(const QPen &pen);
    virtual void handleMinorGridPenChanged(const QPen &pen);
    virtual void handleGridLineColorChanged(const QColor &color);
    virtual void handleMinorGridLineColorChanged(const QColor &color);

private:
    QRectF moveLabelToPosition(qreal angularCoordinate, QPointF labelPoint, QRectF labelRect) const;
    void updateMinorTickGeometry();
};

QT_CHARTS_END_NAMESPACE

#endif // POLARCHARTAXISANGULAR_P_H
