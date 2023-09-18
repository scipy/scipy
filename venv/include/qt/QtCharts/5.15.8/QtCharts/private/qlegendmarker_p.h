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

#ifndef QLEGENDMARKERPRIVATE_H
#define QLEGENDMARKERPRIVATE_H

#include <QtCharts/QChartGlobal>
#include <QGraphicsObject>
#include <QtGui/QBrush>
#include <QtGui/QPen>
#include <QtWidgets/QGraphicsLayoutItem>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class QAbstractSeries;
class QLegend;

class QLegendMarker;
class LegendMarkerItem;

class Q_CHARTS_PRIVATE_EXPORT QLegendMarkerPrivate : public QObject
{
    Q_OBJECT
public:
    explicit QLegendMarkerPrivate(QLegendMarker *q, QLegend *legend);
    virtual ~QLegendMarkerPrivate();

    // Helper for now. (or declare LegendLayout as friend)
    LegendMarkerItem* item() const { return m_item; }

    virtual QAbstractSeries* series() = 0;
    virtual QObject* relatedObject() = 0;

    void invalidateLegend();
    void invalidateAllItems();

public Q_SLOTS:
    virtual void updated() = 0;
    void handleShapeChange();

protected:
    LegendMarkerItem *m_item;
    QLegend *m_legend;
    bool m_customLabel;
    bool m_customBrush;
    bool m_customPen;

private:
    QLegendMarker *q_ptr;

    friend class QLegendPrivate;
    friend class LegendMarkerItem;
    Q_DECLARE_PUBLIC(QLegendMarker)
};

QT_CHARTS_END_NAMESPACE

#endif // QLEGENDMARKERPRIVATE_H
