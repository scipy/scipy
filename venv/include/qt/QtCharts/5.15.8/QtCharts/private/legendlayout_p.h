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

#ifndef LEGENDLAYOUT_H
#define LEGENDLAYOUT_H
#include <QtWidgets/QGraphicsLayout>
#include <QtCharts/QChartGlobal>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class QLegend;
class LegendMarkerItem;

class Q_CHARTS_PRIVATE_EXPORT LegendLayout : public QGraphicsLayout
{
public:

    LegendLayout(QLegend *legend);
    virtual ~LegendLayout();

    void setGeometry(const QRectF &rect);

    void setOffset(qreal x, qreal y);
    QPointF offset() const;

    void invalidate();
protected:
    QSizeF sizeHint(Qt::SizeHint which, const QSizeF &constraint = QSizeF()) const;
    int count() const { return 0; }
    QGraphicsLayoutItem *itemAt(int) const { return 0; };
    void removeAt(int) {};

private:
    void setAttachedGeometry(const QRectF &rect);
    void setDettachedGeometry(const QRectF &rect);

    struct LegendWidthStruct {
        LegendMarkerItem *item;
        qreal width;
    };
    static bool widthLongerThan(const LegendWidthStruct *item1,
                                const LegendWidthStruct *item2);

private:
    QLegend *m_legend;
    qreal m_offsetX;
    qreal m_offsetY;
    qreal m_minOffsetX;
    qreal m_minOffsetY;
    qreal m_maxOffsetX;
    qreal m_maxOffsetY;
    qreal m_width;
    qreal m_height;
};

QT_CHARTS_END_NAMESPACE

#endif
