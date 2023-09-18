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

#ifndef QCHARTVIEW_P_H
#define QCHARTVIEW_P_H

#include <QtCharts/QChartView>
#include <QtCharts/private/qchartglobal_p.h>

QT_BEGIN_NAMESPACE
class QGraphicsScene;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE

class QChart;
class ChartPresenter;
class QChartView;

class Q_CHARTS_PRIVATE_EXPORT QChartViewPrivate
{
public:
    explicit QChartViewPrivate(QChartView *q, QChart *chart = 0);
    ~QChartViewPrivate();
    void setChart(QChart *chart);
    void resize();

protected:
    QChartView *q_ptr;

public:
    QGraphicsScene *m_scene;
    QChart *m_chart;
    QPoint m_rubberBandOrigin;
#ifndef QT_NO_RUBBERBAND
    QRubberBand *m_rubberBand;
#endif
    QChartView::RubberBands m_rubberBandFlags;
};

QT_CHARTS_END_NAMESPACE
#endif
