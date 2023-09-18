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

#ifndef QSPLINESERIES_H
#define QSPLINESERIES_H

#include <QtCharts/QChartGlobal>
#include <QtCharts/qlineseries.h>
#include <QtCore/QList>
#include <QtCore/QPointF>
#include <QtCore/QtGlobal>

QT_CHARTS_BEGIN_NAMESPACE

class QSplineSeriesPrivate;

class Q_CHARTS_EXPORT QSplineSeries : public QLineSeries
{
    Q_OBJECT

public:
    explicit QSplineSeries(QObject *parent = nullptr);
    ~QSplineSeries();
    QAbstractSeries::SeriesType type() const;

private:
    Q_DECLARE_PRIVATE(QSplineSeries)
    Q_DISABLE_COPY(QSplineSeries)
    friend class SplineChartItem;
};

QT_CHARTS_END_NAMESPACE

#endif // QSPLINESERIES_H
