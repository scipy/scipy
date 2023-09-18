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

#ifndef QXYLEGENDMARKER_H
#define QXYLEGENDMARKER_H

#include <QtCharts/QChartGlobal>
#include <QtCharts/QLegendMarker>
#include <QtCharts/QXYSeries>

QT_CHARTS_BEGIN_NAMESPACE

class QXYLegendMarkerPrivate;

class Q_CHARTS_EXPORT QXYLegendMarker : public QLegendMarker
{
    Q_OBJECT
public:
    explicit QXYLegendMarker(QXYSeries *series, QLegend *legend, QObject *parent = nullptr);
    virtual ~QXYLegendMarker();

    virtual LegendMarkerType type() { return LegendMarkerTypeXY; }

    // Related series
    virtual QXYSeries* series();

protected:
    QXYLegendMarker(QXYLegendMarkerPrivate &d, QObject *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QXYLegendMarker)
    Q_DISABLE_COPY(QXYLegendMarker)

};

QT_CHARTS_END_NAMESPACE

#endif // QXYLEGENDMARKER_H
