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

#ifndef QHORIZONTALBARSERIES_H
#define QHORIZONTALBARSERIES_H

#include <QtCharts/QAbstractBarSeries>

QT_CHARTS_BEGIN_NAMESPACE

class QHorizontalBarSeriesPrivate;

class Q_CHARTS_EXPORT QHorizontalBarSeries : public QAbstractBarSeries
{
    Q_OBJECT
public:
    explicit QHorizontalBarSeries(QObject *parent = nullptr);
    ~QHorizontalBarSeries();
    QAbstractSeries::SeriesType type() const;

private:
    Q_DECLARE_PRIVATE(QHorizontalBarSeries)
    Q_DISABLE_COPY(QHorizontalBarSeries)
};

QT_CHARTS_END_NAMESPACE

#endif // QHORIZONTALBARSERIES_H
