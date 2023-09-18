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

#ifndef QSCATTERSERIES_H
#define QSCATTERSERIES_H

#include <QtCharts/QChartGlobal>
#include <QtCharts/qxyseries.h>

QT_CHARTS_BEGIN_NAMESPACE

class QScatterSeriesPrivate;

class Q_CHARTS_EXPORT QScatterSeries : public QXYSeries
{
    Q_OBJECT
    Q_PROPERTY(QColor color READ color WRITE setColor NOTIFY colorChanged)
    Q_PROPERTY(QColor borderColor READ borderColor WRITE setBorderColor NOTIFY borderColorChanged)
    Q_PROPERTY(MarkerShape markerShape READ markerShape WRITE setMarkerShape NOTIFY markerShapeChanged)
    Q_PROPERTY(qreal markerSize READ markerSize WRITE setMarkerSize NOTIFY markerSizeChanged)
    Q_PROPERTY(QBrush brush READ brush WRITE setBrush)
    Q_ENUMS(MarkerShape)

public:
    enum MarkerShape {
        MarkerShapeCircle,
        MarkerShapeRectangle
    };

public:
    explicit QScatterSeries(QObject *parent = nullptr);
    ~QScatterSeries();
    QAbstractSeries::SeriesType type() const;
    void setPen(const QPen &pen);
    void setBrush(const QBrush &brush);
    QBrush brush() const;
    void setColor(const QColor &color);
    QColor color() const;
    void setBorderColor(const QColor &color);
    QColor borderColor() const;
    MarkerShape markerShape() const;
    void setMarkerShape(MarkerShape shape);
    qreal markerSize() const;
    void setMarkerSize(qreal size);

Q_SIGNALS:
    void colorChanged(QColor color);
    void borderColorChanged(QColor color);
    void markerShapeChanged(MarkerShape shape);
    void markerSizeChanged(qreal size);

private:
    Q_DECLARE_PRIVATE(QScatterSeries)
    Q_DISABLE_COPY(QScatterSeries)
    friend class ScatterChartItem;
};

QT_CHARTS_END_NAMESPACE

#endif // QSCATTERSERIES_H
