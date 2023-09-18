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

#ifndef QLEGEND_H
#define QLEGEND_H

#include <QtCharts/QChartGlobal>
#include <QtWidgets/QGraphicsWidget>
#include <QtGui/QPen>
#include <QtGui/QBrush>

QT_CHARTS_BEGIN_NAMESPACE

class QChart;
class QLegendPrivate;
class QLegendMarker;
class QAbstractSeries;

class Q_CHARTS_EXPORT QLegend : public QGraphicsWidget
{
    Q_OBJECT
    Q_PROPERTY(Qt::Alignment alignment READ alignment WRITE setAlignment)
    Q_PROPERTY(bool backgroundVisible READ isBackgroundVisible WRITE setBackgroundVisible NOTIFY backgroundVisibleChanged)
    Q_PROPERTY(QColor color READ color WRITE setColor NOTIFY colorChanged)
    Q_PROPERTY(QColor borderColor READ borderColor WRITE setBorderColor NOTIFY borderColorChanged)
    Q_PROPERTY(QFont font READ font WRITE setFont NOTIFY fontChanged)
    Q_PROPERTY(QColor labelColor READ labelColor WRITE setLabelColor NOTIFY labelColorChanged)
    Q_PROPERTY(bool reverseMarkers READ reverseMarkers WRITE setReverseMarkers NOTIFY reverseMarkersChanged)
    Q_PROPERTY(bool showToolTips READ showToolTips WRITE setShowToolTips NOTIFY showToolTipsChanged)
    Q_PROPERTY(MarkerShape markerShape READ markerShape WRITE setMarkerShape NOTIFY markerShapeChanged)

private:
    explicit QLegend(QChart *chart);

public:
    enum MarkerShape {
        MarkerShapeDefault,
        MarkerShapeRectangle,
        MarkerShapeCircle,
        MarkerShapeFromSeries
    };
    Q_ENUMS(MarkerShape)

    ~QLegend();

    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget = nullptr);

    void setBrush(const QBrush &brush);
    QBrush brush() const;
    void setColor(QColor color);
    QColor color();

    void setPen(const QPen &pen);
    QPen pen() const;
    void setBorderColor(QColor color);
    QColor borderColor();

    void setFont(const QFont &font);
    QFont font() const;
    void setLabelBrush(const QBrush &brush);
    QBrush labelBrush() const;

    void setLabelColor(QColor color);
    QColor labelColor() const;

    void setAlignment(Qt::Alignment alignment);
    Qt::Alignment alignment() const;

    void detachFromChart();
    void attachToChart();
    bool isAttachedToChart();

    void setBackgroundVisible(bool visible = true);
    bool isBackgroundVisible() const;

    QList <QLegendMarker*> markers(QAbstractSeries *series = nullptr) const;

    bool reverseMarkers();
    void setReverseMarkers(bool reverseMarkers = true);

    bool showToolTips() const;
    void setShowToolTips(bool show);

    MarkerShape markerShape() const;
    void setMarkerShape(MarkerShape shape);

protected:
    void hideEvent(QHideEvent *event);
    void showEvent(QShowEvent *event);

Q_SIGNALS:
    void backgroundVisibleChanged(bool visible);
    void colorChanged(QColor color);
    void borderColorChanged(QColor color);
    void fontChanged(QFont font);
    void labelColorChanged(QColor color);
    void reverseMarkersChanged(bool reverseMarkers);
    void showToolTipsChanged(bool showToolTips);
    void markerShapeChanged(MarkerShape shape);

private:
    QScopedPointer<QLegendPrivate> d_ptr;
    Q_DISABLE_COPY(QLegend)
    friend class LegendScroller;
    friend class LegendLayout;
    friend class ChartLayout;
    friend class LegendMarkerItem;
    friend class QLegendMarkerPrivate;
};

QT_CHARTS_END_NAMESPACE

#endif // QLEGEND_H
