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

#ifndef QCHART_H
#define QCHART_H

#include <QtCharts/QAbstractSeries>
#include <QtCharts/QLegend>
#include <QtWidgets/QGraphicsWidget>
#include <QtCore/QMargins>

QT_BEGIN_NAMESPACE
class QGraphicsSceneResizeEvent;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE

class QAbstractSeries;
class QAbstractAxis;
class QLegend;
class QChartPrivate;
class QBoxPlotSeries;

class Q_CHARTS_EXPORT QChart : public QGraphicsWidget
{
    Q_OBJECT
    Q_PROPERTY(QChart::ChartTheme theme READ theme WRITE setTheme)
    Q_PROPERTY(QString title READ title WRITE setTitle)
    Q_PROPERTY(bool backgroundVisible READ isBackgroundVisible WRITE setBackgroundVisible)
    Q_PROPERTY(bool dropShadowEnabled READ isDropShadowEnabled WRITE setDropShadowEnabled)
    Q_PROPERTY(qreal backgroundRoundness READ backgroundRoundness WRITE setBackgroundRoundness)
    Q_PROPERTY(QChart::AnimationOptions animationOptions READ animationOptions WRITE setAnimationOptions)
    Q_PROPERTY(int animationDuration READ animationDuration WRITE setAnimationDuration)
    Q_PROPERTY(QEasingCurve animationEasingCurve READ animationEasingCurve WRITE setAnimationEasingCurve)
    Q_PROPERTY(QMargins margins READ margins WRITE setMargins)
    Q_PROPERTY(QChart::ChartType chartType READ chartType)
    Q_PROPERTY(bool plotAreaBackgroundVisible READ isPlotAreaBackgroundVisible WRITE setPlotAreaBackgroundVisible)
    Q_PROPERTY(bool localizeNumbers READ localizeNumbers WRITE setLocalizeNumbers)
    Q_PROPERTY(QLocale locale READ locale WRITE setLocale)
    Q_PROPERTY(QRectF plotArea READ plotArea WRITE setPlotArea NOTIFY plotAreaChanged)
    Q_ENUMS(ChartTheme)
    Q_ENUMS(AnimationOption)
    Q_ENUMS(ChartType)

public:
    enum ChartType {
        ChartTypeUndefined = 0,
        ChartTypeCartesian,
        ChartTypePolar
    };

    enum ChartTheme {
        ChartThemeLight = 0,
        ChartThemeBlueCerulean,
        ChartThemeDark,
        ChartThemeBrownSand,
        ChartThemeBlueNcs,
        ChartThemeHighContrast,
        ChartThemeBlueIcy,
        ChartThemeQt
    };

    enum AnimationOption {
        NoAnimation = 0x0,
        GridAxisAnimations = 0x1,
        SeriesAnimations = 0x2,
        AllAnimations = 0x3
    };

    Q_DECLARE_FLAGS(AnimationOptions, AnimationOption)

public:
    explicit QChart(QGraphicsItem *parent = nullptr, Qt::WindowFlags wFlags = Qt::WindowFlags());
    ~QChart();

    void addSeries(QAbstractSeries *series);
    void removeSeries(QAbstractSeries *series);
    void removeAllSeries();
    QList<QAbstractSeries *> series() const;

    Q_DECL_DEPRECATED void setAxisX(QAbstractAxis *axis, QAbstractSeries *series = nullptr);
    Q_DECL_DEPRECATED void setAxisY(QAbstractAxis *axis, QAbstractSeries *series = nullptr);
    Q_DECL_DEPRECATED QAbstractAxis *axisX(QAbstractSeries *series = nullptr) const;
    Q_DECL_DEPRECATED QAbstractAxis *axisY(QAbstractSeries *series = nullptr) const;

    void addAxis(QAbstractAxis *axis, Qt::Alignment alignment);
    void removeAxis(QAbstractAxis *axis);
    QList<QAbstractAxis*> axes(Qt::Orientations orientation = Qt::Horizontal|Qt::Vertical, QAbstractSeries *series = nullptr) const;

    void createDefaultAxes();

    void setTheme(QChart::ChartTheme theme);
    QChart::ChartTheme theme() const;

    void setTitle(const QString &title);
    QString title() const;
    void setTitleFont(const QFont &font);
    QFont titleFont() const;
    void setTitleBrush(const QBrush &brush);
    QBrush titleBrush() const;

    void setBackgroundBrush(const QBrush &brush);
    QBrush backgroundBrush() const;
    void setBackgroundPen(const QPen &pen);
    QPen backgroundPen() const;
    void setBackgroundVisible(bool visible = true);
    bool isBackgroundVisible() const;

    void setDropShadowEnabled(bool enabled = true);
    bool isDropShadowEnabled() const;
    void setBackgroundRoundness(qreal diameter);
    qreal backgroundRoundness() const;

    void setAnimationOptions(AnimationOptions options);
    AnimationOptions animationOptions() const;
    void setAnimationDuration(int msecs);
    int animationDuration() const;
    void setAnimationEasingCurve(const QEasingCurve &curve);
    QEasingCurve animationEasingCurve() const;

    void zoomIn();
    void zoomOut();

    void zoomIn(const QRectF &rect);
    void zoom(qreal factor);
    void zoomReset();
    bool isZoomed();

    void scroll(qreal dx, qreal dy);

    QLegend *legend() const;

    void setMargins(const QMargins &margins);
    QMargins margins() const;

    QRectF plotArea() const;
    void setPlotArea(const QRectF &rect);
    void setPlotAreaBackgroundBrush(const QBrush &brush);
    QBrush plotAreaBackgroundBrush() const;
    void setPlotAreaBackgroundPen(const QPen &pen);
    QPen plotAreaBackgroundPen() const;
    void setPlotAreaBackgroundVisible(bool visible = true);
    bool isPlotAreaBackgroundVisible() const;
    void setLocalizeNumbers(bool localize);
    bool localizeNumbers() const;
    void setLocale(const QLocale &locale);
    QLocale locale() const;

    QPointF mapToValue(const QPointF &position, QAbstractSeries *series = nullptr);
    QPointF mapToPosition(const QPointF &value, QAbstractSeries *series = nullptr);

    ChartType chartType() const;

Q_SIGNALS:
    void plotAreaChanged(const QRectF &plotArea);

protected:
    explicit QChart(QChart::ChartType type, QGraphicsItem *parent, Qt::WindowFlags wFlags);
    QScopedPointer<QChartPrivate> d_ptr;
    friend class QLegend;
    friend class DeclarativeChart;
    friend class ChartDataSet;
    friend class ChartPresenter;
    friend class ChartThemeManager;
    friend class QAbstractSeries;
    friend class QBoxPlotSeriesPrivate;
    friend class QCandlestickSeriesPrivate;
    friend class AbstractBarChartItem;

private:
    Q_DISABLE_COPY(QChart)
};

QT_CHARTS_END_NAMESPACE

#ifndef Q_CLANG_QDOC
Q_DECLARE_OPERATORS_FOR_FLAGS(QT_CHARTS_NAMESPACE::QChart::AnimationOptions)
#endif

#endif // QCHART_H
