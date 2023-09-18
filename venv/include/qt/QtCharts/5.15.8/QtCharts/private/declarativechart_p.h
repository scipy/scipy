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

#ifndef DECLARATIVECHART_H
#define DECLARATIVECHART_H

#include <private/glxyseriesdata_p.h>
#include <private/declarativechartglobal_p.h>
#include <private/declarativeabstractrendernode_p.h>

#include <QtCore/QtGlobal>
#include <QtQuick/QQuickItem>
#include <QtWidgets/QGraphicsScene>

#include <QtCharts/QChart>
#include <QtCore/QLocale>
#include <QQmlComponent>

QT_CHARTS_BEGIN_NAMESPACE

class DeclarativeMargins;
class Domain;
class DeclarativeAxes;

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeChart : public QQuickItem
{
    Q_OBJECT
    Q_PROPERTY(Theme theme READ theme WRITE setTheme)
    Q_PROPERTY(Animation animationOptions READ animationOptions WRITE setAnimationOptions)
    Q_PROPERTY(int animationDuration READ animationDuration WRITE setAnimationDuration NOTIFY animationDurationChanged REVISION 5)
    Q_PROPERTY(QEasingCurve animationEasingCurve READ animationEasingCurve WRITE setAnimationEasingCurve NOTIFY animationEasingCurveChanged REVISION 5)
    Q_PROPERTY(QString title READ title WRITE setTitle)
    Q_PROPERTY(QFont titleFont READ titleFont WRITE setTitleFont)
    Q_PROPERTY(QColor titleColor READ titleColor WRITE setTitleColor NOTIFY titleColorChanged)
    Q_PROPERTY(QtCharts::QLegend *legend READ legend CONSTANT)
    Q_PROPERTY(int count READ count)
    Q_PROPERTY(QColor backgroundColor READ backgroundColor WRITE setBackgroundColor NOTIFY backgroundColorChanged)
    Q_PROPERTY(bool dropShadowEnabled READ dropShadowEnabled WRITE setDropShadowEnabled NOTIFY dropShadowEnabledChanged)
    Q_PROPERTY(qreal backgroundRoundness READ backgroundRoundness WRITE setBackgroundRoundness NOTIFY backgroundRoundnessChanged REVISION 3)
    Q_PROPERTY(QtCharts::DeclarativeMargins *margins READ margins NOTIFY marginsChanged REVISION 2)
    Q_PROPERTY(QRectF plotArea READ plotArea WRITE setPlotArea NOTIFY plotAreaChanged REVISION 1)
    Q_PROPERTY(QColor plotAreaColor READ plotAreaColor WRITE setPlotAreaColor NOTIFY plotAreaColorChanged REVISION 3)
    Q_PROPERTY(QQmlListProperty<QAbstractAxis> axes READ axes REVISION 2)
    Q_PROPERTY(bool localizeNumbers READ localizeNumbers WRITE setLocalizeNumbers NOTIFY localizeNumbersChanged REVISION 4)
    Q_PROPERTY(QLocale locale READ locale WRITE setLocale NOTIFY localeChanged REVISION 4)
    Q_ENUMS(Animation)
    Q_ENUMS(Theme)
    Q_ENUMS(SeriesType)

public:
    // duplicating enums from QChart to make the QML api namings 1-to-1 with the C++ api
    enum Theme {
        ChartThemeLight = 0,
        ChartThemeBlueCerulean,
        ChartThemeDark,
        ChartThemeBrownSand,
        ChartThemeBlueNcs,
        ChartThemeHighContrast,
        ChartThemeBlueIcy,
        ChartThemeQt
    };

    enum Animation {
        NoAnimation = 0x0,
        GridAxisAnimations = 0x1,
        SeriesAnimations = 0x2,
        AllAnimations = 0x3
    };

    enum SeriesType {
        SeriesTypeLine,
        SeriesTypeArea,
        SeriesTypeBar,
        SeriesTypeStackedBar,
        SeriesTypePercentBar,
        SeriesTypePie,
        SeriesTypeScatter,
        SeriesTypeSpline,
        SeriesTypeHorizontalBar,
        SeriesTypeHorizontalStackedBar,
        SeriesTypeHorizontalPercentBar,
        SeriesTypeBoxPlot,
        SeriesTypeCandlestick
    };

public:
    DeclarativeChart(QQuickItem *parent = 0);
    ~DeclarativeChart();

public: // From parent classes
    void childEvent(QChildEvent *event);
    void componentComplete();
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry);
    QSGNode *updatePaintNode(QSGNode *oldNode, UpdatePaintNodeData *);
protected:
    void mousePressEvent(QMouseEvent *event);
    void mouseReleaseEvent(QMouseEvent *event);
    void hoverMoveEvent(QHoverEvent *event);
    void mouseDoubleClickEvent(QMouseEvent *event);
private Q_SLOTS:
    void handleAntialiasingChanged(bool enable);
    void sceneChanged(QList<QRectF> region);
    void renderScene();

public:
    void setTheme(DeclarativeChart::Theme theme);
    DeclarativeChart::Theme theme();
    void setAnimationOptions(DeclarativeChart::Animation animations);
    DeclarativeChart::Animation animationOptions();
    void setAnimationDuration(int msecs);
    int animationDuration() const;
    void setAnimationEasingCurve(const QEasingCurve &curve);
    QEasingCurve animationEasingCurve() const;
    void setTitle(QString title);
    QString title();
    QLegend *legend();
    QFont titleFont() const;
    void setTitleFont(const QFont &font);
    void setTitleColor(QColor color);
    QColor titleColor();
    void setBackgroundColor(QColor color);
    QColor backgroundColor();
    void setPlotAreaColor(QColor color);
    QColor plotAreaColor();
    void setLocalizeNumbers(bool localize);
    bool localizeNumbers() const;
    void setLocale(const QLocale &locale);
    QLocale locale() const;

    int count();
    void setDropShadowEnabled(bool enabled);
    bool dropShadowEnabled();
    qreal backgroundRoundness() const;
    void setBackgroundRoundness(qreal diameter);

    // Margins & plotArea
    DeclarativeMargins *margins() { return m_margins; }
    QRectF plotArea() { return m_chart->plotArea(); }
    void setPlotArea(const QRectF &rect);

    // Axis handling
    QAbstractAxis *defaultAxis(Qt::Orientation orientation, QAbstractSeries *series);
    void initializeAxes(QAbstractSeries *series);
    void doInitializeAxes(QAbstractSeries *series, DeclarativeAxes *axes);
    QQmlListProperty<QAbstractAxis> axes();
    static void axesAppendFunc(QQmlListProperty<QAbstractAxis> *list, QAbstractAxis *element);
    static int axesCountFunc(QQmlListProperty<QAbstractAxis> *list);
    static QAbstractAxis *axesAtFunc(QQmlListProperty<QAbstractAxis> *list, int index);
    static void axesClearFunc(QQmlListProperty<QAbstractAxis> *list);

public:
    Q_INVOKABLE QAbstractSeries *series(int index);
    Q_INVOKABLE QAbstractSeries *series(QString seriesName);
    Q_INVOKABLE QAbstractSeries *createSeries(int type, QString name = "", QAbstractAxis *axisX = 0, QAbstractAxis *axisY = 0);
    Q_INVOKABLE void removeSeries(QAbstractSeries *series);
    Q_INVOKABLE void removeAllSeries() { m_chart->removeAllSeries(); }
    Q_INVOKABLE void setAxisX(QAbstractAxis *axis, QAbstractSeries *series = 0);
    Q_INVOKABLE void setAxisY(QAbstractAxis *axis, QAbstractSeries *series = 0);
    Q_INVOKABLE QAbstractAxis *axisX(QAbstractSeries *series = 0);
    Q_INVOKABLE QAbstractAxis *axisY(QAbstractSeries *series = 0);
    Q_INVOKABLE void zoom(qreal factor);
    Q_REVISION(5) Q_INVOKABLE void zoomIn();
    Q_REVISION(5) Q_INVOKABLE void zoomIn(const QRectF &rectangle);
    Q_REVISION(5) Q_INVOKABLE void zoomOut();
    Q_REVISION(5) Q_INVOKABLE void zoomReset();
    Q_REVISION(5) Q_INVOKABLE bool isZoomed();
    Q_INVOKABLE void scrollLeft(qreal pixels);
    Q_INVOKABLE void scrollRight(qreal pixels);
    Q_INVOKABLE void scrollUp(qreal pixels);
    Q_INVOKABLE void scrollDown(qreal pixels);
    Q_REVISION(5) Q_INVOKABLE QPointF mapToValue(const QPointF &position,
                                                 QAbstractSeries *series = 0);
    Q_REVISION(5) Q_INVOKABLE QPointF mapToPosition(const QPointF &value,
                                                    QAbstractSeries *series = 0);


Q_SIGNALS:
    void axisLabelsChanged();
    void titleColorChanged(QColor color);
    void backgroundColorChanged();
    void dropShadowEnabledChanged(bool enabled);
    Q_REVISION(2) void marginsChanged();
    void plotAreaChanged(QRectF plotArea);
    void seriesAdded(QAbstractSeries *series);
    void seriesRemoved(QAbstractSeries *series);
    Q_REVISION(3) void plotAreaColorChanged();
    Q_REVISION(3) void backgroundRoundnessChanged(qreal diameter);
    Q_REVISION(4) void localizeNumbersChanged();
    Q_REVISION(4) void localeChanged();
    Q_REVISION(5) void animationDurationChanged(int msecs);
    Q_REVISION(5) void animationEasingCurveChanged(QEasingCurve curve);
    void needRender();
    void pendingRenderNodeMouseEventResponses();

private Q_SLOTS:
    void changeMargins(int top, int bottom, int left, int right);
    void handleAxisXSet(QAbstractAxis *axis);
    void handleAxisYSet(QAbstractAxis *axis);
    void handleAxisXTopSet(QAbstractAxis *axis);
    void handleAxisYRightSet(QAbstractAxis *axis);
    void handleSeriesAdded(QAbstractSeries *series);
    void handlePendingRenderNodeMouseEventResponses();

protected:
    explicit DeclarativeChart(QChart::ChartType type, QQuickItem *parent);

private:
    void initChart(QChart::ChartType type);
    void seriesAxisAttachHelper(QAbstractSeries *series, QAbstractAxis *axis,
                                Qt::Orientations orientation, Qt::Alignment alignment);
    void findMinMaxForSeries(QAbstractSeries *series,Qt::Orientations orientation,
                             qreal &min, qreal &max);
    void queueRendererMouseEvent(QMouseEvent *event);

    // Extending QChart with DeclarativeChart is not possible because QObject does not support
    // multi inheritance, so we now have a QChart as a member instead
    QChart *m_chart;
    QGraphicsScene *m_scene;
    QPointF m_mousePressScenePoint;
    QPoint m_mousePressScreenPoint;
    QPointF m_lastMouseMoveScenePoint;
    QPoint m_lastMouseMoveScreenPoint;
    Qt::MouseButton m_mousePressButton;
    Qt::MouseButtons m_mousePressButtons;
    QImage *m_sceneImage;
    bool m_sceneImageDirty;
    bool m_updatePending;
    Qt::HANDLE m_paintThreadId;
    Qt::HANDLE m_guiThreadId;
    DeclarativeMargins *m_margins;
    GLXYSeriesDataManager *m_glXYDataManager;
    bool m_sceneImageNeedsClear;
    QVector<QMouseEvent *> m_pendingRenderNodeMouseEvents;
    QVector<MouseEventResponse> m_pendingRenderNodeMouseEventResponses;
    QRectF m_adjustedPlotArea;
};

QT_CHARTS_END_NAMESPACE

#endif // DECLARATIVECHART_H
