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

#ifndef QPIESERIES_P_H
#define QPIESERIES_P_H

#include <QtCharts/QPieSeries>
#include <private/qabstractseries_p.h>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE
class QLegendPrivate;

class Q_CHARTS_PRIVATE_EXPORT QPieSeriesPrivate : public QAbstractSeriesPrivate
{
    Q_OBJECT

public:
    QPieSeriesPrivate(QPieSeries *parent);
    ~QPieSeriesPrivate();

    void initializeDomain();
    void initializeAxes();
    void initializeGraphics(QGraphicsItem* parent);
    void initializeAnimations(QtCharts::QChart::AnimationOptions options, int duration,
                              QEasingCurve &curve);
    void initializeTheme(int index, ChartTheme* theme, bool forced = false);

    QList<QLegendMarker *> createLegendMarkers(QLegend *legend);

    QAbstractAxis::AxisType defaultAxisType(Qt::Orientation orientation) const;
    QAbstractAxis* createDefaultAxis(Qt::Orientation orientation) const;

    void updateDerivativeData();
    void setSizes(qreal innerSize, qreal outerSize);

    static QPieSeriesPrivate *fromSeries(QPieSeries *series);

Q_SIGNALS:
    void calculatedDataChanged();
    void pieSizeChanged();
    void pieStartAngleChanged();
    void pieEndAngleChanged();
    void horizontalPositionChanged();
    void verticalPositionChanged();

public Q_SLOTS:
    void sliceValueChanged();
    void sliceClicked();
    void sliceHovered(bool state);
    void slicePressed();
    void sliceReleased();
    void sliceDoubleClicked();

private:
    QList<QPieSlice *> m_slices;
    qreal m_pieRelativeHorPos;
    qreal m_pieRelativeVerPos;
    qreal m_pieRelativeSize;
    qreal m_pieStartAngle;
    qreal m_pieEndAngle;
    qreal m_sum;
    qreal m_holeRelativeSize;

public:
    friend class QLegendPrivate;
    Q_DECLARE_PUBLIC(QPieSeries)
};

QT_CHARTS_END_NAMESPACE

#endif // QPIESERIES_P_H
