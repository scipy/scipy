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

#ifndef QCANDLESTICKSERIES_P_H
#define QCANDLESTICKSERIES_P_H

#include <private/qabstractseries_p.h>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class CandlestickAnimation;
class QBarCategoryAxis;
class QCandlestickSeries;
class QCandlestickSet;
class QDateTimeAxis;

class Q_CHARTS_PRIVATE_EXPORT QCandlestickSeriesPrivate : public QAbstractSeriesPrivate
{
    Q_OBJECT

public:
    QCandlestickSeriesPrivate(QCandlestickSeries *q);
    ~QCandlestickSeriesPrivate();

    void initializeDomain();
    void initializeAxes();
    void initializeTheme(int index, ChartTheme* theme, bool forced = false);
    void initializeGraphics(QGraphicsItem* parent);
    void initializeAnimations(QChart::AnimationOptions options, int duration, QEasingCurve &curve);

    QList<QLegendMarker *> createLegendMarkers(QLegend *legend);

    virtual QAbstractAxis::AxisType defaultAxisType(Qt::Orientation orientation) const;
    QAbstractAxis *createDefaultAxis(Qt::Orientation orientation) const;

    bool append(const QList<QCandlestickSet *> &sets);
    bool remove(const QList<QCandlestickSet *> &sets);
    bool insert(int index, QCandlestickSet *set);

Q_SIGNALS:
    void clicked(int index, QCandlestickSet *set);
    void pressed(int index, QCandlestickSet *set);
    void released(int index, QCandlestickSet *set);
    void doubleClicked(int index, QCandlestickSet *set);
    void updated();
    void updatedLayout();
    void updatedCandlesticks();

private Q_SLOTS:
    void handleSeriesChange(QAbstractSeries *series);
    void handleSeriesRemove(QAbstractSeries *series);

private:
    void populateBarCategories(QBarCategoryAxis *axis);

protected:
    QList<QCandlestickSet *> m_sets;
    qreal m_maximumColumnWidth;
    qreal m_minimumColumnWidth;
    qreal m_bodyWidth;
    bool m_bodyOutlineVisible;
    qreal m_capsWidth;
    bool m_capsVisible;
    QColor m_increasingColor;
    QColor m_decreasingColor;
    bool m_customIncreasingColor;
    bool m_customDecreasingColor;
    QBrush m_brush;
    QPen m_pen;
    CandlestickAnimation *m_animation;

private:
    Q_DECLARE_PUBLIC(QCandlestickSeries)
};

QT_CHARTS_END_NAMESPACE

#endif // QCANDLESTICKSERIES_P_H
