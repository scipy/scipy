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

#ifndef QBOXPLOTSERIES_P_H
#define QBOXPLOTSERIES_P_H

#include <QtCharts/QBoxPlotSeries>
#include <private/qabstractbarseries_p.h>
#include <private/abstractdomain_p.h>
#include <QtCharts/QBarSet>
#include <QtCharts/private/qchartglobal_p.h>

QT_CHARTS_BEGIN_NAMESPACE

class BoxPlotAnimation;

class Q_CHARTS_PRIVATE_EXPORT QBoxPlotSeriesPrivate : public QAbstractSeriesPrivate
{
    Q_OBJECT

public:
    QBoxPlotSeriesPrivate(QBoxPlotSeries *q);
    ~QBoxPlotSeriesPrivate();

    void initializeGraphics(QGraphicsItem *parent);
    void initializeDomain();
    void initializeAxes();
    void initializeAnimations(QChart::AnimationOptions options, int duration, QEasingCurve &curve);
    void initializeTheme(int index, ChartTheme *theme, bool forced = false);

    QList<QLegendMarker*> createLegendMarkers(QLegend *legend);

    virtual QAbstractAxis::AxisType defaultAxisType(Qt::Orientation orientation) const;
    QAbstractAxis *createDefaultAxis(Qt::Orientation orientation) const;

    bool append(QBoxSet *set);
    bool remove(QBoxSet *set);
    bool append(QList<QBoxSet *> sets);
    bool remove(QList<QBoxSet *> sets);
    bool insert(int index, QBoxSet *set);
    QBoxSet *boxSetAt(int index);

    qreal max();
    qreal min();

private:
    void populateCategories(QBarCategoryAxis *axis);

Q_SIGNALS:
    void updated();
    void clicked(int index, QBoxSet *barset);
    void pressed(int index, QBoxSet *barset);
    void released(int index, QBoxSet *barset);
    void doubleClicked(int index, QBoxSet *barset);
    void updatedBoxes();
    void updatedLayout();
    void restructuredBoxes();

private Q_SLOTS:
    void handleSeriesChange(QAbstractSeries *series);
    void handleSeriesRemove(QAbstractSeries *series);

protected:
    QList<QBoxSet *> m_boxSets;
    QPen m_pen;
    QBrush m_brush;
    bool m_boxOutlineVisible;
    int m_index;
    BoxPlotAnimation *m_animation;
    qreal m_boxWidth;

private:
    Q_DECLARE_PUBLIC(QBoxPlotSeries)
};

QT_CHARTS_END_NAMESPACE

#endif
