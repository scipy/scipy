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

#ifndef QBOXPLOTSERIES_H
#define QBOXPLOTSERIES_H

#include <QtCharts/QChartGlobal>
#include <QtCharts/QBoxSet>
#include <QtCharts/QAbstractSeries>

QT_CHARTS_BEGIN_NAMESPACE

class QBoxPlotSeriesPrivate;

class Q_CHARTS_EXPORT QBoxPlotSeries : public QAbstractSeries
{
    Q_OBJECT
    Q_PROPERTY(bool boxOutlineVisible READ boxOutlineVisible WRITE setBoxOutlineVisible NOTIFY boxOutlineVisibilityChanged)
    Q_PROPERTY(qreal boxWidth READ boxWidth WRITE setBoxWidth NOTIFY boxWidthChanged)
    Q_PROPERTY(QPen pen READ pen WRITE setPen NOTIFY penChanged)
    Q_PROPERTY(QBrush brush READ brush WRITE setBrush NOTIFY brushChanged)
    Q_PROPERTY(int count READ count NOTIFY countChanged REVISION 1)
public:
    explicit QBoxPlotSeries(QObject *parent = nullptr);
    ~QBoxPlotSeries();

    bool append(QBoxSet *box);
    bool remove(QBoxSet *box);
    bool take(QBoxSet *box);
    bool append(QList<QBoxSet *> boxes);
    bool insert(int index, QBoxSet *box);
    int count() const;
    QList<QBoxSet *> boxSets() const;
    void clear();

    QAbstractSeries::SeriesType type() const;

    void setBoxOutlineVisible(bool visible);
    bool boxOutlineVisible();
    void setBoxWidth(qreal width);
    qreal boxWidth();
    void setBrush(const QBrush &brush);
    QBrush brush() const;
    void setPen(const QPen &pen);
    QPen pen() const;

Q_SIGNALS:
    void clicked(QBoxSet *boxset);
    void hovered(bool status, QBoxSet *boxset);
    void pressed(QBoxSet *boxset);
    void released(QBoxSet *boxset);
    void doubleClicked(QBoxSet *boxset);
    void countChanged();
    void penChanged();
    void brushChanged();
    void boxOutlineVisibilityChanged();
    void boxWidthChanged();

    void boxsetsAdded(QList<QBoxSet *> sets);
    void boxsetsRemoved(QList<QBoxSet *> sets);

private:
    Q_DECLARE_PRIVATE(QBoxPlotSeries)
    Q_DISABLE_COPY(QBoxPlotSeries)
    friend class BoxPlotChartItem;
    friend class QBoxPlotLegendMarkerPrivate;
};

QT_CHARTS_END_NAMESPACE

#endif // QBOXPLOTSERIES_H
