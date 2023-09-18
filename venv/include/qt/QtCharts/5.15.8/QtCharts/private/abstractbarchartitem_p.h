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


#ifndef ABSTRACTBARCHARTITEM_H
#define ABSTRACTBARCHARTITEM_H

#include <private/chartitem_p.h>
#include <QtCharts/QAbstractBarSeries>
#include <QtCharts/private/qchartglobal_p.h>
#include <QtGui/QPen>
#include <QtGui/QBrush>

QT_CHARTS_BEGIN_NAMESPACE

class Bar;
class QAxisCategories;
class QChart;
class BarAnimation;
class QBarSetPrivate;
class QAbstractAxis;

class Q_CHARTS_PRIVATE_EXPORT AbstractBarChartItem : public ChartItem
{
    Q_OBJECT
public:
    AbstractBarChartItem(QAbstractBarSeries *series, QGraphicsItem* item = 0);
    virtual ~AbstractBarChartItem();

public:
    void paint(QPainter *painter, const QStyleOptionGraphicsItem *option, QWidget *widget);
    QRectF boundingRect() const;

    virtual QVector<QRectF> calculateLayout() = 0;
    void initializeFullLayout();
    virtual void initializeLayout(int set, int category, int layoutIndex, bool resetAnimation) = 0;
    virtual void applyLayout(const QVector<QRectF> &layout);
    virtual void setAnimation(BarAnimation *animation);
    virtual ChartAnimation *animation() const;
    void setLayout(const QVector<QRectF> &layout);
    QRectF geometry() const { return m_rect;}
    void resetAnimation();

public Q_SLOTS:
    void handleDomainUpdated();
    void handleLayoutChanged();
    void handleLabelsVisibleChanged(bool visible);
    void handleDataStructureChanged();     // structure of of series has changed, recreate graphic items
    void handleVisibleChanged();
    void handleOpacityChanged();
    void handleUpdatedBars();
    void handleLabelsPositionChanged();
    virtual void positionLabels();
    void handleBarValueChange(int index, QBarSet *barset);
    void handleBarValueAdd(int index, int count, QBarSet *barset);
    void handleBarValueRemove(int index, int count, QBarSet *barset);
    void handleSeriesAdded(QAbstractSeries *series);
    void handleSeriesRemoved(QAbstractSeries *series);

protected:
    void positionLabelsVertical();
    void createLabelItems();
    void handleSetStructureChange();
    virtual QString generateLabelText(int set, int category, qreal value);
    void updateBarItems();
    virtual void markLabelsDirty(QBarSet *barset, int index, int count);
    void calculateSeriesPositionAdjustmentAndWidth();

    QRectF m_rect;
    QVector<QRectF> m_layout;

    BarAnimation *m_animation;

    QAbstractBarSeries *m_series; // Not owned.
    QMap<QBarSet *, QList<Bar *> > m_barMap;
    QMap<QBarSet *, QHash<int, Bar *> > m_indexForBarMap;
    int m_firstCategory;
    int m_lastCategory;
    int m_categoryCount;
    QSizeF m_oldSize;
    bool m_labelItemsMissing;
    Qt::Orientation m_orientation;
    bool m_resetAnimation;
    qreal m_seriesPosAdjustment;
    qreal m_seriesWidth;
};

QT_CHARTS_END_NAMESPACE

#endif // ABSTRACTBARCHARTITEM_H
