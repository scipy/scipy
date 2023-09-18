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

#ifndef CHARTDATASET_P_H
#define CHARTDATASET_P_H

#include <QtCharts/QAbstractSeries>
#include <private/abstractdomain_p.h>
#include <private/qabstractaxis_p.h>
#include <QtCharts/private/qchartglobal_p.h>
#include <QtCore/QVector>

QT_CHARTS_BEGIN_NAMESPACE

class QAbstractAxis;
class ChartPresenter;
class GLXYSeriesDataManager;

class Q_CHARTS_PRIVATE_EXPORT ChartDataSet : public QObject
{
    Q_OBJECT
public:
    ChartDataSet(QChart *chart);
    virtual ~ChartDataSet();

    void addSeries(QAbstractSeries *series);
    void removeSeries(QAbstractSeries *series);
    QList<QAbstractSeries *> series() const;

    void addAxis(QAbstractAxis *axis,Qt::Alignment aligment);
    void removeAxis(QAbstractAxis *axis);
    QList<QAbstractAxis*> axes() const;

    bool attachAxis(QAbstractSeries* series,QAbstractAxis *axis);
    bool detachAxis(QAbstractSeries* series,QAbstractAxis *axis);

    void createDefaultAxes();

    void zoomInDomain(const QRectF &rect);
    void zoomOutDomain(const QRectF &rect);
    void zoomResetDomain();
    bool isZoomedDomain();
    void scrollDomain(qreal dx, qreal dy);

    QPointF mapToValue(const QPointF &position, QAbstractSeries *series = 0);
    QPointF mapToPosition(const QPointF &value, QAbstractSeries *series = 0);

    GLXYSeriesDataManager *glXYSeriesDataManager() { return m_glXYSeriesDataManager; }

    AbstractDomain* createDomain(AbstractDomain::DomainType type);
    AbstractDomain* domainForSeries(QAbstractSeries *series) const;

Q_SIGNALS:
    void axisAdded(QAbstractAxis* axis);
    void axisRemoved(QAbstractAxis* axis);
    void seriesAdded(QAbstractSeries* series);
    void seriesRemoved(QAbstractSeries* series);
public Q_SLOTS:
    void reverseChanged();
private:
    void createAxes(QAbstractAxis::AxisTypes type, Qt::Orientation orientation);
    QAbstractAxis *createAxis(QAbstractAxis::AxisType type, Qt::Orientation orientation);
    AbstractDomain::DomainType selectDomain(QList<QAbstractAxis* > axes);
    void deleteAllAxes();
    void deleteAllSeries();
    void findMinMaxForSeries(QList<QAbstractSeries *> series,Qt::Orientations orientation, qreal &min, qreal &max);
private:
    QList<QAbstractSeries *> m_seriesList;
    QList<QAbstractAxis *> m_axisList;
    QChart* m_chart;
    GLXYSeriesDataManager *m_glXYSeriesDataManager;
};

QT_CHARTS_END_NAMESPACE

#endif /* CHARTENGINE_P_H */
