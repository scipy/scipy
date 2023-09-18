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

#ifndef QBARMODELMAPPER_P_H
#define QBARMODELMAPPER_P_H

#include <QtCore/QObject>
#include <QtCharts/QBarModelMapper>
#include <QtCharts/private/qchartglobal_p.h>

QT_BEGIN_NAMESPACE
class QModelIndex;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE

class QBarSet;

class Q_CHARTS_PRIVATE_EXPORT QBarModelMapperPrivate : public QObject
{
    Q_OBJECT
public:
    explicit QBarModelMapperPrivate(QBarModelMapper *q);

public Q_SLOTS:
    // for the model
    void modelUpdated(QModelIndex topLeft, QModelIndex bottomRight);
    void modelHeaderDataUpdated(Qt::Orientation orientation, int first, int last);
    void modelRowsAdded(QModelIndex parent, int start, int end);
    void modelRowsRemoved(QModelIndex parent, int start, int end);
    void modelColumnsAdded(QModelIndex parent, int start, int end);
    void modelColumnsRemoved(QModelIndex parent, int start, int end);
    void handleModelDestroyed();

    // for the series
    void barSetsAdded(QList<QBarSet *> sets);
    void barSetsRemoved(QList<QBarSet *> sets);
    void valuesAdded(int index, int count);
    void valuesRemoved(int index, int count);
    void barLabelChanged();
    void barValueChanged(int index);
    void handleSeriesDestroyed();

    void initializeBarFromModel();

private:
    QBarSet *barSet(QModelIndex index);
    QModelIndex barModelIndex(int barSection, int posInBar);
    void insertData(int start, int end);
    void removeData(int start, int end);
    void blockModelSignals(bool block = true);
    void blockSeriesSignals(bool block = true);

private:
    QAbstractBarSeries *m_series;
    QList<QBarSet *> m_barSets;
    QAbstractItemModel *m_model;
    int m_first;
    int m_count;
    Qt::Orientation m_orientation;
    int m_firstBarSetSection;
    int m_lastBarSetSection;
    bool m_seriesSignalsBlock;
    bool m_modelSignalsBlock;

private:
    QBarModelMapper *q_ptr;
    Q_DECLARE_PUBLIC(QBarModelMapper)
};

QT_CHARTS_END_NAMESPACE

#endif // QBARMODELMAPPER_P_H
