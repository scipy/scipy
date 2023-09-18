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

#ifndef QCANDLESTICKMODELMAPPER_P_H
#define QCANDLESTICKMODELMAPPER_P_H

#include <QtCharts/QCandlestickModelMapper>
#include <QtCharts/private/qchartglobal_p.h>
#include <QtCore/QObject>

QT_BEGIN_NAMESPACE
class QModelIndex;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE

class QCandlestickSet;

class Q_CHARTS_PRIVATE_EXPORT QCandlestickModelMapperPrivate : public QObject
{
    Q_OBJECT

public:
    explicit QCandlestickModelMapperPrivate(QCandlestickModelMapper *q);

Q_SIGNALS:
    void timestampChanged();
    void openChanged();
    void highChanged();
    void lowChanged();
    void closeChanged();
    void firstSetSectionChanged();
    void lastSetSectionChanged();

private Q_SLOTS:
    void initializeCandlestickFromModel();

    // for the model
    void modelDataUpdated(QModelIndex topLeft, QModelIndex bottomRight);
    void modelHeaderDataUpdated(Qt::Orientation orientation, int first, int last);
    void modelRowsInserted(QModelIndex parent, int start, int end);
    void modelRowsRemoved(QModelIndex parent, int start, int end);
    void modelColumnsInserted(QModelIndex parent, int start, int end);
    void modelColumnsRemoved(QModelIndex parent, int start, int end);
    void modelDestroyed();

    // for the series
    void candlestickSetsAdded(const QList<QCandlestickSet *> &sets);
    void candlestickSetsRemoved(const QList<QCandlestickSet *> &sets);
    void candlestickSetChanged();
    void seriesDestroyed();

private:
    QCandlestickSet *candlestickSet(QModelIndex index);
    QModelIndex candlestickModelIndex(int section, int pos);
    void insertData(int start, int end);
    void removeData(int start, int end);
    void blockModelSignals(bool block = true);
    void blockSeriesSignals(bool block = true);

private:
    QAbstractItemModel *m_model;
    QCandlestickSeries *m_series;
    int m_timestamp;
    int m_open;
    int m_high;
    int m_low;
    int m_close;
    int m_firstSetSection;
    int m_lastSetSection;
    QList<QCandlestickSet *> m_sets;
    bool m_modelSignalsBlock;
    bool m_seriesSignalsBlock;

private:
    QCandlestickModelMapper *q_ptr;
    Q_DECLARE_PUBLIC(QCandlestickModelMapper)
};

QT_CHARTS_END_NAMESPACE

#endif // QCANDLESTICKMODELMAPPER_P_H
