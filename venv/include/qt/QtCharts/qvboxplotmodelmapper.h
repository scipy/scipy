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

#ifndef QVBOXPLOTMODELMAPPER_H
#define QVBOXPLOTMODELMAPPER_H

#include <QtCharts/QBoxPlotModelMapper>

QT_CHARTS_BEGIN_NAMESPACE
/* Comment line for syncqt to generate the fwd-include correctly, due to QTBUG-22432 */
class Q_CHARTS_EXPORT QVBoxPlotModelMapper : public QBoxPlotModelMapper
{
    Q_OBJECT
    Q_PROPERTY(QtCharts::QBoxPlotSeries *series READ series WRITE setSeries NOTIFY seriesReplaced)
    Q_PROPERTY(QAbstractItemModel *model READ model WRITE setModel NOTIFY modelReplaced)
    Q_PROPERTY(int firstBoxSetColumn READ firstBoxSetColumn WRITE setFirstBoxSetColumn NOTIFY firstBoxSetColumnChanged)
    Q_PROPERTY(int lastBoxSetColumn READ lastBoxSetColumn WRITE setLastBoxSetColumn NOTIFY lastBoxSetColumnChanged)
    Q_PROPERTY(int firstRow READ firstRow WRITE setFirstRow NOTIFY firstRowChanged)
    Q_PROPERTY(int rowCount READ rowCount WRITE setRowCount NOTIFY rowCountChanged)

public:
    explicit QVBoxPlotModelMapper(QObject *parent = nullptr);

    QAbstractItemModel *model() const;
    void setModel(QAbstractItemModel *model);

    QBoxPlotSeries *series() const;
    void setSeries(QBoxPlotSeries *series);

    int firstBoxSetColumn() const;
    void setFirstBoxSetColumn(int firstBoxSetColumn);

    int lastBoxSetColumn() const;
    void setLastBoxSetColumn(int lastBoxSetColumn);

    int firstRow() const;
    void setFirstRow(int firstRow);

    int rowCount() const;
    void setRowCount(int rowCount);

Q_SIGNALS:
    void seriesReplaced();
    void modelReplaced();
    void firstBoxSetColumnChanged();
    void lastBoxSetColumnChanged();
    void firstRowChanged();
    void rowCountChanged();
};

QT_CHARTS_END_NAMESPACE

#endif // QVBOXPLOTMODELMAPPER_H
