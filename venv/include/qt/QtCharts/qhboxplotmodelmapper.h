/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QHBOXPLOTMODELMAPPER_H
#define QHBOXPLOTMODELMAPPER_H

#include <QtCharts/QBoxPlotModelMapper>

QT_CHARTS_BEGIN_NAMESPACE
/* Comment line for syncqt to generate the fwd-include correctly, due to QTBUG-22432 */
class Q_CHARTS_EXPORT QHBoxPlotModelMapper : public QBoxPlotModelMapper
{
    Q_OBJECT
    Q_PROPERTY(QtCharts::QBoxPlotSeries *series READ series WRITE setSeries NOTIFY seriesReplaced)
    Q_PROPERTY(QAbstractItemModel *model READ model WRITE setModel NOTIFY modelReplaced)
    Q_PROPERTY(int firstBoxSetRow READ firstBoxSetRow WRITE setFirstBoxSetRow NOTIFY firstBoxSetRowChanged)
    Q_PROPERTY(int lastBoxSetRow READ lastBoxSetRow WRITE setLastBoxSetRow NOTIFY lastBoxSetRowChanged)
    Q_PROPERTY(int firstColumn READ firstColumn WRITE setFirstColumn NOTIFY firstColumnChanged)
    Q_PROPERTY(int columnCount READ columnCount WRITE setColumnCount NOTIFY columnCountChanged)

public:
    explicit QHBoxPlotModelMapper(QObject *parent = nullptr);

    QAbstractItemModel *model() const;
    void setModel(QAbstractItemModel *model);

    QBoxPlotSeries *series() const;
    void setSeries(QBoxPlotSeries *series);

    int firstBoxSetRow() const;
    void setFirstBoxSetRow(int firstBoxSetRow);

    int lastBoxSetRow() const;
    void setLastBoxSetRow(int lastBoxSetRow);

    int firstColumn() const;
    void setFirstColumn(int firstColumn);

    int columnCount() const;
    void setColumnCount(int rowCount);

Q_SIGNALS:
    void seriesReplaced();
    void modelReplaced();
    void firstBoxSetRowChanged();
    void lastBoxSetRowChanged();
    void firstColumnChanged();
    void columnCountChanged();
};

QT_CHARTS_END_NAMESPACE

#endif // QHBOXPLOTMODELMAPPER_H
