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

#ifndef QVCANDLESTICKMODELMAPPER_H
#define QVCANDLESTICKMODELMAPPER_H

#include <QtCharts/QCandlestickModelMapper>

QT_CHARTS_BEGIN_NAMESPACE
/* Comment line for syncqt to generate the fwd-include correctly, due to QTBUG-22432 */
class Q_CHARTS_EXPORT QVCandlestickModelMapper : public QCandlestickModelMapper
{
    Q_OBJECT
    Q_PROPERTY(int timestampRow READ timestampRow WRITE setTimestampRow NOTIFY timestampRowChanged)
    Q_PROPERTY(int openRow READ openRow WRITE setOpenRow NOTIFY openRowChanged)
    Q_PROPERTY(int highRow READ highRow WRITE setHighRow NOTIFY highRowChanged)
    Q_PROPERTY(int lowRow READ lowRow WRITE setLowRow NOTIFY lowRowChanged)
    Q_PROPERTY(int closeRow READ closeRow WRITE setCloseRow NOTIFY closeRowChanged)
    Q_PROPERTY(int firstSetColumn READ firstSetColumn WRITE setFirstSetColumn NOTIFY firstSetColumnChanged)
    Q_PROPERTY(int lastSetColumn READ lastSetColumn WRITE setLastSetColumn NOTIFY lastSetColumnChanged)

public:
    explicit QVCandlestickModelMapper(QObject *parent = nullptr);

    Qt::Orientation orientation() const;

    void setTimestampRow(int timestampRow);
    int timestampRow() const;

    void setOpenRow(int openRow);
    int openRow() const;

    void setHighRow(int highRow);
    int highRow() const;

    void setLowRow(int lowRow);
    int lowRow() const;

    void setCloseRow(int closeRow);
    int closeRow() const;

    void setFirstSetColumn(int firstSetColumn);
    int firstSetColumn() const;

    void setLastSetColumn(int lastSetColumn);
    int lastSetColumn() const;

Q_SIGNALS:
    void timestampRowChanged();
    void openRowChanged();
    void highRowChanged();
    void lowRowChanged();
    void closeRowChanged();
    void firstSetColumnChanged();
    void lastSetColumnChanged();
};

QT_CHARTS_END_NAMESPACE

#endif // QVCANDLESTICKMODELMAPPER_H
