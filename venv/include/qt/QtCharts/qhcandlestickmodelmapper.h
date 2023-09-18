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

#ifndef QHCANDLESTICKMODELMAPPER_H
#define QHCANDLESTICKMODELMAPPER_H

#include <QtCharts/QCandlestickModelMapper>

QT_CHARTS_BEGIN_NAMESPACE
/* Comment line for syncqt to generate the fwd-include correctly, due to QTBUG-22432 */
class Q_CHARTS_EXPORT QHCandlestickModelMapper : public QCandlestickModelMapper
{
    Q_OBJECT
    Q_PROPERTY(int timestampColumn READ timestampColumn WRITE setTimestampColumn NOTIFY timestampColumnChanged)
    Q_PROPERTY(int openColumn READ openColumn WRITE setOpenColumn NOTIFY openColumnChanged)
    Q_PROPERTY(int highColumn READ highColumn WRITE setHighColumn NOTIFY highColumnChanged)
    Q_PROPERTY(int lowColumn READ lowColumn WRITE setLowColumn NOTIFY lowColumnChanged)
    Q_PROPERTY(int closeColumn READ closeColumn WRITE setCloseColumn NOTIFY closeColumnChanged)
    Q_PROPERTY(int firstSetRow READ firstSetRow WRITE setFirstSetRow NOTIFY firstSetRowChanged)
    Q_PROPERTY(int lastSetRow READ lastSetRow WRITE setLastSetRow NOTIFY lastSetRowChanged)

public:
    explicit QHCandlestickModelMapper(QObject *parent = nullptr);

    Qt::Orientation orientation() const;

    void setTimestampColumn(int timestampColumn);
    int timestampColumn() const;

    void setOpenColumn(int openColumn);
    int openColumn() const;

    void setHighColumn(int highColumn);
    int highColumn() const;

    void setLowColumn(int lowColumn);
    int lowColumn() const;

    void setCloseColumn(int closeColumn);
    int closeColumn() const;

    void setFirstSetRow(int firstSetRow);
    int firstSetRow() const;

    void setLastSetRow(int lastSetRow);
    int lastSetRow() const;

Q_SIGNALS:
    void timestampColumnChanged();
    void openColumnChanged();
    void highColumnChanged();
    void lowColumnChanged();
    void closeColumnChanged();
    void firstSetRowChanged();
    void lastSetRowChanged();
};

QT_CHARTS_END_NAMESPACE

#endif // QHCANDLESTICKMODELMAPPER_H
