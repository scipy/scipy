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

#ifndef QCANDLESTICKMODELMAPPER_H
#define QCANDLESTICKMODELMAPPER_H

#include <QtCharts/QChartGlobal>
#include <QtCore/QObject>

QT_BEGIN_NAMESPACE
class QAbstractItemModel;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE

class QCandlestickModelMapperPrivate;
class QCandlestickSeries;

class Q_CHARTS_EXPORT QCandlestickModelMapper : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QAbstractItemModel *model READ model WRITE setModel NOTIFY modelReplaced)
    Q_PROPERTY(QtCharts::QCandlestickSeries *series READ series WRITE setSeries NOTIFY seriesReplaced)

public:
    explicit QCandlestickModelMapper(QObject *parent = nullptr);

    void setModel(QAbstractItemModel *model);
    QAbstractItemModel *model() const;

    void setSeries(QCandlestickSeries *series);
    QCandlestickSeries *series() const;

    virtual Qt::Orientation orientation() const = 0;

Q_SIGNALS:
    void modelReplaced();
    void seriesReplaced();

protected:
    void setTimestamp(int timestamp);
    int timestamp() const;

    void setOpen(int open);
    int open() const;

    void setHigh(int high);
    int high() const;

    void setLow(int low);
    int low() const;

    void setClose(int close);
    int close() const;

    void setFirstSetSection(int firstSetSection);
    int firstSetSection() const;

    void setLastSetSection(int lastSetSection);
    int lastSetSection() const;

protected:
    QCandlestickModelMapperPrivate * const d_ptr;
    Q_DECLARE_PRIVATE(QCandlestickModelMapper)
};

QT_CHARTS_END_NAMESPACE

#endif // QCANDLESTICKMODELMAPPER_H
