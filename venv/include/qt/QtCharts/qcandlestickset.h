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

#ifndef QCANDLESTICKSET_H
#define QCANDLESTICKSET_H

#include <QtCharts/QChartGlobal>
#include <QtGui/QBrush>
#include <QtGui/QPen>

QT_CHARTS_BEGIN_NAMESPACE

class QCandlestickSetPrivate;

class Q_CHARTS_EXPORT QCandlestickSet : public QObject
{
    Q_OBJECT
    Q_PROPERTY(qreal timestamp READ timestamp WRITE setTimestamp NOTIFY timestampChanged)
    Q_PROPERTY(qreal open READ open WRITE setOpen NOTIFY openChanged)
    Q_PROPERTY(qreal high READ high WRITE setHigh NOTIFY highChanged)
    Q_PROPERTY(qreal low READ low WRITE setLow NOTIFY lowChanged)
    Q_PROPERTY(qreal close READ close WRITE setClose NOTIFY closeChanged)
    Q_PROPERTY(QBrush brush READ brush WRITE setBrush NOTIFY brushChanged)
    Q_PROPERTY(QPen pen READ pen WRITE setPen NOTIFY penChanged)

public:
    explicit QCandlestickSet(qreal timestamp = 0.0, QObject *parent = nullptr);
    explicit QCandlestickSet(qreal open, qreal high, qreal low, qreal close, qreal timestamp = 0.0,
                             QObject *parent = nullptr);
    virtual ~QCandlestickSet();

    void setTimestamp(qreal timestamp);
    qreal timestamp() const;

    void setOpen(qreal open);
    qreal open() const;

    void setHigh(qreal high);
    qreal high() const;

    void setLow(qreal low);
    qreal low() const;

    void setClose(qreal close);
    qreal close() const;

    void setBrush(const QBrush &brush);
    QBrush brush() const;

    void setPen(const QPen &pen);
    QPen pen() const;

Q_SIGNALS:
    void clicked();
    void hovered(bool status);
    void pressed();
    void released();
    void doubleClicked();
    void timestampChanged();
    void openChanged();
    void highChanged();
    void lowChanged();
    void closeChanged();
    void brushChanged();
    void penChanged();

private:
    QScopedPointer<QCandlestickSetPrivate> d_ptr;
    Q_DECLARE_PRIVATE(QCandlestickSet)
    Q_DISABLE_COPY(QCandlestickSet)
    friend class QCandlestickSeriesPrivate;
};

QT_CHARTS_END_NAMESPACE

#endif // QCANDLESTICKSET_H
