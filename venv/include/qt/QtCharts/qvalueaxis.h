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

#ifndef QVALUEAXIS_H
#define QVALUEAXIS_H

#include <QtCharts/QAbstractAxis>

QT_CHARTS_BEGIN_NAMESPACE

class QValueAxisPrivate;

class Q_CHARTS_EXPORT QValueAxis : public QAbstractAxis
{
    Q_OBJECT
    Q_PROPERTY(int tickCount READ tickCount WRITE setTickCount NOTIFY tickCountChanged)
    Q_PROPERTY(qreal min READ min WRITE setMin NOTIFY minChanged)
    Q_PROPERTY(qreal max READ max WRITE setMax NOTIFY maxChanged)
    Q_PROPERTY(QString labelFormat READ labelFormat WRITE setLabelFormat NOTIFY labelFormatChanged)
    Q_PROPERTY(int minorTickCount READ minorTickCount WRITE setMinorTickCount NOTIFY minorTickCountChanged)
    Q_PROPERTY(qreal tickAnchor READ tickAnchor WRITE setTickAnchor NOTIFY tickAnchorChanged REVISION 1)
    Q_PROPERTY(qreal tickInterval READ tickInterval WRITE setTickInterval NOTIFY tickIntervalChanged REVISION 1)
    Q_PROPERTY(TickType tickType READ tickType WRITE setTickType NOTIFY tickTypeChanged REVISION 1)
    Q_ENUMS(TickType)

public:
    enum TickType {
        TicksDynamic = 0,
        TicksFixed
    };

public:
    explicit QValueAxis(QObject *parent = nullptr);
    ~QValueAxis();

protected:
    QValueAxis(QValueAxisPrivate &d, QObject *parent = nullptr);

public:
    AxisType type() const;

    //range handling
    void setMin(qreal min);
    qreal min() const;
    void setMax(qreal max);
    qreal max() const;
    void setRange(qreal min, qreal max);

    //ticks handling
    void setTickCount(int count);
    int tickCount() const;
    void setMinorTickCount(int count);
    int minorTickCount() const;
    void setTickAnchor(qreal anchor);
    qreal tickAnchor() const;
    void setTickInterval(qreal insterval);
    qreal tickInterval() const;
    void setTickType(QValueAxis::TickType type);
    QValueAxis::TickType tickType() const;

    void setLabelFormat(const QString &format);
    QString labelFormat() const;

public Q_SLOTS:
    void applyNiceNumbers();

Q_SIGNALS:
    void minChanged(qreal min);
    void maxChanged(qreal max);
    void rangeChanged(qreal min, qreal max);
    void tickCountChanged(int tickCount);
    void minorTickCountChanged(int tickCount);
    void labelFormatChanged(const QString &format);
    Q_REVISION(1) void tickIntervalChanged(qreal interval);
    Q_REVISION(1) void tickAnchorChanged(qreal anchor);
    Q_REVISION(1) void tickTypeChanged(QValueAxis::TickType type);

private:
    Q_DECLARE_PRIVATE(QValueAxis)
    Q_DISABLE_COPY(QValueAxis)
};

QT_CHARTS_END_NAMESPACE

#endif // QVALUEAXIS_H
