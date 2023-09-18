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

#ifndef QLOGVALUEAXIS_H
#define QLOGVALUEAXIS_H

#include <QtCharts/QAbstractAxis>

QT_BEGIN_NAMESPACE
class QDateTime;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE

class QLogValueAxisPrivate;

class Q_CHARTS_EXPORT QLogValueAxis : public QAbstractAxis
{
    Q_OBJECT
    Q_PROPERTY(qreal min READ min WRITE setMin NOTIFY minChanged)
    Q_PROPERTY(qreal max READ max WRITE setMax NOTIFY maxChanged)
    Q_PROPERTY(QString labelFormat READ labelFormat WRITE setLabelFormat NOTIFY labelFormatChanged)
    Q_PROPERTY(qreal base READ base WRITE setBase NOTIFY baseChanged)
    Q_PROPERTY(int tickCount READ tickCount NOTIFY tickCountChanged)
    Q_PROPERTY(int minorTickCount READ minorTickCount WRITE setMinorTickCount NOTIFY minorTickCountChanged)

public:
    explicit QLogValueAxis(QObject *parent = nullptr);
    ~QLogValueAxis();

protected:
    QLogValueAxis(QLogValueAxisPrivate &d, QObject *parent = nullptr);

public:
    AxisType type() const;

    //range handling
    void setMin(qreal min);
    qreal min() const;
    void setMax(qreal max);
    qreal max() const;
    void setRange(qreal min, qreal max);

    void setLabelFormat(const QString &format);
    QString labelFormat() const;

    void setBase(qreal base);
    qreal base() const;

    int tickCount() const;

    void setMinorTickCount(int minorTickCount);
    int minorTickCount() const;

Q_SIGNALS:
    void minChanged(qreal min);
    void maxChanged(qreal max);
    void rangeChanged(qreal min, qreal max);
    void labelFormatChanged(const QString &format);
    void baseChanged(qreal base);
    void tickCountChanged(int tickCount);
    void minorTickCountChanged(int minorTickCount);

private:
    Q_DECLARE_PRIVATE(QLogValueAxis)
    Q_DISABLE_COPY(QLogValueAxis)
};

QT_CHARTS_END_NAMESPACE

#endif // QLOGVALUEAXIS_H
