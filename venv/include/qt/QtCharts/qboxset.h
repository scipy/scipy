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

#ifndef QBOXSET_H
#define QBOXSET_H

#include <QtCharts/QChartGlobal>
#include <QtGui/QPen>
#include <QtGui/QBrush>
#include <QtGui/QFont>

QT_CHARTS_BEGIN_NAMESPACE
class QBoxSetPrivate;

class Q_CHARTS_EXPORT QBoxSet : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QPen pen READ pen WRITE setPen NOTIFY penChanged)
    Q_PROPERTY(QBrush brush READ brush WRITE setBrush NOTIFY brushChanged)

public:
    enum ValuePositions {
        LowerExtreme,
        LowerQuartile,
        Median,
        UpperQuartile,
        UpperExtreme
    };

public:
    explicit QBoxSet(const QString label = QString(), QObject *parent = nullptr);
    explicit QBoxSet(const qreal le, const qreal lq, const qreal m, const qreal uq, const qreal ue, const QString label = QString(), QObject *parent = nullptr);
    virtual ~QBoxSet();

    void append(const qreal value);
    void append(const QList<qreal> &values);

    void clear();

    void setLabel(const QString label);
    QString label() const;

    QBoxSet &operator << (const qreal &value);

    void setValue(const int index, const qreal value);
    qreal at(const int index) const;
    qreal operator [](const int index) const;
    int count() const;

    void setPen(const QPen &pen);
    QPen pen() const;

    void setBrush(const QBrush &brush);
    QBrush brush() const;

Q_SIGNALS:
    void clicked();
    void hovered(bool status);
    void pressed();
    void released();
    void doubleClicked();
    void penChanged();
    void brushChanged();

    void valuesChanged();
    void valueChanged(int index);
    void cleared();

private:
    QScopedPointer<QBoxSetPrivate> d_ptr;
    Q_DISABLE_COPY(QBoxSet)
    friend class BarLegendMarker;
    friend class BarChartItem;
    friend class BoxPlotChartItem;
    friend class QBoxPlotSeriesPrivate;
};

QT_CHARTS_END_NAMESPACE

#endif // QBOXSET_H
