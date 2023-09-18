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

#ifndef QBARMODELMAPPER_H
#define QBARMODELMAPPER_H

#include <QtCharts/QChartGlobal>
#include <QtCore/QObject>

QT_BEGIN_NAMESPACE
class QAbstractItemModel;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE

class QBarModelMapperPrivate;
class QAbstractBarSeries;

class Q_CHARTS_EXPORT QBarModelMapper : public QObject
{
    Q_OBJECT

protected:
    explicit QBarModelMapper(QObject *parent = nullptr);

    QAbstractItemModel *model() const;
    void setModel(QAbstractItemModel *model);

    QAbstractBarSeries *series() const;
    void setSeries(QAbstractBarSeries *series);

    int first() const;
    void setFirst(int first);

    int count() const;
    void setCount(int count);

    int firstBarSetSection() const;
    void setFirstBarSetSection(int firstBarSetSection);

    int lastBarSetSection() const;
    void setLastBarSetSection(int lastBarSetSection);

    Qt::Orientation orientation() const;
    void setOrientation(Qt::Orientation orientation);

protected:
    QBarModelMapperPrivate * const d_ptr;
    Q_DECLARE_PRIVATE(QBarModelMapper)
};

QT_CHARTS_END_NAMESPACE

#endif // QBARMODELMAPPER_H
