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

#ifndef QBARCATEGORYAXIS_H
#define QBARCATEGORYAXIS_H

#include <QtCharts/QAbstractAxis>

QT_CHARTS_BEGIN_NAMESPACE

class QBarCategoryAxisPrivate;

class Q_CHARTS_EXPORT QBarCategoryAxis : public QAbstractAxis
{
    Q_OBJECT
    Q_PROPERTY(QStringList categories READ categories WRITE setCategories NOTIFY categoriesChanged)
    Q_PROPERTY(QString min READ min WRITE setMin NOTIFY minChanged)
    Q_PROPERTY(QString max READ max WRITE setMax NOTIFY maxChanged)
    Q_PROPERTY(int count READ count NOTIFY countChanged)

public:
    explicit QBarCategoryAxis(QObject *parent = nullptr);
    ~QBarCategoryAxis();

protected:
    QBarCategoryAxis(QBarCategoryAxisPrivate &d, QObject *parent = nullptr);

public:
    AxisType type() const;
    void append(const QStringList &categories);
    void append(const QString &category);
    void remove(const QString &category);
    void insert(int index, const QString &category);
    void replace(const QString &oldCategory, const QString &newCategory);
    Q_INVOKABLE void clear();
    void setCategories(const QStringList &categories);
    QStringList categories();
    int count() const;
    QString at(int index) const;

    //range handling
    void setMin(const QString &minCategory);
    QString min() const;
    void setMax(const QString &maxCategory);
    QString max() const;
    void setRange(const QString &minCategory, const QString &maxCategory);

Q_SIGNALS:
    void categoriesChanged();
    void minChanged(const QString &min);
    void maxChanged(const QString &max);
    void rangeChanged(const QString &min, const QString &max);
    void countChanged();

private:
    Q_DECLARE_PRIVATE(QBarCategoryAxis)
    Q_DISABLE_COPY(QBarCategoryAxis)
    friend class ChartBarCategoryAxisX;
    friend class ChartBarCategoryAxisY;
};

QT_CHARTS_END_NAMESPACE

#endif // QBARCATEGORYAXIS_H
