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

#ifndef QCATEGORYAXIS_H
#define QCATEGORYAXIS_H

#include <QtCharts/QAbstractAxis>
#include <QtCharts/QValueAxis>

QT_CHARTS_BEGIN_NAMESPACE

class QCategoryAxisPrivate;

class Q_CHARTS_EXPORT QCategoryAxis : public QValueAxis
{
    Q_OBJECT
    Q_PROPERTY(qreal startValue READ startValue WRITE setStartValue)
    Q_PROPERTY(int count READ count)
    Q_PROPERTY(QStringList categoriesLabels READ categoriesLabels)
    Q_PROPERTY(AxisLabelsPosition labelsPosition READ labelsPosition WRITE setLabelsPosition NOTIFY labelsPositionChanged)
    Q_ENUMS(AxisLabelsPosition)

public:

    enum AxisLabelsPosition {
        AxisLabelsPositionCenter = 0x0,
        AxisLabelsPositionOnValue = 0x1
    };

    explicit QCategoryAxis(QObject *parent = nullptr);
    ~QCategoryAxis();

protected:
    QCategoryAxis(QCategoryAxisPrivate &d, QObject *parent = nullptr);

public:
    AxisType type() const;

    void append(const QString &label, qreal categoryEndValue);
    void remove(const QString &label);
    void replaceLabel(const QString &oldLabel, const QString &newLabel);

    qreal startValue(const QString &categoryLabel = QString()) const;
    void setStartValue(qreal min);

    qreal endValue(const QString &categoryLabel) const;

    QStringList categoriesLabels();
    int count() const;

    QCategoryAxis::AxisLabelsPosition labelsPosition() const;
    void setLabelsPosition(QCategoryAxis::AxisLabelsPosition position);

Q_SIGNALS:
    void categoriesChanged();
    void labelsPositionChanged(QCategoryAxis::AxisLabelsPosition position);

private:
    Q_DECLARE_PRIVATE(QCategoryAxis)
    Q_DISABLE_COPY(QCategoryAxis)
};

QT_CHARTS_END_NAMESPACE

#endif // QCATEGORYAXIS_H
