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

//  W A R N I N G
//  -------------
//
// This file is not part of the Qt Chart API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#ifndef DECLARATIVECATEGORYAXIS_H
#define DECLARATIVECATEGORYAXIS_H

#include <QtCharts/QCategoryAxis>
#include <private/declarativechartglobal_p.h>

#include <QtQml/QQmlListProperty>
#include <QtQml/QQmlParserStatus>

QT_CHARTS_BEGIN_NAMESPACE

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeCategoryRange : public QObject
{
    Q_OBJECT
    Q_PROPERTY(qreal endValue READ endValue WRITE setEndValue)
    Q_PROPERTY(QString label READ label WRITE setLabel)

public:
    explicit DeclarativeCategoryRange(QObject *parent = 0);
    qreal endValue() { return m_endValue; }
    void setEndValue(qreal endValue) { m_endValue = endValue; }
    QString label() { return m_label; }
    void setLabel(const QString &label);

private:
    qreal m_endValue;
    QString m_label;
};

class DeclarativeCategoryAxis : public QCategoryAxis, public QQmlParserStatus
{
    Q_OBJECT
    Q_INTERFACES(QQmlParserStatus)
    Q_PROPERTY(QQmlListProperty<QObject> axisChildren READ axisChildren)
    Q_CLASSINFO("DefaultProperty", "axisChildren")
    Q_PROPERTY(AxisLabelsPosition labelsPosition READ labelsPosition WRITE setLabelsPosition NOTIFY labelsPositionChanged REVISION 1)
    Q_ENUMS(AxisLabelsPosition)

public:
    // duplicating enums from QChart to make the QML api namings 1-to-1 with the C++ api
    enum AxisLabelsPosition {
        AxisLabelsPositionCenter = 0x0,
        AxisLabelsPositionOnValue = 0x1
    };

    explicit DeclarativeCategoryAxis(QObject *parent = 0);
    QQmlListProperty<QObject> axisChildren();


public: // from QDeclarativeParserStatus
    void classBegin();
    void componentComplete();

public:
    AxisLabelsPosition labelsPosition() const;
    void setLabelsPosition(AxisLabelsPosition position);

Q_SIGNALS:
    Q_REVISION(1) void labelsPositionChanged(AxisLabelsPosition position);

public Q_SLOTS:
    Q_INVOKABLE void append(const QString &label, qreal categoryEndValue);
    Q_INVOKABLE void remove(const QString &label);
    Q_INVOKABLE void replace(const QString &oldLabel, const QString &newLabel);
    static void appendAxisChildren(QQmlListProperty<QObject> *list, QObject *element);

private:
    static bool endValueLessThan(const QPair<QString, qreal> &value1, const QPair<QString, qreal> &value2);

private:
    AxisLabelsPosition m_labelsPosition;
};

QT_CHARTS_END_NAMESPACE

#endif // DECLARATIVECATEGORYAXIS_H
