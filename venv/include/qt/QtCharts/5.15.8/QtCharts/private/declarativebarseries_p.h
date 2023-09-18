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

#ifndef DECLARATIVEBARSERIES_H
#define DECLARATIVEBARSERIES_H

#include <QtCharts/QBarSeries>
#include <QtCharts/QStackedBarSeries>
#include <QtCharts/QPercentBarSeries>
#include <QtCharts/QHorizontalBarSeries>
#include <QtCharts/QHorizontalStackedBarSeries>
#include <QtCharts/QHorizontalPercentBarSeries>
#include <QtCharts/QBarSet>
#include <private/declarativechartglobal_p.h>
#include <private/declarativeaxes_p.h>

#include <QtQuick/QQuickItem>
#include <QtQml/QQmlParserStatus>

QT_CHARTS_BEGIN_NAMESPACE

class QChart;

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeBarSet : public QBarSet
{
    Q_OBJECT
    Q_PROPERTY(QVariantList values READ values WRITE setValues)
    Q_PROPERTY(qreal borderWidth READ borderWidth WRITE setBorderWidth NOTIFY borderWidthChanged REVISION 1)
    Q_PROPERTY(int count READ count NOTIFY countChanged)
    Q_PROPERTY(QString brushFilename READ brushFilename WRITE setBrushFilename NOTIFY brushFilenameChanged REVISION 2)

public:
    explicit DeclarativeBarSet(QObject *parent = 0);
    QVariantList values();
    void setValues(QVariantList values);
    qreal borderWidth() const;
    void setBorderWidth(qreal borderWidth);
    QString brushFilename() const;
    void setBrushFilename(const QString &brushFilename);

public: // From QBarSet
    Q_INVOKABLE void append(qreal value) { QBarSet::append(value); }
    Q_INVOKABLE void remove(const int index, const int count = 1) { QBarSet::remove(index, count); }
    Q_INVOKABLE void replace(int index, qreal value) { QBarSet::replace(index, value); }
    Q_INVOKABLE qreal at(int index) { return QBarSet::at(index); }

Q_SIGNALS:
    void countChanged(int count);
    Q_REVISION(1) void borderWidthChanged(qreal width);
    Q_REVISION(2) void brushFilenameChanged(const QString &brushFilename);

private Q_SLOTS:
    void handleCountChanged(int index, int count);
    void handleBrushChanged();

private:
    QString m_brushFilename;
    QImage m_brushImage;
};

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeBarSeries : public QBarSeries, public QQmlParserStatus
{
    Q_OBJECT
    Q_INTERFACES(QQmlParserStatus)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisX READ axisX WRITE setAxisX NOTIFY axisXChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisY READ axisY WRITE setAxisY NOTIFY axisYChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisXTop READ axisXTop WRITE setAxisXTop NOTIFY axisXTopChanged REVISION 2)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisYRight READ axisYRight WRITE setAxisYRight NOTIFY axisYRightChanged REVISION 2)
    Q_PROPERTY(QQmlListProperty<QObject> seriesChildren READ seriesChildren)
    Q_CLASSINFO("DefaultProperty", "seriesChildren")

public:
    explicit DeclarativeBarSeries(QQuickItem *parent = 0);
    QAbstractAxis *axisX() { return m_axes->axisX(); }
    void setAxisX(QAbstractAxis *axis) { m_axes->setAxisX(axis); }
    QAbstractAxis *axisY() { return m_axes->axisY(); }
    void setAxisY(QAbstractAxis *axis) { m_axes->setAxisY(axis); }
    QAbstractAxis *axisXTop() { return m_axes->axisXTop(); }
    void setAxisXTop(QAbstractAxis *axis) { m_axes->setAxisXTop(axis); }
    QAbstractAxis *axisYRight() { return m_axes->axisYRight(); }
    void setAxisYRight(QAbstractAxis *axis) { m_axes->setAxisYRight(axis); }
    QQmlListProperty<QObject> seriesChildren();

public:
    Q_INVOKABLE DeclarativeBarSet *at(int index);
    Q_INVOKABLE DeclarativeBarSet *append(QString label, QVariantList values) { return insert(count(), label, values); }
    Q_INVOKABLE DeclarativeBarSet *insert(int index, QString label, QVariantList values);
    Q_INVOKABLE bool remove(QBarSet *barset) { return QBarSeries::remove(barset); }
    Q_INVOKABLE void clear() { return QBarSeries::clear(); }

public: // from QDeclarativeParserStatus
    void classBegin();
    void componentComplete();

Q_SIGNALS:
    Q_REVISION(1) void axisXChanged(QAbstractAxis *axis);
    Q_REVISION(1) void axisYChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisXTopChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisYRightChanged(QAbstractAxis *axis);

public Q_SLOTS:
    static void appendSeriesChildren(QQmlListProperty<QObject> *list, QObject *element);

public:
    DeclarativeAxes *m_axes;
};

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeStackedBarSeries : public QStackedBarSeries, public QQmlParserStatus
{
    Q_OBJECT
    Q_INTERFACES(QQmlParserStatus)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisX READ axisX WRITE setAxisX NOTIFY axisXChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisY READ axisY WRITE setAxisY NOTIFY axisYChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisXTop READ axisXTop WRITE setAxisXTop NOTIFY axisXTopChanged REVISION 2)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisYRight READ axisYRight WRITE setAxisYRight NOTIFY axisYRightChanged REVISION 2)
    Q_PROPERTY(QQmlListProperty<QObject> seriesChildren READ seriesChildren)
    Q_CLASSINFO("DefaultProperty", "seriesChildren")

public:
    explicit DeclarativeStackedBarSeries(QQuickItem *parent = 0);
    QAbstractAxis *axisX() { return m_axes->axisX(); }
    void setAxisX(QAbstractAxis *axis) { m_axes->setAxisX(axis); }
    QAbstractAxis *axisY() { return m_axes->axisY(); }
    void setAxisY(QAbstractAxis *axis) { m_axes->setAxisY(axis); }
    QAbstractAxis *axisXTop() { return m_axes->axisXTop(); }
    void setAxisXTop(QAbstractAxis *axis) { m_axes->setAxisXTop(axis); }
    QAbstractAxis *axisYRight() { return m_axes->axisYRight(); }
    void setAxisYRight(QAbstractAxis *axis) { m_axes->setAxisYRight(axis); }
    QQmlListProperty<QObject> seriesChildren();

public:
    Q_INVOKABLE DeclarativeBarSet *at(int index);
    Q_INVOKABLE DeclarativeBarSet *append(QString label, QVariantList values) { return insert(count(), label, values); }
    Q_INVOKABLE DeclarativeBarSet *insert(int index, QString label, QVariantList values);
    Q_INVOKABLE bool remove(QBarSet *barset) { return QStackedBarSeries::remove(barset); }
    Q_INVOKABLE void clear() { return QStackedBarSeries::clear(); }

public: // from QDeclarativeParserStatus
    void classBegin();
    void componentComplete();

Q_SIGNALS:
    Q_REVISION(1) void axisXChanged(QAbstractAxis *axis);
    Q_REVISION(1) void axisYChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisXTopChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisYRightChanged(QAbstractAxis *axis);

public Q_SLOTS:
    static void appendSeriesChildren(QQmlListProperty<QObject> *list, QObject *element);

public:
    DeclarativeAxes *m_axes;
};

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativePercentBarSeries : public QPercentBarSeries, public QQmlParserStatus
{
    Q_OBJECT
    Q_INTERFACES(QQmlParserStatus)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisX READ axisX WRITE setAxisX NOTIFY axisXChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisY READ axisY WRITE setAxisY NOTIFY axisYChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisXTop READ axisXTop WRITE setAxisXTop NOTIFY axisXTopChanged REVISION 2)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisYRight READ axisYRight WRITE setAxisYRight NOTIFY axisYRightChanged REVISION 2)
    Q_PROPERTY(QQmlListProperty<QObject> seriesChildren READ seriesChildren)
    Q_CLASSINFO("DefaultProperty", "seriesChildren")

public:
    explicit DeclarativePercentBarSeries(QQuickItem *parent = 0);
    QAbstractAxis *axisX() { return m_axes->axisX(); }
    void setAxisX(QAbstractAxis *axis) { m_axes->setAxisX(axis); }
    QAbstractAxis *axisY() { return m_axes->axisY(); }
    void setAxisY(QAbstractAxis *axis) { m_axes->setAxisY(axis); }
    QAbstractAxis *axisXTop() { return m_axes->axisXTop(); }
    void setAxisXTop(QAbstractAxis *axis) { m_axes->setAxisXTop(axis); }
    QAbstractAxis *axisYRight() { return m_axes->axisYRight(); }
    void setAxisYRight(QAbstractAxis *axis) { m_axes->setAxisYRight(axis); }
    QQmlListProperty<QObject> seriesChildren();

public:
    Q_INVOKABLE DeclarativeBarSet *at(int index);
    Q_INVOKABLE DeclarativeBarSet *append(QString label, QVariantList values) { return insert(count(), label, values); }
    Q_INVOKABLE DeclarativeBarSet *insert(int index, QString label, QVariantList values);
    Q_INVOKABLE bool remove(QBarSet *barset) { return QPercentBarSeries::remove(barset); }
    Q_INVOKABLE void clear() { return QPercentBarSeries::clear(); }

public: // from QDeclarativeParserStatus
    void classBegin();
    void componentComplete();

Q_SIGNALS:
    Q_REVISION(1) void axisXChanged(QAbstractAxis *axis);
    Q_REVISION(1) void axisYChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisXTopChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisYRightChanged(QAbstractAxis *axis);

public Q_SLOTS:
    static void appendSeriesChildren(QQmlListProperty<QObject> *list, QObject *element);

public:
    DeclarativeAxes *m_axes;
};

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeHorizontalBarSeries : public QHorizontalBarSeries, public QQmlParserStatus
{
    Q_OBJECT
    Q_INTERFACES(QQmlParserStatus)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisX READ axisX WRITE setAxisX NOTIFY axisXChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisY READ axisY WRITE setAxisY NOTIFY axisYChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisXTop READ axisXTop WRITE setAxisXTop NOTIFY axisXTopChanged REVISION 2)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisYRight READ axisYRight WRITE setAxisYRight NOTIFY axisYRightChanged REVISION 2)
    Q_PROPERTY(QQmlListProperty<QObject> seriesChildren READ seriesChildren)
    Q_CLASSINFO("DefaultProperty", "seriesChildren")

public:
    explicit DeclarativeHorizontalBarSeries(QQuickItem *parent = 0);
    QAbstractAxis *axisX() { return m_axes->axisX(); }
    void setAxisX(QAbstractAxis *axis) { m_axes->setAxisX(axis); }
    QAbstractAxis *axisY() { return m_axes->axisY(); }
    void setAxisY(QAbstractAxis *axis) { m_axes->setAxisY(axis); }
    QAbstractAxis *axisXTop() { return m_axes->axisXTop(); }
    void setAxisXTop(QAbstractAxis *axis) { m_axes->setAxisXTop(axis); }
    QAbstractAxis *axisYRight() { return m_axes->axisYRight(); }
    void setAxisYRight(QAbstractAxis *axis) { m_axes->setAxisYRight(axis); }
    QQmlListProperty<QObject> seriesChildren();

public:
    Q_INVOKABLE DeclarativeBarSet *at(int index);
    Q_INVOKABLE DeclarativeBarSet *append(QString label, QVariantList values) { return insert(count(), label, values); }
    Q_INVOKABLE DeclarativeBarSet *insert(int index, QString label, QVariantList values);
    Q_INVOKABLE bool remove(QBarSet *barset) { return QHorizontalBarSeries::remove(barset); }
    Q_INVOKABLE void clear() { return QHorizontalBarSeries::clear(); }

public: // from QDeclarativeParserStatus
    void classBegin();
    void componentComplete();

Q_SIGNALS:
    Q_REVISION(1) void axisXChanged(QAbstractAxis *axis);
    Q_REVISION(1) void axisYChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisXTopChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisYRightChanged(QAbstractAxis *axis);

public Q_SLOTS:
    static void appendSeriesChildren(QQmlListProperty<QObject> *list, QObject *element);

public:
    DeclarativeAxes *m_axes;
};

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeHorizontalStackedBarSeries : public QHorizontalStackedBarSeries, public QQmlParserStatus
{
    Q_OBJECT
    Q_INTERFACES(QQmlParserStatus)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisX READ axisX WRITE setAxisX NOTIFY axisXChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisY READ axisY WRITE setAxisY NOTIFY axisYChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisXTop READ axisXTop WRITE setAxisXTop NOTIFY axisXTopChanged REVISION 2)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisYRight READ axisYRight WRITE setAxisYRight NOTIFY axisYRightChanged REVISION 2)
    Q_PROPERTY(QQmlListProperty<QObject> seriesChildren READ seriesChildren)
    Q_CLASSINFO("DefaultProperty", "seriesChildren")

public:
    explicit DeclarativeHorizontalStackedBarSeries(QQuickItem *parent = 0);
    QAbstractAxis *axisX() { return m_axes->axisX(); }
    void setAxisX(QAbstractAxis *axis) { m_axes->setAxisX(axis); }
    QAbstractAxis *axisY() { return m_axes->axisY(); }
    void setAxisY(QAbstractAxis *axis) { m_axes->setAxisY(axis); }
    QAbstractAxis *axisXTop() { return m_axes->axisXTop(); }
    void setAxisXTop(QAbstractAxis *axis) { m_axes->setAxisXTop(axis); }
    QAbstractAxis *axisYRight() { return m_axes->axisYRight(); }
    void setAxisYRight(QAbstractAxis *axis) { m_axes->setAxisYRight(axis); }
    QQmlListProperty<QObject> seriesChildren();

public:
    Q_INVOKABLE DeclarativeBarSet *at(int index);
    Q_INVOKABLE DeclarativeBarSet *append(QString label, QVariantList values) { return insert(count(), label, values); }
    Q_INVOKABLE DeclarativeBarSet *insert(int index, QString label, QVariantList values);
    Q_INVOKABLE bool remove(QBarSet *barset) { return QHorizontalStackedBarSeries::remove(barset); }
    Q_INVOKABLE void clear() { return QHorizontalStackedBarSeries::clear(); }

public: // from QDeclarativeParserStatus
    void classBegin();
    void componentComplete();

Q_SIGNALS:
    Q_REVISION(1) void axisXChanged(QAbstractAxis *axis);
    Q_REVISION(1) void axisYChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisXTopChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisYRightChanged(QAbstractAxis *axis);

public Q_SLOTS:
    static void appendSeriesChildren(QQmlListProperty<QObject> *list, QObject *element);

public:
    DeclarativeAxes *m_axes;
};

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeHorizontalPercentBarSeries : public QHorizontalPercentBarSeries, public QQmlParserStatus
{
    Q_OBJECT
    Q_INTERFACES(QQmlParserStatus)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisX READ axisX WRITE setAxisX NOTIFY axisXChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisY READ axisY WRITE setAxisY NOTIFY axisYChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisXTop READ axisXTop WRITE setAxisXTop NOTIFY axisXTopChanged REVISION 2)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisYRight READ axisYRight WRITE setAxisYRight NOTIFY axisYRightChanged REVISION 2)
    Q_PROPERTY(QQmlListProperty<QObject> seriesChildren READ seriesChildren)
    Q_CLASSINFO("DefaultProperty", "seriesChildren")

public:
    explicit DeclarativeHorizontalPercentBarSeries(QQuickItem *parent = 0);
    QAbstractAxis *axisX() { return m_axes->axisX(); }
    void setAxisX(QAbstractAxis *axis) { m_axes->setAxisX(axis); }
    QAbstractAxis *axisY() { return m_axes->axisY(); }
    void setAxisY(QAbstractAxis *axis) { m_axes->setAxisY(axis); }
    QAbstractAxis *axisXTop() { return m_axes->axisXTop(); }
    void setAxisXTop(QAbstractAxis *axis) { m_axes->setAxisXTop(axis); }
    QAbstractAxis *axisYRight() { return m_axes->axisYRight(); }
    void setAxisYRight(QAbstractAxis *axis) { m_axes->setAxisYRight(axis); }
    QQmlListProperty<QObject> seriesChildren();

public:
    Q_INVOKABLE DeclarativeBarSet *at(int index);
    Q_INVOKABLE DeclarativeBarSet *append(QString label, QVariantList values) { return insert(count(), label, values); }
    Q_INVOKABLE DeclarativeBarSet *insert(int index, QString label, QVariantList values);
    Q_INVOKABLE bool remove(QBarSet *barset) { return QHorizontalPercentBarSeries::remove(barset); }
    Q_INVOKABLE void clear() { return QHorizontalPercentBarSeries::clear(); }

public: // from QDeclarativeParserStatus
    void classBegin();
    void componentComplete();

Q_SIGNALS:
    Q_REVISION(1) void axisXChanged(QAbstractAxis *axis);
    Q_REVISION(1) void axisYChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisXTopChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisYRightChanged(QAbstractAxis *axis);

public Q_SLOTS:
    static void appendSeriesChildren(QQmlListProperty<QObject> *list, QObject *element);

public:
    DeclarativeAxes *m_axes;
};

QT_CHARTS_END_NAMESPACE

#endif // DECLARATIVEBARSERIES_H
