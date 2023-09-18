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

#ifndef DECLARATIVEAREASERIES_H
#define DECLARATIVEAREASERIES_H

#include <QtCharts/QAreaSeries>
#include <private/declarativechartglobal_p.h>
#include <private/declarativeaxes_p.h>

QT_CHARTS_BEGIN_NAMESPACE
class DeclarativeLineSeries;

class Q_QMLCHARTS_PRIVATE_EXPORT DeclarativeAreaSeries : public QAreaSeries
{
    Q_OBJECT
    Q_PROPERTY(QtCharts::DeclarativeLineSeries *upperSeries READ upperSeries WRITE setUpperSeries)
    Q_PROPERTY(QtCharts::DeclarativeLineSeries *lowerSeries READ lowerSeries WRITE setLowerSeries)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisX READ axisX WRITE setAxisX NOTIFY axisXChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisY READ axisY WRITE setAxisY NOTIFY axisYChanged REVISION 1)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisXTop READ axisXTop WRITE setAxisXTop NOTIFY axisXTopChanged REVISION 2)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisYRight READ axisYRight WRITE setAxisYRight NOTIFY axisYRightChanged REVISION 2)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisAngular READ axisAngular WRITE setAxisAngular NOTIFY axisAngularChanged REVISION 3)
    Q_PROPERTY(QtCharts::QAbstractAxis *axisRadial READ axisRadial WRITE setAxisRadial NOTIFY axisRadialChanged REVISION 3)
    Q_PROPERTY(qreal borderWidth READ borderWidth WRITE setBorderWidth NOTIFY borderWidthChanged REVISION 1)
    Q_PROPERTY(QString brushFilename READ brushFilename WRITE setBrushFilename NOTIFY brushFilenameChanged REVISION 4)
    Q_PROPERTY(QBrush brush READ brush WRITE setBrush NOTIFY brushChanged REVISION 4)

public:
    explicit DeclarativeAreaSeries(QObject *parent = 0);
    void setUpperSeries(DeclarativeLineSeries *series);
    DeclarativeLineSeries *upperSeries() const;
    void setLowerSeries(DeclarativeLineSeries *series);
    DeclarativeLineSeries *lowerSeries() const;
    QAbstractAxis *axisX() { return m_axes->axisX(); }
    void setAxisX(QAbstractAxis *axis) { m_axes->setAxisX(axis); }
    QAbstractAxis *axisY() { return m_axes->axisY(); }
    void setAxisY(QAbstractAxis *axis) { m_axes->setAxisY(axis); }
    QAbstractAxis *axisXTop() { return m_axes->axisXTop(); }
    void setAxisXTop(QAbstractAxis *axis) { m_axes->setAxisXTop(axis); }
    QAbstractAxis *axisYRight() { return m_axes->axisYRight(); }
    void setAxisYRight(QAbstractAxis *axis) { m_axes->setAxisYRight(axis); }
    QAbstractAxis *axisAngular() { return m_axes->axisX(); }
    void setAxisAngular(QAbstractAxis *axis) { m_axes->setAxisX(axis); }
    QAbstractAxis *axisRadial() { return m_axes->axisY(); }
    void setAxisRadial(QAbstractAxis *axis) { m_axes->setAxisY(axis); }
    qreal borderWidth() const;
    void setBorderWidth(qreal borderWidth);
    QString brushFilename() const;
    void setBrushFilename(const QString &brushFilename);
    void setBrush(const QBrush &brush);
    QBrush brush() const;

Q_SIGNALS:
    Q_REVISION(1) void axisXChanged(QAbstractAxis *axis);
    Q_REVISION(1) void axisYChanged(QAbstractAxis *axis);
    Q_REVISION(1) void borderWidthChanged(qreal width);
    Q_REVISION(2) void axisXTopChanged(QAbstractAxis *axis);
    Q_REVISION(2) void axisYRightChanged(QAbstractAxis *axis);
    Q_REVISION(3) void axisAngularChanged(QAbstractAxis *axis);
    Q_REVISION(3) void axisRadialChanged(QAbstractAxis *axis);
    Q_REVISION(4) void brushChanged();
    Q_REVISION(4) void brushFilenameChanged(const QString &brushFilename);

private Q_SLOTS:
    void handleBrushChanged();

public:
    DeclarativeAxes *m_axes;

private:
    QString m_brushFilename;
    QImage m_brushImage;
};

QT_CHARTS_END_NAMESPACE

#endif // DECLARATIVEAREASERIES_H
