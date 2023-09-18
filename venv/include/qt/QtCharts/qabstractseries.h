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

#ifndef QABSTRACTSERIES_H
#define QABSTRACTSERIES_H

#include <QtCharts/QChartGlobal>
#include <QtCharts/QAbstractAxis>
#include <QtCore/QObject>
#include <QtGui/QPen>

QT_CHARTS_BEGIN_NAMESPACE

class QAbstractSeriesPrivate;
class QChart;

class Q_CHARTS_EXPORT QAbstractSeries : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QString name READ name WRITE setName NOTIFY nameChanged)
    Q_PROPERTY(bool visible READ isVisible WRITE setVisible NOTIFY visibleChanged)
    Q_PROPERTY(qreal opacity READ opacity WRITE setOpacity NOTIFY opacityChanged)
    Q_PROPERTY(SeriesType type READ type)
    Q_PROPERTY(bool useOpenGL READ useOpenGL WRITE setUseOpenGL NOTIFY useOpenGLChanged)
    Q_ENUMS(SeriesType)

public:
    enum SeriesType {
        SeriesTypeLine,
        SeriesTypeArea,
        SeriesTypeBar,
        SeriesTypeStackedBar,
        SeriesTypePercentBar,
        SeriesTypePie,
        SeriesTypeScatter,
        SeriesTypeSpline,
        SeriesTypeHorizontalBar,
        SeriesTypeHorizontalStackedBar,
        SeriesTypeHorizontalPercentBar,
        SeriesTypeBoxPlot,
        SeriesTypeCandlestick
    };

protected:
    QAbstractSeries(QAbstractSeriesPrivate &d, QObject *parent = nullptr);

public:
    ~QAbstractSeries();
    virtual SeriesType type() const = 0;

    void setName(const QString &name);
    QString name() const;
    void setVisible(bool visible = true);
    bool isVisible() const;
    qreal opacity() const;
    void setOpacity(qreal opacity);
    void setUseOpenGL(bool enable = true);
    bool useOpenGL() const;

    QChart *chart() const;

    bool attachAxis(QAbstractAxis *axis);
    bool detachAxis(QAbstractAxis *axis);
    QList<QAbstractAxis*> attachedAxes();

    void show();
    void hide();

Q_SIGNALS:
    void nameChanged();
    void visibleChanged();
    void opacityChanged();
    void useOpenGLChanged();

protected:
    QScopedPointer<QAbstractSeriesPrivate> d_ptr;
    friend class ChartDataSet;
    friend class ChartPresenter;
    friend class ChartThemeManager;
    friend class QLegendPrivate;
    friend class DeclarativeChart;
    friend class QAreaSeries;
    friend class GLWidget;
};

QT_CHARTS_END_NAMESPACE

#endif // QABSTRACTSERIES_H
