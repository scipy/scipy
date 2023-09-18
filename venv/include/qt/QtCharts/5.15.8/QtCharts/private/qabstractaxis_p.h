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

#ifndef QABSTRACTAXIS_P_H
#define QABSTRACTAXIS_P_H

#include <QtCharts/QAbstractAxis>
#include <private/chartaxiselement_p.h>
#include <QtCharts/QChart>
#include <QtCharts/private/qchartglobal_p.h>
#include <QtCore/QDebug>

QT_BEGIN_NAMESPACE
class QGraphicsItem;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE

class ChartPresenter;
class AbstractDomain;
class QChart;
class QAbstractSeries;
class ChartTheme;
class ChartElement;

class Q_CHARTS_PRIVATE_EXPORT QAbstractAxisPrivate : public QObject
{
    Q_OBJECT
public:
    QAbstractAxisPrivate(QAbstractAxis *q);
    ~QAbstractAxisPrivate();

public:
    Qt::Alignment alignment() const { return m_alignment; }
    Qt::Orientation orientation() const { return m_orientation; }
    void setAlignment( Qt::Alignment alignment);

    virtual void initializeDomain(AbstractDomain *domain) = 0;
    virtual void initializeGraphics(QGraphicsItem *parent) = 0;
    virtual void initializeTheme(ChartTheme* theme, bool forced = false);
    virtual void initializeAnimations(QChart::AnimationOptions options, int duration,
                                      QEasingCurve &curve);

    //interface for manipulating range form base class
    virtual void setMin(const QVariant &min) = 0;
    virtual void setMax(const QVariant &max) = 0;
    virtual void setRange(const QVariant &min, const QVariant &max) = 0;

    //interface manipulating range form domain
    virtual void setRange(qreal min, qreal max) = 0;
    virtual qreal min() = 0;
    virtual qreal max() = 0;

    ChartAxisElement *axisItem() { return m_item.data(); }

public Q_SLOTS:
    void handleRangeChanged(qreal min, qreal max);

Q_SIGNALS:
    void rangeChanged(qreal min, qreal max);

protected:
    QAbstractAxis *q_ptr;
    QChart *m_chart = nullptr;
    QScopedPointer<ChartAxisElement> m_item;

private:
    QList<QAbstractSeries*> m_series;

    Qt::Alignment m_alignment;
    Qt::Orientation m_orientation = Qt::Orientation(0);

    bool m_visible = true;

    bool m_arrowVisible = true;
    QPen m_axisPen;
    QBrush m_axisBrush;

    bool m_gridLineVisible = true;
    QPen m_gridLinePen;
    bool m_minorGridLineVisible = true;
    QPen m_minorGridLinePen;

    bool m_labelsVisible = true;
    bool m_labelsEditable = false;
    QBrush m_labelsBrush;
    QFont m_labelsFont;
    int m_labelsAngle = 0;

    bool m_titleVisible = true;
    QBrush m_titleBrush;
    QFont m_titleFont;
    QString m_title;

    bool m_shadesVisible = false;
    QPen m_shadesPen;
    QBrush m_shadesBrush;
    qreal m_shadesOpacity = 1;

    bool m_dirty = false;

    bool m_reverse = false;

    friend class QAbstractAxis;
    friend class ChartDataSet;
    friend class ChartPresenter;
};

QT_CHARTS_END_NAMESPACE

#endif
