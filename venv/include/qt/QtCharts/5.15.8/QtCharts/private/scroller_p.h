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

#ifndef SCROLLER_P_H
#define SCROLLER_P_H

#include <QtCharts/QChartGlobal>
#include <QtCore/QObject>
#include <QtCore/QBasicTimer>
#include <QtCore/QElapsedTimer>
#include <QtCore/QPointF>
#include <QtCharts/private/qchartglobal_p.h>

QT_BEGIN_NAMESPACE
class QGraphicsSceneMouseEvent;
QT_END_NAMESPACE

QT_CHARTS_BEGIN_NAMESPACE

class Scroller;
class QLegend;

class Q_CHARTS_PRIVATE_EXPORT ScrollTicker : public QObject
{
    Q_OBJECT
public:
    explicit ScrollTicker(Scroller *scroller, QObject *parent = 0);
    void start(int interval);
    void stop();
protected:
    void timerEvent(QTimerEvent *event);

private:
    QBasicTimer m_timer;
    Scroller *m_scroller;
};

class Q_CHARTS_PRIVATE_EXPORT Scroller
{
public:
    enum State {
        Idle,
        Pressed,
        Move,
        Scroll
    };

    Scroller();
    virtual ~Scroller();

    virtual void setOffset(const QPointF &point) = 0;
    virtual QPointF offset() const = 0;

    void move(const QPointF &delta);
    void scrollTo(const QPointF &delta);

    void handleMousePressEvent(QGraphicsSceneMouseEvent *event);
    void handleMouseMoveEvent(QGraphicsSceneMouseEvent *event);
    void handleMouseReleaseEvent(QGraphicsSceneMouseEvent *event);

    void scrollTick();

private:
    void startTicker(int interval);
    void stopTicker();

private:
    void calculateSpeed(const QPointF &position);
    void lowerSpeed(QPointF &speed, qreal maxSpeed = 100);

private:
    ScrollTicker m_ticker;
    QElapsedTimer m_timeStamp;
    QPointF m_speed;
    QPointF m_fraction;
    int m_timeTresholdMin;
    int m_timeTresholdMax;

    State m_state;
    QPointF m_pressPos;
    QPointF m_lastPos;
    qreal m_treshold;
};

QT_CHARTS_END_NAMESPACE

#endif /* SCROLLER_P_H */
