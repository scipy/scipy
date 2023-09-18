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

#ifndef XYANIMATION_P_H
#define XYANIMATION_P_H

#include <private/chartanimation_p.h>
#include <QtCharts/private/qchartglobal_p.h>
#include <QtCore/QPointF>

QT_CHARTS_BEGIN_NAMESPACE

class XYChart;

class Q_CHARTS_PRIVATE_EXPORT XYAnimation : public ChartAnimation
{
protected:
    enum Animation { AddPointAnimation, RemovePointAnimation, ReplacePointAnimation, NewAnimation };
public:
    XYAnimation(XYChart *item, int duration, QEasingCurve &curve);
    ~XYAnimation();
    void setup(const QVector<QPointF> &oldPoints, const QVector<QPointF> &newPoints, int index = -1);
    Animation animationType() const { return m_type; };

protected:
    QVariant interpolated(const QVariant &start, const QVariant &end, qreal progress) const;
    void updateCurrentValue(const QVariant &value);
    void updateState(QAbstractAnimation::State newState, QAbstractAnimation::State oldState);
    XYChart *chartItem() { return m_item; }
protected:
    Animation m_type;
    bool m_dirty;
    int m_index;
private:
    XYChart *m_item;
    QVector<QPointF> m_oldPoints;
    QVector<QPointF> m_newPoints;
};

QT_CHARTS_END_NAMESPACE

#endif
