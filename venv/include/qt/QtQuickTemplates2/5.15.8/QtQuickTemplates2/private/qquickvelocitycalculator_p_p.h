/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Templates 2 module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QQUICKVELOCITYCALCULATOR_P_P_H
#define QQUICKVELOCITYCALCULATOR_P_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtCore/qpoint.h>
#include <QtCore/qelapsedtimer.h>

QT_BEGIN_NAMESPACE

class QQuickVelocityCalculator
{
public:
    void startMeasuring(const QPointF &point1, qint64 timestamp = 0);
    void stopMeasuring(const QPointF &m_point2, qint64 timestamp = 0);
    void reset();
    QPointF velocity() const;

private:
    QPointF m_point1;
    QPointF m_point2;
    qint64 m_point1Timestamp = 0;
    qint64 m_point2Timestamp = 0;
    // When a timestamp isn't available, we must use a timer.
    // When stopMeasuring() has been called, we store the elapsed time in point2timestamp.
    QElapsedTimer m_timer;
};

QT_END_NAMESPACE

#endif // QQUICKVELOCITYCALCULATOR_P_P_H
