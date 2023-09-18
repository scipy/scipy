/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQuick module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QQuickANGLEDDIRECTION_H
#define QQuickANGLEDDIRECTION_H

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
#include "qquickdirection_p.h"
#include <QtQml/qqml.h>

QT_BEGIN_NAMESPACE

class QQuickAngleDirection : public QQuickDirection
{
    Q_OBJECT
    Q_PROPERTY(qreal angle READ angle WRITE setAngle NOTIFY angleChanged)
    Q_PROPERTY(qreal magnitude READ magnitude WRITE setMagnitude NOTIFY magnitudeChanged)
    Q_PROPERTY(qreal angleVariation READ angleVariation WRITE setAngleVariation NOTIFY angleVariationChanged)
    Q_PROPERTY(qreal magnitudeVariation READ magnitudeVariation WRITE setMagnitudeVariation NOTIFY magnitudeVariationChanged)
    QML_NAMED_ELEMENT(AngleDirection)
public:
    explicit QQuickAngleDirection(QObject *parent = 0);
    QPointF sample(const QPointF &from) override;
    qreal angle() const
    {
        return m_angle;
    }

    qreal magnitude() const
    {
        return m_magnitude;
    }

    qreal angleVariation() const
    {
        return m_angleVariation;
    }

    qreal magnitudeVariation() const
    {
        return m_magnitudeVariation;
    }

Q_SIGNALS:

    void angleChanged(qreal arg);

    void magnitudeChanged(qreal arg);

    void angleVariationChanged(qreal arg);

    void magnitudeVariationChanged(qreal arg);

public Q_SLOTS:
void setAngle(qreal arg)
{
    if (m_angle != arg) {
        m_angle = arg;
        Q_EMIT angleChanged(arg);
    }
}

void setMagnitude(qreal arg)
{
    if (m_magnitude != arg) {
        m_magnitude = arg;
        Q_EMIT magnitudeChanged(arg);
    }
}

void setAngleVariation(qreal arg)
{
    if (m_angleVariation != arg) {
        m_angleVariation = arg;
        Q_EMIT angleVariationChanged(arg);
    }
}

void setMagnitudeVariation(qreal arg)
{
    if (m_magnitudeVariation != arg) {
        m_magnitudeVariation = arg;
        Q_EMIT magnitudeVariationChanged(arg);
    }
}

private:
qreal m_angle;
qreal m_magnitude;
qreal m_angleVariation;
qreal m_magnitudeVariation;
};

QT_END_NAMESPACE
#endif // QQuickANGLEDDIRECTION_H
