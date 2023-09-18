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

#ifndef ATTRACTORAFFECTOR_H
#define ATTRACTORAFFECTOR_H

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
#include "qquickparticleaffector_p.h"

QT_BEGIN_NAMESPACE

class QQuickAttractorAffector : public QQuickParticleAffector
{
    Q_OBJECT
    Q_PROPERTY(qreal strength READ strength WRITE setStrength NOTIFY strengthChanged)
    Q_PROPERTY(qreal pointX READ pointX WRITE setPointX NOTIFY pointXChanged)
    Q_PROPERTY(qreal pointY READ pointY WRITE setPointY NOTIFY pointYChanged)
    Q_PROPERTY(AffectableParameters affectedParameter READ affectedParameter WRITE setAffectedParameter NOTIFY affectedParameterChanged)
    Q_PROPERTY(Proportion proportionalToDistance READ proportionalToDistance WRITE setProportionalToDistance NOTIFY proportionalToDistanceChanged)
    QML_NAMED_ELEMENT(Attractor)

public:
    enum Proportion{
        Constant,
        Linear,
        Quadratic,
        InverseLinear,
        InverseQuadratic
    };
    Q_ENUM(Proportion)

    enum AffectableParameters {
        Position,
        Velocity,
        Acceleration
    };
    Q_ENUM(AffectableParameters)

    explicit QQuickAttractorAffector(QQuickItem *parent = 0);

    qreal strength() const
    {
        return m_strength;
    }

    qreal pointX() const
    {
        return m_x;
    }

    qreal pointY() const
    {
        return m_y;
    }

    AffectableParameters affectedParameter() const
    {
        return m_physics;
    }

    Proportion proportionalToDistance() const
    {
        return m_proportionalToDistance;
    }

Q_SIGNALS:

    void strengthChanged(qreal arg);

    void pointXChanged(qreal arg);

    void pointYChanged(qreal arg);

    void affectedParameterChanged(AffectableParameters arg);

    void proportionalToDistanceChanged(Proportion arg);

public Q_SLOTS:
void setStrength(qreal arg)
{
    if (m_strength != arg) {
        m_strength = arg;
        Q_EMIT strengthChanged(arg);
    }
}

void setPointX(qreal arg)
{
    if (m_x != arg) {
        m_x = arg;
        Q_EMIT pointXChanged(arg);
    }
}

void setPointY(qreal arg)
{
    if (m_y != arg) {
        m_y = arg;
        Q_EMIT pointYChanged(arg);
    }
}
void setAffectedParameter(AffectableParameters arg)
{
    if (m_physics != arg) {
        m_physics = arg;
        Q_EMIT affectedParameterChanged(arg);
    }
}

void setProportionalToDistance(Proportion arg)
{
    if (m_proportionalToDistance != arg) {
        m_proportionalToDistance = arg;
        Q_EMIT proportionalToDistanceChanged(arg);
    }
}

protected:
    bool affectParticle(QQuickParticleData *d, qreal dt) override;

private:
qreal m_strength;
qreal m_x;
qreal m_y;
AffectableParameters m_physics;
Proportion m_proportionalToDistance;
};

QT_END_NAMESPACE
#endif // ATTRACTORAFFECTOR_H
