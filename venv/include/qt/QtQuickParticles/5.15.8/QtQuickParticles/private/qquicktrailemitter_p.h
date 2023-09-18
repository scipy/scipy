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

#ifndef FOLLOWEMITTER_H
#define FOLLOWEMITTER_H

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
#include "qquickparticleemitter_p.h"
#include "qquickparticleaffector_p.h"

QT_BEGIN_NAMESPACE

class QQuickTrailEmitter : public QQuickParticleEmitter
{
    Q_OBJECT
    Q_PROPERTY(QString follow READ follow WRITE setFollow NOTIFY followChanged)
    Q_PROPERTY(int emitRatePerParticle READ particlesPerParticlePerSecond WRITE setParticlesPerParticlePerSecond NOTIFY particlesPerParticlePerSecondChanged)

    Q_PROPERTY(QQuickParticleExtruder* emitShape READ emissonShape WRITE setEmissionShape NOTIFY emissionShapeChanged)
    Q_PROPERTY(qreal emitHeight READ emitterYVariation WRITE setEmitterYVariation NOTIFY emitterYVariationChanged)
    Q_PROPERTY(qreal emitWidth READ emitterXVariation WRITE setEmitterXVariation NOTIFY emitterXVariationChanged)
    QML_NAMED_ELEMENT(TrailEmitter)

public:
    enum EmitSize {
        ParticleSize = -2//Anything less than 0 will do
    };
    Q_ENUM(EmitSize)
    explicit QQuickTrailEmitter(QQuickItem *parent = 0);
    void emitWindow(int timeStamp) override;
    void reset() override;

    int particlesPerParticlePerSecond() const
    {
        return m_particlesPerParticlePerSecond;
    }

    qreal emitterXVariation() const
    {
        return m_emitterXVariation;
    }

    qreal emitterYVariation() const
    {
        return m_emitterYVariation;
    }

    QString follow() const
    {
        return m_follow;
    }

    QQuickParticleExtruder* emissonShape() const
    {
        return m_emissionExtruder;
    }

Q_SIGNALS:
    void emitFollowParticles(const QJSValue &particles, const QJSValue &followed);

    void particlesPerParticlePerSecondChanged(int arg);

    void emitterXVariationChanged(qreal arg);

    void emitterYVariationChanged(qreal arg);

    void followChanged(const QString &arg);

    void emissionShapeChanged(QQuickParticleExtruder* arg);

public Q_SLOTS:

    void setParticlesPerParticlePerSecond(int arg)
    {
        if (m_particlesPerParticlePerSecond != arg) {
            m_particlesPerParticlePerSecond = arg;
            Q_EMIT particlesPerParticlePerSecondChanged(arg);
        }
    }
    void setEmitterXVariation(qreal arg)
    {
        if (m_emitterXVariation != arg) {
            m_emitterXVariation = arg;
            Q_EMIT emitterXVariationChanged(arg);
        }
    }

    void setEmitterYVariation(qreal arg)
    {
        if (m_emitterYVariation != arg) {
            m_emitterYVariation = arg;
            Q_EMIT emitterYVariationChanged(arg);
        }
    }

    void setFollow(const QString &arg)
    {
        if (m_follow != arg) {
            m_follow = arg;
            Q_EMIT followChanged(arg);
        }
    }

    void setEmissionShape(QQuickParticleExtruder* arg)
    {
        if (m_emissionExtruder != arg) {
            m_emissionExtruder = arg;
            Q_EMIT emissionShapeChanged(arg);
        }
    }

private Q_SLOTS:
    void recalcParticlesPerSecond();

private:
    QSet<QQuickParticleData*> m_pending;
    QVector<qreal> m_lastEmission;
    int m_particlesPerParticlePerSecond;
    qreal m_lastTimeStamp;
    qreal m_emitterXVariation;
    qreal m_emitterYVariation;
    QString m_follow;
    int m_followCount;
    QQuickParticleExtruder* m_emissionExtruder;
    QQuickParticleExtruder* m_defaultEmissionExtruder;
    bool isEmitFollowConnected();
};

QT_END_NAMESPACE
#endif // FOLLOWEMITTER_H
