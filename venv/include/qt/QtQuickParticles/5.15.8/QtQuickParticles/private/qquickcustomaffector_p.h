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

#ifndef CUSTOMAFFECTOR_H
#define CUSTOMAFFECTOR_H

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

#include <QObject>
#include <QtQml/qqml.h>
#include "qquickparticlesystem_p.h"
#include "qquickparticleextruder_p.h"
#include "qquickparticleaffector_p.h"
#include "qquickdirection_p.h"

QT_BEGIN_NAMESPACE

class QQuickCustomAffector : public QQuickParticleAffector
{
    Q_OBJECT
    Q_PROPERTY(bool relative READ relative WRITE setRelative NOTIFY relativeChanged)
    Q_PROPERTY(QQuickDirection *position READ position WRITE setPosition NOTIFY positionChanged RESET positionReset)
    Q_PROPERTY(QQuickDirection *velocity READ velocity WRITE setVelocity NOTIFY velocityChanged RESET velocityReset)
    Q_PROPERTY(QQuickDirection *acceleration READ acceleration WRITE setAcceleration NOTIFY accelerationChanged RESET accelerationReset)
    QML_NAMED_ELEMENT(Affector)

public:
    explicit QQuickCustomAffector(QQuickItem *parent = 0);
    void affectSystem(qreal dt) override;

    QQuickDirection * position() const
    {
        return m_position;
    }

    QQuickDirection * velocity() const
    {
        return m_velocity;
    }

    QQuickDirection * acceleration() const
    {
        return m_acceleration;
    }

    void positionReset()
    {
        m_position = &m_nullVector;
    }

    void velocityReset()
    {
        m_velocity = &m_nullVector;
    }

    void accelerationReset()
    {
        m_acceleration = &m_nullVector;
    }

    bool relative() const
    {
        return m_relative;
    }


Q_SIGNALS:
    void affectParticles(const QJSValue &particles, qreal dt);

    void positionChanged(QQuickDirection * arg);

    void velocityChanged(QQuickDirection * arg);

    void accelerationChanged(QQuickDirection * arg);

    void relativeChanged(bool arg);

public Q_SLOTS:
    void setPosition(QQuickDirection * arg)
    {
        if (m_position != arg) {
            m_position = arg;
            Q_EMIT positionChanged(arg);
        }
    }

    void setVelocity(QQuickDirection * arg)
    {
        if (m_velocity != arg) {
            m_velocity = arg;
            Q_EMIT velocityChanged(arg);
        }
    }

    void setAcceleration(QQuickDirection * arg)
    {
        if (m_acceleration != arg) {
            m_acceleration = arg;
            Q_EMIT accelerationChanged(arg);
        }
    }

    void setRelative(bool arg)
    {
        if (m_relative != arg) {
            m_relative = arg;
            Q_EMIT relativeChanged(arg);
        }
    }

protected:
    bool isAffectConnected();
    bool affectParticle(QQuickParticleData *d, qreal dt) override;

private:
    void affectProperties(const QList<QQuickParticleData*> &particles, qreal dt);
    QQuickDirection * m_position;
    QQuickDirection * m_velocity;
    QQuickDirection * m_acceleration;

    QQuickDirection m_nullVector;
    bool m_relative;
};

QT_END_NAMESPACE
#endif // CUSTOMAFFECTOR_H
