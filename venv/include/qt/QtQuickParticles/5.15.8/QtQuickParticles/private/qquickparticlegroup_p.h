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
#ifndef QQuickPARTICLEGROUP
#define QQuickPARTICLEGROUP

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
#include <private/qquickspriteengine_p.h>
#include "qquickparticlesystem_p.h"
#include "qqmlparserstatus.h"

QT_BEGIN_NAMESPACE

class QQuickParticleGroup : public QQuickStochasticState, public QQmlParserStatus
{
    Q_OBJECT
    //### Would setting limits per group be useful? Or clutter the API?
    //Q_PROPERTY(int maximumAlive READ maximumAlive WRITE setMaximumAlive NOTIFY maximumAliveChanged)

    Q_PROPERTY(QQuickParticleSystem* system READ system WRITE setSystem NOTIFY systemChanged)

    //Intercept children requests and assign to the group & system
    Q_PROPERTY(QQmlListProperty<QObject> particleChildren READ particleChildren DESIGNABLE false)//### Hidden property for in-state system definitions - ought not to be used in actual "Sprite" states
    Q_CLASSINFO("DefaultProperty", "particleChildren")
    QML_NAMED_ELEMENT(ParticleGroup)
    Q_INTERFACES(QQmlParserStatus)

public:
    explicit QQuickParticleGroup(QObject* parent = 0);

    QQmlListProperty<QObject> particleChildren();

    int maximumAlive() const
    {
        return m_maximumAlive;
    }

    QQuickParticleSystem* system() const
    {
        return m_system;
    }

public Q_SLOTS:

    void setMaximumAlive(int arg)
    {
        if (m_maximumAlive != arg) {
            m_maximumAlive = arg;
            Q_EMIT maximumAliveChanged(arg);
        }
    }

    void setSystem(QQuickParticleSystem* arg);

    void delayRedirect(QObject* obj);

Q_SIGNALS:

    void maximumAliveChanged(int arg);

    void systemChanged(QQuickParticleSystem* arg);

protected:
    void componentComplete() override;
    void classBegin() override {}

private:

    void performDelayedRedirects();

    int m_maximumAlive;
    QQuickParticleSystem* m_system;
    QList<QObject*> m_delayedRedirects;
};

QT_END_NAMESPACE

#endif
