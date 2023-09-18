/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
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

#ifndef QQUICK3DQUATERNIONANIMATION_H
#define QQUICK3DQUATERNIONANIMATION_H

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

#include <QtGui/qquaternion.h>

#include <QtQuick3D/private/qquick3dobject_p.h>
#include <QtQuick/private/qquickanimation_p.h>

QT_BEGIN_NAMESPACE

class QQuick3DQuaternionAnimationPrivate;

class Q_QUICK3D_EXPORT QQuick3DQuaternionAnimation : public QQuickPropertyAnimation
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QQuick3DQuaternionAnimation)
    Q_PROPERTY(QQuaternion from READ from WRITE setFrom)
    Q_PROPERTY(QQuaternion to READ to WRITE setTo)
    Q_PROPERTY(Type type READ type WRITE setType NOTIFY typeChanged)

    Q_PROPERTY(float fromXRotation READ fromXRotation WRITE setFromXRotation NOTIFY fromXRotationChanged)
    Q_PROPERTY(float fromYRotation READ fromYRotation WRITE setFromYRotation NOTIFY fromYRotationChanged)
    Q_PROPERTY(float fromZRotation READ fromZRotation WRITE setFromZRotation NOTIFY fromZRotationChanged)
    Q_PROPERTY(float toXRotation READ toXRotation WRITE setToXRotation NOTIFY toXRotationChanged)
    Q_PROPERTY(float toYRotation READ toYRotation WRITE setToYRotation NOTIFY toYRotationChanged)
    Q_PROPERTY(float toZRotation READ toZRotation WRITE setToZRotation NOTIFY toZRotationChanged)

public:
    enum Type {
        Slerp = 0,
        Nlerp
    };
    Q_ENUM(Type)

    QQuick3DQuaternionAnimation(QObject *parent = nullptr);

    QQuaternion from() const;
    void setFrom(const QQuaternion &f);

    QQuaternion to() const;
    void setTo(const QQuaternion &t);

    Type type() const;
    void setType(Type type);

    float fromXRotation() const;
    void setFromXRotation(float f);

    float fromYRotation() const;
    void setFromYRotation(float f);

    float fromZRotation() const;
    void setFromZRotation(float f);

    float toXRotation() const;
    void setToXRotation(float f);

    float toYRotation() const;
    void setToYRotation(float f);

    float toZRotation() const;
    void setToZRotation(float f);

Q_SIGNALS:
    void typeChanged(Type type);
    void fromXRotationChanged(float value);
    void fromYRotationChanged(float value);
    void fromZRotationChanged(float value);
    void toXRotationChanged(float value);
    void toYRotationChanged(float value);
    void toZRotationChanged(float value);
};

Q_QUICK3D_EXPORT QVariant q_quaternionInterpolator(const QQuaternion &from, const QQuaternion &to, qreal progress);

QT_END_NAMESPACE

#endif // QQUICK3DQUATERNIONANIMATION_H
