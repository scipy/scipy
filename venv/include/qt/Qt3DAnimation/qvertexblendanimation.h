/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DANIMATION_QVERTEXBLENDANIMATION_H
#define QT3DANIMATION_QVERTEXBLENDANIMATION_H

#include <Qt3DRender/qgeometryrenderer.h>
#include <Qt3DAnimation/qabstractanimation.h>
#include <Qt3DAnimation/qmorphtarget.h>

#include <Qt3DAnimation/qt3danimation_global.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

class QVertexBlendAnimationPrivate;

class Q_3DANIMATIONSHARED_EXPORT QVertexBlendAnimation : public QAbstractAnimation
{
    Q_OBJECT
    Q_PROPERTY(QVector<float> targetPositions READ targetPositions WRITE setTargetPositions NOTIFY targetPositionsChanged)
    Q_PROPERTY(float interpolator READ interpolator NOTIFY interpolatorChanged)
    Q_PROPERTY(Qt3DRender::QGeometryRenderer *target READ target WRITE setTarget NOTIFY targetChanged)
    Q_PROPERTY(QString targetName READ targetName WRITE setTargetName NOTIFY targetNameChanged)

public:
    explicit QVertexBlendAnimation(QObject *parent = nullptr);

    QVector<float> targetPositions() const;
    float interpolator() const;
    Qt3DRender::QGeometryRenderer *target() const;
    QString targetName() const;

    void setMorphTargets(const QVector<Qt3DAnimation::QMorphTarget *> &targets);
    void addMorphTarget(Qt3DAnimation::QMorphTarget *target);
    void removeMorphTarget(Qt3DAnimation::QMorphTarget *target);

    QVector<Qt3DAnimation::QMorphTarget *> morphTargetList();

public Q_SLOTS:
    void setTargetPositions(const QVector<float> &targetPositions);
    void setTarget(Qt3DRender::QGeometryRenderer *target);
    void setTargetName(const QString name);

Q_SIGNALS:
    void targetPositionsChanged(const QVector<float> &targetPositions);
    void interpolatorChanged(float interpolator);
    void targetChanged(Qt3DRender::QGeometryRenderer *target);
    void targetNameChanged(const QString &name);

private:

    void updateAnimation(float position);

    Q_DECLARE_PRIVATE(QVertexBlendAnimation)
};

} // Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_QVERTEXBLENDANIMATION_H
