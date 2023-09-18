/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
** Contact: http://www.qt-project.org/legal
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

#ifndef QT3DANIMATION_QADDITIVECLIPBLEND_H
#define QT3DANIMATION_QADDITIVECLIPBLEND_H

#include <Qt3DAnimation/qabstractclipblendnode.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

class QAdditiveClipBlendPrivate;

class Q_3DANIMATIONSHARED_EXPORT QAdditiveClipBlend : public QAbstractClipBlendNode
{
    Q_OBJECT
    Q_PROPERTY(Qt3DAnimation::QAbstractClipBlendNode *baseClip READ baseClip WRITE setBaseClip NOTIFY baseClipChanged)
    Q_PROPERTY(Qt3DAnimation::QAbstractClipBlendNode *additiveClip READ additiveClip WRITE setAdditiveClip NOTIFY additiveClipChanged)
    Q_PROPERTY(float additiveFactor READ additiveFactor WRITE setAdditiveFactor NOTIFY additiveFactorChanged)
public:
    explicit QAdditiveClipBlend(Qt3DCore::QNode *parent = nullptr);
    ~QAdditiveClipBlend();

    float additiveFactor() const;
    Qt3DAnimation::QAbstractClipBlendNode *baseClip() const;
    Qt3DAnimation::QAbstractClipBlendNode *additiveClip() const;

public Q_SLOTS:
    void setAdditiveFactor(float additiveFactor);
    void setBaseClip(Qt3DAnimation::QAbstractClipBlendNode *baseClip);
    void setAdditiveClip(Qt3DAnimation::QAbstractClipBlendNode *additiveClip);

Q_SIGNALS:
    void additiveFactorChanged(float additiveFactor);
    void baseClipChanged(Qt3DAnimation::QAbstractClipBlendNode *baseClip);
    void additiveClipChanged(Qt3DAnimation::QAbstractClipBlendNode *additiveClip);

protected:
    explicit QAdditiveClipBlend(QAdditiveClipBlendPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QAdditiveClipBlend)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_QADDITIVECLIPBLEND_H
