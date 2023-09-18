/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DANIMATION_QCLIPANIMATOR_H
#define QT3DANIMATION_QCLIPANIMATOR_H

#include <Qt3DAnimation/qt3danimation_global.h>
#include <Qt3DAnimation/qabstractclipanimator.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {

class QAbstractAnimationClip;
class QChannelMapper;
class QClipAnimatorPrivate;

class Q_3DANIMATIONSHARED_EXPORT QClipAnimator : public Qt3DAnimation::QAbstractClipAnimator
{
    Q_OBJECT
    Q_PROPERTY(Qt3DAnimation::QAbstractAnimationClip *clip READ clip WRITE setClip NOTIFY clipChanged)

public:
    explicit QClipAnimator(Qt3DCore::QNode *parent = nullptr);
    ~QClipAnimator();

    Qt3DAnimation::QAbstractAnimationClip *clip() const;

public Q_SLOTS:
    void setClip(Qt3DAnimation::QAbstractAnimationClip *clip);

Q_SIGNALS:
    void clipChanged(Qt3DAnimation::QAbstractAnimationClip *clip);

protected:
    QClipAnimator(QClipAnimatorPrivate &dd, Qt3DCore::QNode *parent = nullptr);

private:
    Q_DECLARE_PRIVATE(QClipAnimator)
    Qt3DCore::QNodeCreatedChangeBasePtr createNodeCreationChange() const override;
};

} // namespace Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_QCLIPANIMATOR_H
