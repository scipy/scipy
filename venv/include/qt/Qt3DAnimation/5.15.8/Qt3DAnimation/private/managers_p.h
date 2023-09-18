/****************************************************************************
**
** Copyright (C) 2017 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DANIMATION_ANIMATION_MANAGERS_P_H
#define QT3DANIMATION_ANIMATION_MANAGERS_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtGlobal>
#include <Qt3DAnimation/private/handle_types_p.h>
#include <Qt3DAnimation/private/animationclip_p.h>
#include <Qt3DAnimation/private/clock_p.h>
#include <Qt3DAnimation/private/blendedclipanimator_p.h>
#include <Qt3DAnimation/private/clipanimator_p.h>
#include <Qt3DAnimation/private/channelmapping_p.h>
#include <Qt3DAnimation/private/channelmapper_p.h>
#include <Qt3DAnimation/private/skeleton_p.h>
#include <Qt3DCore/private/qresourcemanager_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DAnimation {
namespace Animation {

class ClipBlendNode;

class AnimationClipLoaderManager : public Qt3DCore::QResourceManager<
        AnimationClip,
        Qt3DCore::QNodeId>
{
public:
    AnimationClipLoaderManager() {}
};

class ClockManager : public Qt3DCore::QResourceManager<
    Clock,
    Qt3DCore::QNodeId>
{
public:
    ClockManager() {}
};

class ClipAnimatorManager : public Qt3DCore::QResourceManager<
        ClipAnimator,
        Qt3DCore::QNodeId>
{
public:
    ClipAnimatorManager() {}
};

class BlendedClipAnimatorManager : public Qt3DCore::QResourceManager<
        BlendedClipAnimator,
        Qt3DCore::QNodeId>
{
public:
    BlendedClipAnimatorManager() {}
};

class ChannelMappingManager : public Qt3DCore::QResourceManager<
        ChannelMapping,
        Qt3DCore::QNodeId>
{
public:
    ChannelMappingManager() {}
};

class ChannelMapperManager : public Qt3DCore::QResourceManager<
        ChannelMapper,
        Qt3DCore::QNodeId>
{
public:
    ChannelMapperManager() {}
};

class Q_AUTOTEST_EXPORT ClipBlendNodeManager
{
public:
    ClipBlendNodeManager();
    ~ClipBlendNodeManager();

    bool containsNode(Qt3DCore::QNodeId id) const;
    void appendNode(Qt3DCore::QNodeId id, ClipBlendNode *node);
    ClipBlendNode *lookupNode(Qt3DCore::QNodeId id) const;
    void releaseNode(Qt3DCore::QNodeId id);

private:
    QHash<Qt3DCore::QNodeId, ClipBlendNode *> m_nodes;
};

class SkeletonManager : public Qt3DCore::QResourceManager<
        Skeleton,
        Qt3DCore::QNodeId>
{
public:
    SkeletonManager() {}
};

} // namespace Animation
} // namespace Qt3DAnimation

Q_DECLARE_RESOURCE_INFO(Qt3DAnimation::Animation::AnimationClip, Q_REQUIRES_CLEANUP)
Q_DECLARE_RESOURCE_INFO(Qt3DAnimation::Animation::ClipAnimator, Q_REQUIRES_CLEANUP)
Q_DECLARE_RESOURCE_INFO(Qt3DAnimation::Animation::BlendedClipAnimator, Q_REQUIRES_CLEANUP)
Q_DECLARE_RESOURCE_INFO(Qt3DAnimation::Animation::ChannelMapping, Q_REQUIRES_CLEANUP)
Q_DECLARE_RESOURCE_INFO(Qt3DAnimation::Animation::ChannelMapper, Q_REQUIRES_CLEANUP)
Q_DECLARE_RESOURCE_INFO(Qt3DAnimation::Animation::Skeleton, Q_REQUIRES_CLEANUP)

QT_END_NAMESPACE

#endif // QT3DANIMATION_ANIMATION_MANAGERS_P_H
