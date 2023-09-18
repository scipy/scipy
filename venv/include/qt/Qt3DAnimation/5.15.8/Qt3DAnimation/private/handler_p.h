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

#ifndef QT3DANIMATION_ANIMATION_HANDLER_H
#define QT3DANIMATION_ANIMATION_HANDLER_H

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
#include <Qt3DCore/qaspectjob.h>
#include <Qt3DCore/qnodeid.h>
#include <QtCore/qscopedpointer.h>
#include <QtCore/qmutex.h>

QT_BEGIN_NAMESPACE

#if defined(QT_BUILD_INTERNAL)
class tst_Handler;
#endif

namespace Qt3DAnimation {
namespace Animation {

class AnimationClip;
class AnimationClipLoaderManager;
class ClockManager;
class ClipAnimator;
class ClipAnimatorManager;
class BlendedClipAnimator;
class BlendedClipAnimatorManager;
class ChannelMapping;
class ChannelMappingManager;
class ChannelMapper;
class ChannelMapperManager;
class ClipBlendNodeManager;
class SkeletonManager;

class FindRunningClipAnimatorsJob;
class LoadAnimationClipJob;
class EvaluateClipAnimatorJob;
class BuildBlendTreesJob;
class EvaluateBlendClipAnimatorJob;

using BuildBlendTreesJobPtr = QSharedPointer<BuildBlendTreesJob>;
using EvaluateBlendClipAnimatorJobPtr = QSharedPointer<EvaluateBlendClipAnimatorJob>;

class Q_AUTOTEST_EXPORT Handler
{
public:
    Handler();
    ~Handler();

    enum DirtyFlag {
        AnimationClipDirty,
        ChannelMappingsDirty,
        ClipAnimatorDirty,
        BlendedClipAnimatorDirty
    };

    qint64 simulationTime() const { return m_simulationTime; }

    void setDirty(DirtyFlag flag, Qt3DCore::QNodeId nodeId);

    void setClipAnimatorRunning(const HClipAnimator &handle, bool running);
    QVector<HClipAnimator> runningClipAnimators() const { return m_runningClipAnimators; }

    void setBlendedClipAnimatorRunning(const HBlendedClipAnimator &handle, bool running);
    QVector<HBlendedClipAnimator> runningBlenndedClipAnimators() const { return m_runningBlendedClipAnimators; }

    AnimationClipLoaderManager *animationClipLoaderManager() const Q_DECL_NOTHROW { return m_animationClipLoaderManager.data(); }
    ClockManager *clockManager() const Q_DECL_NOTHROW { return m_clockManager.data(); }
    ClipAnimatorManager *clipAnimatorManager() const Q_DECL_NOTHROW { return m_clipAnimatorManager.data(); }
    BlendedClipAnimatorManager *blendedClipAnimatorManager() const Q_DECL_NOTHROW { return m_blendedClipAnimatorManager.data(); }
    ChannelMappingManager *channelMappingManager() const Q_DECL_NOTHROW { return m_channelMappingManager.data(); }
    ChannelMapperManager *channelMapperManager() const Q_DECL_NOTHROW { return m_channelMapperManager.data(); }
    ClipBlendNodeManager *clipBlendNodeManager() const Q_DECL_NOTHROW { return m_clipBlendNodeManager.data(); }
    SkeletonManager *skeletonManager() const Q_DECL_NOTHROW { return m_skeletonManager.data(); }

    QVector<Qt3DCore::QAspectJobPtr> jobsToExecute(qint64 time);

    void cleanupHandleList(QVector<HAnimationClip> *clips);
    void cleanupHandleList(QVector<HClipAnimator> *animators);
    void cleanupHandleList(QVector<HBlendedClipAnimator> *animators);

private:
    QMutex m_mutex;
    QScopedPointer<AnimationClipLoaderManager> m_animationClipLoaderManager;
    QScopedPointer<ClockManager> m_clockManager;
    QScopedPointer<ClipAnimatorManager> m_clipAnimatorManager;
    QScopedPointer<BlendedClipAnimatorManager> m_blendedClipAnimatorManager;
    QScopedPointer<ChannelMappingManager> m_channelMappingManager;
    QScopedPointer<ChannelMapperManager> m_channelMapperManager;
    QScopedPointer<ClipBlendNodeManager> m_clipBlendNodeManager;
    QScopedPointer<SkeletonManager> m_skeletonManager;

    QVector<HAnimationClip> m_dirtyAnimationClips;
    QVector<HClipAnimator> m_dirtyClipAnimators;
    QVector<HBlendedClipAnimator> m_dirtyBlendedAnimators;

    QVector<HClipAnimator> m_runningClipAnimators;
    QVector<HBlendedClipAnimator> m_runningBlendedClipAnimators;

    QSharedPointer<LoadAnimationClipJob> m_loadAnimationClipJob;
    QSharedPointer<FindRunningClipAnimatorsJob> m_findRunningClipAnimatorsJob;
    QVector<QSharedPointer<EvaluateClipAnimatorJob>> m_evaluateClipAnimatorJobs;
    QVector<EvaluateBlendClipAnimatorJobPtr> m_evaluateBlendClipAnimatorJobs;
    BuildBlendTreesJobPtr m_buildBlendTreesJob;

    qint64 m_simulationTime;

#if defined(QT_BUILD_INTERNAL)
    friend class QT_PREPEND_NAMESPACE(tst_Handler);
#endif
};

} // namespace Animation
} // namespace Qt3DAnimation

QT_END_NAMESPACE

#endif // QT3DANIMATION_ANIMATION_HANDLER_H
