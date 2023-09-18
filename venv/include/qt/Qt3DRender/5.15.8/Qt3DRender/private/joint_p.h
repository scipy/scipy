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

#ifndef QT3DRENDER_RENDER_JOINT_H
#define QT3DRENDER_RENDER_JOINT_H

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

#include <Qt3DRender/private/backendnode_p.h>
#include <Qt3DCore/private/sqt_p.h>
#include <Qt3DRender/private/handle_types_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
namespace Render {

class JointManager;
class SkeletonManager;

class Q_3DRENDERSHARED_PRIVATE_EXPORT Joint : public BackendNode
{
public:
    Joint();

    void cleanup();
    void syncFromFrontEnd(const Qt3DCore::QNode *frontEnd, bool firstTime) override;

    Qt3DCore::Sqt localPose() const { return m_localPose; }
    QMatrix4x4 inverseBindMatrix() const { return m_inverseBindMatrix; }
    QString name() const { return m_name; }
    QVector<Qt3DCore::QNodeId> childJointIds() const { return m_childJointIds; }

    QVector3D translation() const { return m_localPose.translation; }
    QQuaternion rotation() const { return m_localPose.rotation; }
    QVector3D scale() const { return m_localPose.scale; }

    void setOwningSkeleton(HSkeleton skeletonHandle) { m_owningSkeleton = skeletonHandle; }
    HSkeleton owningSkeleton() const { return m_owningSkeleton; }

    void setJointManager(JointManager *jointManager) { m_jointManager = jointManager; }
    JointManager *jointManager() const { return m_jointManager; }

    void setSkeletonManager(SkeletonManager *skeletonManager) { m_skeletonManager = skeletonManager; }
    SkeletonManager *skeletonManager() const { return m_skeletonManager; }

private:
    QMatrix4x4 m_inverseBindMatrix;
    Qt3DCore::Sqt m_localPose;
    QVector<Qt3DCore::QNodeId> m_childJointIds;
    QString m_name;
    JointManager *m_jointManager;
    SkeletonManager *m_skeletonManager;
    HSkeleton m_owningSkeleton;
};

class JointFunctor : public Qt3DCore::QBackendNodeMapper
{
public:
    explicit JointFunctor(AbstractRenderer *renderer,
                          JointManager *jointManager,
                          SkeletonManager *skeletonManager);
    Qt3DCore::QBackendNode *create(const Qt3DCore::QNodeCreatedChangeBasePtr &change) const final;
    Qt3DCore::QBackendNode *get(Qt3DCore::QNodeId id) const final;
    void destroy(Qt3DCore::QNodeId id) const final;

private:
    AbstractRenderer *m_renderer;
    JointManager *m_jointManager;
    SkeletonManager *m_skeletonManager;
};

} // namespace Render
} // namespace Qt3DRender


QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_JOINT_H
