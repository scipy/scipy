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

#ifndef QT3DRENDER_RENDER_LOADSKELETONJOB_P_H
#define QT3DRENDER_RENDER_LOADSKELETONJOB_P_H

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

#include <Qt3DCore/qaspectjob.h>

#include <QtCore/qsharedpointer.h>
#include <Qt3DRender/private/skeletondata_p.h>
#include <Qt3DRender/private/handle_types_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {
class QJoint;
}

namespace Qt3DRender {
namespace Render {

class NodeManagers;
class LoadSkeletonJobPrivate;

class LoadSkeletonJob : public Qt3DCore::QAspectJob
{
public:
    explicit LoadSkeletonJob(const HSkeleton &handle);

    void setNodeManagers(NodeManagers *nodeManagers) { m_nodeManagers = nodeManagers; }

protected:
    void run() override;
    void loadSkeleton(Skeleton *skeleton);
    void loadSkeletonFromUrl(Skeleton *skeleton);
    void loadSkeletonFromData(Skeleton *skeleton);
    Qt3DCore::QJoint *createFrontendJoints(const SkeletonData &skeletonData) const;
    Qt3DCore::QJoint *createFrontendJoint(const QString &jointName,
                                          const Qt3DCore::Sqt &localPose,
                                          const QMatrix4x4 &inverseBindMatrix) const;
    void processJointHierarchy(Qt3DCore::QNodeId jointId,
                               int parentJointIndex,
                               SkeletonData &skeletonData);

    HSkeleton m_handle;
    NodeManagers *m_nodeManagers;

private:
    Q_DECLARE_PRIVATE(LoadSkeletonJob)
};

typedef QSharedPointer<LoadSkeletonJob> LoadSkeletonJobPtr;

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_LOADSKELETONJOB_P_H
