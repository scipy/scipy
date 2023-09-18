/****************************************************************************
**
** Copyright (C) 2018 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_RENDER_ABSTRACTPICKINGJOB_P_H
#define QT3DRENDER_RENDER_ABSTRACTPICKINGJOB_P_H

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
#include <Qt3DRender/private/qray3d_p.h>
#include <Qt3DRender/private/handle_types_p.h>
#include <Qt3DRender/private/qboundingvolumeprovider_p.h>
#include <Qt3DRender/private/qcollisionqueryresult_p.h>
#include <Qt3DRender/private/pickboundingvolumeutils_p.h>
#include <Qt3DRender/private/qt3drender_global_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {
class QNodeId;
}

namespace Qt3DRender {
namespace Render {

class Entity;
class Renderer;
class NodeManagers;
class RenderSettings;

class Q_3DRENDERSHARED_PRIVATE_EXPORT AbstractPickingJob : public Qt3DCore::QAspectJob
{
public:
    AbstractPickingJob();

    void setRoot(Entity *root);
    void setFrameGraphRoot(FrameGraphNode *frameGraphRoot);
    void setRenderSettings(RenderSettings *settings);
    void setManagers(NodeManagers *manager);

    // public for unit tests
    virtual bool runHelper() = 0;
    static RayCasting::QRay3D intersectionRay(const QPoint &pos,
                                              const Matrix4x4 &viewMatrix,
                                              const Matrix4x4 &projectionMatrix,
                                              const QRect &viewport);

protected:
    AbstractPickingJob(Qt3DCore::QAspectJobPrivate &dd);

    void run() final;

    NodeManagers *m_manager;
    Entity *m_node;
    FrameGraphNode *m_frameGraphRoot;
    RenderSettings *m_renderSettings;

    bool m_oneEnabledAtLeast;

    QRect windowViewport(const QSize &area, const QRectF &relativeViewport) const;
    RayCasting::QRay3D rayForViewportAndCamera(const PickingUtils::ViewportCameraAreaDetails &vca,
                                               QObject *eventSource,
                                               const QPoint &pos) const;
};

} // namespace Render

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_ABSTRACTPICKINGJOB_P_H
