/****************************************************************************
**
** Copyright (C) 2016 Paul Lemire
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

#ifndef QT3DRENDER_RENDER_FRUSTUMCULLINGJOB_P_H
#define QT3DRENDER_RENDER_FRUSTUMCULLINGJOB_P_H

#include <Qt3DCore/qaspectjob.h>
#include <Qt3DCore/private/matrix4x4_p.h>
#include <Qt3DCore/private/vector3d_p.h>
#include <Qt3DCore/private/vector4d_p.h>
#include <Qt3DRender/private/aligned_malloc_p.h>
#include <Qt3DRender/private/qt3drender_global_p.h>

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

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

namespace Render {

class Entity;
class EntityManager;
class NodeManagers;
struct Plane;

class Q_3DRENDERSHARED_PRIVATE_EXPORT FrustumCullingJob : public Qt3DCore::QAspectJob
{
public:
    FrustumCullingJob();
    ~FrustumCullingJob();

    QT3D_ALIGNED_MALLOC_AND_FREE()

    inline void setRoot(Entity *root) Q_DECL_NOTHROW { m_root = root; }
    inline void setManagers(NodeManagers *manager) Q_DECL_NOTHROW { m_manager = manager; }
    inline void setActive(bool active) Q_DECL_NOTHROW { m_active = active; }
    inline bool isActive() const Q_DECL_NOTHROW { return m_active; }
    inline void setViewProjection(const Matrix4x4 &viewProjection) Q_DECL_NOTHROW { m_viewProjection = viewProjection; }
    inline Matrix4x4 viewProjection() const Q_DECL_NOTHROW { return m_viewProjection; }

    QVector<Entity *> visibleEntities() const Q_DECL_NOTHROW { return m_visibleEntities; }

    void run() final;

private:
    struct Q_AUTOTEST_EXPORT Plane
    {
        explicit Plane(const Vector4D &planeEquation)
            : planeEquation(planeEquation)
            , normal(Vector3D(planeEquation).normalized())
            , d(planeEquation.w() / Vector3D(planeEquation).length())
        {}

        const Vector4D planeEquation;
        const Vector3D normal;
        const float d;
    };

    void cullScene(Entity *e, const Plane *planes);
    Matrix4x4 m_viewProjection;
    Entity *m_root;
    NodeManagers *m_manager;
    QVector<Entity *> m_visibleEntities;
    bool m_active;
};

typedef QSharedPointer<FrustumCullingJob> FrustumCullingJobPtr;

} // Render

} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_FRUSTUMCULLINGJOB_P_H
