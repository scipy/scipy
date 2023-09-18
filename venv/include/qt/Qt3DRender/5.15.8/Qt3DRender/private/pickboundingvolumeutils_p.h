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

#ifndef QT3DRENDER_RENDER_PICKBOUNDINGVOLUMEUTILS_H
#define QT3DRENDER_RENDER_PICKBOUNDINGVOLUMEUTILS_H

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

#include <Qt3DCore/QNodeId>
#include <Qt3DRender/QAbstractRayCaster>
#include <Qt3DRender/private/qray3d_p.h>
#include <Qt3DRender/private/qraycastingservice_p.h>
#include <Qt3DRender/qpickingsettings.h>


QT_BEGIN_NAMESPACE

class QSurface;

namespace Qt3DRender {
namespace RayCasting {
class QAbstractCollisionQueryService;
}

namespace Render {

class Entity;
class Renderer;
class FrameGraphNode;
class NodeManagers;

namespace PickingUtils {

struct Q_AUTOTEST_EXPORT ViewportCameraAreaDetails
{
    Qt3DCore::QNodeId cameraId;
    Qt3DCore::QNodeId viewportNodeId;
    QRectF viewport;
    QSize area;
    QSurface *surface = nullptr;
    Qt3DCore::QNodeIdVector layersFilters;
};
QT3D_DECLARE_TYPEINFO_3(Qt3DRender, Render, PickingUtils, ViewportCameraAreaDetails, Q_COMPLEX_TYPE)

class Q_AUTOTEST_EXPORT ViewportCameraAreaGatherer
{
public:
    ViewportCameraAreaGatherer(const Qt3DCore::QNodeId &nodeId = Qt3DCore::QNodeId()) : m_targetCamera(nodeId) { }
    QVector<ViewportCameraAreaDetails> gather(FrameGraphNode *root);

private:
    Qt3DCore::QNodeId m_targetCamera;
    QVector<FrameGraphNode *> m_leaves;

    void visit(FrameGraphNode *node);
    ViewportCameraAreaDetails gatherUpViewportCameraAreas(Render::FrameGraphNode *node) const;
    bool isUnique(const QVector<ViewportCameraAreaDetails> &vcaList, const ViewportCameraAreaDetails &vca) const;
};

typedef QVector<RayCasting::QCollisionQueryResult::Hit> HitList;

class Q_AUTOTEST_EXPORT HierarchicalEntityPicker
{
public:
    explicit HierarchicalEntityPicker(const RayCasting::QRay3D &ray, bool requireObjectPicker = true);

    void setLayerFilterIds(const Qt3DCore::QNodeIdVector &layerFilterIds);
    void setLayerIds(const Qt3DCore::QNodeIdVector &layerIds, QAbstractRayCaster::FilterMode mode);

    bool collectHits(NodeManagers *manager, Entity *root);
    inline HitList hits() const { return m_hits; }
    inline QVector<Entity *> entities() const { return m_entities; }
    inline QHash<Qt3DCore::QNodeId, int> entityToPriorityTable() const { return m_entityToPriorityTable; }

private:
    RayCasting::QRay3D m_ray;
    HitList m_hits;
    QVector<Entity *> m_entities;
    bool m_objectPickersRequired;
    Qt3DCore::QNodeIdVector m_layerFilterIds;
    Qt3DCore::QNodeIdVector m_layerIds;
    QAbstractRayCaster::FilterMode m_layerFilterMode = QAbstractRayCaster::AcceptAnyMatchingLayers;
    QHash<Qt3DCore::QNodeId, int> m_entityToPriorityTable;
};

struct Q_AUTOTEST_EXPORT AbstractCollisionGathererFunctor
{
    AbstractCollisionGathererFunctor();
    virtual ~AbstractCollisionGathererFunctor();

    bool m_objectPickersRequired = true;
    NodeManagers *m_manager = nullptr;
    RayCasting::QRay3D m_ray;
    QHash<Qt3DCore::QNodeId, int> m_entityToPriorityTable;

    virtual HitList computeHits(const QVector<Entity *> &entities, Qt3DRender::QPickingSettings::PickResultMode mode) = 0;

    // This define is required to work with QtConcurrent
    typedef HitList result_type;
    HitList operator ()(const Entity *entity) const;
    virtual HitList pick(const Entity *entity) const = 0;

    bool rayHitsEntity(const Entity *entity) const;
    static void sortHits(HitList &results);
};

struct Q_AUTOTEST_EXPORT EntityCollisionGathererFunctor : public AbstractCollisionGathererFunctor
{
    HitList computeHits(const QVector<Entity *> &entities, Qt3DRender::QPickingSettings::PickResultMode mode) override;
    HitList pick(const Entity *entity) const override;
};

struct Q_AUTOTEST_EXPORT TriangleCollisionGathererFunctor : public AbstractCollisionGathererFunctor
{
    bool m_frontFaceRequested;
    bool m_backFaceRequested;

    HitList computeHits(const QVector<Entity *> &entities, Qt3DRender::QPickingSettings::PickResultMode mode) override;
    HitList pick(const Entity *entity) const override;
};

struct Q_AUTOTEST_EXPORT LineCollisionGathererFunctor : public AbstractCollisionGathererFunctor
{
    float m_pickWorldSpaceTolerance;

    HitList computeHits(const QVector<Entity *> &entities, Qt3DRender::QPickingSettings::PickResultMode mode) override;
    HitList pick(const Entity *entity) const override;
};

struct Q_AUTOTEST_EXPORT PointCollisionGathererFunctor : public AbstractCollisionGathererFunctor
{
    float m_pickWorldSpaceTolerance;

    HitList computeHits(const QVector<Entity *> &entities, Qt3DRender::QPickingSettings::PickResultMode mode) override;
    HitList pick(const Entity *entity) const override;
};

} // PickingUtils

} // Render

} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_PICKBOUNDINGVOLUMEUTILS_H
