/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Copyright (C) 2016 The Qt Company Ltd and/or its subsidiary(-ies).
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

#ifndef QT3DRENDER_RENDER_ENTITY_H
#define QT3DRENDER_RENDER_ENTITY_H

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
#include <Qt3DRender/private/abstractrenderer_p.h>
#include <Qt3DRender/private/handle_types_p.h>
#include <Qt3DCore/qnodecreatedchange.h>
#include <Qt3DCore/private/qentity_p.h>
#include <Qt3DCore/private/qhandle_p.h>
#include <QVector>

QT_BEGIN_NAMESPACE

class QMatrix4x4;

namespace Qt3DCore {
class QNode;
class QEntity;
class QComponent;
}

namespace Qt3DRender {

class QRenderAspect;

namespace Render {

class Sphere;
class Renderer;
class NodeManagers;
class EntityPrivate;

class Q_3DRENDERSHARED_PRIVATE_EXPORT Entity : public BackendNode
{
public:
    Entity();
    ~Entity();
    void cleanup();

    void setParentHandle(HEntity parentHandle);
    void setNodeManagers(NodeManagers *manager);
    void syncFromFrontEnd(const Qt3DCore::QNode *frontEnd, bool firstTime) override;

    void dump() const;

    void  setHandle(HEntity handle);
    HEntity handle() const { return m_handle; }
    Entity *parent() const;
    HEntity parentHandle() const { return m_parentHandle; }

    void removeFromParentChildHandles();
    void appendChildHandle(HEntity childHandle);
    void removeChildHandle(HEntity childHandle) { m_childrenHandles.removeOne(childHandle); }
    QVector<HEntity> childrenHandles() const { return m_childrenHandles; }
    QVector<Entity *> children() const;
    bool hasChildren() const { return !m_childrenHandles.empty(); }
    void traverse(const std::function<void(Entity *)> &operation);
    void traverse(const std::function<void(const Entity *)> &operation) const;

    Matrix4x4 *worldTransform();
    const Matrix4x4 *worldTransform() const;
    Sphere *localBoundingVolume() const { return m_localBoundingVolume.data(); }
    Sphere *worldBoundingVolume() const { return m_worldBoundingVolume.data(); }
    Sphere *worldBoundingVolumeWithChildren() const { return m_worldBoundingVolumeWithChildren.data(); }

    void addComponent(Qt3DCore::QNodeIdTypePair idAndType);
    void removeComponent(Qt3DCore::QNodeId nodeId);

    bool isBoundingVolumeDirty() const;
    void unsetBoundingVolumeDirty();

    void setTreeEnabled(bool enabled) { m_treeEnabled = enabled; }
    bool isTreeEnabled() const { return m_treeEnabled; }

    Qt3DCore::QNodeIdVector layerIds() const { return m_layerComponents + m_recursiveLayerComponents; }
    void addRecursiveLayerId(const Qt3DCore::QNodeId layerId);
    void removeRecursiveLayerId(const Qt3DCore::QNodeId layerId);
    void clearRecursiveLayerIds() { m_recursiveLayerComponents.clear(); }

    template<class Backend>
    Qt3DCore::QHandle<Backend> componentHandle() const
    {
        return Qt3DCore::QHandle<Backend>();
    }

    template<class Backend>
    QVector<Qt3DCore::QHandle<Backend> > componentsHandle() const
    {
        return QVector<Qt3DCore::QHandle<Backend> >();
    }

    template<class Backend>
    Backend *renderComponent() const
    {
        return nullptr;
    }

    template<class Backend>
    QVector<Backend *> renderComponents() const
    {
        return QVector<Backend *>();
    }

    template<class Backend>
    Qt3DCore::QNodeId componentUuid() const
    {
        return Qt3DCore::QNodeId();
    }

    template<class Backend>
    QVector<Qt3DCore::QNodeId> componentsUuid() const
    {
        return QVector<Qt3DCore::QNodeId>();
    }

    template<typename T>
    bool containsComponentsOfType() const
    {
        return !componentUuid<T>().isNull();
    }

    template<typename T, typename Ts, typename ... Ts2>
    bool containsComponentsOfType() const
    {
        return containsComponentsOfType<T>() && containsComponentsOfType<Ts, Ts2...>();
    }

protected:
    Q_DECLARE_PRIVATE(Entity)

private:
    NodeManagers *m_nodeManagers;
    HEntity m_handle;
    HEntity m_parentHandle;
    QVector<HEntity > m_childrenHandles;

    HMatrix m_worldTransform;
    QSharedPointer<Sphere> m_localBoundingVolume;
    QSharedPointer<Sphere> m_worldBoundingVolume;
    QSharedPointer<Sphere> m_worldBoundingVolumeWithChildren;

    // Handles to Components
    Qt3DCore::QNodeId m_transformComponent;
    Qt3DCore::QNodeId m_materialComponent;
    Qt3DCore::QNodeId m_cameraComponent;
    QVector<Qt3DCore::QNodeId> m_layerComponents;
    QVector<Qt3DCore::QNodeId> m_levelOfDetailComponents;
    QVector<Qt3DCore::QNodeId> m_rayCasterComponents;
    QVector<Qt3DCore::QNodeId> m_shaderDataComponents;
    QVector<Qt3DCore::QNodeId> m_lightComponents;
    QVector<Qt3DCore::QNodeId> m_environmentLightComponents;
    Qt3DCore::QNodeId m_geometryRendererComponent;
    Qt3DCore::QNodeId m_objectPickerComponent;
    Qt3DCore::QNodeId m_boundingVolumeDebugComponent;
    Qt3DCore::QNodeId m_computeComponent;
    Qt3DCore::QNodeId m_armatureComponent;

    // Includes recursive layers
    Qt3DCore::QNodeIdVector m_recursiveLayerComponents;

    QString m_objectName;
    bool m_boundingDirty;
    // true only if this and all parent nodes are enabled
    bool m_treeEnabled;
};

#define ENTITY_COMPONENT_TEMPLATE_SPECIALIZATION(Type, Handle) \
    /* Handle */ \
    template<> \
    Q_3DRENDERSHARED_PRIVATE_EXPORT Handle Entity::componentHandle<Type>() const; \
    /* Component */ \
    template<> \
    Q_3DRENDERSHARED_PRIVATE_EXPORT Type *Entity::renderComponent<Type>() const; \
    /* Uuid */ \
    template<> \
    Q_3DRENDERSHARED_PRIVATE_EXPORT Qt3DCore::QNodeId Entity::componentUuid<Type>() const;


#define ENTITY_COMPONENT_LIST_TEMPLATE_SPECIALIZATION(Type, Handle) \
    /* Handle */ \
    template<> \
    Q_3DRENDERSHARED_PRIVATE_EXPORT QVector<Handle> Entity::componentsHandle<Type>() const; \
    /* Component */ \
    template<> \
    Q_3DRENDERSHARED_PRIVATE_EXPORT QVector<Type *> Entity::renderComponents<Type>() const; \
    /* Uuid */ \
    template<> \
    Q_3DRENDERSHARED_PRIVATE_EXPORT Qt3DCore::QNodeIdVector Entity::componentsUuid<Type>() const;

#define ENTITY_COMPONENT_TEMPLATE_IMPL(Type, Handle, Manager, variable) \
    /* Handle */ \
    template<> \
    Handle Entity::componentHandle<Type>() const \
    { \
        return m_nodeManagers->lookupHandle<Type, Manager, Handle>(variable); \
    } \
    /* Component */ \
    template<> \
    Type *Entity::renderComponent<Type>() const \
    { \
        return m_nodeManagers->lookupResource<Type, Manager>(variable); \
    } \
    /* Uuid */ \
    template<> \
    Qt3DCore::QNodeId Entity::componentUuid<Type>() const \
    { \
        return variable; \
    }

#define ENTITY_COMPONENT_LIST_TEMPLATE_IMPL(Type, Handle, Manager, variable) \
    /* Handle */ \
    template<> \
    QVector<Handle> Entity::componentsHandle<Type>() const \
    { \
        Manager *manager = m_nodeManagers->manager<Type, Manager>(); \
        QVector<Handle> entries; \
        entries.reserve(variable.size()); \
        for (const QNodeId id : variable) \
            entries.push_back(manager->lookupHandle(id)); \
        return entries; \
        } \
    /* Component */ \
    template<> \
    QVector<Type *> Entity::renderComponents<Type>() const \
    { \
        Manager *manager = m_nodeManagers->manager<Type, Manager>(); \
        QVector<Type *> entries; \
        entries.reserve(variable.size()); \
        for (const QNodeId id : variable) \
            entries.push_back(manager->lookupResource(id)); \
        return entries; \
    } \
    /* Uuid */ \
    template<> \
    Qt3DCore::QNodeIdVector Entity::componentsUuid<Type>() const \
    { \
        return variable; \
    }

ENTITY_COMPONENT_TEMPLATE_SPECIALIZATION(Material, HMaterial)
ENTITY_COMPONENT_TEMPLATE_SPECIALIZATION(CameraLens, HCamera)
ENTITY_COMPONENT_TEMPLATE_SPECIALIZATION(Transform, HTransform)
ENTITY_COMPONENT_TEMPLATE_SPECIALIZATION(GeometryRenderer, HGeometryRenderer)
ENTITY_COMPONENT_TEMPLATE_SPECIALIZATION(ObjectPicker, HObjectPicker)
ENTITY_COMPONENT_TEMPLATE_SPECIALIZATION(ComputeCommand, HComputeCommand)
ENTITY_COMPONENT_TEMPLATE_SPECIALIZATION(Armature, HArmature)
ENTITY_COMPONENT_LIST_TEMPLATE_SPECIALIZATION(Layer, HLayer)
ENTITY_COMPONENT_LIST_TEMPLATE_SPECIALIZATION(LevelOfDetail, HLevelOfDetail)
ENTITY_COMPONENT_LIST_TEMPLATE_SPECIALIZATION(RayCaster, HRayCaster)
ENTITY_COMPONENT_LIST_TEMPLATE_SPECIALIZATION(ShaderData, HShaderData)
ENTITY_COMPONENT_LIST_TEMPLATE_SPECIALIZATION(Light, HLight)
ENTITY_COMPONENT_LIST_TEMPLATE_SPECIALIZATION(EnvironmentLight, HEnvironmentLight)

class Q_AUTOTEST_EXPORT RenderEntityFunctor : public Qt3DCore::QBackendNodeMapper
{
public:
    explicit RenderEntityFunctor(AbstractRenderer *renderer, NodeManagers *manager);
    Qt3DCore::QBackendNode *create(const Qt3DCore::QNodeCreatedChangeBasePtr &change) const override;
    Qt3DCore::QBackendNode *get(Qt3DCore::QNodeId id) const override;
    void destroy(Qt3DCore::QNodeId id) const override;

private:
    NodeManagers *m_nodeManagers;
    AbstractRenderer *m_renderer;
};

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_ENTITY_H
