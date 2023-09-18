/****************************************************************************
**
** Copyright (C) 2015 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_RENDER_NODEMANAGERS_H
#define QT3DRENDER_RENDER_NODEMANAGERS_H

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

#include <Qt3DCore/private/qresourcemanager_p.h>
#include <Qt3DRender/private/qt3drender_global_p.h>

QT_BEGIN_NAMESPACE

class QMatrix4x4;

namespace Qt3DRender {

class QTextureImageData;

namespace Render {

class CameraManager;
class EntityManager;
class BufferManager;
class AttributeManager;
class GeometryManager;
class GeometryRendererManager;
class ObjectPickerManager;
class RayCasterManager;
class BoundingVolumeDebugManager;
class MaterialManager;
class MatrixManager;
class ShaderManager;
class ShaderBuilderManager;
class TechniqueManager;
class EffectManager;
class RenderPassManager;
class RenderTargetManager;
class SceneManager;
class AttachmentManager;
class ParameterManager;
class ShaderDataManager;
class TextureImageManager;
class FilterKeyManager;
class FrameGraphManager;
class TransformManager;
class TextureManager;
class TextureDataManager;
class TextureImageDataManager;
class LayerManager;
class LevelOfDetailManager;
class LightManager;
class EnvironmentLightManager;
class ComputeCommandManager;
class RenderStateManager;
class ArmatureManager;
class SkeletonManager;
class JointManager;
class ShaderImageManager;

class FrameGraphNode;
class Entity;
class CameraLens;
class Material;
class Shader;
class ShaderBuilder;
class Technique;
class Effect;
class RenderPass;
class Texture;
class Layer;
class LevelOfDetail;
class FilterKey;
class FrameGraphNode;
class Transform;
class Scene;
class RenderTargetOutput;
class RenderTarget;
class ShaderData;
class Parameter;
class TextureImage;
class Buffer;
class Attribute;
class Geometry;
class GeometryRenderer;
class ObjectPicker;
class RayCaster;
//class BoundingVolumeDebug;
class Light;
class EnvironmentLight;
class ComputeCommand;
class RenderStateNode;
class Armature;
class Skeleton;
class Joint;
class ShaderImage;


class Q_3DRENDERSHARED_PRIVATE_EXPORT NodeManagers
{
public:
    NodeManagers();
    ~NodeManagers();

    template<class Backend, typename Manager>
    Manager *manager() const noexcept
    {
        return nullptr;
    }

    template<class Backend, typename Manager, typename Key>
    Backend *lookupResource(const Key &id) const
    {
        Manager *mgr = manager<Backend, Manager>();
        if (mgr != nullptr)
            return mgr->lookupResource(id);
        return nullptr;
    }

    template<class Backend, typename Manager, typename Handle, typename Key>
    Handle lookupHandle(const Key &id) const
    {
        Manager *mgr = manager<Backend, Manager>();
        if (mgr != nullptr)
            return mgr->lookupHandle(id);
        return Handle();
    }

    template<class Backend, typename Manager, typename Handle>
    Backend *data(Handle handle)
    {
        Manager *mgr = manager<Backend, Manager>();
        if (mgr != nullptr)
            return mgr->data(handle);
        return nullptr;
    }


    inline CameraManager *cameraManager() const noexcept { return m_cameraManager; }
    inline EntityManager *renderNodesManager() const noexcept { return m_renderNodesManager; }
    inline MaterialManager *materialManager() const noexcept { return m_materialManager; }
    inline MatrixManager *worldMatrixManager() const noexcept { return m_worldMatrixManager; }
    inline ShaderManager *shaderManager() const noexcept { return m_shaderManager; }
    inline ShaderBuilderManager *shaderBuilderManager() const noexcept { return m_shaderBuilderManager; }
    inline TechniqueManager *techniqueManager() const noexcept { return m_techniqueManager; }
    inline EffectManager *effectManager() const noexcept { return m_effectManager; }
    inline RenderPassManager *renderPassManager() const noexcept { return m_renderPassManager; }
    inline TextureManager *textureManager() const noexcept { return m_textureManager; }
    inline LayerManager *layerManager() const noexcept { return m_layerManager; }
    inline LevelOfDetailManager *levelOfDetailManager() const noexcept { return m_levelOfDetailManager; }
    inline FilterKeyManager *filterKeyManager() const noexcept { return m_filterKeyManager; }
    inline FrameGraphManager *frameGraphManager() const noexcept { return m_frameGraphManager; }
    inline TransformManager *transformManager() const noexcept { return m_transformManager; }
    inline RenderTargetManager *renderTargetManager() const noexcept { return m_renderTargetManager; }
    inline SceneManager *sceneManager() const noexcept { return m_sceneManager; }
    inline AttachmentManager *attachmentManager() const noexcept { return m_attachmentManager; }
    inline ParameterManager *parameterManager() const noexcept { return m_parameterManager; }
    inline ShaderDataManager *shaderDataManager() const noexcept { return m_shaderDataManager; }
    inline TextureImageManager *textureImageManager() const noexcept { return m_textureImageManager; }
    inline BufferManager *bufferManager() const noexcept { return m_bufferManager; }
    inline AttributeManager *attributeManager() const noexcept { return m_attributeManager; }
    inline GeometryManager *geometryManager() const noexcept { return m_geometryManager; }
    inline GeometryRendererManager *geometryRendererManager() const noexcept { return m_geometryRendererManager; }
    inline ObjectPickerManager *objectPickerManager() const noexcept { return m_objectPickerManager; }
    inline RayCasterManager *rayCasterManager() const noexcept { return m_rayCasterManager; }
    //    inline BoundingVolumeDebugManager *boundingVolumeDebugManager() const noexcept { return m_boundingVolumeDebugManager; }
    inline LightManager *lightManager() const noexcept { return m_lightManager; }
    inline EnvironmentLightManager *environmentLightManager() const noexcept { return m_environmentLightManager; }
    inline ComputeCommandManager *computeJobManager() const noexcept { return m_computeJobManager; }
    inline RenderStateManager *renderStateManager() const noexcept { return m_renderStateManager; }
    inline ArmatureManager *armatureManager() const noexcept { return m_armatureManager; }
    inline SkeletonManager *skeletonManager() const noexcept { return m_skeletonManager; }
    inline JointManager *jointManager() const noexcept { return m_jointManager; }
    inline ShaderImageManager *shaderImageManager() const noexcept { return m_shaderImageManager; }

private:
    CameraManager *m_cameraManager;
    EntityManager *m_renderNodesManager;
    MaterialManager *m_materialManager;
    MatrixManager *m_worldMatrixManager;
    ShaderManager *m_shaderManager;
    ShaderBuilderManager *m_shaderBuilderManager;
    TechniqueManager *m_techniqueManager;
    EffectManager *m_effectManager;
    RenderPassManager *m_renderPassManager;
    TextureManager *m_textureManager;
    TextureImageManager *m_textureImageManager;
    LayerManager *m_layerManager;
    LevelOfDetailManager *m_levelOfDetailManager;
    FilterKeyManager *m_filterKeyManager;
    FrameGraphManager *m_frameGraphManager;
    TransformManager *m_transformManager;
    RenderTargetManager *m_renderTargetManager;
    SceneManager *m_sceneManager;
    AttachmentManager *m_attachmentManager;
    ParameterManager *m_parameterManager;
    ShaderDataManager *m_shaderDataManager;
    BufferManager *m_bufferManager;
    AttributeManager *m_attributeManager;
    GeometryManager *m_geometryManager;
    GeometryRendererManager *m_geometryRendererManager;
    ObjectPickerManager *m_objectPickerManager;
    RayCasterManager *m_rayCasterManager;
    //    BoundingVolumeDebugManager *m_boundingVolumeDebugManager;
    LightManager *m_lightManager;
    EnvironmentLightManager *m_environmentLightManager;
    ComputeCommandManager *m_computeJobManager;
    RenderStateManager *m_renderStateManager;
    ArmatureManager *m_armatureManager;
    SkeletonManager *m_skeletonManager;
    JointManager *m_jointManager;
    ShaderImageManager *m_shaderImageManager;
};

// Specializations

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT CameraManager *NodeManagers::manager<CameraLens>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT EntityManager *NodeManagers::manager<Entity>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT MaterialManager *NodeManagers::manager<Material>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT MatrixManager *NodeManagers::manager<QMatrix4x4*>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT ShaderManager *NodeManagers::manager<Shader>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT ShaderBuilderManager *NodeManagers::manager<ShaderBuilder>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT TechniqueManager *NodeManagers::manager<Technique>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT EffectManager *NodeManagers::manager<Effect>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT RenderPassManager *NodeManagers::manager<RenderPass>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT TextureManager *NodeManagers::manager<Texture>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT LayerManager *NodeManagers::manager<Layer>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT LevelOfDetailManager *NodeManagers::manager<LevelOfDetail>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT FilterKeyManager *NodeManagers::manager<FilterKey>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT FrameGraphManager *NodeManagers::manager<FrameGraphNode*>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT TransformManager *NodeManagers::manager<Transform>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT RenderTargetManager *NodeManagers::manager<RenderTarget>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT SceneManager *NodeManagers::manager<Scene>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT AttachmentManager *NodeManagers::manager<RenderTargetOutput>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT ParameterManager *NodeManagers::manager<Parameter>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT ShaderDataManager *NodeManagers::manager<ShaderData>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT BufferManager *NodeManagers::manager<Buffer>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT TextureImageManager *NodeManagers::manager<TextureImage>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT AttributeManager *NodeManagers::manager<Attribute>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT GeometryManager *NodeManagers::manager<Geometry>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT GeometryRendererManager *NodeManagers::manager<GeometryRenderer>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT ObjectPickerManager *NodeManagers::manager<ObjectPicker>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT RayCasterManager *NodeManagers::manager<RayCaster>() const noexcept;

//template<>
//Q_3DRENDERSHARED_PRIVATE_EXPORT BoundingVolumeDebugManager *NodeManagers::manager<BoundingVolumeDebug>() const noexcept;

template<>
LightManager *NodeManagers::manager<Light>() const noexcept;

template<>
EnvironmentLightManager *NodeManagers::manager<EnvironmentLight>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT ComputeCommandManager *NodeManagers::manager<ComputeCommand>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT RenderStateManager *NodeManagers::manager<RenderStateNode>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT ArmatureManager *NodeManagers::manager<Armature>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT SkeletonManager *NodeManagers::manager<Skeleton>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT JointManager *NodeManagers::manager<Joint>() const noexcept;

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT ShaderImageManager *NodeManagers::manager<ShaderImage>() const noexcept;

} // Render

} // Qt3DRender

QT_END_NAMESPACE


#endif // NODEMANAGERS_H
