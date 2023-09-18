/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_RENDER_HANDLE_TYPES_P_H
#define QT3DRENDER_RENDER_HANDLE_TYPES_P_H

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

#include <Qt3DRender/qt3drender_global.h>
#include <Qt3DCore/private/qhandle_p.h>
#include <Qt3DCore/private/matrix4x4_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QTextureImageData;

namespace Render {

class RenderTargetOutput;
class CameraLens;
class FilterKey;
class Effect;
class Entity;
class Shader;
class ShaderBuilder;
class FrameGraphNode;
class Layer;
class LevelOfDetail;
class Material;
class Technique;
class Texture;
class Transform;
class RenderTarget;
class RenderPass;
class Parameter;
class ShaderData;
class TextureImage;
class Buffer;
class Attribute;
class Geometry;
class GeometryRenderer;
class ObjectPicker;
class RayCaster;
class BoundingVolumeDebug;
class Light;
class EnvironmentLight;
class ComputeCommand;
class RenderStateNode;
class Armature;
class Skeleton;
class Joint;
class ShaderImage;

typedef Qt3DCore::QHandle<RenderTargetOutput> HAttachment;
typedef Qt3DCore::QHandle<CameraLens> HCamera;
typedef Qt3DCore::QHandle<FilterKey> HFilterKey;
typedef Qt3DCore::QHandle<Effect> HEffect;
typedef Qt3DCore::QHandle<Entity> HEntity;
typedef Qt3DCore::QHandle<FrameGraphNode *> HFrameGraphNode;
typedef Qt3DCore::QHandle<Layer> HLayer;
typedef Qt3DCore::QHandle<LevelOfDetail> HLevelOfDetail;
typedef Qt3DCore::QHandle<Material> HMaterial;
typedef Qt3DCore::QHandle<Matrix4x4> HMatrix;
typedef Qt3DCore::QHandle<Shader> HShader;
typedef Qt3DCore::QHandle<ShaderBuilder> HShaderBuilder;
typedef Qt3DCore::QHandle<Technique> HTechnique;
typedef Qt3DCore::QHandle<Texture> HTexture;
typedef Qt3DCore::QHandle<Transform> HTransform;
typedef Qt3DCore::QHandle<RenderTarget> HTarget;
typedef Qt3DCore::QHandle<RenderPass> HRenderPass;
typedef Qt3DCore::QHandle<QTextureImageData> HTextureData;
typedef Qt3DCore::QHandle<Parameter> HParameter;
typedef Qt3DCore::QHandle<ShaderData> HShaderData;
typedef Qt3DCore::QHandle<TextureImage> HTextureImage;
typedef Qt3DCore::QHandle<Buffer> HBuffer;
typedef Qt3DCore::QHandle<Attribute> HAttribute;
typedef Qt3DCore::QHandle<Geometry> HGeometry;
typedef Qt3DCore::QHandle<GeometryRenderer> HGeometryRenderer;
typedef Qt3DCore::QHandle<ObjectPicker> HObjectPicker;
typedef Qt3DCore::QHandle<RayCaster> HRayCaster;
typedef Qt3DCore::QHandle<BoundingVolumeDebug> HBoundingVolumeDebug;
typedef Qt3DCore::QHandle<Light> HLight;
typedef Qt3DCore::QHandle<EnvironmentLight> HEnvironmentLight;
typedef Qt3DCore::QHandle<ComputeCommand> HComputeCommand;
typedef Qt3DCore::QHandle<RenderStateNode> HRenderState;
typedef Qt3DCore::QHandle<Armature> HArmature;
typedef Qt3DCore::QHandle<Skeleton> HSkeleton;
typedef Qt3DCore::QHandle<Joint> HJoint;
typedef Qt3DCore::QHandle<ShaderImage> HShaderImage;

} // namespace Render

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_HANDLE_TYPES_P_H
