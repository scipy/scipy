/****************************************************************************
**
** Copyright (C) 2008-2012 NVIDIA Corporation.
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QSSG_RENDER_IMPL_RENDERABLE_OBJECTS_H
#define QSSG_RENDER_IMPL_RENDERABLE_OBJECTS_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <QtQuick3DRuntimeRender/private/qssgrendermodel_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderdefaultmaterial_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendercustommaterial_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendermesh_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendershaderkeys_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendershadercache_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderableimage_p.h>

#include <QtQuick3DUtils/private/qssginvasivelinkedlist_p.h>

QT_BEGIN_NAMESPACE

enum class QSSGRenderableObjectFlag
{
    HasTransparency = 1 << 0,
    CompletelyTransparent = 1 << 1,
    Dirty = 1 << 2,
    Pickable = 1 << 3,
    DefaultMaterialMeshSubset = 1 << 4,
    Custom = 1 << 6,
    CustomMaterialMeshSubset = 1 << 7,
    HasRefraction = 1 << 8,
    CastsShadows = 1 << 9,
    ReceivesShadows = 1 << 10,
    HasAttributePosition = 1 << 11,
    HasAttributeNormal = 1 << 12,
    HasAttributeTexCoord0 = 1 << 13,
    HasAttributeTexCoord1 = 1 << 14,
    HasAttributeTangent = 1 << 15,
    HasAttributeBinormal = 1 << 16,
    HasAttributeColor = 1 << 17,
    HasSkeletalAnimation = 1 << 18
};

struct QSSGRenderableObjectFlags : public QFlags<QSSGRenderableObjectFlag>
{
    void setHasTransparency(bool inHasTransparency)
    {
        setFlag(QSSGRenderableObjectFlag::HasTransparency, inHasTransparency);
    }
    bool hasTransparency() const { return this->operator&(QSSGRenderableObjectFlag::HasTransparency); }
    bool hasRefraction() const { return this->operator&(QSSGRenderableObjectFlag::HasRefraction); }
    void setCompletelyTransparent(bool inTransparent)
    {
        setFlag(QSSGRenderableObjectFlag::CompletelyTransparent, inTransparent);
    }
    bool isCompletelyTransparent() const
    {
        return this->operator&(QSSGRenderableObjectFlag::CompletelyTransparent);
    }
    void setDirty(bool inDirty) { setFlag(QSSGRenderableObjectFlag::Dirty, inDirty); }
    bool isDirty() const { return this->operator&(QSSGRenderableObjectFlag::Dirty); }
    void setPickable(bool inPickable) { setFlag(QSSGRenderableObjectFlag::Pickable, inPickable); }
    bool isPickable() const { return this->operator&(QSSGRenderableObjectFlag::Pickable); }

    void setCastsShadows(bool inCastsShadows) { setFlag(QSSGRenderableObjectFlag::CastsShadows, inCastsShadows); }
    bool castsShadows() const { return this->operator&(QSSGRenderableObjectFlag::CastsShadows); }

    void setReceivesShadows(bool inReceivesShadows) { setFlag(QSSGRenderableObjectFlag::ReceivesShadows, inReceivesShadows); }
    bool receivesShadows() const { return this->operator&(QSSGRenderableObjectFlag::ReceivesShadows); }

    void setHasAttributePosition(bool b) { setFlag(QSSGRenderableObjectFlag::HasAttributePosition, b); }
    bool hasAttributePosition() const { return this->operator&(QSSGRenderableObjectFlag::HasAttributePosition); }

    void setHasAttributeNormal(bool b) { setFlag(QSSGRenderableObjectFlag::HasAttributeNormal, b); }
    bool hasAttributeNormal() const { return this->operator&(QSSGRenderableObjectFlag::HasAttributeNormal); }

    void setHasAttributeTexCoord0(bool b) { setFlag(QSSGRenderableObjectFlag::HasAttributeTexCoord0, b); }
    bool hasAttributeTexCoord0() const { return this->operator&(QSSGRenderableObjectFlag::HasAttributeTexCoord0); }

    void setHasAttributeTexCoord1(bool b) { setFlag(QSSGRenderableObjectFlag::HasAttributeTexCoord1, b); }
    bool hasAttributeTexCoord1() const { return this->operator&(QSSGRenderableObjectFlag::HasAttributeTexCoord1); }

    void setHasAttributeTangent(bool b) { setFlag(QSSGRenderableObjectFlag::HasAttributeTangent, b); }
    bool hasAttributeTangent() const { return this->operator&(QSSGRenderableObjectFlag::HasAttributeTangent); }

    void setHasAttributeBinormal(bool b) { setFlag(QSSGRenderableObjectFlag::HasAttributeBinormal, b); }
    bool hasAttributeBinormal() const { return this->operator&(QSSGRenderableObjectFlag::HasAttributeBinormal); }

    void setHasAttributeColor(bool b) { setFlag(QSSGRenderableObjectFlag::HasAttributeColor, b); }
    bool hasAttributeColor() const { return this->operator&(QSSGRenderableObjectFlag::HasAttributeColor); }

    // Mutually exclusive values
    void setDefaultMaterialMeshSubset(bool inMeshSubset)
    {
        setFlag(QSSGRenderableObjectFlag::DefaultMaterialMeshSubset, inMeshSubset);
    }
    bool isDefaultMaterialMeshSubset() const
    {
        return this->operator&(QSSGRenderableObjectFlag::DefaultMaterialMeshSubset);
    }

    void setCustomMaterialMeshSubset(bool inMeshSubset)
    {
        setFlag(QSSGRenderableObjectFlag::CustomMaterialMeshSubset, inMeshSubset);
    }
    bool isCustomMaterialMeshSubset() const
    {
        return this->operator&(QSSGRenderableObjectFlag::CustomMaterialMeshSubset);
    }

    void setCustom(bool inCustom) { setFlag(QSSGRenderableObjectFlag::Custom, inCustom); }
    bool isCustom() const { return this->operator&(QSSGRenderableObjectFlag::Custom); }
};

struct QSSGNodeLightEntry
{
    QSSGRenderLight *light = nullptr;
    qint32 lightIndex;
    QSSGNodeLightEntry *nextNode = nullptr;
    QSSGNodeLightEntry() = default;
    QSSGNodeLightEntry(QSSGRenderLight *inLight, qint32 inLightIndex)
        : light(inLight), lightIndex(inLightIndex), nextNode(nullptr)
    {
    }
};

using QSSGNodeLightEntryList = QSSGInvasiveSingleLinkedList<QSSGNodeLightEntry, &QSSGNodeLightEntry::nextNode>;

struct QSSGRenderableObject;

typedef void (*TRenderFunction)(QSSGRenderableObject &inObject, const QVector2D &inCameraProperties);

struct QSSGRenderableObject;

// Used for sorting
struct QSSGRenderableObjectHandle
{
    QSSGRenderableObject *obj;
    float cameraDistanceSq;
    static inline QSSGRenderableObjectHandle create(QSSGRenderableObject *o, float camDistSq = 0.0f) { return {o, camDistSq};}
};
Q_DECLARE_TYPEINFO(QSSGRenderableObjectHandle, Q_PRIMITIVE_TYPE);

struct QSSGRenderableObject
{
    // Variables used for picking
    const QMatrix4x4 &globalTransform;
    const QSSGBounds3 &bounds;
    QSSGRenderableObjectFlags renderableFlags;
    // For rough sorting for transparency and for depth
    QVector3D worldCenterPoint;
    TessellationModeValues tessellationMode;
    // For custom renderable objects the render function must be defined
    TRenderFunction renderFunction;
    QSSGNodeLightEntryList scopedLights;
    QSSGRenderableObject(QSSGRenderableObjectFlags inFlags,
                           const QVector3D &inWorldCenterPt,
                           const QMatrix4x4 &inGlobalTransform,
                           const QSSGBounds3 &inBounds,
                           TessellationModeValues inTessMode = TessellationModeValues::NoTessellation,
                           TRenderFunction inFunction = nullptr)

        : globalTransform(inGlobalTransform)
        , bounds(inBounds)
        , renderableFlags(inFlags)
        , worldCenterPoint(inWorldCenterPt)
        , tessellationMode(inTessMode)
        , renderFunction(inFunction)
    {
    }
};

Q_STATIC_ASSERT(std::is_trivially_destructible<QSSGRenderableObject>::value);

// Different subsets from the same model will get the same
// model context so we can generate the MVP and normal matrix once
// and only once per subset.
struct QSSGModelContext
{
    const QSSGRenderModel &model;
    QMatrix4x4 modelViewProjection;
    QMatrix3x3 normalMatrix;

    QSSGModelContext(const QSSGRenderModel &inModel, const QMatrix4x4 &inViewProjection) : model(inModel)
    {
        model.calculateMVPAndNormalMatrix(inViewProjection, modelViewProjection, normalMatrix);
    }
};

Q_STATIC_ASSERT(std::is_trivially_destructible<QSSGModelContext>::value);

class QSSGRendererImpl;
struct QSSGLayerRenderData;
struct QSSGShadowMapEntry;

struct QSSGSubsetRenderableBase : public QSSGRenderableObject
{
    const QSSGRef<QSSGRendererImpl> &generator;
    const QSSGModelContext &modelContext;
    const QSSGRenderSubset &subset;
    float opacity;

    QSSGSubsetRenderableBase(QSSGRenderableObjectFlags inFlags,
                               const QVector3D &inWorldCenterPt,
                               const QSSGRef<QSSGRendererImpl> &gen,
                               const QSSGRenderSubset &inSubset,
                               const QSSGModelContext &inModelContext,
                               float inOpacity);
    void renderShadowMapPass(const QVector2D &inCameraVec,
                             const QSSGRenderLight *inLight,
                             const QSSGRenderCamera &inCamera,
                             QSSGShadowMapEntry *inShadowMapEntry) const;

    void renderDepthPass(const QVector2D &inCameraVec, QSSGRenderableImage *inDisplacementImage, float inDisplacementAmount, QSSGCullFaceMode cullFaceMode = QSSGCullFaceMode::Back);
};

Q_STATIC_ASSERT(std::is_trivially_destructible<QSSGSubsetRenderableBase>::value);

/**
 *	A renderable that corresponds to a subset (a part of a model).
 *	These are created per subset per layer and are responsible for actually
 *	rendering this type of object.
 */
struct QSSGSubsetRenderable : public QSSGSubsetRenderableBase
{
    const QSSGRenderDefaultMaterial &material;
    QSSGRenderableImage *firstImage;
    QSSGShaderDefaultMaterialKey shaderDescription;
    QSSGDataView<QMatrix4x4> bones;

    QSSGSubsetRenderable(QSSGRenderableObjectFlags inFlags,
                           const QVector3D &inWorldCenterPt,
                           const QSSGRef<QSSGRendererImpl> &gen,
                           const QSSGRenderSubset &inSubset,
                           const QSSGRenderDefaultMaterial &mat,
                           const QSSGModelContext &inModelContext,
                           float inOpacity,
                           QSSGRenderableImage *inFirstImage,
                           QSSGShaderDefaultMaterialKey inShaderKey,
                           const QSSGDataView<QMatrix4x4> &inBoneGlobals);

    void render(const QVector2D &inCameraVec, const ShaderFeatureSetList &inFeatureSet);

    void renderDepthPass(const QVector2D &inCameraVec);

    QSSGRenderDefaultMaterial::MaterialBlendMode getBlendingMode() { return material.blendMode; }
};

Q_STATIC_ASSERT(std::is_trivially_destructible<QSSGSubsetRenderable>::value);

struct QSSGCustomMaterialRenderable : public QSSGSubsetRenderableBase
{
    const QSSGRenderCustomMaterial &material;
    QSSGRenderableImage *firstImage;
    QSSGShaderDefaultMaterialKey shaderDescription;

    QSSGCustomMaterialRenderable(QSSGRenderableObjectFlags inFlags,
                                   const QVector3D &inWorldCenterPt,
                                   const QSSGRef<QSSGRendererImpl> &gen,
                                   const QSSGRenderSubset &inSubset,
                                   const QSSGRenderCustomMaterial &mat,
                                   const QSSGModelContext &inModelContext,
                                   float inOpacity,
                                   QSSGRenderableImage *inFirstImage,
                                   QSSGShaderDefaultMaterialKey inShaderKey);

    void render(const QVector2D &inCameraVec,
                const QSSGLayerRenderData &inLayerData,
                const QSSGRenderLayer &inLayer,
                const QVector<QSSGRenderLight *> &inLights,
                const QSSGRenderCamera &inCamera,
                const QSSGRef<QSSGRenderTexture2D> &inDepthTexture,
                const QSSGRef<QSSGRenderTexture2D> &inSsaoTexture,
                const ShaderFeatureSetList &inFeatureSet);

    void renderDepthPass(const QVector2D &inCameraVec,
                         const QSSGRenderLayer &inLayer,
                         const QVector<QSSGRenderLight *> &inLights,
                         const QSSGRenderCamera &inCamera,
                         const QSSGRenderTexture2D *inDepthTexture);
};

Q_STATIC_ASSERT(std::is_trivially_destructible<QSSGCustomMaterialRenderable>::value);

QT_END_NAMESPACE

#endif
