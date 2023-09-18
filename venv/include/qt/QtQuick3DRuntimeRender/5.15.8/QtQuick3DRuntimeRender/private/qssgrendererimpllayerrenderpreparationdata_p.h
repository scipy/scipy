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

#ifndef QSSG_RENDERER_IMPL_LAYER_RENDER_PREPARATION_DATA_H
#define QSSG_RENDERER_IMPL_LAYER_RENDER_PREPARATION_DATA_H

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

#include <QtQuick3DRuntimeRender/private/qssgrendererimpllayerrenderhelper_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendershadercache_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderableobjects_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderclippingfrustum_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderresourcetexture2d_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendergpuprofiler_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendershadowmap_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderableobjects_p.h>

#include <QtQuick3DRuntimeRender/private/qssgrenderitem2d_p.h>

#define QSSG_RENDER_MINIMUM_RENDER_OPACITY .01f

QT_BEGIN_NAMESPACE
struct QSSGLayerRenderData;
class QSSGRendererImpl;
struct QSSGRenderableObject;

enum class QSSGLayerRenderPreparationResultFlag
{
    // Was the data in this layer dirty (meaning re-render to texture, possibly)
    WasLayerDataDirty = 1,
    // Was the data in this layer dirty *or* this layer *or* any effect dirty.
    WasDirty = 1 << 1,

    // Should create independent viewport
    // If we aren't rendering to texture we still may have width/height manipulations
    // that require our own viewport.
    ShouldCreateIndependentViewport = 1 << 2,

    RequiresDepthTexture = 1 << 3,

    // SSAO should be done in a separate pass
    // Note that having an AO pass necessitates a DepthTexture so this flag should
    // never be set without the RequiresDepthTexture flag as well.
    RequiresSsaoPass = 1 << 4,

    // if some light cause shadow
    // we need a separate per light shadow map pass
    RequiresShadowMapPass = 1 << 5,

    // This is the case when direct rendering, and need to clear an FBO, but not a Window
    RequiresTransparentClear = 1 << 6
};

struct QSSGLayerRenderPreparationResultFlags : public QFlags<QSSGLayerRenderPreparationResultFlag>
{
    bool wasLayerDataDirty() const
    {
        return this->operator&(QSSGLayerRenderPreparationResultFlag::WasLayerDataDirty);
    }
    void setLayerDataDirty(bool inValue)
    {
        setFlag(QSSGLayerRenderPreparationResultFlag::WasLayerDataDirty, inValue);
    }

    bool wasDirty() const { return this->operator&(QSSGLayerRenderPreparationResultFlag::WasDirty); }
    void setWasDirty(bool inValue) { setFlag(QSSGLayerRenderPreparationResultFlag::WasDirty, inValue); }

    bool shouldCreateIndependentViewport() const
    {
        return this->operator&(QSSGLayerRenderPreparationResultFlag::ShouldCreateIndependentViewport);
    }
    void setShouldCreateIndependentViewport(bool inValue)
    {
        setFlag(QSSGLayerRenderPreparationResultFlag::ShouldCreateIndependentViewport, inValue);
    }

    bool requiresDepthTexture() const
    {
        return this->operator&(QSSGLayerRenderPreparationResultFlag::RequiresDepthTexture);
    }
    void setRequiresDepthTexture(bool inValue)
    {
        setFlag(QSSGLayerRenderPreparationResultFlag::RequiresDepthTexture, inValue);
    }

    bool requiresSsaoPass() const { return this->operator&(QSSGLayerRenderPreparationResultFlag::RequiresSsaoPass); }
    void setRequiresSsaoPass(bool inValue)
    {
        setFlag(QSSGLayerRenderPreparationResultFlag::RequiresSsaoPass, inValue);
    }

    bool requiresShadowMapPass() const
    {
        return this->operator&(QSSGLayerRenderPreparationResultFlag::RequiresShadowMapPass);
    }
    void setRequiresShadowMapPass(bool inValue)
    {
        setFlag(QSSGLayerRenderPreparationResultFlag::RequiresShadowMapPass, inValue);
    }

    bool requiresTransparentClear() const
    {
        return this->operator&(QSSGLayerRenderPreparationResultFlag::RequiresTransparentClear);
    }

    void setRequiresTransparentClear(bool inValue)
    {
        setFlag(QSSGLayerRenderPreparationResultFlag::RequiresTransparentClear, inValue);
    }

};

struct QSSGLayerRenderPreparationResult : public QSSGLayerRenderHelper
{
    QSSGRenderEffect *lastEffect = nullptr;
    QSSGLayerRenderPreparationResultFlags flags;
    quint32 maxAAPassIndex = 0;
    QSSGLayerRenderPreparationResult() = default;
    QSSGLayerRenderPreparationResult(const QSSGLayerRenderHelper &inHelper)
        : QSSGLayerRenderHelper(inHelper), lastEffect(nullptr), maxAAPassIndex(0)
    {
    }
};

struct QSSGRenderableNodeEntry
{
    QSSGRenderNode *node = nullptr;
    QSSGNodeLightEntryList lights;
    QSSGRenderableNodeEntry() = default;
    QSSGRenderableNodeEntry(QSSGRenderNode &inNode) : node(&inNode) {}
};

struct QSSGScopedLightsListScope
{
    QVector<QSSGRenderLight *> &lightsList;
    QVector<QVector3D> &lightDirList;
    qint32 listOriginalSize;
    QSSGScopedLightsListScope(QVector<QSSGRenderLight *> &inLights,
                                QVector<QVector3D> &inDestLightDirList,
                                QVector<QVector3D> &inSrcLightDirList,
                                QSSGNodeLightEntryList &inScopedLights)
        : lightsList(inLights), lightDirList(inDestLightDirList), listOriginalSize(lightsList.size())
    {
        auto iter = inScopedLights.begin();
        const auto end = inScopedLights.end();
        while (iter != end) {
            lightsList.push_back(iter->light);
            lightDirList.push_back(inSrcLightDirList.at(iter->lightIndex));
            ++iter;
        }
    }
    ~QSSGScopedLightsListScope()
    {
        lightsList.resize(listOriginalSize);
        lightDirList.resize(listOriginalSize);
    }
};

struct QSSGDefaultMaterialPreparationResult
{
    QSSGRenderableImage *firstImage;
    float opacity;
    QSSGRenderableObjectFlags renderableFlags;
    QSSGShaderDefaultMaterialKey materialKey;
    bool dirty;

    QSSGDefaultMaterialPreparationResult(QSSGShaderDefaultMaterialKey inMaterialKey);
};

// Data used strictly in the render preparation step.
struct QSSGLayerRenderPreparationData
{
    typedef void (*TRenderRenderableFunction)(QSSGLayerRenderData &inData,
                                              QSSGRenderableObject &inObject,
                                              const QVector2D &inCameraProps,
                                              const ShaderFeatureSetList &inShaderFeatures,
                                              quint32 lightIndex,
                                              const QSSGRenderCamera &inCamera);
    typedef QHash<QSSGRenderLight *, QSSGRenderNode *> TLightToNodeMap;
    typedef QVector<QSSGModelContext *> TModelContextPtrList;
    typedef QVector<QSSGRenderableObjectHandle> TRenderableObjectList;

    // typedef Pool<SNodeLightEntry, ForwardingAllocator> TNodeLightEntryPoolType;

    enum Enum {
        MAX_AA_LEVELS = 8,
        MAX_TEMPORAL_AA_LEVELS = 2,
    };

    QSSGRenderLayer &layer;
    QSSGRef<QSSGRendererImpl> renderer;
    // List of nodes we can render, not all may be active.  Found by doing a depth-first
    // search through m_FirstChild if length is zero.

    // TNodeLightEntryPoolType m_RenderableNodeLightEntryPool;
    QVector<QSSGRenderableNodeEntry> renderableNodes;
    QVector<QSSGRenderableNodeEntry> renderableItem2Ds;
    QVector<QSSGRenderableNodeEntry> renderedItem2Ds;
    TLightToNodeMap lightToNodeMap; // map of lights to nodes to cache if we have looked up a
    // given scoped light yet.
    // Built at the same time as the renderable nodes map.
    // these are processed so they are available when the shaders for the models
    // are being generated.
    QVector<QSSGRenderCamera *> cameras;
    QVector<QSSGRenderLight *> lights;

    // Results of prepare for render.
    QSSGRenderCamera *camera;
    QVector<QSSGRenderLight *> globalLights; // Only contains lights that are global.
    TRenderableObjectList opaqueObjects;
    TRenderableObjectList transparentObjects;
    // Sorted lists of the rendered objects.  There may be other transforms applied so
    // it is simplest to duplicate the lists.
    TRenderableObjectList renderedOpaqueObjects;
    TRenderableObjectList renderedTransparentObjects;
    QMatrix4x4 viewProjection;
    QSSGOption<QSSGClippingFrustum> clippingFrustum;
    QSSGOption<QSSGLayerRenderPreparationResult> layerPrepResult;
    QSSGOption<QVector3D> cameraDirection;
    // Scoped lights need a level of indirection into a light direction list.  The source light
    // directions list is as long as there are lights on the layer.  It holds invalid
    // information for
    // any lights that are not both active and scoped; but the relative position for a given
    // light
    // in this list is completely constant and immutable; this relative position is saved on a
    // structure
    // and used when looking up the light direction for a given light.
    QVector<QVector3D> sourceLightDirections;
    QVector<QVector3D> lightDirections;
    TModelContextPtrList modelContexts;

    ShaderFeatureSetList features;
    bool featuresDirty;
    size_t featureSetHash;
    bool tooManyLightsError;

    // shadow mapps
    QSSGRef<QSSGRenderShadowMap> shadowMapManager;

    QSSGLayerRenderPreparationData(QSSGRenderLayer &inLayer, const QSSGRef<QSSGRendererImpl> &inRenderer);
    virtual ~QSSGLayerRenderPreparationData();
    void createShadowMapManager();

    static QByteArray cgLightingFeatureName();

    QSSGShaderDefaultMaterialKey generateLightingKey(QSSGRenderDefaultMaterial::MaterialLighting inLightingType, bool receivesShadows = true);

    void prepareImageForRender(QSSGRenderImage &inImage,
                               QSSGImageMapTypes inMapType,
                               QSSGRenderableImage *&ioFirstImage,
                               QSSGRenderableImage *&ioNextImage,
                               QSSGRenderableObjectFlags &ioFlags,
                               QSSGShaderDefaultMaterialKey &ioGeneratedShaderKey,
                               quint32 inImageIndex, QSSGRenderDefaultMaterial *inMaterial = nullptr);

    void setVertexInputPresence(const QSSGRenderableObjectFlags &renderableFlags,
                                QSSGShaderDefaultMaterialKey &key);

    QSSGDefaultMaterialPreparationResult prepareDefaultMaterialForRender(QSSGRenderDefaultMaterial &inMaterial,
                                                                           QSSGRenderableObjectFlags &inExistingFlags,
                                                                           float inOpacity);

    QSSGDefaultMaterialPreparationResult prepareCustomMaterialForRender(QSSGRenderCustomMaterial &inMaterial,
                                                                          QSSGRenderableObjectFlags &inExistingFlags,
                                                                          float inOpacity, bool alreadyDirty);

    bool prepareModelForRender(QSSGRenderModel &inModel,
                               const QMatrix4x4 &inViewProjection,
                               const QSSGOption<QSSGClippingFrustum> &inClipFrustum,
                               QSSGNodeLightEntryList &inScopedLights);

    // Helper function used during PRepareForRender and PrepareAndRender
    bool prepareRenderablesForRender(const QMatrix4x4 &inViewProjection,
                                     const QSSGOption<QSSGClippingFrustum> &inClipFrustum,
                                     QSSGLayerRenderPreparationResultFlags &ioFlags);

    // returns true if this object will render something different than it rendered the last
    // time.
    virtual void prepareForRender(const QSize &inViewportDimensions);
    bool checkLightProbeDirty(QSSGRenderImage &inLightProbe);
    void setShaderFeature(const char *inName, bool inValue);
    ShaderFeatureSetList getShaderFeatureSet();
    size_t getShaderFeatureSetHash();

    QVector3D getCameraDirection();
    // Per-frame cache of renderable objects post-sort.
    const QVector<QSSGRenderableObjectHandle> &getOpaqueRenderableObjects(bool performSort = true);
    // If layer depth test is false, this may also contain opaque objects.
    const QVector<QSSGRenderableObjectHandle> &getTransparentRenderableObjects();
    const QVector<QSSGRenderableNodeEntry> &getRenderableItem2Ds();

    virtual void resetForFrame();
};
QT_END_NAMESPACE
#endif
