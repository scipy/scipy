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

#ifndef QSSG_RENDERER_IMPL_LAYER_RENDER_DATA_H
#define QSSG_RENDERER_IMPL_LAYER_RENDER_DATA_H

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

#include <QtQuick3DRuntimeRender/private/qssgrendererimpllayerrenderpreparationdata_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderresourcebufferobjects_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderresourcetexture2d_p.h>

QT_BEGIN_NAMESPACE

struct QSSGLayerRenderData : public QSSGLayerRenderPreparationData
{
    QAtomicInt ref;

    // Layers can be rendered offscreen for many reasons; effects, progressive aa,
    // or just because a flag forces it.  If they are rendered offscreen we can then
    // cache the result so we don't render the layer again if it isn't dirty.
    QSSGResourceTexture2D m_layerTexture;
    QSSGResourceTexture2D m_temporalAATexture;
    QSSGResourceTexture2D m_prevTemporalAATexture;
    // Sometimes we need to render our depth buffer to a depth texture.
    QSSGResourceTexture2D m_layerDepthTexture;
    QSSGResourceTexture2D m_layerPrepassDepthTexture;
    QSSGResourceTexture2D m_layerSsaoTexture;
    // if we render multisampled we need resolve buffers
    QSSGResourceTexture2D m_layerMultisampleTexture;
    QSSGResourceTexture2D m_layerMultisamplePrepassDepthTexture;
    QSSGResourceTexture2D m_layerMultisampleWidgetTexture;

    // GPU profiler per layer
    QScopedPointer<QSSGRenderGPUProfiler> m_layerProfilerGpu;

    QSSGRenderCamera m_sceneCamera;
    QVector2D m_sceneDimensions;

    // ProgressiveAA algorithm details.
    quint32 m_progressiveAAPassIndex;
    // Increments every frame regardless to provide appropriate jittering
    quint32 m_temporalAAPassIndex;
    // Ensures we don't stop on an in-between frame; we will run two frames after the dirty flag
    // is clear.
    quint32 m_nonDirtyTemporalAAPassIndex;
    float m_textScale;

    QSSGOption<QVector3D> m_boundingRectColor;
    QSSGRenderTextureFormat m_depthBufferFormat;

    QSize m_previousDimensions;

    QSSGLayerRenderData(QSSGRenderLayer &inLayer, const QSSGRef<QSSGRendererImpl> &inRenderer);

    virtual ~QSSGLayerRenderData() override;

    void prepareForRender();

    // Internal Call
    void prepareForRender(const QSize &inViewportDimensions) override;

    QSSGRenderTextureFormat getDepthBufferFormat();
    QSSGRenderFrameBufferAttachment getFramebufferDepthAttachmentFormat(QSSGRenderTextureFormat depthFormat);

    // Render this layer assuming viewport and RT are setup.  Just renders exactly this item
    // no effects.
    void renderClearPass();
    void renderDepthPass(bool inEnableTransparentDepthWrite = false);
    void renderAoPass();
#ifdef QT_QUICK3D_DEBUG_SHADOWS
    void renderDebugDepthMap(QSSGRenderTexture2D *theDepthTex, QSSGRenderTextureCube *theDepthCube);
#endif
    void renderShadowMapPass(QSSGResourceFrameBuffer *theFB);
    void renderShadowCubeBlurPass(QSSGResourceFrameBuffer *theFB,
                                  const QSSGRef<QSSGRenderTextureCube> &target0,
                                  const QSSGRef<QSSGRenderTextureCube> &target1,
                                  float filterSz,
                                  float clipFar);
    void renderShadowMapBlurPass(QSSGResourceFrameBuffer *theFB,
                                 const QSSGRef<QSSGRenderTexture2D> &target0,
                                 const QSSGRef<QSSGRenderTexture2D> &target1,
                                 float filterSz,
                                 float clipFar);

    void render(QSSGResourceFrameBuffer *theFB = nullptr);
    void resetForFrame() override;

    void createGpuProfiler();
    void startProfiling(QString &nameID, bool sync);
    void endProfiling(QString &nameID);
    void startProfiling(const char *nameID, bool sync);
    void endProfiling(const char *nameID);
    void addVertexCount(quint32 count);

    void applyLayerPostEffects(const QSSGRef<QSSGRenderFrameBuffer> &theFB);
    void runnableRenderToViewport(const QSSGRef<QSSGRenderFrameBuffer> &theFB);

    // test method to render this layer to a given view projection without running the entire
    // layer setup system.  This assumes the client has setup the viewport, scissor, and render
    // target
    // the way they want them.
    void prepareAndRender(const QMatrix4x4 &inViewProjection);

    bool progressiveAARenderRequest() const;

protected:
    // Used for both the normal passes and the depth pass.
    // When doing the depth pass, we disable blending completely because it does not really make
    // sense
    // to write blend equations into
    void runRenderPass(TRenderRenderableFunction renderFn,
                       bool inEnableBlending,
                       bool inEnableDepthWrite,
                       bool inEnableTransparentDepthWrite,
                       bool inSortOpaqueRenderables,
                       quint32 indexLight,
                       const QSSGRenderCamera &inCamera,
                       QSSGResourceFrameBuffer *theFB = nullptr);
};
QT_END_NAMESPACE
#endif
