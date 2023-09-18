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

#ifndef QSSG_RENDER_CONTEXT_CORE_H
#define QSSG_RENDER_CONTEXT_CORE_H

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

#include <QtQuick3DRuntimeRender/private/qssgrenderinputstreamfactory_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderthreadpool_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderdynamicobjectsystem_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendercustommaterialsystem_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendereffectsystem_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderimagebatchloader_p.h>
#include <QtQuick3DRuntimeRender/private/qtquick3druntimerenderglobal_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderinputstreamfactory_p.h>
#include <QtQuick3DRuntimeRender/private/qssgperframeallocator_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendershadercache_p.h>

#include <QtQuick3DUtils/private/qssgperftimer_p.h>

#include <QtCore/QPair>
#include <QtCore/QSize>

#include <QtGui/qcolor.h>

QT_BEGIN_NAMESPACE

enum class ScaleModes
{
    ExactSize = 0, // Ensure the viewport is exactly same size as application
    ScaleToFit = 1, // Resize viewport keeping aspect ratio
    ScaleToFill = 2, // Resize viewport to entire window
    FitSelected = 3, // Resize presentation to fit into viewport
};

class QSSGMaterialSystem;
class QSSGRendererInterface;

class Q_QUICK3DRUNTIMERENDER_EXPORT QSSGRenderContextInterface
{
    Q_DISABLE_COPY(QSSGRenderContextInterface)
public:
    QAtomicInt ref;
private:
    const QSSGRef<QSSGRenderContext> m_renderContext;
    QSSGPerfTimer m_perfTimer;

    const QSSGRef<QSSGInputStreamFactory> m_inputStreamFactory;
    const QSSGRef<QSSGBufferManager> m_bufferManager;
    const QSSGRef<QSSGResourceManager> m_resourceManager;
    const QSSGRef<QSSGRendererInterface> m_renderer;
    const QSSGRef<QSSGDynamicObjectSystem> m_dynamicObjectSystem;
    const QSSGRef<QSSGEffectSystem> m_effectSystem;
    const QSSGRef<QSSGShaderCache> m_shaderCache;
    const QSSGRef<QSSGAbstractThreadPool> m_threadPool;
    const QSSGRef<IImageBatchLoader> m_imageBatchLoader;
    const QSSGRef<QSSGMaterialSystem> m_customMaterialSystem;
    const QSSGRef<QSSGShaderProgramGeneratorInterface> m_shaderProgramGenerator;
    const QSSGRef<QSSGDefaultMaterialShaderGeneratorInterface> m_defaultMaterialShaderGenerator;
    const QSSGRef<QSSGMaterialShaderGeneratorInterface> m_customMaterialShaderGenerator;
    QSSGPerFrameAllocator m_perFrameAllocator;
    quint32 m_frameCount = 0;
    // Viewport that this render context should use
    QRect m_viewport;
    QRect m_scissorRect;
    QSize m_windowDimensions {800, 480};
    bool m_wireframeMode = false;
    QColor m_sceneColor;
    QSSGRef<QSSGRenderFrameBuffer> m_contextRenderTarget;
    QPair<float, int> m_fps = qMakePair(0.0f, 0);

    QSSGRenderContextInterface(const QSSGRef<QSSGRenderContext> &ctx, const QString &inApplicationDirectory);

public:
    static QSSGRef<QSSGRenderContextInterface> getRenderContextInterface(const QSSGRef<QSSGRenderContext> &ctx,
                                                                         const QString &inApplicationDirectory,
                                                                         quintptr wid);
    static QSSGRef<QSSGRenderContextInterface> getRenderContextInterface(quintptr wid);

    ~QSSGRenderContextInterface();
    const QSSGRef<QSSGRendererInterface> &renderer() const;
    const QSSGRef<QSSGBufferManager> &bufferManager() const;
    const QSSGRef<QSSGResourceManager> &resourceManager() const;
    const QSSGRef<QSSGRenderContext> &renderContext() const;
    const QSSGRef<QSSGInputStreamFactory> &inputStreamFactory() const;
    const QSSGRef<QSSGEffectSystem> &effectSystem() const;
    const QSSGRef<QSSGShaderCache> &shaderCache() const;
    const QSSGRef<QSSGAbstractThreadPool> &threadPool() const;
    const QSSGRef<IImageBatchLoader> &imageBatchLoader() const;
    const QSSGRef<QSSGDynamicObjectSystem> &dynamicObjectSystem() const;
    const QSSGRef<QSSGMaterialSystem> &customMaterialSystem() const;
    QSSGPerfTimer *performanceTimer() { return &m_perfTimer; }
    const QSSGRef<QSSGShaderProgramGeneratorInterface> &shaderProgramGenerator() const;
    const QSSGRef<QSSGDefaultMaterialShaderGeneratorInterface> &defaultMaterialShaderGenerator() const;
    const QSSGRef<QSSGMaterialShaderGeneratorInterface> &customMaterialShaderGenerator() const;
    // The memory used for the per frame allocator is released as the first step in BeginFrame.
    // This is useful for short lived objects and datastructures.
    QSSGPerFrameAllocator &perFrameAllocator() { return m_perFrameAllocator; }

    // Get the number of times EndFrame has been called
    quint32 frameCount() { return m_frameCount; }

    // Get fps
    QPair<float, int> getFPS() { return m_fps; }
    // Set fps by higher level, etc application
    void setFPS(QPair<float, int> inFPS) { m_fps = inFPS; }

    void setSceneColor(const QColor &inSceneColor) { m_sceneColor = inSceneColor; }

    // render Gpu profiler values
    void dumpGpuProfilerStats();

    // The reason you can set both window dimensions and an overall viewport is that the mouse
    // needs to be inverted
    // which requires the window height, and then the rest of the system really requires the
    // viewport.
    void setWindowDimensions(const QSize &inWindowDimensions) { m_windowDimensions = inWindowDimensions; }
    QSize windowDimensions() { return m_windowDimensions; }

    // In addition to the window dimensions which really have to be set, you can optionally
    // set the viewport which will force the entire viewer to render specifically to this
    // viewport.
    void setViewport(QRect inViewport) { m_viewport = inViewport; }
    QRect viewport() const { return m_viewport; }
    QRect contextViewport() const;

    void setScissorRect(QRect inScissorRect) { m_scissorRect = inScissorRect; }
    QRect scissorRect() const { return m_scissorRect; }

    void setWireframeMode(bool inEnable) { m_wireframeMode = inEnable; }
    bool wireframeMode() { return m_wireframeMode; }


    QVector2D mousePickViewport() const;
    QVector2D mousePickMouseCoords(const QVector2D &inMouseCoords) const;

    // Steps needed to render:
    // 1.  BeginFrame - sets up new target in render graph
    // 2.  Add everything you need to the render graph
    // 3.  RunRenderGraph - runs the graph, rendering things to main render target
    // 4.  Render any additional stuff to main render target on top of previously rendered
    // information
    // 5.  EndFrame

    // Clients need to call this every frame in order for various subsystems to release
    // temporary per-frame allocated objects.
    // Also sets up the viewport according to SetViewportInfo
    // and the topmost presentation dimensions.  Expects there to be exactly one presentation
    // dimension pushed at this point.
    // This also starts a render target in the render graph.
    void beginFrame();

    bool prepareLayerForRender(QSSGRenderLayer &inLayer);

    void renderLayer(QSSGRenderLayer &inLayer, bool needsClear);

    // Now you can render to the main render target if you want to render over the top
    // of everything.
    // Next call end frame.
    void endFrame();
};
QT_END_NAMESPACE

#endif
