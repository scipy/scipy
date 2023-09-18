/****************************************************************************
**
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

#ifndef QSSGSCENERENDERER_H
#define QSSGSCENERENDERER_H

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

#include <QtQuick3DRender/private/qssgrendercontext_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrendercontextcore_p.h>

#include <qsgtextureprovider.h>
#include <qsgrendernode.h>
#include <QSGSimpleTextureNode>

#include <QtQuick3D/private/qquick3dviewport_p.h>

QT_BEGIN_NAMESPACE


class QQuick3DSceneManager;
class QQuick3DViewport;
struct QSSGRenderLayer;

class QQuick3DSceneRenderer
{
public:
    struct FramebufferObject {
        FramebufferObject(const QSize &s, const QSSGRef<QSSGRenderContext> &context,
                          int msaaSamples = 1);
        ~FramebufferObject();
        QSize size;
        QSSGRef<QSSGRenderContext> renderContext;
        QSSGRef<QSSGRenderFrameBuffer> fbo;
        QSSGRef<QSSGRenderTexture2D> color0;
        QSSGRef<QSSGRenderTexture2D> depthStencil;
        int samples;
    };

    QQuick3DSceneRenderer(QWindow *window);
    ~QQuick3DSceneRenderer();
protected:
    GLuint render();
    void render(const QRect &viewport, bool clearFirst = false);
    void synchronize(QQuick3DViewport *item, const QSize &size, bool useFBO = true);
    void update();
    void invalidateFramebufferObject();
    void onRenderModeChanged();
    QSize surfaceSize() const { return m_surfaceSize; }
    QSSGRenderPickResult pick(const QPointF &pos);
    QSSGRenderPickResult syncPick(const QPointF &pos);
    QQuick3DRenderStats *renderStats();

private:
    void updateLayerNode(QQuick3DViewport *view3D);
    void addNodeToLayer(QSSGRenderNode *node);
    void removeNodeFromLayer(QSSGRenderNode *node);
    void deleteResources();
    QSSGRef<QSSGRenderContextInterface> m_sgContext;
    QSharedPointer<QQuick3DSceneManager> m_sceneManager;
    QSSGRenderLayer *m_layer = nullptr;
    QSize m_surfaceSize;
    void *data = nullptr;
    bool m_layerSizeIsDirty = true;
    bool m_aaIsDirty = true;
    QWindow *m_window = nullptr;
    FramebufferObject *m_antialiasingFbo = nullptr;
    FramebufferObject *m_fbo = nullptr;
    QQuick3DRenderStats *m_renderStats = nullptr;

    QSSGRenderNode *m_sceneRootNode = nullptr;
    QSSGRenderNode *m_importRootNode = nullptr;

    float m_ssaaMultiplier = 1.5f;

    friend class SGFramebufferObjectNode;
    friend class QQuick3DSGRenderNode;
    friend class QQuick3DSGDirectRenderer;
    friend class QQuick3DViewport;
};

class QOpenGLVertexArrayObjectHelper;

class SGFramebufferObjectNode final : public QSGTextureProvider, public QSGSimpleTextureNode
{
    Q_OBJECT

public:
    SGFramebufferObjectNode();
    ~SGFramebufferObjectNode() override;

    void scheduleRender();

    QSGTexture *texture() const override;

    void preprocess() override;

public Q_SLOTS:
    void render();

    void handleScreenChange();

public:
    QQuickWindow *window;
    QQuick3DSceneRenderer *renderer;
    QQuick3DViewport *quickFbo;

    bool renderPending;
    bool invalidatePending;

    qreal devicePixelRatio;
    int requestedFramesCount;
};

class QQuick3DSGRenderNode final : public QSGRenderNode
{
public:

    StateFlags changedStates() const override;
    void render(const RenderState *state) override;
    void releaseResources() override;
    RenderingFlags flags() const override;
public:
    QQuickWindow *window = nullptr;
    QQuick3DSceneRenderer *renderer = nullptr;
};

class QQuick3DSGDirectRenderer : public QObject
{
    Q_OBJECT
public:
    enum QQuick3DSGDirectRendererMode {
        Underlay,
        Overlay
    };
    QQuick3DSGDirectRenderer(QQuick3DSceneRenderer *renderer, QQuickWindow *window, QQuick3DSGDirectRendererMode mode = Underlay);
    ~QQuick3DSGDirectRenderer();

    QQuick3DSceneRenderer *renderer() { return m_renderer; }
    void onRenderModeChanged() { m_renderer = nullptr; }
    void setViewport(const QRectF &viewport);

    void requestRender();
    void setVisibility(bool visible);

private Q_SLOTS:
    void render();

private:
    QQuick3DSceneRenderer *m_renderer = nullptr;
    QQuickWindow *m_window = nullptr;
    QQuick3DSGDirectRendererMode m_mode;
    QRectF m_viewport;
    bool m_isVisible = true;
};

QT_END_NAMESPACE

#endif // QSSGSCENERENDERER_H
