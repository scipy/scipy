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
#ifndef QT3DRENDER_RENDER_ABSTRACTRENDERER_P_H
#define QT3DRENDER_RENDER_ABSTRACTRENDERER_P_H

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

#include <QtCore/qflags.h>
#include <QtCore/qmutex.h>
#include <Qt3DRender/private/qt3drender_global_p.h>
#include <Qt3DRender/private/handle_types_p.h>
#include <Qt3DRender/qrenderapi.h>
#include <Qt3DCore/qaspectjob.h>
#include <Qt3DCore/qnodeid.h>
#include <QtGui/qsurfaceformat.h>

#include <qopenglcontext.h>

QT_BEGIN_NAMESPACE

class QSurface;
class QSize;
class QScreen;
class QOpenGLTexture;
class QMouseEvent;
class QKeyEvent;

namespace Qt3DCore {
class QAbstractFrameAdvanceService;
class QEventFilterService;
class QAbstractAspectJobManager;
class QServiceLocator;
class QAspectManager;
}

namespace Qt3DRender {

class QRenderAspect;
struct GraphicsApiFilterData;

namespace Render {

class NodeManagers;
class Entity;
class FrameGraphNode;
class RenderSettings;
class BackendNode;
class OffscreenSurfaceHelper;
class Shader;
class RenderBackendResourceAccessor;
class Q_3DRENDERSHARED_PRIVATE_EXPORT AbstractRenderer
{
public:
    virtual ~AbstractRenderer() {}

    // Changes made to backend nodes are reported to the Renderer
    enum BackendNodeDirtyFlag {
        TransformDirty      = 1 << 0,
        MaterialDirty       = 1 << 1,
        GeometryDirty       = 1 << 2,
        ComputeDirty        = 1 << 3,
        ParameterDirty      = 1 << 4,
        FrameGraphDirty     = 1 << 5,
        EntityEnabledDirty  = 1 << 6,
        BuffersDirty        = 1 << 7,
        TexturesDirty       = 1 << 8,
        ShadersDirty        = 1 << 9,
        SkeletonDataDirty   = 1 << 10,
        JointDirty          = 1 << 11,
        LayersDirty         = 1 << 12,
        TechniquesDirty     = 1 << 13,
        LightsDirty         = 1 << 15,
        AllDirty            = 0xffffff
    };
    Q_DECLARE_FLAGS(BackendNodeDirtySet, BackendNodeDirtyFlag)

    virtual void dumpInfo() const = 0;

    virtual API api() const = 0;

    virtual qint64 time() const = 0;
    virtual void setTime(qint64 time) = 0;
    virtual void setJobsInLastFrame(int jobsInLastFrame) = 0;

    virtual void setAspect(QRenderAspect *aspect) = 0;
    virtual void setNodeManagers(NodeManagers *managers) = 0;
    virtual void setServices(Qt3DCore::QServiceLocator *services) = 0;
    virtual void setSurfaceExposed(bool exposed) = 0;

    virtual QRenderAspect *aspect() const = 0;
    virtual NodeManagers *nodeManagers() const = 0;
    virtual Qt3DCore::QServiceLocator *services() const = 0;

    virtual void initialize() = 0;
    virtual void shutdown() = 0;
    virtual void releaseGraphicsResources() = 0;

    // Threaded renderer
    virtual void render() = 0;
    // Synchronous renderer
    virtual void doRender(bool swapBuffers) = 0;

    virtual void cleanGraphicsResources() = 0;

    virtual bool isRunning() const = 0;

    virtual void markDirty(BackendNodeDirtySet changes, BackendNode *node) = 0;
    virtual BackendNodeDirtySet dirtyBits() = 0;
#if defined(QT_BUILD_INTERNAL)
    virtual void clearDirtyBits(BackendNodeDirtySet changes) = 0;
#endif
    virtual bool shouldRender() const = 0;
    virtual void skipNextFrame() = 0;
    virtual void jobsDone(Qt3DCore::QAspectManager *manager) = 0;

    virtual void setPendingEvents(const QList<QPair<QObject *, QMouseEvent>> &mouseEvents, const QList<QKeyEvent> &keyEvents) = 0;

    virtual QVector<Qt3DCore::QAspectJobPtr> preRenderingJobs() = 0;
    virtual QVector<Qt3DCore::QAspectJobPtr> renderBinJobs() = 0;

    virtual void setSceneRoot(Entity *root) = 0;

    virtual Entity *sceneRoot() const = 0;
    virtual FrameGraphNode *frameGraphRoot() const = 0;

    virtual Qt3DCore::QAbstractFrameAdvanceService *frameAdvanceService() const = 0;

    virtual void setSettings(RenderSettings *settings) = 0;
    virtual RenderSettings *settings() const = 0;

    virtual QVariant executeCommand(const QStringList &args) = 0;

    // For QtQuick rendering (Scene2D)
    virtual void setOpenGLContext(QOpenGLContext *ctx) = 0;
    virtual void setScreen(QScreen *) {}
    virtual QScreen *screen() const { return nullptr; }
    virtual bool accessOpenGLTexture(Qt3DCore::QNodeId nodeId, QOpenGLTexture **texture, QMutex **lock, bool readonly) = 0;
    virtual QSharedPointer<RenderBackendResourceAccessor> resourceAccessor() const = 0;

    virtual void setOffscreenSurfaceHelper(OffscreenSurfaceHelper *helper) = 0;
    virtual QSurfaceFormat format() = 0;
    virtual QOpenGLContext *shareContext() const = 0;
    virtual const GraphicsApiFilterData *contextInfo() const = 0;

    // These commands are executed in a dedicated command thread
    // More will be added later
    virtual void loadShader(Shader *shader, Qt3DRender::Render::HShader shaderHandle) = 0;

    // Runtime Cache for Generated Shader Graph
    bool containsGeneratedShaderGraph(const QByteArray &key) const { return m_cachedGeneratedShaderGraphes.contains(key); };
    QByteArray cachedGeneratedShaderGraph(const QByteArray &key) const { return m_cachedGeneratedShaderGraphes.value(key); };
    void insertGeneratedShaderGraph(const QByteArray &key, const QByteArray shaderCode) { m_cachedGeneratedShaderGraphes.insert(key, shaderCode); }
    void removeGeneratedShaderGraph(const QByteArray &key) { m_cachedGeneratedShaderGraphes.remove(key); };

private:
    QHash<QByteArray, QByteArray> m_cachedGeneratedShaderGraphes;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(AbstractRenderer::BackendNodeDirtySet)

} // Render

} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_ABSTRACTRENDERER_P_H

