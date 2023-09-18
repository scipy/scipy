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

#ifndef QT3DRENDER_QRENDERASPECT_P_H
#define QT3DRENDER_QRENDERASPECT_P_H

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

#include <Qt3DRender/qrenderaspect.h>
#include <Qt3DCore/private/qabstractaspect_p.h>
#include <Qt3DRender/private/qt3drender_global_p.h>
#include <Qt3DRender/private/expandboundingvolumejob_p.h>
#include <Qt3DRender/private/updateworldtransformjob_p.h>
#include <Qt3DRender/private/updateworldboundingvolumejob_p.h>
#include <Qt3DRender/private/calcboundingvolumejob_p.h>
#include <Qt3DRender/private/updateskinningpalettejob_p.h>
#include <Qt3DRender/private/updateentitylayersjob_p.h>
#include <Qt3DRender/private/updatetreeenabledjob_p.h>
#include <Qt3DRender/private/genericlambdajob_p.h>
#include <Qt3DRender/private/pickboundingvolumejob_p.h>
#include <Qt3DRender/private/raycastingjob_p.h>

#include <QtCore/qmutex.h>

QT_BEGIN_NAMESPACE

class QSurface;
class QScreen;

namespace Qt3DRender {

class QSceneImporter;

namespace Render {
class AbstractRenderer;
class NodeManagers;
class QRenderPlugin;
}

namespace Render {
class OffscreenSurfaceHelper;
class PickEventFilter;

using SynchronizerJobPtr = GenericLambdaJobPtr<std::function<void()>>;

class UpdateLevelOfDetailJob;
typedef QSharedPointer<UpdateLevelOfDetailJob> UpdateLevelOfDetailJobPtr;
}

class Q_3DRENDERSHARED_PRIVATE_EXPORT QRenderAspectPrivate : public Qt3DCore::QAbstractAspectPrivate
{
public:
    QRenderAspectPrivate(QRenderAspect::RenderType type);
    ~QRenderAspectPrivate();

    Q_DECLARE_PUBLIC(QRenderAspect)

    static QRenderAspectPrivate* findPrivate(Qt3DCore::QAspectEngine *engine);
    static QRenderAspectPrivate *get(QRenderAspect *q);

    void syncDirtyFrontEndNode(Qt3DCore::QNode *node, Qt3DCore::QBackendNode *backend, bool firstTime) const override;
    void jobsDone() override;
    void frameDone() override;

    void createNodeManagers();
    void onEngineStartup();

    void registerBackendTypes();
    void unregisterBackendTypes();
    void loadSceneParsers();
    void loadRenderPlugin(const QString &pluginName);
    void renderInitialize(QOpenGLContext *context);
    void renderSynchronous(bool swapBuffers = true);
    void renderShutdown();
    void registerBackendType(const QMetaObject &, const Qt3DCore::QBackendNodeMapperPtr &functor);
    QVector<Qt3DCore::QAspectJobPtr> createGeometryRendererJobs() const;
    QVector<Qt3DCore::QAspectJobPtr> createPreRendererJobs() const;
    QVector<Qt3DCore::QAspectJobPtr> createRenderBufferJobs() const;
    Render::AbstractRenderer *loadRendererPlugin();

    Render::NodeManagers *m_nodeManagers;
    Render::AbstractRenderer *m_renderer;

    bool m_initialized;
    bool m_renderAfterJobs;
    QList<QSceneImporter *> m_sceneImporter;
    QVector<QString> m_loadedPlugins;
    QVector<Render::QRenderPlugin *> m_renderPlugins;
    QRenderAspect::RenderType m_renderType;
    Render::OffscreenSurfaceHelper *m_offscreenHelper;
    QScreen *m_screen = nullptr;

    Render::UpdateTreeEnabledJobPtr m_updateTreeEnabledJob;
    Render::UpdateWorldTransformJobPtr m_worldTransformJob;
    Render::ExpandBoundingVolumeJobPtr m_expandBoundingVolumeJob;
    Render::CalculateBoundingVolumeJobPtr m_calculateBoundingVolumeJob;
    Render::UpdateWorldBoundingVolumeJobPtr m_updateWorldBoundingVolumeJob;
    Render::UpdateSkinningPaletteJobPtr m_updateSkinningPaletteJob;
    Render::UpdateLevelOfDetailJobPtr m_updateLevelOfDetailJob;
    Render::UpdateEntityLayersJobPtr m_updateEntityLayersJob;
    Render::SynchronizerJobPtr m_syncLoadingJobs;
    Render::PickBoundingVolumeJobPtr m_pickBoundingVolumeJob;
    Render::RayCastingJobPtr m_rayCastingJob;

    QScopedPointer<Render::PickEventFilter> m_pickEventFilter;

    static QMutex m_pluginLock;
    static QVector<QString> m_pluginConfig;
    static QVector<QRenderAspectPrivate *> m_instances;
    static void configurePlugin(const QString &plugin);
};

}

QT_END_NAMESPACE

#endif // QT3DRENDER_QRENDERASPECT_P_H
