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

#ifndef QSSGVIEW3D_H
#define QSSGVIEW3D_H

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

#include <QtGui/QOpenGLFramebufferObject>
#include <QtQuick/QQuickItem>
#include <QtCore/qurl.h>

#include <QtQuick3D/qtquick3dglobal.h>
#include <QtQuick3D/private/qquick3dpickresult_p.h>

#include <QtQuick3DRender/private/qssgrenderframebuffer_p.h>

QT_BEGIN_NAMESPACE

class QSSGView3DPrivate;
class QQuick3DCamera;
class QQuick3DSceneEnvironment;
class QQuick3DNode;
class QQuick3DSceneRootNode;
class QQuick3DSceneRenderer;
class QQuick3DRenderStats;

class SGFramebufferObjectNode;
class QQuick3DSGRenderNode;
class QQuick3DSGDirectRenderer;

class Q_QUICK3D_EXPORT QQuick3DViewport : public QQuickItem
{
    Q_OBJECT
    Q_PROPERTY(QQmlListProperty<QObject> data READ data DESIGNABLE false FINAL)
    Q_PROPERTY(QQuick3DCamera *camera READ camera WRITE setCamera NOTIFY cameraChanged FINAL)
    Q_PROPERTY(QQuick3DSceneEnvironment *environment READ environment WRITE setEnvironment NOTIFY environmentChanged FINAL)
    Q_PROPERTY(QQuick3DNode *scene READ scene NOTIFY sceneChanged)
    Q_PROPERTY(QQuick3DNode *importScene READ importScene WRITE setImportScene NOTIFY importSceneChanged FINAL)
    Q_PROPERTY(RenderMode renderMode READ renderMode WRITE setRenderMode NOTIFY renderModeChanged FINAL)
    Q_PROPERTY(QQuick3DRenderStats *renderStats READ renderStats CONSTANT)
    Q_CLASSINFO("DefaultProperty", "data")
public:
    enum RenderMode {
        Offscreen,
        Underlay,
        Overlay,
        Inline
    };
    Q_ENUM(RenderMode)

    explicit QQuick3DViewport(QQuickItem *parent = nullptr);
    ~QQuick3DViewport() override;

    QQmlListProperty<QObject> data();

    QQuick3DCamera *camera() const;
    QQuick3DSceneEnvironment *environment() const;
    QQuick3DNode *scene() const;
    QQuick3DNode *importScene() const;
    RenderMode renderMode() const;
    QQuick3DRenderStats *renderStats() const;

    QQuick3DSceneRenderer *createRenderer() const;

    bool isTextureProvider() const override;
    QSGTextureProvider *textureProvider() const override;
    void releaseResources() override;

    Q_INVOKABLE QVector3D mapFrom3DScene(const QVector3D &scenePos) const;
    Q_INVOKABLE QVector3D mapTo3DScene(const QVector3D &viewPos) const;

    Q_INVOKABLE QQuick3DPickResult pick(float x, float y) const;

    void setShaderCacheFile(const QUrl &shaderCacheFile);
    QUrl shaderCacheFile();
    void setShaderCache(const QByteArray &shaderCache);
    void exportShaderCache(const QUrl &shaderCacheFile, bool binaryShaders, int compressionLevel = -1);
    void exportShaderCache(bool binaryShaders, int compressionLevel = -1);
    QByteArray shaderCacheData() const;

protected:
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;
    QSGNode *updatePaintNode(QSGNode *, UpdatePaintNodeData *) override;
    void itemChange(QQuickItem::ItemChange change, const QQuickItem::ItemChangeData &value) override;

public Q_SLOTS:
    void setCamera(QQuick3DCamera *camera);
    void setEnvironment(QQuick3DSceneEnvironment * environment);
    void setImportScene(QQuick3DNode *inScene);
    void setRenderMode(RenderMode renderMode);
    void cleanupDirectRenderer();

private Q_SLOTS:
    void invalidateSceneGraph();

Q_SIGNALS:
    void cameraChanged();
    void environmentChanged();
    void sceneChanged();
    void importSceneChanged();
    void renderModeChanged();
    void shaderCacheLoadErrors(const QByteArray &errors);
    void shaderCacheExported(bool success);

private:
    Q_DISABLE_COPY(QQuick3DViewport)
    QQuick3DSceneRenderer *getRenderer() const;
    void updateDynamicTextures();
    void setupDirectRenderer(RenderMode mode);
    void updateClearBeforeRendering();
    bool checkIsVisible() const;

    void readShaderCache();
    void writeShaderCache(const QUrl &shaderCacheFile);
    void doExportShaderCache();
    void doImportShaderCache();

    QQuick3DCamera *m_camera = nullptr;
    QQuick3DSceneEnvironment *m_environment = nullptr;
    QQuick3DSceneRootNode *m_sceneRoot = nullptr;
    QQuick3DNode *m_importScene = nullptr;
    mutable SGFramebufferObjectNode *m_node = nullptr;
    mutable QQuick3DSGRenderNode *m_renderNode = nullptr;
    mutable QQuick3DSGDirectRenderer *m_directRenderer = nullptr;
    bool m_renderModeDirty = false;
    RenderMode m_renderMode = Offscreen;
    QQuick3DRenderStats *m_renderStats = nullptr;
    QUrl m_shaderCacheFile;
    QByteArray m_shaderCacheData;
    QByteArray m_shaderCacheIO;
    QUrl m_exportShaderCacheFile;
    bool m_exportShaderCacheRequested = false;
    bool m_binaryShaders = false;
    bool m_fileExport = false;
    int m_compressionLevel = -1;
    QHash<QObject*, QMetaObject::Connection> m_connections;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QQuick3DViewport)

#endif // QSSGVIEW3D_H
