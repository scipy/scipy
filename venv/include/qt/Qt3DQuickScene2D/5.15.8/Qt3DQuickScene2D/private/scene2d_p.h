/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DRENDER_RENDER_QUICK3DSCENE2D_SCENE2D_P_H
#define QT3DRENDER_RENDER_QUICK3DSCENE2D_SCENE2D_P_H

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

#include <Qt3DCore/qnodeid.h>
#include <Qt3DCore/qpropertyupdatedchange.h>
#include <Qt3DRender/qpickevent.h>
#include <Qt3DQuickScene2D/qscene2d.h>

#include <private/qscene2d_p.h>
#include <private/qrendertargetoutput_p.h>
#include <private/backendnode_p.h>
#include <private/attachmentpack_p.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

namespace Render {

class GraphicsContext;

namespace Quick {

class Scene2D;

class RenderQmlEventHandler : public QObject
{
    Q_OBJECT
public:
    RenderQmlEventHandler(Scene2D *node);
    bool event(QEvent *e) override;

private:
    Scene2D *m_node;
};

class Q_3DQUICKSCENE2DSHARED_EXPORT Scene2D : public Qt3DRender::Render::BackendNode
{
public:
    Scene2D();
    ~Scene2D();

    void render();
    void initializeRender();
    void setSharedObject(Qt3DRender::Quick::Scene2DSharedObjectPtr sharedObject);
    void cleanup();
    void setOutput(Qt3DCore::QNodeId outputId);
    void initializeSharedObject();

    void syncFromFrontEnd(const Qt3DCore::QNode *frontEnd, bool firstTime) override;

    bool updateFbo(QOpenGLTexture *texture);
    void syncRenderControl();
    bool registerObjectPickerEvents(Qt3DCore::QEntity *qentity);
    void unregisterObjectPickerEvents(Qt3DCore::QNodeId entityId);
    void handlePickEvent(int type, const QPickEvent *ev);

    QOpenGLContext *m_context;
    QOpenGLContext *m_shareContext;
    QThread *m_renderThread;
    Qt3DCore::QNodeId m_outputId;
    QSharedPointer<Qt3DRender::Quick::Scene2DSharedObject> m_sharedObject;
    Qt3DCore::QNodeId m_peerId;
    Qt3DRender::Render::Attachment m_attachmentData;

    GLuint m_fbo;
    GLuint m_rbo;
    QSize m_textureSize;

    bool m_initialized;
    bool m_renderInitialized;
    bool m_mouseEnabled;
    Qt3DRender::Quick::QScene2D::RenderPolicy m_renderPolicy;
    QVector<Qt3DCore::QNodeId> m_entities;
    Qt3DRender::QPickEventPtr m_cachedPickEvent;
#ifdef QT_OPENGL_ES_2_ANGLE
    bool m_usingAngle;
#endif
    QVector<QMetaObject::Connection> m_connections;
};

} // Quick
} // Render
} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_QUICK3DSCENE2D_SCENE2D_P_H
