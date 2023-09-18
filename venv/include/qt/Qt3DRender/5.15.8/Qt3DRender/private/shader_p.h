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

#ifndef QT3DRENDER_RENDER_SHADER_H
#define QT3DRENDER_RENDER_SHADER_H

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

#include <Qt3DRender/private/backendnode_p.h>
#include <Qt3DRender/qshaderprogram.h>
#include <QVector>

QT_BEGIN_NAMESPACE

class QOpenGLShaderProgram;

namespace Qt3DRender {

namespace Render {

class ShaderManager;
class AttachmentPack;

class Q_3DRENDERSHARED_PRIVATE_EXPORT Shader : public BackendNode
{
public:
    static const int modelMatrixNameId;
    static const int viewMatrixNameId;
    static const int projectionMatrixNameId;
    static const int modelViewMatrixNameId;
    static const int viewProjectionMatrixNameId;
    static const int modelViewProjectionNameId;
    static const int mvpNameId;
    static const int inverseModelMatrixNameId;
    static const int inverseViewMatrixNameId;
    static const int inverseProjectionMatrixNameId;
    static const int inverseModelViewNameId;
    static const int inverseViewProjectionMatrixNameId;
    static const int inverseModelViewProjectionNameId;
    static const int modelNormalMatrixNameId;
    static const int modelViewNormalNameId;
    static const int viewportMatrixNameId;
    static const int inverseViewportMatrixNameId;
    static const int textureTransformMatrixNameId;
    static const int aspectRatioNameId;
    static const int exposureNameId;
    static const int gammaNameId;
    static const int timeNameId;
    static const int eyePositionNameId;
    static const int skinningPaletteNameId;

    Shader();
    ~Shader();

    void cleanup();

    QVector<QByteArray> shaderCode() const;
    void setShaderCode(QShaderProgram::ShaderType type, const QByteArray &code);

    void syncFromFrontEnd(const Qt3DCore::QNode *frontEnd, bool firstTime) override;

    inline QString log() const { return m_log; }
    inline QShaderProgram::Status status() const { return m_status; }

    void setFormat(QShaderProgram::Format format);
    QShaderProgram::Format format() const { return m_format; }
    bool isDirty() const { return m_dirty; }
    void unsetDirty() { m_dirty = false; }

    inline bool requiresFrontendSync() const { return m_requiresFrontendSync; }
    inline void unsetRequiresFrontendSync() { m_requiresFrontendSync = false; }

    // Set by Renderer plugin
    void setLog(const QString &log);
    void setStatus(QShaderProgram::Status status);
    void initializeFromReference(const Shader &other);

    void requestCacheRebuild();

private:
    QVector<QByteArray> m_shaderCode;

    QString m_log;
    bool m_requiresFrontendSync;
    QShaderProgram::Status m_status;
    QShaderProgram::Format m_format;
    bool m_dirty;

};

#ifndef QT_NO_DEBUG_STREAM
inline QDebug operator<<(QDebug dbg, const Shader &shader)
{
    QDebugStateSaver saver(dbg);
    dbg << "QNodeId =" << shader.peerId() << Qt::endl;
    return dbg;
}
#endif

class Q_AUTOTEST_EXPORT ShaderFunctor : public Qt3DCore::QBackendNodeMapper
{
public:
    explicit ShaderFunctor(AbstractRenderer *renderer,
                           ShaderManager *manager);
    Qt3DCore::QBackendNode *create(const Qt3DCore::QNodeCreatedChangeBasePtr &change) const final;
    Qt3DCore::QBackendNode *get(Qt3DCore::QNodeId id) const final;
    void destroy(Qt3DCore::QNodeId id) const final;

private:
    AbstractRenderer *m_renderer;
    ShaderManager *m_shaderManager;
};

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_RENDER_SHADER_H
