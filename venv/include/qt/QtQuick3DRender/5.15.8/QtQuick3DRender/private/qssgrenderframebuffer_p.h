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

#ifndef QSSG_RENDER_RENDER_FRAME_BUFFER_H
#define QSSG_RENDER_RENDER_FRAME_BUFFER_H

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

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DRender/private/qssgrenderbackend_p.h>

QT_BEGIN_NAMESPACE

class QSSGRenderContext;
class QSSGRenderTexture2D;
class QSSGRenderRenderBuffer;
class QSSGRenderTextureCube;

class Q_QUICK3DRENDER_EXPORT QSSGRenderTextureOrRenderBuffer
{
    // ### this could be a union
    QSSGRef<QSSGRenderTexture2D> m_texture2D;
    QSSGRef<QSSGRenderTextureCube> m_textureCube;
    QSSGRef<QSSGRenderRenderBuffer> m_renderBuffer;

public:
    QSSGRenderTextureOrRenderBuffer(const QSSGRef<QSSGRenderTexture2D> &texture);
    QSSGRenderTextureOrRenderBuffer(const QSSGRef<QSSGRenderRenderBuffer> &render);
    QSSGRenderTextureOrRenderBuffer(const QSSGRef<QSSGRenderTextureCube> &textureCube);
    QSSGRenderTextureOrRenderBuffer();
    QSSGRenderTextureOrRenderBuffer(const QSSGRenderTextureOrRenderBuffer &other);
    ~QSSGRenderTextureOrRenderBuffer();

    QSSGRenderTextureOrRenderBuffer &operator=(const QSSGRenderTextureOrRenderBuffer &other);

    bool hasTexture2D() const { return m_texture2D != nullptr; }
    bool hasTextureCube() const { return m_textureCube != nullptr; }
    bool hasRenderBuffer() const { return m_renderBuffer != nullptr; }

    QSSGRef<QSSGRenderTexture2D> texture2D() const;
    QSSGRef<QSSGRenderTextureCube> textureCube() const;
    QSSGRef<QSSGRenderRenderBuffer> renderBuffer() const;
};

class Q_QUICK3DRENDER_EXPORT QSSGRenderFrameBuffer
{
    Q_DISABLE_COPY(QSSGRenderFrameBuffer)
public:
    QAtomicInt ref;

private:
    QSSGRef<QSSGRenderContext> m_context; ///< pointer to context
    QSSGRef<QSSGRenderBackend> m_backend; ///< pointer to backend

    QSSGRenderTextureOrRenderBuffer m_attachments[static_cast<int>(QSSGRenderFrameBufferAttachment::LastAttachment)]; ///< attachments array
    QSSGRenderBackend::QSSGRenderBackendRenderTargetObject m_bufferHandle; ///< opaque backend handle

public:
    /**
     * @brief constructor
     *
     * @param[in] context		Pointer to context
     * @param[in] fnd			Pointer to foundation
     *
     * @return No return.
     */
    QSSGRenderFrameBuffer(const QSSGRef<QSSGRenderContext> &context);

    /// destructor
    ~QSSGRenderFrameBuffer();

    /**
     * @brief query attachment
     *
     *
     * @return buffer format
     */
    QSSGRenderTextureOrRenderBuffer attachment(QSSGRenderFrameBufferAttachment attachment);

    /**
     * @brief Attach a render or texture buffer to a render target
     *		  For texture attachments we use always level 0
     *
     * @param[in] attachment		Attachment point (e.g. COLOR0, DEPTH...)
     * @param[in] buffer			Contains a pointer to the attachment
     * @param[in] target			Attachment texture target
     *
     * @return no return
     */
    void attach(QSSGRenderFrameBufferAttachment attachment,
                        const QSSGRenderTextureOrRenderBuffer &buffer,
                        QSSGRenderTextureTargetType target = QSSGRenderTextureTargetType::Texture2D);

    /**
     * @brief Attach a particular face of the texture cubemap to a render target
     *
     * @param[in] attachment		Attachment point (e.g. COLOR0, DEPTH...)
     * @param[in] buffer			Pointer to the Texture Array which contains the
     * layers
     * @param[in] face				The face of the cubemap that will be attached to the
     * target
     * @param[in] level				Mip level of the texture that will be attached
     * (default 0)
     *
     * @return no return
     */
    // ### currently unused
    void attachFace(QSSGRenderFrameBufferAttachment attachment,
                            const QSSGRenderTextureOrRenderBuffer &buffer,
                            QSSGRenderTextureCubeFace face);

    /**
     * @brief Check that this framebuffer is complete and can be rendered to.
     *
     *
     * @return true if complete
     */
    bool isComplete();

    /**
     * @brief query if framebuffer has any attachment
     *
     * @return true if any attachment
     */
    bool hasAnyAttachment() { return (m_attachmentBits != 0); }

    /**
     * @brief get the backend object handle
     *
     * @return the backend object handle.
     */
    QSSGRenderBackend::QSSGRenderBackendRenderTargetObject handle()
    {
        return m_bufferHandle;
    }

private:
    /**
     * @brief releaes an attached object
     *
     * @return which target we released
     */
    QSSGRenderTextureTargetType releaseAttachment(QSSGRenderFrameBufferAttachment idx);

    quint32 m_attachmentBits; ///< holds flags for current attached buffers
};

QT_END_NAMESPACE

#endif
