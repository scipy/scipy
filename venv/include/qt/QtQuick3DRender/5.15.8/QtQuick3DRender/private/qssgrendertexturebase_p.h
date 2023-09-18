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

#ifndef QSSG_RENDER_TEXTURE_BUFFER_H
#define QSSG_RENDER_TEXTURE_BUFFER_H

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
class QSSGRenderTextureSampler;

struct QSSGTextureDetails
{
    qint32 width = 0;
    qint32 height = 0;
    qint32 depth = 0;
    qint32 sampleCount = 1;
    QSSGRenderTextureFormat format = QSSGRenderTextureFormat::Unknown;

    QSSGTextureDetails(qint32 w, qint32 h, qint32 d, qint32 samples, QSSGRenderTextureFormat f)
        : width(w), height(h), depth(d), sampleCount(samples), format(f)
    {
    }
    QSSGTextureDetails() = default;
};

class Q_QUICK3DRENDER_EXPORT QSSGRenderTextureBase
{
    Q_DISABLE_COPY(QSSGRenderTextureBase)
public:
    QAtomicInt ref;

protected:
    QSSGRef<QSSGRenderContext> m_context; ///< pointer to context
    QSSGRef<QSSGRenderBackend> m_backend; ///< pointer to backend
    QSSGRenderBackend::QSSGRenderBackendTextureObject m_handle; ///< opaque backend handle
    qint32 m_textureUnit; ///< texture unit this texture should use
    bool m_samplerParamsDirty; ///< true if sampler state is dirty
    bool m_texStateDirty; ///< true if texture object state is dirty
    qint32 m_sampleCount; ///< texture height
    QSSGRenderTextureFormat m_format; ///< texture format
    QSSGRenderTextureTargetType m_texTarget; ///< texture target
    QSSGRenderTextureSampler *m_sampler; ///< current texture sampler state
    qint32 m_baseLevel; ///< minimum lod specified
    qint32 m_maxLevel; ///< maximum lod specified
    qint32 m_maxMipLevel; ///< highest mip level
    bool m_immutable; ///< true if this is a immutable texture ( size and format )
    bool m_ownsTexture = true;

public:
    /**
     * @brief constructor
     *
     * @param[in] context		Pointer to context
     * @param[in] fnd			Pointer to foundation
     * @param[in] texTarget		Texture target
     *
     * @return No return.
     */
    QSSGRenderTextureBase(const QSSGRef<QSSGRenderContext> &context, QSSGRenderTextureTargetType texTarget, bool ownsTexture = true);

    virtual ~QSSGRenderTextureBase();

    virtual void setMinFilter(QSSGRenderTextureMinifyingOp value);
    virtual void setMagFilter(QSSGRenderTextureMagnifyingOp value);

    virtual void setBaseLevel(qint32 value);
    virtual void setMaxLevel(qint32 value);

    virtual void setTextureWrapS(QSSGRenderTextureCoordOp value);
    virtual void setTextureWrapT(QSSGRenderTextureCoordOp value);

    virtual void setTextureCompareMode(QSSGRenderTextureCompareMode value);
    virtual void setTextureCompareFunc(QSSGRenderTextureCompareOp value);

    virtual void setTextureUnit(qint32 unit) { m_textureUnit = unit; }
    virtual qint32 textureUnit() const { return m_textureUnit; }

    // Get the texture details for mipmap level 0 if it was set.
    virtual QSSGTextureDetails textureDetails() const = 0;

    bool isMultisampleTexture() const { return (m_texTarget == QSSGRenderTextureTargetType::Texture2D_MS); }
    qint32 sampleCount() const { return m_sampleCount; }
    bool isImmutable() const { return m_immutable; }
    QSSGRenderTextureTargetType target() const { return m_texTarget; }
    void setsamplerParamsDirty() { m_samplerParamsDirty = true; }

    /**
     * @brief Bind a texture for shader access
     *
     *
     * @return No return.
     */
    virtual void bind() = 0;

    virtual quint32 numMipmaps() const { return m_maxMipLevel; }

    /**
     * @brief Query if texture needs coordinate swizzle
     *
     * @return texture swizzle mode
     */
    virtual QSSGRenderTextureSwizzleMode textureSwizzleMode()
    {
        // if our backend supports hardware texture swizzle then there is no need for a shader
        // swizzle
        return (m_backend->getRenderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::TexSwizzle))
                ? QSSGRenderTextureSwizzleMode::NoSwizzle
                : m_backend->getTextureSwizzleMode(m_format);
    }

    /**
     * @brief get the backend object handle
     *
     * @return the backend object handle.
     */
    virtual QSSGRenderBackend::QSSGRenderBackendTextureObject handle() { return m_handle; }

protected:
    void applyTexParams();
    void applyTexSwizzle();
};

QT_END_NAMESPACE

#endif
