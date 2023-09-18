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

#ifndef QSSG_RENDER_RENDER_TEXTURE_2D_H
#define QSSG_RENDER_RENDER_TEXTURE_2D_H

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
#include <QtQuick3DRender/private/qssgrendertexturebase_p.h>

QT_BEGIN_NAMESPACE

class QSSGRenderContext;
class QSGTexture;

class Q_QUICK3DRENDER_EXPORT QSSGRenderTexture2D : public QSSGRenderTextureBase
{
    Q_DISABLE_COPY(QSSGRenderTexture2D)
private:
    qint32 m_width; ///< texture width
    qint32 m_height; ///< texture height

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
    QSSGRenderTexture2D(const QSSGRef<QSSGRenderContext> &context);

    QSSGRenderTexture2D(const QSSGRef<QSSGRenderContext> &context, QSGTexture *qsgTexture);

    virtual ~QSSGRenderTexture2D() override;

    // Get the texture details for mipmap level 0 if it was set.
    QSSGTextureDetails textureDetails() const override;

    /**
     * @brief Create GL texture object and upload data
     *
     * @param[in] newBuffer			Texture data for level 0
     * @param[in] inMipLevel		Texture level count
     * @param[in] width				Texture width
     * @param[in] height			Texture height
     * @param[in] format			Texture data format
     * @param[in] formaInternal		Texture internal format
     *
     * @return No return.
     */
    void setTextureData(QSSGByteView newBuffer,
                                quint8 inMipLevel,
                                qint32 width,
                                qint32 height,
                                QSSGRenderTextureFormat format,
                                QSSGRenderTextureFormat formaInternal = QSSGRenderTextureFormat::Unknown);

    /**
     * @brief Create memory storage for a texture object
     *		  This create a texture storage which is immutable in size and format
     *		  Use this for textures used within compute shaders
     *
     * @param[in] inLevels			Texture level count
     * @param[in] width				Texture width
     * @param[in] height			Texture height
     * @param[in] formaInternal		Texture internal format
     * @param[in] format			Texture data format of dataBuffer
     * @param[in] dataBuffer		Texture data for level 0
     *
     * @return No return.
     */
    void setTextureStorage(qint32 inLevels,
                                   qint32 width,
                                   qint32 height,
                                   QSSGRenderTextureFormat formaInternal,
                                   QSSGRenderTextureFormat format = QSSGRenderTextureFormat::Unknown,
                                   QSSGByteView dataBuffer = QSSGByteView());

    void setTextureDataMultisample(qint32 sampleCount, qint32 width, qint32 height, QSSGRenderTextureFormat format);

    // Update a sub-rect of the image.  newBuffer is expected to be a continguous subrect of the
    // image.
    void setTextureSubData(QSSGByteView newBuffer,
                                   quint8 inMipLevel,
                                   qint32 inXOffset,
                                   qint32 inYOffset,
                                   qint32 inSubImageWidth,
                                   qint32 inSubImageHeight,
                                   QSSGRenderTextureFormat format);
    // Generate a set of mipmaps from mipLevel( 0 ).  Uses the graphis layer to do this if
    // possible
    // glGenerateMipmap
    void generateMipmaps(QSSGRenderHint genType = QSSGRenderHint::Nicest);

    /**
     * @brief Bind a texture for shader access
     *
     *
     * @return No return.
     */
    void bind() override;

    /**
     * @brief Query if texture needs coordinate swizzle
     *
     * @return texture swizzle mode
     */
    QSSGRenderTextureSwizzleMode textureSwizzleMode() override
    {
        // if our backend supports hardware texture swizzle then there is no need for a shader
        // swizzle
        return (m_backend->getRenderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::TexSwizzle))
                ? QSSGRenderTextureSwizzleMode::NoSwizzle
                : m_backend->getTextureSwizzleMode(m_format);
    }

};

QT_END_NAMESPACE

#endif
