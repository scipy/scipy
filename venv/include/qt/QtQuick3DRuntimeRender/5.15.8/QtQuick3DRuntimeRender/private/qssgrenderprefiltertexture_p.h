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

#ifndef QSSG_RENDER_PREFILTER_TEXTURE_H
#define QSSG_RENDER_PREFILTER_TEXTURE_H

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

#include <QtQuick3DRender/private/qssgrendertexture2d_p.h>
#include <QtQuick3DRender/private/qssgrendercontext_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderloadedtexture_p.h>

QT_BEGIN_NAMESPACE

class QSSGRenderPrefilterTexture
{
public:
    QAtomicInt ref;
    QSSGRenderPrefilterTexture(const QSSGRef<QSSGRenderContext> &inQSSGRenderContext,
                                 qint32 inWidth,
                                 qint32 inHeight,
                                 const QSSGRef<QSSGRenderTexture2D> &inTexture,
                                 QSSGRenderTextureFormat inDestFormat);
    virtual ~QSSGRenderPrefilterTexture();

    virtual void build(void *inTextureData, qint32 inTextureDataSize, QSSGRenderTextureFormat inFormat) = 0;

    static QSSGRef<QSSGRenderPrefilterTexture> create(const QSSGRef<QSSGRenderContext> &inQSSGRenderContext,
                                                          qint32 inWidth,
                                                          qint32 inHeight,
                                                          const QSSGRef<QSSGRenderTexture2D> &inTexture,
                                                          QSSGRenderTextureFormat inDestFormat);

protected:
    QSSGRef<QSSGRenderTexture2D> m_texture2D;
    QSSGRenderTextureFormat m_destinationFormat;

    qint32 m_width;
    qint32 m_height;
    qint32 m_maxMipMapLevel;
    qint32 m_sizeOfFormat;
    qint32 m_sizeOfInternalFormat;
    qint32 m_internalNoOfComponent;
    qint32 m_noOfComponent;
    QSSGRef<QSSGRenderContext> m_renderContext;
};

class QSSGRenderPrefilterTextureCPU : public QSSGRenderPrefilterTexture
{
public:
    QSSGRenderPrefilterTextureCPU(const QSSGRef<QSSGRenderContext> &inQSSGRenderContext,
                                    qint32 inWidth,
                                    qint32 inHeight,
                                    const QSSGRef<QSSGRenderTexture2D> &inTexture,
                                    QSSGRenderTextureFormat inDestFormat);

    void build(void *inTextureData, qint32 inTextureDataSize, QSSGRenderTextureFormat inFormat) override;

    QSSGTextureData createBsdfMipLevel(QSSGTextureData &inCurMipLevel, QSSGTextureData &inPrevMipLevel, qint32 width, qint32 height);

    int wrapMod(int a, int base);
    void getWrappedCoords(int &sX, int &sY, int width, int height);
};

class QSSGRenderPrefilterTextureCompute : public QSSGRenderPrefilterTexture
{
public:
    QSSGRenderPrefilterTextureCompute(const QSSGRef<QSSGRenderContext> &inQSSGRenderContext,
                                        qint32 inWidth,
                                        qint32 inHeight,
                                        const QSSGRef<QSSGRenderTexture2D> &inTexture,
                                        QSSGRenderTextureFormat inDestFormat);
    ~QSSGRenderPrefilterTextureCompute() override;

    void build(void *inTextureData, qint32 inTextureDataSize, QSSGRenderTextureFormat inFormat) override;

private:
    void createLevel0Tex(void *inTextureData, qint32 inTextureDataSize, QSSGRenderTextureFormat inFormat);

    QSSGRef<QSSGRenderShaderProgram> m_bsdfProgram;
    QSSGRef<QSSGRenderShaderProgram> m_bsdfRGBEProgram;
    QSSGRef<QSSGRenderShaderProgram> m_uploadProgram_RGBA8;
    QSSGRef<QSSGRenderShaderProgram> m_uploadProgram_RGB8;
    QSSGRef<QSSGRenderTexture2D> m_level0Tex;
    bool m_textureCreated = false;

    QSSGRenderShaderProgram *createComputeProgram(const QSSGRef<QSSGRenderContext> &context,
                                                  QSSGRenderTextureFormat inFormat);
    QSSGRef<QSSGRenderShaderProgram> getOrCreateUploadComputeProgram(const QSSGRef<QSSGRenderContext> &context,
                                                                         QSSGRenderTextureFormat inFormat);
};
QT_END_NAMESPACE

#endif
