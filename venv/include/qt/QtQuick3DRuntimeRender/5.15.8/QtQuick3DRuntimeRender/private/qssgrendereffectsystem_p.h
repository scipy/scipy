/****************************************************************************
**
** Copyright (C) 2020 The Qt Company Ltd.
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

#ifndef QSSG_RENDER_EFFECT_SYSTEM_H
#define QSSG_RENDER_EFFECT_SYSTEM_H

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
#include <QtQuick3DUtils/private/qssgoption_p.h>
#include <QtQuick3DRuntimeRender/private/qssgrenderdynamicobjectsystem_p.h>
#include <QtGui/QVector2D>

QT_BEGIN_NAMESPACE
struct QSSGRenderEffect;
struct QSSGEffectContext;
struct QSSGEffectShader;
class QSSGEffectSystemInterface;
class QSSGResourceManager;
class QSSGRenderDepthStencilState;

namespace dynamic {
struct QSSGCommand;
struct QSSGApplyInstanceValue;
struct QSSGBindBuffer;
struct QSSGApplyValue;
struct QSSGApplyBlending;
struct QSSGAllocateBuffer;
struct QSSGBindShader;
struct QSSGAllocateImage;
struct QSSGAllocateDataBuffer;
struct QSSGApplyBufferValue;
struct QSSGApplyDepthValue;
struct QSSGApplyDataBufferValue;
struct QSSGApplyRenderState;
struct QSSGApplyImageValue;
struct QSSGDepthStencil;
}

struct QSSGEffectRenderArgument
{
    QSSGRenderEffect *m_effect;
    QSSGRef<QSSGRenderTexture2D> m_colorBuffer;
    // Some effects need the camera near and far ranges.
    QVector2D m_cameraClipRange;
    // Some effects require the depth buffer from the rendering of the layer
    // most do not.
    QSSGRef<QSSGRenderTexture2D> m_depthTexture;
    // this is a depth pre-pass texture we need for some effects like bloom
    // actually we need the stencil values
    QSSGRef<QSSGRenderTexture2D> m_depthStencilBuffer;

    QSSGEffectRenderArgument(QSSGRenderEffect *inEffect,
                             const QSSGRef<QSSGRenderTexture2D> &inColorBuffer,
                             const QVector2D &inCameraClipRange,
                             const QSSGRef<QSSGRenderTexture2D> &inDepthTexture = nullptr,
                             const QSSGRef<QSSGRenderTexture2D> &inDepthBuffer = nullptr);
};

struct QSSGEffectTextureData
{
    QSSGRef<QSSGRenderTexture2D> texture;
    bool needsAlphaMultiply = false;
    QSSGEffectTextureData(const QSSGRef<QSSGRenderTexture2D> &inTexture, bool inNeedsMultiply);
    QSSGEffectTextureData() = default;
};

/**
 * An effect is essentially a function that takes a image and produces a new image.
 * The source and dest images aren't guaranteed to be the same size, the effect may enlarge or
 * shrink the result.
 * A specialization is when you want the effect to render to the final render target instead of
 * to a separate image. In this case the effect cannot enlarge or shrink the final target and it
 * will render to the destination buffer using the given MVP.
 */
class Q_QUICK3DRUNTIMERENDER_EXPORT QSSGEffectSystem
{
    Q_DISABLE_COPY(QSSGEffectSystem)

    typedef QHash<QString, char *> TPathDataMap;
    typedef QSet<QString> TPathSet;
    typedef QHash<TStrStrPair, QSSGRef<QSSGEffectShader>> TShaderMap;
    typedef QVector<QSSGRef<QSSGEffectContext>> TContextList;

    QSSGRenderContextInterface *m_context;
    QSSGRef<QSSGResourceManager> m_resourceManager;
    QVector<QByteArray> m_effectList;
    TContextList m_contexts;
    QString m_textureStringBuilder;
    QString m_textureStringBuilder2;
    TShaderMap m_shaderMap;
    QSSGRef<QSSGRenderDepthStencilState> m_defaultStencilState;
    QVector<QSSGRef<QSSGRenderDepthStencilState>> m_depthStencilStates;

public:
    QAtomicInt ref;
    QSSGEffectSystem(QSSGRenderContextInterface *inContext);

    ~QSSGEffectSystem();

    QSSGEffectContext &getEffectContext(QSSGRenderEffect &inEffect);

    void allocateBuffer(QSSGRenderEffect &inEffect,
                        const dynamic::QSSGAllocateBuffer &inCommand,
                        qint32 inFinalWidth,
                        qint32 inFinalHeight,
                        QSSGRenderTextureFormat inSourceTextureFormat);

    void allocateImage(QSSGRenderEffect &inEffect, const dynamic::QSSGAllocateImage &inCommand, qint32 inFinalWidth, qint32 inFinalHeight);

    void allocateDataBuffer(QSSGRenderEffect &inEffect, const dynamic::QSSGAllocateDataBuffer &inCommand);

    QSSGRef<QSSGRenderTexture2D> findTexture(QSSGRenderEffect *inEffect, const QByteArray &inName);

    QSSGRef<QSSGRenderFrameBuffer> bindBuffer(QSSGRenderEffect &inEffect,
                                              const dynamic::QSSGBindBuffer &inCommand,
                                              QMatrix4x4 &outMVP,
                                              QVector2D &outDestSize);

    QSSGRef<QSSGEffectShader> bindShader(const QSSGRenderEffect &effect, const dynamic::QSSGBindShader &inCommand);

    void doApplyInstanceValue(QSSGRenderEffect *inEffect,
                              const QByteArray &inPropertyName,
                              const QVariant &propertyValue,
                              QSSGRenderShaderDataType inPropertyType,
                              const QSSGRef<QSSGRenderShaderProgram> &inShader);

    void applyInstanceValue(QSSGRenderEffect &inEffect,
                            const QSSGRef<QSSGRenderShaderProgram> &inShader,
                            const dynamic::QSSGApplyInstanceValue &inCommand);

    void applyValue(QSSGRenderEffect &inEffect, const QSSGRef<QSSGRenderShaderProgram> &inShader, const dynamic::QSSGApplyValue &inCommand);

    bool applyBlending(const dynamic::QSSGApplyBlending &inCommand);

    // This has the potential to change the source texture for the current render pass
    QSSGEffectTextureData applyBufferValue(QSSGRenderEffect *inEffect,
                                           const QSSGRef<QSSGRenderShaderProgram> &inShader,
                                           const dynamic::QSSGApplyBufferValue &inCommand,
                                           const QSSGRef<QSSGRenderTexture2D> &inSourceTexture,
                                           const QSSGEffectTextureData &inCurrentSourceTexture);

    void applyDepthValue(QSSGRenderEffect *inEffect,
                         const QSSGRef<QSSGRenderShaderProgram> &inShader,
                         const dynamic::QSSGApplyDepthValue &inCommand,
                         const QSSGRef<QSSGRenderTexture2D> &inTexture);

    void applyImageValue(QSSGRenderEffect *inEffect,
                         const QSSGRef<QSSGRenderShaderProgram> &inShader,
                         const dynamic::QSSGApplyImageValue &inCommand);

    void applyDataBufferValue(QSSGRenderEffect *inEffect,
                              const QSSGRef<QSSGRenderShaderProgram> &inShader,
                              const dynamic::QSSGApplyDataBufferValue &inCommand);

    void applyRenderStateValue(const QSSGRef<QSSGRenderFrameBuffer> &inTarget,
                               const QSSGRef<QSSGRenderTexture2D> &inDepthStencilTexture,
                               const dynamic::QSSGApplyRenderState &theCommand);

    static bool compareDepthStencilState(QSSGRenderDepthStencilState &inState, dynamic::QSSGDepthStencil &inStencil);

    void renderPass(QSSGEffectShader &inShader,
                    const QMatrix4x4 &inMVP,
                    const QSSGEffectTextureData &inSourceTexture,
                    const QSSGRef<QSSGRenderFrameBuffer> &inFrameBuffer,
                    QVector2D &inDestSize,
                    const QVector2D &inCameraClipRange,
                    const QSSGRef<QSSGRenderTexture2D> &inDepthStencil,
                    QSSGOption<dynamic::QSSGDepthStencil> inDepthStencilCommand);

    void doRenderEffect(QSSGRenderEffect *inEffect,
                        const QSSGRef<QSSGRenderTexture2D> &inSourceTexture,
                        QMatrix4x4 &inMVP,
                        const QSSGRef<QSSGRenderFrameBuffer> inTarget,
                        bool inEnableBlendWhenRenderToTarget,
                        const QSSGRef<QSSGRenderTexture2D> &inDepthTexture,
                        const QSSGRef<QSSGRenderTexture2D> &inDepthStencilTexture,
                        const QVector2D &inCameraClipRange);

    // Render this effect.  Returns false in the case the effect wasn't rendered and
    // the render state is guaranteed to be the same as before.
    // The texture returned is allocated using the resource manager, and it is up to the caller
    // to deallocate it or return it to the temporary pool if items when necessary.
    // Pass in true if you want the result image premultiplied.  Most of the functions in the
    // system assume non-premultiplied color for images so probably this is false.
    QSSGRef<QSSGRenderTexture2D> renderEffect(const QSSGEffectRenderArgument &inRenderArgument);

    // Render the effect to the currently bound render target using this MVP and optionally
    // enabling blending when rendering to the target
    bool renderEffect(const QSSGEffectRenderArgument &inRenderArgument, QMatrix4x4 &inMVP, bool inEnableBlendWhenRenderToTarget);

    // Calling release effect context with no context results in no problems.
    void releaseEffectContext(QSSGEffectContext *inContext);

    // If the effect has a context you can call this to clear persistent buffers back to their
    // original value.
    void resetEffectFrameData(QSSGEffectContext &inContext);

    // Set the shader data for a given path.  Used when a path doesn't correspond to a file but
    // the data has been auto-generated. The system will look for data under this path key
    // during the BindShader effect command.
    void setShaderData(const QByteArray &path, const char *data, const char *inShaderType, const char *inShaderVersion, bool inHasGeomShader, bool inIsComputeShader);

    void init();

    QSSGRef<QSSGResourceManager> getResourceManager();
};

QT_END_NAMESPACE
#endif
