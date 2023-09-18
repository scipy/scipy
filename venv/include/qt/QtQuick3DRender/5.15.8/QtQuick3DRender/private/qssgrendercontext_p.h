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

#ifndef QDEMON_RENDER_CONTEXT_H
#define QDEMON_RENDER_CONTEXT_H

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

#include <private/qssgdataref_p.h>

#include <QtQuick3DRender/private/qtquick3drenderglobal_p.h>
#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DRender/private/qssgglimplobjects_p.h>
#include <QtQuick3DRender/private/qssgrenderbackendgles2_p.h>
#include <QtQuick3DRender/private/qssgrenderbackendgl3_p.h>
#include <QtQuick3DRender/private/qssgrenderbackendgl4_p.h>
#include <QtQuick3DRender/private/qssgrenderbackendnull_p.h>
#include <QtQuick3DRender/private/qssgrendervertexbuffer_p.h>
#include <QtQuick3DRender/private/qssgrenderindexbuffer_p.h>
#include <QtQuick3DRender/private/qssgrenderconstantbuffer_p.h>
#include <QtQuick3DRender/private/qssgrenderframebuffer_p.h>
#include <QtQuick3DRender/private/qssgrenderrenderbuffer_p.h>
#include <QtQuick3DRender/private/qssgrenderdepthstencilstate_p.h>
#include <QtQuick3DRender/private/qssgrenderrasterizerstate_p.h>
#include <QtQuick3DRender/private/qssgrenderinputassembler_p.h>
#include <QtQuick3DRender/private/qssgrenderattriblayout_p.h>
#include <QtQuick3DRender/private/qssgrenderimagetexture_p.h>
#include <QtQuick3DRender/private/qssgrendertimerquery_p.h>
#include <QtQuick3DRender/private/qssgrendersync_p.h>
#include <QtQuick3DRender/private/qssgrendertexturecube_p.h>
#include <QtQuick3DRender/private/qssgrenderstoragebuffer_p.h>

#include <QtCore/QString>
#include <QtCore/QVector>
#include <QtCore/QSharedPointer>

#include <QtGui/QSurfaceFormat>

QT_BEGIN_NAMESPACE

enum class QSSGRenderShaderProgramBinaryType
{
    Unknown = 0,
    NVBinary = 1,
};

typedef QHash<QByteArray, QSSGRef<QSSGRenderConstantBuffer>> TContextConstantBufferMap;
typedef QHash<QByteArray, QSSGRef<QSSGRenderStorageBuffer>> TContextStorageBufferMap;
typedef QHash<QSSGRenderBackend::QSSGRenderBackendRasterizerStateObject, QSSGRenderRasterizerState *> TContextRasterizerStateMap;

class QSSGRenderProgramPipeline;

// Now for scoped property access.
template<typename TDataType>
struct QSSGRenderContextScopedProperty : public QSSGRenderGenericScopedProperty<QSSGRenderContext, TDataType>
{
    // the static assert does not compile with MSVC2015
    //Q_STATIC_ASSERT_X(!std::is_reference<TDataType>::value, "Changing the same data!!!");
    typedef typename QSSGRenderGenericScopedProperty<QSSGRenderContext, TDataType>::TGetter TGetter;
    typedef typename QSSGRenderGenericScopedProperty<QSSGRenderContext, TDataType>::TSetter TSetter;
    QSSGRenderContextScopedProperty(QSSGRenderContext &ctx, TGetter getter, TSetter setter)
        : QSSGRenderGenericScopedProperty<QSSGRenderContext, TDataType>(ctx, getter, setter)
    {
    }
    QSSGRenderContextScopedProperty(QSSGRenderContext &ctx, TGetter getter, TSetter setter, const TDataType &inNewValue)
        : QSSGRenderGenericScopedProperty<QSSGRenderContext, TDataType>(ctx, getter, setter, inNewValue)
    {
    }
};

class Q_QUICK3DRENDER_EXPORT QSSGRenderContext
{
    Q_DISABLE_COPY(QSSGRenderContext)
public:
    QAtomicInt ref;
    // these variables represent the current hardware state of the render context.
    QSSGGLHardPropertyContext m_hardwarePropertyContext;

private:
    friend class QSSGRenderContextInterface;
    void releaseResources();
    const QSSGRef<QSSGRenderBackend> m_backend; ///< pointer to our render backend

    QSSGRenderBackend::QSSGRenderBackendRenderTargetObject m_defaultOffscreenRenderTarget; ///< this is a special target set from outside if we
    /// never render to a window directly (GL only)
    qint32 m_depthBits; ///< this is the depth bits count of the default window render target
    qint32 m_stencilBits; ///< this is the stencil bits count of the default window render target
    qint32 m_maxSamples; ///< this is the max samples of the default window render target

protected:
    TContextConstantBufferMap m_constantToImpMap;
    TContextStorageBufferMap m_storageToImpMap;

    qint32 m_maxTextureUnits;
    qint32 m_nextTextureUnit;
    qint32 m_maxConstantBufferUnits;
    qint32 m_nextConstantBufferUnit;

    QVarLengthArray<QSSGGLHardPropertyContext, 4> m_propertyStack;

    bool bindShaderToInputAssembler(const QSSGRef<QSSGRenderInputAssembler> &inputAssembler,
                                    const QSSGRef<QSSGRenderShaderProgram> &shader);
    bool applyPreDrawProperties();
    void onPostDraw();

public:
    QSSGRenderContext(const QSSGRef<QSSGRenderBackend> &inBackend);
    ~QSSGRenderContext();

    const QSSGRef<QSSGRenderBackend> &backend() { return m_backend; }

    void maxTextureSize(qint32 &oWidth, qint32 &oHeight);

    QByteArray shadingLanguageVersion() { return m_backend->getShadingLanguageVersion(); }

    QSSGRenderContextType renderContextType() const { return m_backend->getRenderContextType(); }

    qint32 depthBits() const
    {
        // only query this if a framebuffer is bound
        if (m_hardwarePropertyContext.m_frameBuffer)
            return m_backend->getDepthBits();

        return m_depthBits;
    }

    qint32 stencilBits() const
    {
        // only query this if a framebuffer is bound
        if (m_hardwarePropertyContext.m_frameBuffer)
            return m_backend->getStencilBits();

        return m_stencilBits;
    }

    bool renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps inCap) const
    {
        return m_backend->getRenderBackendCap(inCap);
    }

    bool supportsMultisampleTextures() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::MsTexture);
    }

    qint32 maxSamples() const
    {
        // only query this if a framebuffer is bound
        if (m_hardwarePropertyContext.m_frameBuffer)
            return m_backend->getMaxSamples();

        return m_maxSamples;
    }

    bool supportsConstantBuffer() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::ConstantBuffer);
    }

    bool supportsDXTImages() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::DxtImages);
    }

    bool supportsDepthStencil() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::DepthStencilTexture);
    }

    bool supportsFpRenderTarget() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::FpRenderTarget);
    }

    bool supportsTessellation() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::Tessellation);
    }

    bool supportsGeometryStage() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::Geometry);
    }

    bool supportsCompute() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::Compute);
    }

    bool supportsSampleQuery() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::SampleQuery);
    }

    bool supportsTimerQuery() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::TimerQuery);
    }

    bool supportsCommandSync() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::CommandSync);
    }
    bool supportsTextureArray() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::TextureArray);
    }
    bool supportsStorageBuffer() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::StorageBuffer);
    }
    bool supportsShaderImageLoadStore() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::ShaderImageLoadStore);
    }
    bool supportsProgramPipeline() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::ProgramPipeline);
    }
    // Are blend modes really supported in HW?
    bool supportsAdvancedBlendHW() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::AdvancedBlend);
    }
    bool supportsAdvancedBlendHwKHR() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::AdvancedBlendKHR);
    }
    bool supportsBlendCoherency() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::BlendCoherency);
    }
    bool supportsStandardDerivatives() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::StandardDerivatives);
    }
    bool supportsTextureLod() const
    {
        return renderBackendCap(QSSGRenderBackend::QSSGRenderBackendCaps::TextureLod);
    }

    void setDefaultRenderTarget(quint64 targetID)
    {
        m_defaultOffscreenRenderTarget = reinterpret_cast<QSSGRenderBackend::QSSGRenderBackendRenderTargetObject>(targetID);
    }

    void setDefaultDepthBufferBitCount(qint32 depthBits) { m_depthBits = depthBits; }

    void setDepthStencilState(const QSSGRef<QSSGRenderDepthStencilState> &inDepthStencilState);

    void setRasterizerState(const QSSGRef<QSSGRenderRasterizerState> &inRasterizerState);

    void registerConstantBuffer(QSSGRenderConstantBuffer *buffer);
    QSSGRef<QSSGRenderConstantBuffer> getConstantBuffer(const QByteArray &bufferName) const;
    void bufferDestroyed(QSSGRenderConstantBuffer *buffer);

    qint32 nextConstantBufferUnit();

    void registerStorageBuffer(QSSGRenderStorageBuffer *buffer);
    QSSGRef<QSSGRenderStorageBuffer> getStorageBuffer(const QByteArray &bufferName);
    void bufferDestroyed(QSSGRenderStorageBuffer *buffer);

    void setMemoryBarrier(QSSGRenderBufferBarrierFlags barriers);

    qint32 nextTextureUnit();

    QSSGRef<QSSGRenderAttribLayout> createAttributeLayout(QSSGDataView<QSSGRenderVertexBufferEntry> attribs);
    QSSGRef<QSSGRenderInputAssembler> createInputAssembler(const QSSGRef<QSSGRenderAttribLayout> &attribLayout,
                                                               QSSGDataView<QSSGRef<QSSGRenderVertexBuffer>> buffers,
                                                               const QSSGRef<QSSGRenderIndexBuffer> &indexBuffer,
                                                               QSSGDataView<quint32> strides,
                                                               QSSGDataView<quint32> offsets,
                                                               QSSGRenderDrawMode primType = QSSGRenderDrawMode::Triangles,
                                                               quint32 patchVertexCount = 1);
    void setInputAssembler(const QSSGRef<QSSGRenderInputAssembler> &inputAssembler, bool forceSet = false);

    QSSGRenderVertFragCompilationResult compileSource(
            const char *shaderName,
            QSSGByteView vertShader,
            QSSGByteView fragShader,
            QSSGByteView tessControlShaderSource = QSSGByteView(),
            QSSGByteView tessEvaluationShaderSource = QSSGByteView(),
            QSSGByteView geometryShaderSource = QSSGByteView(),
            bool separateProgram = false,
            QSSGRenderShaderProgramBinaryType type = QSSGRenderShaderProgramBinaryType::Unknown,
            bool binaryProgram = false);

    QSSGRenderVertFragCompilationResult compileBinary(const char *shaderName,
            QSSGRenderShaderProgramBinaryType type,
            QSSGByteView vertShader,
            QSSGByteView fragShader,
            QSSGByteView tessControlShaderSource = QSSGByteView(),
            QSSGByteView tessEvaluationShaderSource = QSSGByteView(),
            QSSGByteView geometryShaderSource = QSSGByteView());
    QSSGRenderVertFragCompilationResult compileBinary(const char *shaderName,
                                                      quint32 format,
                                                      const QByteArray &binary);

    QSSGRenderVertFragCompilationResult compileComputeSource(const QByteArray &shaderName,
                                                                       QSSGByteView computeShaderSource);

    void shaderDestroyed(QSSGRenderShaderProgram *shader);

    QSSGRef<QSSGRenderProgramPipeline> createProgramPipeline();

    void setClearColor(QVector4D inClearColor, bool forceSet = false);
    QVector4D clearColor() const { return m_hardwarePropertyContext.m_clearColor; }

    void setBlendFunction(QSSGRenderBlendFunctionArgument inFunctions, bool forceSet = false);
    QSSGRenderBlendFunctionArgument blendFunction() const
    {
        return m_hardwarePropertyContext.m_blendFunction;
    }

    void setBlendEquation(QSSGRenderBlendEquationArgument inEquations, bool forceSet = false);
    QSSGRenderBlendEquationArgument blendEquation() const
    {
        return m_hardwarePropertyContext.m_blendEquation;
    }
    void resetBlendEquation(bool forceSet = false);

    void setCullingEnabled(bool inEnabled, bool forceSet = false);
    bool isCullingEnabled() const { return m_hardwarePropertyContext.m_cullingEnabled; }

    void setCullFaceMode(QSSGCullFaceMode inCullFaceMode, bool forceSet = false);
    QSSGCullFaceMode cullFaceMode() const { return m_hardwarePropertyContext.m_cullFaceMode; }

    void setDepthFunction(QSSGRenderBoolOp inFunction, bool forceSet = false);
    QSSGRenderBoolOp depthFunction() const { return m_hardwarePropertyContext.m_depthFunction; }

    void setBlendingEnabled(bool inEnabled, bool forceSet = false);
    bool isBlendingEnabled() const { return m_hardwarePropertyContext.m_blendingEnabled; }

    void setDepthWriteEnabled(bool inEnabled, bool forceSet = false);
    bool isDepthWriteEnabled() const { return m_hardwarePropertyContext.m_depthWriteEnabled; }
    void setDepthTestEnabled(bool inEnabled, bool forceSet = false);
    bool isDepthTestEnabled() const { return m_hardwarePropertyContext.m_depthTestEnabled; }

    void setStencilTestEnabled(bool inEnabled, bool forceSet = false);
    bool isStencilTestEnabled() const { return m_hardwarePropertyContext.m_stencilTestEnabled; }

    void setScissorTestEnabled(bool inEnabled, bool forceSet = false);
    bool isScissorTestEnabled() const { return m_hardwarePropertyContext.m_scissorTestEnabled; }
    void setScissorRect(QRect inRect, bool forceSet = false);
    QRect scissorRect() const { return m_hardwarePropertyContext.m_scissorRect; }

    void setViewport(QRect inViewport, bool forceSet = false);
    QRect viewport() const { return m_hardwarePropertyContext.m_viewport; }

    void setColorWritesEnabled(bool inEnabled, bool forceSet = false);
    bool isColorWritesEnabled() const { return m_hardwarePropertyContext.m_colorWritesEnabled; }

    void setMultisampleEnabled(bool inEnabled, bool forceSet = false);
    bool isMultisampleEnabled() const { return m_hardwarePropertyContext.m_multisampleEnabled; }

    void setActiveShader(const QSSGRef<QSSGRenderShaderProgram> &inShader, bool forceSet = false);
    QSSGRef<QSSGRenderShaderProgram> activeShader() const;

    void setActiveProgramPipeline(const QSSGRef<QSSGRenderProgramPipeline> &inProgramPipeline, bool forceSet = false);
    QSSGRef<QSSGRenderProgramPipeline> activeProgramPipeline() const;

    void dispatchCompute(const QSSGRef<QSSGRenderShaderProgram> &inShader, quint32 numGroupsX, quint32 numGroupsY, quint32 numGroupsZ);

    void setDrawBuffers(QSSGDataView<qint32> inDrawBufferSet);
    void setReadBuffer(QSSGReadFace inReadFace);

    void readPixels(QRect inRect, QSSGRenderReadPixelFormat inFormat, QSSGByteRef inWriteBuffer);

    void setRenderTarget(const QSSGRef<QSSGRenderFrameBuffer> &inBuffer, bool forceSet = false);
    void setReadTarget(const QSSGRef<QSSGRenderFrameBuffer> &inBuffer, bool forceSet = false);
    const QSSGRef<QSSGRenderFrameBuffer> &renderTarget() const
    {
        return m_hardwarePropertyContext.m_frameBuffer;
    }

    void solveCullingOptions(const QSSGCullFaceMode);

    void resetBlendState();

    // Push the entire set of properties.
    void pushPropertySet();

    // Pop the entire set of properties, potentially forcing the values
    // to opengl.
    void popPropertySet(bool inForceSetProperties);

    // clear current bound render target
    void clear(QSSGRenderClearFlags flags);
    // clear passed in rendertarget
    void clear(const QSSGRef<QSSGRenderFrameBuffer> &fb, QSSGRenderClearFlags flags);

    // copy framebuffer content between read target and render target
    void blitFramebuffer(qint32 srcX0,
                         qint32 srcY0,
                         qint32 srcX1,
                         qint32 srcY1,
                         qint32 dstX0,
                         qint32 dstY0,
                         qint32 dstX1,
                         qint32 dstY1,
                         QSSGRenderClearFlags flags,
                         QSSGRenderTextureMagnifyingOp filter);

    void copyFramebufferTexture(qint32 srcX0,
                                qint32 srcY0,
                                qint32 width,
                                qint32 height,
                                qint32 dstX0,
                                qint32 dstY0,
                                const QSSGRenderTextureOrRenderBuffer &buffer);

    void draw(QSSGRenderDrawMode drawMode, quint32 count, quint32 offset);

    QSurfaceFormat format() const { return m_backend->format(); }
    void resetStates()
    {
        pushPropertySet();
        popPropertySet(true);
        m_backend->resetStates();
    }

    // Used during layer rendering because we can't set the *actual* viewport to what it should
    // be due to hardware problems.
    // Set during begin render.
    static QMatrix4x4 applyVirtualViewportToProjectionMatrix(const QMatrix4x4 &inProjection,
                                                             const QRectF &inViewport,
                                                             const QRectF &inVirtualViewport);

    static QSSGRef<QSSGRenderContext> createGl(const QSurfaceFormat &format);

    static QSSGRef<QSSGRenderContext> createNull();
};

QT_END_NAMESPACE

#endif
