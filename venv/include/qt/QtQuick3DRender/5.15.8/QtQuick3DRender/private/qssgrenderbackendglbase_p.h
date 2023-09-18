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

#ifndef QSSG_RENDER_BACKEND_GL_BASE_H
#define QSSG_RENDER_BACKEND_GL_BASE_H

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

/// @file QSSGrenderbackendglbase.h
///       NVRender OpenGL Core backend definition.

#include <QtQuick3DRender/private/qssgrenderbackend_p.h>
#include <QtQuick3DRender/private/qssgopenglutil_p.h>
#include <QtQuick3DRender/private/qssgrenderbackendshaderprogramgl_p.h>

#include <QtCore/QVector>

#include <QtGui/QSurfaceFormat>
#include <QtGui/QOpenGLFunctions>

QT_BEGIN_NAMESPACE

// Enable this to log opengl errors instead of an assert
//#define RENDER_BACKEND_LOG_GL_ERRORS

///< forward declaration
class QSSGRenderBackendRasterizerStateGL;
class QSSGRenderBackendDepthStencilStateGL;

namespace QSSGGlExtStrings {
QByteArray exts3tc();
QByteArray extsdxt();
QByteArray extsAniso();
QByteArray extsTexSwizzle();
QByteArray extsFPRenderTarget();
QByteArray extsTimerQuery();
QByteArray extsGpuShader5();
}

class QSSGRenderBackendGLBase : public QSSGRenderBackend
{
    QSSGRenderBackendGLBase() = default;

public:
    /// constructor
    QSSGRenderBackendGLBase(const QSurfaceFormat &format);
    /// destructor
    ~QSSGRenderBackendGLBase() override;

public:
    /// API Interface
    QSSGRenderContextType getRenderContextType() const override;
    bool isESCompatible() const;

    QByteArray getShadingLanguageVersion() override;
    /// get implementation depended values
    qint32 getMaxCombinedTextureUnits() override;
    bool getRenderBackendCap(QSSGRenderBackendCaps inCap) const override;
    qint32 getDepthBits() const override;
    qint32 getStencilBits() const override;
    qint32 getMaxSamples() const override;
    void getRenderBackendValue(QSSGRenderBackendQuery inQuery, qint32 *params) const override;

    /// state get/set functions
    void setRenderState(bool bEnable, const QSSGRenderState value) override;
    bool getRenderState(const QSSGRenderState value) override;
    virtual QSSGRenderBackendDepthStencilStateObject createDepthStencilState(
            bool enableDepth,
            bool depthMask,
            QSSGRenderBoolOp depthFunc,
            bool enableStencil,
            QSSGRenderStencilFunction &stencilFuncFront,
            QSSGRenderStencilFunction &stencilFuncBack,
            QSSGRenderStencilOperation &depthStencilOpFront,
            QSSGRenderStencilOperation &depthStencilOpBack) override;
    virtual void releaseDepthStencilState(QSSGRenderBackendDepthStencilStateObject inDepthStencilState) override;
    virtual QSSGRenderBackendRasterizerStateObject createRasterizerState(float depthBias, float depthScale) override;
    void releaseRasterizerState(QSSGRenderBackendRasterizerStateObject rasterizerState) override;
    virtual void setDepthStencilState(QSSGRenderBackendDepthStencilStateObject inDepthStencilState) override;
    void setRasterizerState(QSSGRenderBackendRasterizerStateObject rasterizerState) override;
    QSSGRenderBoolOp getDepthFunc() override;
    void setDepthFunc(const QSSGRenderBoolOp func) override;
    bool getDepthWrite() override;
    void setDepthWrite(bool bEnable) override;
    void setColorWrites(bool bRed, bool bGreen, bool bBlue, bool bAlpha) override;
    void setMultisample(bool bEnable) override;
    void getBlendFunc(QSSGRenderBlendFunctionArgument *pBlendFuncArg) override;
    void setBlendFunc(const QSSGRenderBlendFunctionArgument &blendFuncArg) override;
    void setBlendEquation(const QSSGRenderBlendEquationArgument &pBlendEquArg) override;
    void setBlendBarrier(void) override;
    QSSGCullFaceMode getCullFaceMode() override;
    void setCullFaceMode(const QSSGCullFaceMode cullFaceMode) override;
    void getScissorRect(QRect *pRect) override;
    void setScissorRect(const QRect &rect) override;
    void getViewportRect(QRect *pRect) override;
    void setViewportRect(const QRect &rect) override;

    void setClearColor(const QVector4D *pClearColor) override;
    void clear(QSSGRenderClearFlags flags) override;

    /// resource handling
    QSSGRenderBackendBufferObject createBuffer(QSSGRenderBufferType bindFlags,
                                                 QSSGRenderBufferUsageType usage,
                                                 QSSGByteView hostData) override;
    void bindBuffer(QSSGRenderBackendBufferObject bo, QSSGRenderBufferType bindFlags) override;
    void releaseBuffer(QSSGRenderBackendBufferObject bo) override;
    void updateBuffer(QSSGRenderBackendBufferObject bo,
                      QSSGRenderBufferType bindFlags,
                      QSSGRenderBufferUsageType usage,
                      QSSGByteView data) override;
    void updateBufferRange(QSSGRenderBackendBufferObject bo,
                           QSSGRenderBufferType bindFlags,
                           size_t offset,
                           QSSGByteView data) override;
    void *mapBuffer(QSSGRenderBackendBufferObject bo,
                    QSSGRenderBufferType bindFlags,
                    size_t offset,
                    size_t length,
                    QSSGRenderBufferAccessFlags accessFlags) override;
    bool unmapBuffer(QSSGRenderBackendBufferObject bo, QSSGRenderBufferType bindFlags) override;
    void setMemoryBarrier(QSSGRenderBufferBarrierFlags barriers) override;

    QSSGRenderBackendQueryObject createQuery() override;
    void releaseQuery(QSSGRenderBackendQueryObject qo) override;
    void beginQuery(QSSGRenderBackendQueryObject qo, QSSGRenderQueryType type) override;
    void endQuery(QSSGRenderBackendQueryObject qo, QSSGRenderQueryType type) override;
    void getQueryResult(QSSGRenderBackendQueryObject qo, QSSGRenderQueryResultType resultType, quint32 *params) override;
    void getQueryResult(QSSGRenderBackendQueryObject qo, QSSGRenderQueryResultType resultType, quint64 *params) override;
    void setQueryTimer(QSSGRenderBackendQueryObject qo) override;

    QSSGRenderBackendSyncObject createSync(QSSGRenderSyncType tpye, QSSGRenderSyncFlags syncFlags) override;
    void releaseSync(QSSGRenderBackendSyncObject so) override;
    void waitSync(QSSGRenderBackendSyncObject so, QSSGRenderCommandFlushFlags syncFlags, quint64 timeout) override;

    QSSGRenderBackendRenderTargetObject createRenderTarget() override;
    void releaseRenderTarget(QSSGRenderBackendRenderTargetObject rto) override;

    void renderTargetAttach(QSSGRenderBackendRenderTargetObject rto,
                            QSSGRenderFrameBufferAttachment attachment,
                            QSSGRenderBackendRenderbufferObject rbo) override;

    void renderTargetAttach(QSSGRenderBackendRenderTargetObject rto,
                            QSSGRenderFrameBufferAttachment attachment,
                            QSSGRenderBackendTextureObject to,
                            QSSGRenderTextureTargetType target = QSSGRenderTextureTargetType::Texture2D) override;

    void renderTargetAttach(QSSGRenderBackendRenderTargetObject rto,
                            QSSGRenderFrameBufferAttachment attachment,
                            QSSGRenderBackendTextureObject to,
                            qint32 level,
                            qint32 layer) override;

    void setRenderTarget(QSSGRenderBackendRenderTargetObject rto) override;
    bool renderTargetIsValid(QSSGRenderBackendRenderTargetObject rto) override;

    QSSGRenderBackendRenderbufferObject createRenderbuffer(QSSGRenderRenderBufferFormat storageFormat,
                                                             qint32 width,
                                                             qint32 height) override;
    void releaseRenderbuffer(QSSGRenderBackendRenderbufferObject rbo) override;
    bool resizeRenderbuffer(QSSGRenderBackendRenderbufferObject rbo,
                            QSSGRenderRenderBufferFormat storageFormat,
                            qint32 width,
                            qint32 height) override;

    QSSGRenderBackendTextureObject createTexture() override;
    void bindTexture(QSSGRenderBackendTextureObject to, QSSGRenderTextureTargetType target, qint32 unit) override;
    void setActiveTexture(qint32 unit) override;
    void bindImageTexture(QSSGRenderBackendTextureObject to,
                          quint32 unit,
                          qint32 level,
                          bool layered,
                          qint32 layer,
                          QSSGRenderImageAccessType access,
                          QSSGRenderTextureFormat format) override;
    void releaseTexture(QSSGRenderBackendTextureObject to) override;
    void setTextureData2D(QSSGRenderBackendTextureObject to,
                          QSSGRenderTextureTargetType target,
                          qint32 level,
                          QSSGRenderTextureFormat internalFormat,
                          qint32 width,
                          qint32 height,
                          qint32 border,
                          QSSGRenderTextureFormat format,
                          QSSGByteView hostData) override;
    void setTextureDataCubeFace(QSSGRenderBackendTextureObject to,
                                QSSGRenderTextureTargetType target,
                                qint32 level,
                                QSSGRenderTextureFormat internalFormat,
                                qint32 width,
                                qint32 height,
                                qint32 border,
                                QSSGRenderTextureFormat format,
                                QSSGByteView hostData) override;
    void createTextureStorage2D(QSSGRenderBackendTextureObject to,
                                QSSGRenderTextureTargetType target,
                                qint32 levels,
                                QSSGRenderTextureFormat internalFormat,
                                qint32 width,
                                qint32 height) override;
    void setTextureSubData2D(QSSGRenderBackendTextureObject to,
                             QSSGRenderTextureTargetType target,
                             qint32 level,
                             qint32 xOffset,
                             qint32 yOffset,
                             qint32 width,
                             qint32 height,
                             QSSGRenderTextureFormat format,
                             QSSGByteView hostData) override;
    void setCompressedTextureData2D(QSSGRenderBackendTextureObject to,
                                    QSSGRenderTextureTargetType target,
                                    qint32 level,
                                    QSSGRenderTextureFormat internalFormat,
                                    qint32 width,
                                    qint32 height,
                                    qint32 border,
                                    QSSGByteView hostData) override;
    void setCompressedTextureDataCubeFace(QSSGRenderBackendTextureObject to,
                                          QSSGRenderTextureTargetType target,
                                          qint32 level,
                                          QSSGRenderTextureFormat internalFormat,
                                          qint32 width,
                                          qint32 height,
                                          qint32 border,
                                          QSSGByteView hostData) override;
    void setCompressedTextureSubData2D(QSSGRenderBackendTextureObject to,
                                       QSSGRenderTextureTargetType target,
                                       qint32 level,
                                       qint32 xOffset,
                                       qint32 yOffset,
                                       qint32 width,
                                       qint32 height,
                                       QSSGRenderTextureFormat format,
                                       QSSGByteView hostData) override;
    void setMultisampledTextureData2D(QSSGRenderBackendTextureObject to,
                                      QSSGRenderTextureTargetType target,
                                      qint32 samples,
                                      QSSGRenderTextureFormat internalFormat,
                                      qint32 width,
                                      qint32 height,
                                      bool fixedsamplelocations) override = 0;

    void setTextureData3D(QSSGRenderBackendTextureObject to,
                          QSSGRenderTextureTargetType target,
                          qint32 level,
                          QSSGRenderTextureFormat internalFormat,
                          qint32 width,
                          qint32 height,
                          qint32 depth,
                          qint32 border,
                          QSSGRenderTextureFormat format,
                          QSSGByteView hostData) override;

    void generateMipMaps(QSSGRenderBackendTextureObject to,
                         QSSGRenderTextureTargetType target,
                         QSSGRenderHint genType) override;

    virtual QSSGRenderTextureSwizzleMode getTextureSwizzleMode(const QSSGRenderTextureFormat inFormat) const override;

    QSSGRenderBackendSamplerObject createSampler(
            QSSGRenderTextureMinifyingOp minFilter = QSSGRenderTextureMinifyingOp::Linear,
            QSSGRenderTextureMagnifyingOp magFilter = QSSGRenderTextureMagnifyingOp::Linear,
            QSSGRenderTextureCoordOp wrapS = QSSGRenderTextureCoordOp::ClampToEdge,
            QSSGRenderTextureCoordOp wrapT = QSSGRenderTextureCoordOp::ClampToEdge,
            QSSGRenderTextureCoordOp wrapR = QSSGRenderTextureCoordOp::ClampToEdge,
            qint32 minLod = -1000,
            qint32 maxLod = 1000,
            float lodBias = 0.0,
            QSSGRenderTextureCompareMode compareMode = QSSGRenderTextureCompareMode::NoCompare,
            QSSGRenderTextureCompareOp compareFunc = QSSGRenderTextureCompareOp::LessThanOrEqual,
            float anisotropy = 1.0,
            float *borderColor = nullptr) override;

    void updateSampler(QSSGRenderBackendSamplerObject so,
                       QSSGRenderTextureTargetType target,
                       QSSGRenderTextureMinifyingOp minFilter = QSSGRenderTextureMinifyingOp::Linear,
                       QSSGRenderTextureMagnifyingOp magFilter = QSSGRenderTextureMagnifyingOp::Linear,
                       QSSGRenderTextureCoordOp wrapS = QSSGRenderTextureCoordOp::ClampToEdge,
                       QSSGRenderTextureCoordOp wrapT = QSSGRenderTextureCoordOp::ClampToEdge,
                       QSSGRenderTextureCoordOp wrapR = QSSGRenderTextureCoordOp::ClampToEdge,
                       float minLod = -1000.0,
                       float maxLod = 1000.0,
                       float lodBias = 0.0,
                       QSSGRenderTextureCompareMode compareMode = QSSGRenderTextureCompareMode::NoCompare,
                       QSSGRenderTextureCompareOp compareFunc = QSSGRenderTextureCompareOp::LessThanOrEqual,
                       float anisotropy = 1.0,
                       float *borderColor = nullptr) override;

    void updateTextureObject(QSSGRenderBackendTextureObject to,
                             QSSGRenderTextureTargetType target,
                             qint32 baseLevel,
                             qint32 maxLevel) override;

    void updateTextureSwizzle(QSSGRenderBackendTextureObject to,
                              QSSGRenderTextureTargetType target,
                              QSSGRenderTextureSwizzleMode swizzleMode) override;

    void releaseSampler(QSSGRenderBackendSamplerObject so) override;

    virtual QSSGRenderBackendAttribLayoutObject createAttribLayout(QSSGDataView<QSSGRenderVertexBufferEntry> attribs) override;
    void releaseAttribLayout(QSSGRenderBackendAttribLayoutObject ao) override;

    virtual QSSGRenderBackendInputAssemblerObject createInputAssembler(QSSGRenderBackendAttribLayoutObject attribLayout,
                                                                         QSSGDataView<QSSGRenderBackendBufferObject> buffers,
                                                                         const QSSGRenderBackendBufferObject indexBuffer,
                                                                         QSSGDataView<quint32> strides,
                                                                         QSSGDataView<quint32> offsets,
                                                                         quint32 patchVertexCount) override;
    void releaseInputAssembler(QSSGRenderBackendInputAssemblerObject iao) override;

    bool setInputAssembler(QSSGRenderBackendInputAssemblerObject iao, QSSGRenderBackendShaderProgramObject po) override = 0;
    void resetStates() override;
    void setPatchVertexCount(QSSGRenderBackendInputAssemblerObject, quint32) override { Q_ASSERT(false); }

    // shader
    virtual QSSGRenderBackendVertexShaderObject createVertexShader(QSSGByteView source,
                                                                     QByteArray &errorMessage,
                                                                     bool binary) override;
    virtual QSSGRenderBackendFragmentShaderObject createFragmentShader(QSSGByteView source,
                                                                         QByteArray &errorMessage,
                                                                         bool binary) override;
    virtual QSSGRenderBackendTessControlShaderObject createTessControlShader(QSSGByteView source,
                                                                               QByteArray &errorMessage,
                                                                               bool binary) override;
    virtual QSSGRenderBackendTessEvaluationShaderObject createTessEvaluationShader(QSSGByteView source,
                                                                                     QByteArray &errorMessage,
                                                                                     bool binary) override;
    virtual QSSGRenderBackendGeometryShaderObject createGeometryShader(QSSGByteView source,
                                                                         QByteArray &errorMessage,
                                                                         bool binary) override;
    virtual QSSGRenderBackendComputeShaderObject createComputeShader(QSSGByteView source,
                                                                       QByteArray &errorMessage,
                                                                       bool binary) override;
    void releaseVertexShader(QSSGRenderBackendVertexShaderObject vso) override;
    void releaseFragmentShader(QSSGRenderBackendFragmentShaderObject fso) override;
    void releaseTessControlShader(QSSGRenderBackendTessControlShaderObject tcso) override;
    void releaseTessEvaluationShader(QSSGRenderBackendTessEvaluationShaderObject teso) override;
    void releaseGeometryShader(QSSGRenderBackendGeometryShaderObject gso) override;
    void releaseComputeShader(QSSGRenderBackendComputeShaderObject cso) override;
    void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendVertexShaderObject vso) override;
    void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendFragmentShaderObject fso) override;
    void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendTessControlShaderObject tcso) override;
    void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendTessEvaluationShaderObject teso) override;
    void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendGeometryShaderObject gso) override;
    void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendComputeShaderObject cso) override;
    void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendVertexShaderObject vso) override;
    void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendFragmentShaderObject fso) override;
    void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendTessControlShaderObject tcso) override;
    void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendTessEvaluationShaderObject teso) override;
    void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendGeometryShaderObject gso) override;
    void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendComputeShaderObject cso) override;
    QSSGRenderBackendShaderProgramObject createShaderProgram(bool isSeparable) override;
    void releaseShaderProgram(QSSGRenderBackendShaderProgramObject po) override;
    bool linkProgram(QSSGRenderBackendShaderProgramObject po, QByteArray &errorMessage) override;
    bool linkProgram(QSSGRenderBackendShaderProgramObject po, QByteArray &errorMessage,
                     quint32 format, const QByteArray &binary) override;
    void getProgramBinary(QSSGRenderBackendShaderProgramObject po, quint32 &format, QByteArray &binary) override;
    void setActiveProgram(QSSGRenderBackendShaderProgramObject po) override;
    void dispatchCompute(QSSGRenderBackendShaderProgramObject po, quint32 numGroupsX, quint32 numGroupsY, quint32 numGroupsZ) override;
    QSSGRenderBackendProgramPipeline createProgramPipeline() override;
    void releaseProgramPipeline(QSSGRenderBackendProgramPipeline po) override;
    void setActiveProgramPipeline(QSSGRenderBackendProgramPipeline po) override;
    void setProgramStages(QSSGRenderBackendProgramPipeline ppo,
                          QSSGRenderShaderTypeFlags flags,
                          QSSGRenderBackendShaderProgramObject po) override;

    // uniforms
    qint32 getConstantCount(QSSGRenderBackendShaderProgramObject po) override;
    qint32 getConstantInfoByID(QSSGRenderBackendShaderProgramObject po,
                               quint32 id,
                               quint32 bufSize,
                               qint32 *numElem,
                               QSSGRenderShaderDataType *type,
                               qint32 *binding,
                               char *nameBuf) override;
    void setConstantValue(QSSGRenderBackendShaderProgramObject po,
                          quint32 id,
                          QSSGRenderShaderDataType type,
                          qint32 count,
                          const void *value,
                          bool transpose) override;

    // uniform buffers
    qint32 getConstantBufferCount(QSSGRenderBackendShaderProgramObject po) override;
    qint32 getConstantBufferInfoByID(QSSGRenderBackendShaderProgramObject po,
                                     quint32 id,
                                     quint32 nameBufSize,
                                     qint32 *paramCount,
                                     qint32 *bufferSize,
                                     qint32 *length,
                                     char *nameBuf) override;
    void getConstantBufferParamIndices(QSSGRenderBackendShaderProgramObject po, quint32 id, qint32 *indices) override;
    void getConstantBufferParamInfoByIndices(QSSGRenderBackendShaderProgramObject po,
                                             quint32 count,
                                             quint32 *indices,
                                             QSSGRenderShaderDataType *type,
                                             qint32 *size,
                                             qint32 *offset) override;
    void programSetConstantBlock(QSSGRenderBackendShaderProgramObject po, quint32 blockIndex, quint32 binding) override;
    void programSetConstantBuffer(quint32 index, QSSGRenderBackendBufferObject bo) override;

    // storage buffers
    qint32 getStorageBufferCount(QSSGRenderBackendShaderProgramObject po) override;
    qint32 getStorageBufferInfoByID(QSSGRenderBackendShaderProgramObject po,
                                    quint32 id,
                                    quint32 nameBufSize,
                                    qint32 *paramCount,
                                    qint32 *bufferSize,
                                    qint32 *length,
                                    char *nameBuf) override;
    void programSetStorageBuffer(quint32 index, QSSGRenderBackendBufferObject bo) override;

    /// draw calls
    void draw(QSSGRenderDrawMode drawMode, quint32 start, quint32 count) override;
    void drawIndexed(QSSGRenderDrawMode drawMode, quint32 count, QSSGRenderComponentType type, const void *indices) override;

    // read calls
    void readPixel(QSSGRenderBackendRenderTargetObject rto,
                   qint32 x,
                   qint32 y,
                   qint32 width,
                   qint32 height,
                   QSSGRenderReadPixelFormat inFormat,
                   QSSGByteRef pixels) override;

    QSurfaceFormat format() const override { return m_format; }

private:
    const char *getShadingLanguageVersionString();
    const char *getVersionString();
    const char *getVendorString();
    const char *getRendererString();
    void getAttributes(QSSGRenderBackendShaderProgramGL *pProgram);

protected:
    const char *getExtensionString(); // Used to resolve caps in the different backends

private:
    static const qint32 ACTIVATED_TEXTURE_UNIT_UNKNOWN = -1;

protected:
    virtual bool compileSource(GLuint shaderID, QSSGByteView source, QByteArray &errorMessage, bool binary);
    virtual void setAndInspectHardwareCaps();

    GLConversion m_conversion; ///< Class for conversion from base type to GL types
    QList<QByteArray> m_extensions; ///< contains the OpenGL extension string
    qint32 m_maxAttribCount; ///< Maximum attributes which can be used
    qint32 m_usedAttribCount; ///< Number of attributes which have possibly been used
    qint32 m_activatedTextureUnit = ACTIVATED_TEXTURE_UNIT_UNKNOWN; ///< Activated Texture Unit
    QVector<GLenum> m_drawBuffersArray; ///< Contains the drawbuffer enums
    QSurfaceFormat m_format;

    QSSGRenderBackendRasterizerStateGL *m_currentRasterizerState = nullptr; ///< this holds the current rasterizer state
    QSSGRenderBackendDepthStencilStateGL *m_currentDepthStencilState = nullptr; ///< this holds the current depth stencil state

#ifdef RENDER_BACKEND_LOG_GL_ERRORS
    void checkGLError(const char *function, const char *file, const unsigned int line) const
    {
        GLenum error = m_glFunctions->glGetError();
        if (error != GL_NO_ERROR) {
            qCCritical(RENDER_GL_ERROR) << GLConversion::processGLError(error) << " " << function << " " << file << " " << line;
        }
    }
#else
    void checkGLError() const
    {
#ifdef QT_DEBUG
        const GLenum error = m_glFunctions->glGetError();
        if (error != GL_NO_ERROR)
            qCCritical(RENDER_GL_ERROR, "GL Error: %s", GLConversion::processGLError(error));
#endif
    }
#endif
    QOpenGLFunctions *m_glFunctions = nullptr;
    QOpenGLExtraFunctions *m_glExtraFunctions = nullptr;
};

QT_END_NAMESPACE

#endif
