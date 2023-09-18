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

#ifndef QSSG_RENDER_BACKEND_H
#define QSSG_RENDER_BACKEND_H

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

#include <QtQuick3DRender/private/qtquick3drenderglobal_p.h>
#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>
#include <QtQuick3DUtils/private/qssgbounds3_p.h>

#include <QtGui/qsurfaceformat.h>

QT_BEGIN_NAMESPACE

#define HandleToID_cast(staticType, dynamicType, handle) static_cast<staticType>(reinterpret_cast<dynamicType>(handle))

class Q_QUICK3DRENDER_EXPORT QSSGRenderBackend
{
    QSSGRenderBackend(const QSSGRenderBackend &) = delete;
    QSSGRenderBackend &operator=(const QSSGRenderBackend &) = delete;

public:
    QAtomicInt ref;

    QSSGRenderBackend() = default;
    virtual ~QSSGRenderBackend() {}
    /// opaque buffer object handle
    typedef struct _QSSGRenderBackendBufferObject *QSSGRenderBackendBufferObject;
    /// opaque attribute layout object handle
    typedef struct _QSSGRenderBackendAttribLayoutObject *QSSGRenderBackendAttribLayoutObject;
    /// opaque input assembler object handle
    typedef struct _QSSGRenderBackendInputAssemblerObject *QSSGRenderBackendInputAssemblerObject;
    /// opaque texture object handle
    typedef struct _QSSGRenderBackendTextureObject *QSSGRenderBackendTextureObject;
    /// opaque sampler object handle
    typedef struct _QSSGRenderBackendSamplerObject *QSSGRenderBackendSamplerObject;
    /// opaque renderbuffer object handle
    typedef struct _QSSGRenderBackendRenderbufferObject *QSSGRenderBackendRenderbufferObject;
    /// opaque framebuffer object handle
    typedef struct _QSSGRenderBackendRenderTargetObject *QSSGRenderBackendRenderTargetObject;
    /// opaque vertex shader object handle
    typedef struct _QSSGRenderBackendVertexShaderObject *QSSGRenderBackendVertexShaderObject;
    /// opaque fragment shader object handle
    typedef struct _QSSGRenderBackendFragmentShaderObject *QSSGRenderBackendFragmentShaderObject;
    /// opaque tesselation control shader object handle
    typedef struct _QSSGRenderBackendTessControlShaderObject *QSSGRenderBackendTessControlShaderObject;
    /// opaque tesselation evaluation shader object handle
    typedef struct _QSSGRenderBackendTessEvaluationShaderObject *QSSGRenderBackendTessEvaluationShaderObject;
    /// opaque geometry shader object handle
    typedef struct _QSSGRenderBackendGeometryShaderObject *QSSGRenderBackendGeometryShaderObject;
    /// opaque compute shader object handle
    typedef struct _QSSGRenderBackendComputeShaderObject *QSSGRenderBackendComputeShaderObject;
    /// opaque shader program object handle
    typedef struct _QSSGRenderBackendShaderProgramObject *QSSGRenderBackendShaderProgramObject;
    /// opaque depth stencil state object handle
    typedef struct _QSSGRenderBackendDepthStencilStateObject *QSSGRenderBackendDepthStencilStateObject;
    /// opaque rasterizer state object handle
    typedef struct _QSSGRenderBackendRasterizerStateObject *QSSGRenderBackendRasterizerStateObject;
    /// opaque query object handle
    typedef struct _QSSGRenderBackendQueryObject *QSSGRenderBackendQueryObject;
    /// opaque sync object handle
    typedef struct _QSSGRenderBackendSyncObject *QSSGRenderBackendSyncObject;
    /// opaque sync object handle
    typedef struct _QSSGRenderBackendProgramPipeline *QSSGRenderBackendProgramPipeline;
    /// opaque sync object handle
    typedef struct _QSSGRenderBackendPathObject *QSSGRenderBackendPathObject;

    // backend capability caps
    enum class QSSGRenderBackendCaps
    {
        ConstantBuffer, ///< Constant buffer support query
        DepthStencilTexture, ///< depth stencil texture format suport query
        DxtImages, ///< DXT image support query
        FpRenderTarget, ///< render to floating point target support query
        MsTexture, ///< Multisample texture support query
        TexSwizzle, ///< Texture swizzle support query
        FastBlits, ///< Hardware supports fast blits
        Tessellation, ///< Hardware supports tessellation
        Compute, ///< Hardware supports compute shader
        Geometry, ///< Hardware supports geometry shader
        SampleQuery, ///< Hardware supports query calls of type samples
        TimerQuery, ///< Hardware supports query calls of type timer
        CommandSync, ///< Hardware supports command sync object
        TextureArray, ///< Hardware supports texture arrays
        StorageBuffer, ///< Hardware supports shader storage buffers
        ShaderImageLoadStore, ///< Hardware supports shader image load store operations
        ProgramPipeline, ///< Driver supports separate programs
        AdvancedBlend, ///< Driver supports advanced blend modes
        BlendCoherency, ///< Hardware supports blend coherency
        gpuShader5, // for high precision sampling
        AdvancedBlendKHR, ///< Driver supports advanced blend modes
        VertexArrayObject,
        StandardDerivatives,
        TextureLod
    };

    // backend queries
    enum class QSSGRenderBackendQuery
    {
        MaxTextureSize, ///< Return max supported texture size
        MaxTextureArrayLayers, ///< Return max supported layer count for texture arrays
        MaxConstantBufferSlots, ///< Return max supported constant buffe slots for a single
        /// shader stage
        MaxConstantBufferBlockSize ///< Return max supported size for a single constant
        /// buffer block
    };

    /// backend interface

    /**
     * @brief get the backend type
     *
     * @return true backend type
     */
    virtual QSSGRenderContextType getRenderContextType() const = 0;

    /**
     * @brief get the version of the shading language
     * @return version string.
     */
    virtual QByteArray getShadingLanguageVersion() = 0;

    /**
     * @brief get maximum supported texture image units that
     *	can be used to access texture maps from the vertex shader and the fragment processor
     *combined.
     *
     * @return max texture size
     */
    virtual qint32 getMaxCombinedTextureUnits() = 0;

    /**
     * @brief query Backend capabilities
     *
     * @param[in] inCap		CAPS flag to query
     *						@ConstantBuffer, @DepthStencilTexture, ...
     *
     * @return true if supported
     */
    virtual bool getRenderBackendCap(QSSGRenderBackendCaps inCap) const = 0;

    /**
     * @brief query Backend values
     *
     * @param[in] inQuery		Query flag to get value  for
     *							@MaxTextureSize, @MaxTextureArrayLayers,
     *...
     * @param[in/out] params	the query result is stored here
     *
     * @return no return
     */
    virtual void getRenderBackendValue(QSSGRenderBackendQuery inQuery, qint32 *params) const = 0;

    /**
     * @brief query for bit depth of the depth buffer
     *
     * @return depth buffer bitplanes
     */
    virtual qint32 getDepthBits() const = 0;

    /**
     * @brief query for bit depth of the stencil buffer
     *
     * @return stencil buffer bitplanes
     */
    virtual qint32 getStencilBits() const = 0;

    /**
     * @brief query for maximum supported number of samples for multisampling.
     *
     * @return maximum samples. The value must be at least 4.
     */
    virtual qint32 getMaxSamples() const = 0;

    /*
     * @brief set a backend rende state
     *
     * @param[in] bEnable	enable/disable state
     * @param[in] value		type of state
     *
     * @return no return
     */
    virtual void setRenderState(bool bEnable, const QSSGRenderState value) = 0;

    /**
     * @brief get a backend rende state
     *
     * @param[in] value		type of state
     *
     * @return  true if state enabled otherwise false
     */
    virtual bool getRenderState(const QSSGRenderState value) = 0;

    /**
     * @brief get current depth function
     *
     * @return  active depth function
     */
    virtual QSSGRenderBoolOp getDepthFunc() = 0;

    /**
     * @brief create a depth stencil state object
     *
     * @param[in] enableDepth			enable depth test
     * @param[in] depthMask				enable depth writes
     * @param[in] depthFunc				depth compare function
     * @param[in] enableStencil			enable stencil test
     * @param[in] stencilFuncFront		stencil setup front faces
     * @param[in] stencilFuncBack		stencil setup back faces
     * @param[in] depthStencilOpFront	depth/stencil operations front faces
     * @param[in] depthStencilOpBack	depth/stencil operations back faces
     *
     * @return  opaque handle to state object
     */
    virtual QSSGRenderBackendDepthStencilStateObject createDepthStencilState(
            bool enableDepth,
            bool depthMask,
            QSSGRenderBoolOp depthFunc,
            bool enableStencil,
            QSSGRenderStencilFunction &stencilFuncFront,
            QSSGRenderStencilFunction &stencilFuncBack,
            QSSGRenderStencilOperation &depthStencilOpFront,
            QSSGRenderStencilOperation &depthStencilOpBack) = 0;

    /**
     * @brief release a depth stencil state object
     *
     * @param[in] depthStencilState		pointer to state object
     *
     * @return  none
     */
    virtual void releaseDepthStencilState(QSSGRenderBackendDepthStencilStateObject depthStencilState) = 0;

    /**
     * @brief create a rasterizer state object
     *
     * @param[in] depthBias			any othe value than 0 enables depth bias
     * @param[in] depthScale		any othe value than 0 enables depth scale
     *
     * @return  opaque handle to state object
     */
    virtual QSSGRenderBackendRasterizerStateObject createRasterizerState(float depthBias,
                                                                         float depthScale) = 0;

    /**
     * @brief release a rasterizer state object
     *
     * @param[in] rasterizerState		pointer to state object
     *
     * @return  none
     */
    virtual void releaseRasterizerState(QSSGRenderBackendRasterizerStateObject rasterizerState) = 0;

    /**
     * @brief set depth stencil state
     *
     * @param[in] depthStencilState		pointer to state object
     *
     * @return  none
     */
    virtual void setDepthStencilState(QSSGRenderBackendDepthStencilStateObject depthStencilState) = 0;

    /**
     * @brief set rasterizer state
     *
     * @param[in] rasterizerState		pointer to state object
     *
     * @return  none
     */
    virtual void setRasterizerState(QSSGRenderBackendRasterizerStateObject rasterizerState) = 0;

    /**
     * @brief set current depth function
     *
     * @param[in] func		type of function
     *
     * @return no return
     */
    virtual void setDepthFunc(const QSSGRenderBoolOp func) = 0;

    /**
     * @brief query if depth write is enabled
     *
     * @return true if enabled
     */
    virtual bool getDepthWrite() = 0;

    /**
     * @brief enable / disable depth writes
     *
     * @param[in] bEnable	true for enable
     *
     * @return no return
     */
    virtual void setDepthWrite(bool bEnable) = 0;

    /**
     * @brief enable / disable color channel writes
     *
     * @param[in] bRed		true for enable red channel
     * @param[in] bGreen	true for enable green channel
     * @param[in] bBlue		true for enable blue channel
     * @param[in] bAlpha	true for enable alpha channel
     *
     * @return no return
     */
    virtual void setColorWrites(bool bRed, bool bGreen, bool bBlue, bool bAlpha) = 0;

    /**
     * @brief enable / disable multisample rendering
     *
     * @param[in] bEnable	true for enable
     *
     * @return no return
     */
    virtual void setMultisample(bool bEnable) = 0;

    /**
     * @brief query blend functions
     *
     * @param[out] pBlendFuncArg	blending functions
     *
     * @return no return
     */
    virtual void getBlendFunc(QSSGRenderBlendFunctionArgument *pBlendFuncArg) = 0;

    /**
     * @brief set blend functions
     *
     * @param[in] pBlendFuncArg	blending functions
     *
     * @return no return
     */
    virtual void setBlendFunc(const QSSGRenderBlendFunctionArgument &blendFuncArg) = 0;

    /**
     * @brief set blend equation
     *
     * @param[in] pBlendEquArg	blending equation
     *
     * @return no return
     */
    virtual void setBlendEquation(const QSSGRenderBlendEquationArgument &pBlendEquArg) = 0;

    /**
     * @brief guarantee blend coherency
     *
     *
     * @return no return
     */
    virtual void setBlendBarrier(void) = 0;

    /**
     * @brief query cull face mode
     *
     * @return cull face mode
     */
    virtual QSSGCullFaceMode getCullFaceMode() = 0;

    /**
     * @brief set cull face mode
     *
     * @param[in] cullFaceMode
     *
     * @return no return
     */
    virtual void setCullFaceMode(const QSSGCullFaceMode cullFaceMode) = 0;

    /**
     * @brief query scissor rectangle
     *
     * @param[out] pRect	contains scissor rect
     *
     * @return no return
     */
    virtual void getScissorRect(QRect *pRect) = 0;

    /**
     * @brief set scissor rectangle
     *
     * @param[out] pRect	contains scissor rect
     *
     * @return no return
     */
    virtual void setScissorRect(const QRect &rect) = 0;

    /**
     * @brief query viewport rectangle
     *
     * @param[out] pRect	contains viewport rect
     *
     * @return no return
     */
    virtual void getViewportRect(QRect *pRect) = 0;

    /**
     * @brief set viewport rectangle
     *
     * @param[out] pRect	contains viewport rect
     *
     * @return no return
     */
    virtual void setViewportRect(const QRect &rect) = 0;

    /**
     * @brief query viewport rectangle
     *
     * @param[in] clearColor	clear color
     *
     * @return no return
     */
    virtual void setClearColor(const QVector4D *pClearColor) = 0;

    /**
     * @brief query viewport rectangle
     *
     * @param[in] flags	clear flags
     *
     * @return no return
     */
    virtual void clear(QSSGRenderClearFlags flags) = 0;

    /**
     * @brief create a buffer object
     *
     * @param[in] size			Size of the buffer
     * @param[in] bindFlags		Where to bind this buffer (e.g. vertex, index, ...)
     *							For OpenGL this should be a single
     *value
     * @param[in] usage			Usage of the buffer (e.g. static, dynamic...)
     * @param[in] hostPtr       A pointer to the buffer data that is allocated by the
     *application.
     *
     * @return The created buffer object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendBufferObject createBuffer(QSSGRenderBufferType bindFlags,
                                                         QSSGRenderBufferUsageType usage,
                                                         QSSGByteView hostData) = 0;

    /**
     * @brief bind a buffer object
     *
     * @param[in] bo			Pointer to buffer object
     * @param[in] bindFlags		Where to bind this buffer (e.g. vertex, index, ...)
     *							For OpenGL this should be a single
     *value
     *
     * @return no return.
     */
    virtual void bindBuffer(QSSGRenderBackendBufferObject bo, QSSGRenderBufferType bindFlags) = 0;

    /**
     * @brief Release a single buffer object
     *
     * @param[in] bo			Pointer to buffer object
     *
     * @return no return.
     */
    virtual void releaseBuffer(QSSGRenderBackendBufferObject bo) = 0;

    /**
     * @brief update a whole buffer object
     *
     * @param[in] bo			Pointer to buffer object
     * @param[in] bindFlags		Where to bind this buffer (e.g. vertex, index, ...)
     *							For OpenGL this should be a single
     *value
     * @param[in] size			Size of the data buffer
     * @param[in] usage			Usage of the buffer (e.g. static, dynamic...)
     * @param[in] data			A pointer to the buffer data that is allocated by the
     *application.
     *
     * @return no return.
     */
    virtual void updateBuffer(QSSGRenderBackendBufferObject bo,
                              QSSGRenderBufferType bindFlags,
                              QSSGRenderBufferUsageType usage,
                              QSSGByteView data) = 0;

    /**
     * @brief update a range of a buffer object
     *
     * @param[in] bo			Pointer to buffer object
     * @param[in] bindFlags		Where to bind this buffer (e.g. vertex, index, ...)
     *							For OpenGL this should be a single
     *value
     * @param[in] size			Size of the data buffer
     * @param[in] usage			Usage of the buffer (e.g. static, dynamic...)
     * @param[in] data			A pointer to the buffer data that is allocated by the
     *application.
     *
     * @return no return.
     */
    virtual void updateBufferRange(QSSGRenderBackendBufferObject bo,
                                   QSSGRenderBufferType bindFlags,
                                   size_t offset,
                                   QSSGByteView data) = 0;

    /**
     * @brief Get a pointer to the buffer data ( GL(ES) >= 3 only )
     *
     * @param[in] bo			Pointer to buffer object
     * @param[in] bindFlags		Where to bind this buffer (e.g. vertex, index, ...)
     *							For OpenGL this should be a single
     *value
     * @param[in] offset		Byte offset into the data buffer
     * @param[in] length		Byte length of mapping size
     * @param[in] access		Access of the buffer (e.g. read, write, ...)
     *
     * @return pointer to mapped data or null.
     */
    virtual void *mapBuffer(QSSGRenderBackendBufferObject bo,
                            QSSGRenderBufferType bindFlags,
                            size_t offset,
                            size_t length,
                            QSSGRenderBufferAccessFlags accessFlags) = 0;

    /**
     * @brief Unmap a previously mapped buffer ( GL(ES) >= 3 only )
     *		  This functions transfers the content to the hardware buffer
     *
     * @param[in] bo			Pointer to buffer object
     * @param[in] bindFlags		Where to bind this buffer (e.g. vertex, index, ...)
     *							For OpenGL this should be a single
     *value
     *
     * @return true if successful
     */
    virtual bool unmapBuffer(QSSGRenderBackendBufferObject bo, QSSGRenderBufferType bindFlags) = 0;

    /**
     * @brief Set a memory barrier
     *
     * @param[in] barriers		Flags for barriers
     *
     * @return no return.
     */
    virtual void setMemoryBarrier(QSSGRenderBufferBarrierFlags barriers) = 0;

    /**
     * @brief create a query object
     *
     * @return The created query object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendQueryObject createQuery() = 0;

    /**
     * @brief delete query objects
     *
     * @param[in] qo		Handle to query object
     *
     * @return  no return
     */
    virtual void releaseQuery(QSSGRenderBackendQueryObject qo) = 0;

    /**
     * @brief Start query recording
     *
     * @param[in] qo		Handle to query object
     * @param[in] type		Type of query
     *
     * @return  no return
     */
    virtual void beginQuery(QSSGRenderBackendQueryObject qo, QSSGRenderQueryType type) = 0;

    /**
     * @brief End query recording
     *
     * @param[in] qo		Handle to query object
     * @param[in] type		Type of query
     *
     * @return  no return
     */
    virtual void endQuery(QSSGRenderBackendQueryObject qo, QSSGRenderQueryType type) = 0;

    /**
     * @brief Get a query result
     *
     * @param[in]  qo		Handle to query object
     * @param[in]  type		Type of query
     * @param[out] params	Contains result of query regarding query type
     *
     * @return  no return
     */
    virtual void getQueryResult(QSSGRenderBackendQueryObject qo, QSSGRenderQueryResultType resultType, quint32 *params) = 0;

    /**
     * @brief Get a query result
     *
     * @param[in]  qo		Handle to query object
     * @param[in]  type		Type of query
     * @param[out] params	Contains result of query regarding query type 64 bit returns
     *
     * @return  no return
     */
    virtual void getQueryResult(QSSGRenderBackendQueryObject qo, QSSGRenderQueryResultType resultType, quint64 *params) = 0;

    /**
     * @brief Record the GPU time using the query object
     *
     * @param[in]  qo		Handle to query object
     *
     * @return  no return
     */
    virtual void setQueryTimer(QSSGRenderBackendQueryObject qo) = 0;

    /**
     * @brief create a sync object and place it in the command queue
     *
     * @param[in] tpye			Type to sync
     * @param[in] syncFlags		Currently unused
     *
     * @return The created sync object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendSyncObject createSync(QSSGRenderSyncType tpye, QSSGRenderSyncFlags syncFlags) = 0;

    /**
     * @brief delete sync object
     *
     * @param[in] so		Handle to sync object
     *
     * @return  no return
     */
    virtual void releaseSync(QSSGRenderBackendSyncObject so) = 0;

    /**
     * @brief wait for sync object to be signaled
     *
     * @param[in] so			Handle to sync object
     * @param[in] syncFlags		Currently unused
     * @param[in] timeout		Currently ignored
     *
     * @return  no return
     */
    virtual void waitSync(QSSGRenderBackendSyncObject so, QSSGRenderCommandFlushFlags syncFlags, quint64 timeout) = 0;

    /**
     * @brief create a render target object
     *
     *
     * @return The created render target object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendRenderTargetObject createRenderTarget() = 0;

    /**
     * @brief Release a single render target object
     *
     * @param[in] rto			Pointer to render target object
     *
     * @return no return.
     */
    virtual void releaseRenderTarget(QSSGRenderBackendRenderTargetObject rto) = 0;

    /**
     * @brief Attach a renderbuffer object to the framebuffer
     *
     * @param[in] rto			Pointer to render target object
     * @param[in] attachment	Attachment point (e.g COLOR0, DEPTH)
     * @param[in] rbo			Pointer to renderbuffer object
     *
     * @return no return.
     */
    virtual void renderTargetAttach(QSSGRenderBackendRenderTargetObject rto,
                                    QSSGRenderFrameBufferAttachment attachment,
                                    QSSGRenderBackendRenderbufferObject rbo) = 0;

    /**
     * @brief Attach a texture object to the render target
     *
     * @param[in] rto			Pointer to render target object
     * @param[in] attachment	Attachment point (e.g COLOR0, DEPTH)
     * @param[in] to			Pointer to texture object
     * @param[in] target		Attachment texture target
     *
     * @return no return.
     */
    virtual void renderTargetAttach(QSSGRenderBackendRenderTargetObject rto,
                                    QSSGRenderFrameBufferAttachment attachment,
                                    QSSGRenderBackendTextureObject to,
                                    QSSGRenderTextureTargetType target = QSSGRenderTextureTargetType::Texture2D) = 0;

    /**
     * @brief Attach a texture object to the render target
     *
     * @param[in] rto			Pointer to render target object
     * @param[in] attachment	Attachment point (e.g COLOR0, DEPTH)
     * @param[in] to			Pointer to texture object
     * @param[in] level			Texture mip level
     * @param[in] layer			Texture layer or slice
     *
     * @return no return.
     */
    virtual void renderTargetAttach(QSSGRenderBackendRenderTargetObject rto,
                                    QSSGRenderFrameBufferAttachment attachment,
                                    QSSGRenderBackendTextureObject to,
                                    qint32 level,
                                    qint32 layer) = 0;

    /**
     * @brief Make a render target active
     *
     * @param[in] rto			Pointer to render target object
     *
     * @return no return.
     */
    virtual void setRenderTarget(QSSGRenderBackendRenderTargetObject rto) = 0;

    /**
     * @brief Check if a render target is ready for render
     *
     * @param[in] rto			Pointer to render target object
     *
     * @return true if usable.
     */
    virtual bool renderTargetIsValid(QSSGRenderBackendRenderTargetObject rto) = 0;

    /**
     * @brief Make a render target active for reading
     *
     * @param[in] rto			Pointer to render target object
     *
     * @return no return.
     */
    virtual void setReadTarget(QSSGRenderBackendRenderTargetObject rto) = 0;

    /**
     * @brief Set active buffers for drawing
     *
     * @param[in] rto				Pointer to render target object
     * @param[in] inDrawBufferSet	Pointer to array of enabled render targets
     *
     * @return no return.
     */
    virtual void setDrawBuffers(QSSGRenderBackendRenderTargetObject rto, QSSGDataView<qint32> inDrawBufferSet) = 0;

    /**
     * @brief Set active buffer for reading
     *
     * @param[in] rto				Pointer to render target object
     * @param[in] inReadFace		Buffer to read from
     *
     * @return no return.
     */
    virtual void setReadBuffer(QSSGRenderBackendRenderTargetObject rto, QSSGReadFace inReadFace) = 0;

    /**
     * @brief Copy framebuffer attachments. Source is set with SetReadTarget dest with
     * SetRenderTarget
     *
     * @param[in] srcX0				Lower left X coord of source rectangle
     * @param[in] srcY0				Lower left Y coord of source rectangle
     * @param[in] srcX1				Upper right X coord of source rectangle
     * @param[in] srcY1				Upper right Y coord of source rectangle
     * @param[in] dstX0				Lower left X coord of dest rectangle
     * @param[in] dstY0				Lower left Y coord of dest rectangle
     * @param[in] dstX1				Upper right X coord of dest rectangle
     * @param[in] dstY1				Upper right Y coord of dest rectangle
     * @param[in] inDrawBufferSet	pointer to array of enabled render targets
     * @param[in] filter			Copy filter method (NEAREST or LINEAR)
     *
     * @return no return.
     */
    virtual void blitFramebuffer(qint32 srcX0,
                                 qint32 srcY0,
                                 qint32 srcX1,
                                 qint32 srcY1,
                                 qint32 dstX0,
                                 qint32 dstY0,
                                 qint32 dstX1,
                                 qint32 dstY1,
                                 QSSGRenderClearFlags flags,
                                 QSSGRenderTextureMagnifyingOp filter) = 0;

    /**
     * @brief Copy framebuffer attachments to texture. Source is set with SetReadTarget dest with
     * SetRenderTarget
     *
     * @param[in] srcX0             Lower left X coord of source rectangle
     * @param[in] srcY0             Lower left Y coord of source rectangle
     * @param[in] srcX1             Width source rectangle
     * @param[in] srcY1             Height source rectangle
     * @param[in] dstX0             Lower left X coord of dest rectangle
     * @param[in] dstY0             Lower left Y coord of dest rectangle
     * @param[in] flags             Attachment to copy
     * @param[in] texture           The destination texture
     * @param[in] target            The texture target
     *
     * @return no return.
     */
    virtual void copyFramebufferTexture(qint32 srcX0,
                                        qint32 srcY0,
                                        qint32 width,
                                        qint32 height,
                                        qint32 dstX0,
                                        qint32 dstY0,
                                        QSSGRenderBackendTextureObject texture,
                                        QSSGRenderTextureTargetType target
                                            = QSSGRenderTextureTargetType::Texture2D) = 0;

    /**
     * @brief create a render buffer object
     *
     * @param[in] storageFormat	Format of the buffer
     * @param[in] width			Buffer with
     * @param[in] height		Buffer height
     *
     * @return The created render buffer object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendRenderbufferObject createRenderbuffer(QSSGRenderRenderBufferFormat storageFormat,
                                                                     qint32 width,
                                                                     qint32 height) = 0;

    /**
     * @brief Release a single renderbuffer object
     *
     * @param[in] bo			Pointer to renderbuffer object
     *
     * @return no return.
     */
    virtual void releaseRenderbuffer(QSSGRenderBackendRenderbufferObject rbo) = 0;

    /**
     * @brief resize a render buffer object
     *
     * @param[in] storageFormat	Format of the buffer
     * @param[in] width			Buffer with
     * @param[in] height		Buffer height
     *
     * @return True on success
     */
    virtual bool resizeRenderbuffer(QSSGRenderBackendRenderbufferObject rbo,
                                    QSSGRenderRenderBufferFormat storageFormat,
                                    qint32 width,
                                    qint32 height) = 0;

    /**
     * @brief create a texture object
     *
     * @return The created texture object or nullptr if the creation failed..
     */
    virtual QSSGRenderBackendTextureObject createTexture() = 0;

    /**
     * @brief set texture data for a 2D texture
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target 2D, 3D
     * @param[in] level				Texture mip level
     * @param[in] internalFormat	format of the texture
     * @param[in] width				texture width
     * @param[in] height			texture height
     * @param[in] border			border
     * @param[in] format			format of provided pixel data
     * @param[in] hostPtr			A pointer to the buffer data that is allocated by the
     * application.
     *
     * @return No return
     */
    virtual void setTextureData2D(QSSGRenderBackendTextureObject to,
                                  QSSGRenderTextureTargetType target,
                                  qint32 level,
                                  QSSGRenderTextureFormat internalFormat,
                                  qint32 width,
                                  qint32 height,
                                  qint32 border,
                                  QSSGRenderTextureFormat format,
                                  QSSGByteView hostData) = 0;

    /**
     * @brief set texture data for the face of a Cube map
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target face
     * @param[in] level				Texture mip level
     * @param[in] internalFormat	format of the texture
     * @param[in] width				texture width
     * @param[in] height			texture height
     * @param[in] border			border
     * @param[in] format			format of provided pixel data
     * @param[in] hostPtr			A pointer to the buffer data that is allocated by the
     * application.
     *
     * @return No return
     */
    virtual void setTextureDataCubeFace(QSSGRenderBackendTextureObject to,
                                        QSSGRenderTextureTargetType target,
                                        qint32 level,
                                        QSSGRenderTextureFormat internalFormat,
                                        qint32 width,
                                        qint32 height,
                                        qint32 border,
                                        QSSGRenderTextureFormat format,
                                        QSSGByteView hostData) = 0;

    /**
     * @brief create a storage for a 2D texture including mip levels
     *		  Note that this makes texture immutable in size and format
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target 2D, 3D
     * @param[in] levels			Texture mip level count
     * @param[in] internalFormat	format of the texture
     * @param[in] width				texture width
     * @param[in] height			texture height
     *
     * @return No return
     */
    virtual void createTextureStorage2D(QSSGRenderBackendTextureObject to,
                                        QSSGRenderTextureTargetType target,
                                        qint32 levels,
                                        QSSGRenderTextureFormat internalFormat,
                                        qint32 width,
                                        qint32 height) = 0;

    /**
     * @brief set texture sub data for a 2D texture
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target 2D, 3D
     * @param[in] level				Texture mip level
     * @param[in] xOffset			Texture x offset
     * @param[in] yOffset			Texture y offset
     * @param[in] width				Texture width
     * @param[in] height			Texture height
     * @param[in] border			border
     * @param[in] format			format of texture
     * @param[in] hostPtr			A pointer to the buffer data that is allocated by the
     * application.
     *
     * @return No return
     */
    virtual void setTextureSubData2D(QSSGRenderBackendTextureObject to,
                                     QSSGRenderTextureTargetType target,
                                     qint32 level,
                                     qint32 xOffset,
                                     qint32 yOffset,
                                     qint32 width,
                                     qint32 height,
                                     QSSGRenderTextureFormat format,
                                     QSSGByteView hostData) = 0;

    /**
     * @brief set compressed texture data for a 2D texture
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target 2D, 3D
     * @param[in] level				Texture mip level
     * @param[in] internalFormat	format of the texture
     * @param[in] width				texture width
     * @param[in] height			texture height
     * @param[in] border			border
     * @param[in] imageSize			image size in bytes located at hostPtr
     * @param[in] hostPtr			A pointer to the buffer data that is allocated by the
     * application.
     *
     * @return No return
     */
    virtual void setCompressedTextureData2D(QSSGRenderBackendTextureObject to,
                                            QSSGRenderTextureTargetType target,
                                            qint32 level,
                                            QSSGRenderTextureFormat internalFormat,
                                            qint32 width,
                                            qint32 height,
                                            qint32 border,
                                            QSSGByteView hostData) = 0;

    /**
     * @brief set compressed texture data for a Cubemap face
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target 2D, 3D
     * @param[in] level				Texture mip level
     * @param[in] internalFormat	format of the texture
     * @param[in] width				texture width
     * @param[in] height			texture height
     * @param[in] border			border
     * @param[in] imageSize			image size in bytes located at hostPtr
     * @param[in] hostPtr			A pointer to the buffer data that is allocated by the
     * application.
     *
     * @return No return
     */
    virtual void setCompressedTextureDataCubeFace(QSSGRenderBackendTextureObject to,
                                                  QSSGRenderTextureTargetType target,
                                                  qint32 level,
                                                  QSSGRenderTextureFormat internalFormat,
                                                  qint32 width,
                                                  qint32 height,
                                                  qint32 border,
                                                  QSSGByteView hostData) = 0;

    /**
     * @brief set compressed texture sub data for a 2D texture
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target 2D, 3D
     * @param[in] level				Texture mip level
     * @param[in] xOffset			Texture x offset
     * @param[in] yOffset			Texture y offset
     * @param[in] width				texture width
     * @param[in] height			texture height
     * @param[in] format			format of provided pixel data
     * @param[in] imageSize			image size in bytes located at hostPtr
     * @param[in] hostPtr			A pointer to the buffer data that is allocated by the
     * application.
     *
     * @return No return
     */
    virtual void setCompressedTextureSubData2D(QSSGRenderBackendTextureObject to,
                                               QSSGRenderTextureTargetType target,
                                               qint32 level,
                                               qint32 xOffset,
                                               qint32 yOffset,
                                               qint32 width,
                                               qint32 height,
                                               QSSGRenderTextureFormat format,
                                               QSSGByteView hostData) = 0;

    /**
     * @brief establish a multisampled 2D texture
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target 2D MS
     * @param[in] samples			Textures sample count
     * @param[in] internalFormat	Format of the texture
     * @param[in] width				Texture width
     * @param[in] height			Texture height
     * @param[in] bool				Fixed sample locations
     *
     * @return No return
     */
    virtual void setMultisampledTextureData2D(QSSGRenderBackendTextureObject to,
                                              QSSGRenderTextureTargetType target,
                                              qint32 samples,
                                              QSSGRenderTextureFormat internalFormat,
                                              qint32 width,
                                              qint32 height,
                                              bool fixedsamplelocations) = 0;

    /**
     * @brief set texture data for a 3D texture or 2D texture array
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target 2D, 3D
     * @param[in] level				Texture mip level
     * @param[in] internalFormat	format of the texture
     * @param[in] width				texture width
     * @param[in] height			texture height
     * @param[in] depth				texture depth or slice count
     * @param[in] border			border
     * @param[in] format			format of provided pixel data
     * @param[in] hostPtr			A pointer to the buffer data that is allocated by the
     * application.
     *
     * @return No return
     */
    virtual void setTextureData3D(QSSGRenderBackendTextureObject to,
                                  QSSGRenderTextureTargetType target,
                                  qint32 level,
                                  QSSGRenderTextureFormat internalFormat,
                                  qint32 width,
                                  qint32 height,
                                  qint32 depth,
                                  qint32 border,
                                  QSSGRenderTextureFormat format,
                                  QSSGByteView hostData) = 0;

    /**
     * @brief generate mipmap levels
     *
     * @param[in] to				Pointer to texture object
     * @param[in] target			Texture target 2D,...
     * @param[in] hint				How to generate mips (Nicest)
     *
     * @return No return
     */
    virtual void generateMipMaps(QSSGRenderBackendTextureObject to,
                                 QSSGRenderTextureTargetType target,
                                 QSSGRenderHint genType) = 0;

    /**
     * @brief bind a texture object
     *
     * @param[in] to			Pointer to texture object
     * @param[in] target		Where to bind this texture (e.g. 2D, 3D, ...)
     * @param[in] unit			Which unit to bind this texture
     *
     * @return no return.
     */
    virtual void bindTexture(QSSGRenderBackendTextureObject to,
                             QSSGRenderTextureTargetType target,
                             qint32 unit) = 0;

    /**
     * @brief select active texture unit
     *
     * @param[in] unit          Which unit to bind for texture
     *
     * @return no return.
     */
    virtual void setActiveTexture(qint32 unit) = 0;

    /**
     * @brief bind a image/texture object
     *
     * @param[in] to			Pointer to texture object
     * @param[in] unit			Which unit to bind this texture
     * @param[in] level			Which level to bind
     * @param[in] layered		Bind layered texture (cube map, array,... )
     * @param[in] level			Specify layer. Only valid of layered=false.
     * @param[in] access		Access mode ( read, write, read-write )
     * @param[in] format		Texture format must be compatible with Image format
     *
     * @return no return.
     */
    virtual void bindImageTexture(QSSGRenderBackendTextureObject to,
                                  quint32 unit,
                                  qint32 level,
                                  bool layered,
                                  qint32 layer,
                                  QSSGRenderImageAccessType accessFlags,
                                  QSSGRenderTextureFormat format) = 0;

    /**
     * @brief Release a single texture object
     *
     * @param[in] to			Pointer to buffer object
     *
     * @return no return.
     */
    virtual void releaseTexture(QSSGRenderBackendTextureObject to) = 0;

    /**
     * @brief query texture swizzle mode
     *		  This is mainly for luminance, alpha replacement with R8 formats
     *
     * @param[in]  inFormat			input texture format to check
     *
     * @return texture swizzle mode
     */
    virtual QSSGRenderTextureSwizzleMode getTextureSwizzleMode(const QSSGRenderTextureFormat inFormat) const = 0;

    /**
     * @ brief create a sampler
     *
     * @param[in] minFilter		Texture min filter
     * @param[in] magFilter		Texture mag filter
     * @param[in] wrapS			Texture coord generation for S
     * @param[in] wrapT			Texture coord generation for T
     * @param[in] wrapR			Texture coord generation for R
     * @param[in] minLod		Texture min level of detail
     * @param[in] maxLod		Texture max level of detail
     * @param[in] lodBias		Texture level of detail example
     * @param[in] compareMode	Texture compare mode
     * @param[in] compareFunc	Texture compare function
     * @param[in] anisoFilter	Aniso filter value [1.0, 16.0]
     * @param[in] borderColor	Texture border color float[4]
     *
     * @return The created sampler object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendSamplerObject createSampler(
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
            float *borderColor = nullptr) = 0;

    /**
     * @ brief update a sampler
     *
     * @param[in] so			Pointer to sampler object
     * @param[in] target		Texture target 2D, 3D
     * @param[in] minFilter		Texture min filter
     * @param[in] magFilter		Texture mag filter
     * @param[in] wrapS			Texture coord generation for S
     * @param[in] wrapT			Texture coord generation for T
     * @param[in] wrapR			Texture coord generation for R
     * @param[in] minLod		Texture min level of detail
     * @param[in] maxLod		Texture max level of detail
     * @param[in] lodBias		Texture level of detail bias (unused)
     * @param[in] compareMode	Texture compare mode
     * @param[in] compareFunc	Texture compare function
     * @param[in] anisoFilter	Aniso filter value [1.0, 16.0]
     * @param[in] borderColor	Texture border color float[4] (unused)
     *
     * @return No return
     */
    virtual void updateSampler(QSSGRenderBackendSamplerObject so,
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
                               float *borderColor = nullptr) = 0;

    /**
     * @ brief Update a textures swizzle mode
     *
     * @param[in] so			Pointer to texture object
     * @param[in] target		Texture target 2D, 3D
     * @param[in] swizzleMode	Texture swizzle mode
     *
     * @return No return
     */
    virtual void updateTextureSwizzle(QSSGRenderBackendTextureObject to,
                                      QSSGRenderTextureTargetType target,
                                      QSSGRenderTextureSwizzleMode swizzleMode) = 0;

    /**
     * @ brief Update state belonging to a texture object
     *
     * @param[in] so			Pointer to texture object
     * @param[in] target		Texture target 2D, 3D, Cube
     * @param[in] baseLevel		Texture base level
     * @param[in] maxLevel		Texture max level
     *
     * @return No return
     */
    virtual void updateTextureObject(QSSGRenderBackendTextureObject to,
                                     QSSGRenderTextureTargetType target,
                                     qint32 baseLevel,
                                     qint32 maxLevel) = 0;

    /**
     * @brief Release a single sampler object
     *
     * @param[in] so			Pointer to sampler object
     *
     * @return no return.
     */
    virtual void releaseSampler(QSSGRenderBackendSamplerObject so) = 0;

    /**
     * @brief create a attribute layout object
     *
     * @param[in] attribs	Array off vertex attributes.
     *
     * @return The created attribute layout object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendAttribLayoutObject createAttribLayout(QSSGDataView<QSSGRenderVertexBufferEntry> attribs) = 0;

    /**
     * @brief Release a attribute layoutr object
     *
     * @param[in] ao			Pointer to attribute layout object
     *
     * @return no return.
     */
    virtual void releaseAttribLayout(QSSGRenderBackendAttribLayoutObject ao) = 0;

    /**
     * @brief create a input assembler object
     *
     * @param[in] attribLayout		Pointer to QSSGRenderBackendAttribLayoutObject object
     * @param[in] buffers			list of vertex buffers
     * @param[in] indexBuffer		index buffer object
     * @param[in] strides			list of strides of the buffer
     * @param[in] offsets			list of offsets into the buffer
     * @param[in] patchVertexCount	vertext count for a patch. Only valid for patch primitives
     *
     * @return The created input assembler object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendInputAssemblerObject createInputAssembler(QSSGRenderBackendAttribLayoutObject attribLayout,
                                                                         QSSGDataView<QSSGRenderBackendBufferObject> buffers,
                                                                         const QSSGRenderBackendBufferObject indexBuffer,
                                                                         QSSGDataView<quint32> strides,
                                                                         QSSGDataView<quint32> offsets,
                                                                         quint32 patchVertexCount) = 0;

    /**
     * @brief Release a input assembler object
     *
     * @param[in] iao					Pointer to attribute layout object
     *
     * @return no return.
     */
    virtual void releaseInputAssembler(QSSGRenderBackendInputAssemblerObject iao) = 0;

    /**
     * @brief Set a input assembler object.
     *		  This setup the render engine vertex assmebly
     *
     * @param[in] iao					Pointer to attribute layout object
     * @param[in] po					Pointer program object
     *
     * @return false if it fails.
     */
    virtual bool setInputAssembler(QSSGRenderBackendInputAssemblerObject iao, QSSGRenderBackendShaderProgramObject po) = 0;

    /**
     * @brief Reset all the states cached in the backend
     *        so as to ensure that backend commands to set necessary states are called.
     *        This will be called every frame just before QtQuick3D rendering commands are called
     *        because such states might have been changed outside of QtQuick3D.
     */
    virtual void resetStates() = 0;

    /**
     * @brief Set the per patch vertex count
     *
     * @param[in] iao					Pointer to attribute layout object
     * @param[in] count					Count of vertices per patch
     *
     * @return false if it fails.
     */
    virtual void setPatchVertexCount(QSSGRenderBackendInputAssemblerObject iao, quint32 count) = 0;

    /**
     * @brief create a vertex shader object
     *
     * @param[in] source			Pointer to shader source
     * @param[in/out] errorMessage	Pointer to copy the error message
     * @param[in] binary			True if the source is actually a binary program
     *
     * @return The created vertex shader object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendVertexShaderObject createVertexShader(QSSGByteView source,
                                                                     QByteArray &errorMessage,
                                                                     bool binary) = 0;

    /**
     * @brief release a vertex shader object
     *
     * @param[in] vso		Pointer to vertex shader object
     *
     * @return No Return.
     */
    virtual void releaseVertexShader(QSSGRenderBackendVertexShaderObject vso) = 0;

    /**
     * @brief create a fragment shader object
     *
     * @param[in] source			Pointer to shader source
     * @param[in/out] errorMessage	Pointer to copy the error message
     * @param[in] binary			True if the source is actually a binary program
     *
     * @return The created vertex shader object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendFragmentShaderObject createFragmentShader(QSSGByteView source,
                                                                         QByteArray &errorMessage,
                                                                         bool binary) = 0;

    /**
     * @brief release a fragment shader object
     *
     * @param[in] vso		Pointer to fragment shader object
     *
     * @return No Return.
     */
    virtual void releaseFragmentShader(QSSGRenderBackendFragmentShaderObject fso) = 0;

    /**
     * @brief create a tessellation control shader object
     *
     * @param[in] source			Pointer to shader source
     * @param[in/out] errorMessage	Pointer to copy the error message
     * @param[in] binary			True if the source is actually a binary program
     *
     * @return The created tessellation control shader object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendTessControlShaderObject createTessControlShader(QSSGByteView source,
                                                                               QByteArray &errorMessage,
                                                                               bool binary) = 0;

    /**
     * @brief release a tessellation control shader object
     *
     * @param[in] tcso		Pointer to tessellation control shader object
     *
     * @return No Return.
     */
    virtual void releaseTessControlShader(QSSGRenderBackendTessControlShaderObject tcso) = 0;

    /**
     * @brief create a tessellation evaluation shader object
     *
     * @param[in] source			Pointer to shader source
     * @param[in/out] errorMessage	Pointer to copy the error message
     * @param[in] binary			True if the source is actually a binary program
     *
     * @return The created tessellation evaluation shader object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendTessEvaluationShaderObject createTessEvaluationShader(QSSGByteView source,
                                                                                     QByteArray &errorMessage,
                                                                                     bool binary) = 0;

    /**
     * @brief release a tessellation evaluation shader object
     *
     * @param[in] tcso		Pointer to tessellation evaluation shader object
     *
     * @return No Return.
     */
    virtual void releaseTessEvaluationShader(QSSGRenderBackendTessEvaluationShaderObject teso) = 0;

    /**
     * @brief create a geometry shader object
     *
     * @param[in] source			Pointer to shader source
     * @param[in/out] errorMessage	Pointer to copy the error message
     * @param[in] binary			True if the source is actually a binary program
     *
     * @return The created geometry shader object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendGeometryShaderObject createGeometryShader(QSSGByteView source,
                                                                         QByteArray &errorMessage,
                                                                         bool binary) = 0;

    /**
     * @brief release a geometry shader object
     *
     * @param[in] tcso		Pointer to geometry shader object
     *
     * @return No Return.
     */
    virtual void releaseGeometryShader(QSSGRenderBackendGeometryShaderObject gso) = 0;

    /**
     * @brief create a compute shader object
     *
     * @param[in] source			Pointer to shader source
     * @param[in/out] errorMessage	Pointer to copy the error message
     * @param[in] binary			True if the source is actually a binary program
     *
     * @return The created compute shader object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendComputeShaderObject createComputeShader(QSSGByteView source,
                                                                       QByteArray &errorMessage,
                                                                       bool binary) = 0;

    /**
     * @brief release a compute shader object
     *
     * @param[in] cso		Pointer to compute shader object
     *
     * @return No Return.
     */
    virtual void releaseComputeShader(QSSGRenderBackendComputeShaderObject cso) = 0;

    /**
     * @brief attach a vertex shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] vso		Pointer to vertex shader object
     *
     * @return No Return.
     */
    virtual void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendVertexShaderObject vso) = 0;

    /**
     * @brief detach a vertex shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] vso		Pointer to vertex shader object
     *
     * @return No Return.
     */
    virtual void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendVertexShaderObject vso) = 0;

    /**
     * @brief attach a fragment shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] fso		Pointer to fragment shader object
     *
     * @return No Return.
     */
    virtual void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendFragmentShaderObject fso) = 0;

    /**
     * @brief detach a fragment shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] fso		Pointer to fragment shader object
     *
     * @return No Return.
     */
    virtual void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendFragmentShaderObject fso) = 0;

    /**
     * @brief attach a tessellation control shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] tcso		Pointer to tessellation control shader object
     *
     * @return No Return.
     */
    virtual void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendTessControlShaderObject tcso) = 0;

    /**
     * @brief detach a tessellation control shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] tcso		Pointer to tessellation control shader object
     *
     * @return No Return.
     */
    virtual void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendTessControlShaderObject tcso) = 0;

    /**
     * @brief attach a tessellation evaluation shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] teso		Pointer to tessellation evaluation shader object
     *
     * @return No Return.
     */
    virtual void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendTessEvaluationShaderObject teso) = 0;

    /**
     * @brief detach a tessellation evaluation shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] teso		Pointer to tessellation evaluation shader object
     *
     * @return No Return.
     */
    virtual void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendTessEvaluationShaderObject teso) = 0;

    /**
     * @brief attach a geometry shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] teso		Pointer to geometry shader object
     *
     * @return No Return.
     */
    virtual void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendGeometryShaderObject gso) = 0;

    /**
     * @brief detach a geometry shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] teso		Pointer to geometry shader object
     *
     * @return No Return.
     */
    virtual void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendGeometryShaderObject gso) = 0;

    /**
     * @brief attach a compute shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] cso		Pointer to compute shader object
     *
     * @return No Return.
     */
    virtual void attachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendComputeShaderObject cso) = 0;

    /**
     * @brief detach a compute shader object to a program object
     *
     * @param[in] po		Pointer to program object
     * @param[in] cso		Pointer to compute shader object
     *
     * @return No Return.
     */
    virtual void detachShader(QSSGRenderBackendShaderProgramObject po, QSSGRenderBackendComputeShaderObject cso) = 0;

    /**
     * @brief create a shader program object
     *
     * @param[in] isSeparable   Tell the backend that this program is separable
     *
     * @return The created shader program object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendShaderProgramObject createShaderProgram(bool isSeparable) = 0;

    /**
     * @brief release a shader program object
     *
     * @param[in] po            Pointer to shader program object
     *
     * @return No Return.
     */
    virtual void releaseShaderProgram(QSSGRenderBackendShaderProgramObject po) = 0;

    /**
     * @brief link a shader program object
     *
     * @param[in] po                Pointer to shader program object
     * @param[in/out] errorMessage  Pointer to copy the error message
     *
     * @return True if program is successfully linked.
     */
    virtual bool linkProgram(QSSGRenderBackendShaderProgramObject po, QByteArray &errorMessage) = 0;

    /**
     * @brief link a binary shader program object
     *
     * @param[in] po                Pointer to shader program object
     * @param[in/out] errorMessage  Pointer to copy the error message
     * @param[in] format            Binary format
     * @param[in] binary            Binary data
     *
     * @return True if program is successfully linked.
     */
    virtual bool linkProgram(QSSGRenderBackendShaderProgramObject po,
                             QByteArray &errorMessage,
                             quint32 format, const QByteArray &binary) = 0;

    /**
     * @brief Get shader program binary
     *
     * @param[in] po            Pointer to shader program object
     * @param[out] format       The format of the program binary
     * @param[out] binary       The program binary data
     *
     */
    virtual void getProgramBinary(QSSGRenderBackendShaderProgramObject po, quint32 &format,
                                  QByteArray &binary) = 0;

    /**
     * @brief Make a program current
     *
     * @param[in] po        Pointer to shader program object
     *
     * @return No return
     */
    virtual void setActiveProgram(QSSGRenderBackendShaderProgramObject po) = 0;

    /**
     * @brief create a program pipeline object
     *
     *
     * @return The created program pipeline object or nullptr if the creation failed.
     */
    virtual QSSGRenderBackendProgramPipeline createProgramPipeline() = 0;

    /**
     * @brief release a program pipeline object
     *
     * @param[in] ppo		Pointer to program pipeline object
     *
     * @return No Return.
     */
    virtual void releaseProgramPipeline(QSSGRenderBackendProgramPipeline ppo) = 0;

    /**
     * @brief Make a program pipeline current
     *
     * @param[in] ppo		Pointer to program pipeline object
     *
     * @return No return
     */
    virtual void setActiveProgramPipeline(QSSGRenderBackendProgramPipeline ppo) = 0;

    /**
     * @brief Make a program stage active for this pipeline
     *
     * @param[in] ppo		Pointer to program pipeline object
     * @param[in] flags		Shader stage flags to which this po is bound to
     * @param[in] po		Pointer to shader program object
     *
     * @return No return
     */
    virtual void setProgramStages(QSSGRenderBackendProgramPipeline ppo,
                                  QSSGRenderShaderTypeFlags flags,
                                  QSSGRenderBackendShaderProgramObject po) = 0;

    /**
     * @brief Runs a compute program
     *
     * @param[in] po				Pointer to shader program object
     * @param[in] numGroupsX		The number of work groups to be launched in the X
     * dimension
     * @param[in] numGroupsY		The number of work groups to be launched in the Y
     * dimension
     * @param[in] numGroupsZ		The number of work groups to be launched in the Z
     * dimension
     *
     * @return No return
     */
    virtual void dispatchCompute(QSSGRenderBackendShaderProgramObject po, quint32 numGroupsX, quint32 numGroupsY, quint32 numGroupsZ) = 0;

    /**
     * @brief Query constant count for a program object
     *
     * @param[in] po				Pointer to shader program object
     *
     * @return Return active constant count
     */
    virtual qint32 getConstantCount(QSSGRenderBackendShaderProgramObject po) = 0;

    /**
     * @brief Query constant buffer count for a program object
     *
     * @param[in] po				Pointer to shader program object
     *
     * @return Return active constant buffer count
     */
    virtual qint32 getConstantBufferCount(QSSGRenderBackendShaderProgramObject po) = 0;

    /**
     * @brief Query constant information by ID
     *
     * @param[in] po				Pointer to shader program object
     * @param[in] id				Constant ID
     * @param[in] bufSize			Max char for nameBuf
     * @param[out] numElem			Usually one unless for arrays
     * @param[out] type				Constant data type (QVector4D, QVector3D,...)
     * @param[out] binding			Unit binding point for samplers and images
     * @param[out] nameBuf			Name of the constant
     *
     * @return Return current constant location or -1 if not found
     */
    virtual qint32 getConstantInfoByID(QSSGRenderBackendShaderProgramObject po,
                                       quint32 id,
                                       quint32 bufSize,
                                       qint32 *numElem,
                                       QSSGRenderShaderDataType *type,
                                       qint32 *binding,
                                       char *nameBuf) = 0;

    /**
     * @brief Query constant buffer information by ID
     *
     * @param[in] po				Pointer to shader program object
     * @param[in] id				Constant buffer ID
     * @param[in] nameBufSize		Size of nameBuf
     * @param[out] paramCount		Count ot parameter contained in the buffer
     * @param[out] bufferSize		Data size of the constant buffer
     * @param[out] length			Actual characters written
     * @param[out] nameBuf			Receives the name of the buffer
     *
     * @return Return current constant buffer location or -1 if not found
     */
    virtual qint32 getConstantBufferInfoByID(QSSGRenderBackendShaderProgramObject po,
                                             quint32 id,
                                             quint32 nameBufSize,
                                             qint32 *paramCount,
                                             qint32 *bufferSize,
                                             qint32 *length,
                                             char *nameBuf) = 0;

    /**
     * @brief Query constant buffer param indices
     *
     * @param[in] po				Pointer to shader program object
     * @param[in] id				Constant buffer ID
     * @param[out] indices			Receives the indices of the uniforms within the
     * constant buffer
     *
     * @return no return value
     */
    virtual void getConstantBufferParamIndices(QSSGRenderBackendShaderProgramObject po, quint32 id, qint32 *indices) = 0;

    /**
     * @brief Query constant buffer param info by indices
     *
     * @param[in] po				Pointer to shader program object
     * @param[in] count				Number of indices
     * @param[in] indices			The indices of the uniforms within the constant
     * buffer
     * @param[out] type				Array of param types ( float ,int, ...)
     * @param[out] size				Array of param size
     * @param[out] offset			Array of param offsets within the constant buffer
     *
     * @return no return value
     */
    virtual void getConstantBufferParamInfoByIndices(QSSGRenderBackendShaderProgramObject po,
                                                     quint32 count,
                                                     quint32 *indices,
                                                     QSSGRenderShaderDataType *type,
                                                     qint32 *size,
                                                     qint32 *offset) = 0;

    /**
     * @brief Bind program constant block
     *
     * @param[in] po				Pointer to shader program object
     * @param[in] blockIndex		Constant block index returned by
     * GetConstantBufferInfoByID
     * @param[in] binding			Block binding location which should be the same as index
     * in ProgramSetConstantBlock
     *
     * @return No return
     */
    virtual void programSetConstantBlock(QSSGRenderBackendShaderProgramObject po, quint32 blockIndex, quint32 binding) = 0;

    /**
     * @brief Bind constant buffer for usage in the current active shader program
     *
     * @param[in] index				Constant ID
     * @param[in] bo				Pointer to constant buffer object
     *
     * @return No return
     */
    virtual void programSetConstantBuffer(quint32 index, QSSGRenderBackendBufferObject bo) = 0;

    /**
     * @brief Query storage buffer count for a program object
     *
     * @param[in] po				Pointer to shader program object
     *
     * @return Return active storage buffer count
     */
    virtual qint32 getStorageBufferCount(QSSGRenderBackendShaderProgramObject po) = 0;

    /**
     * @brief Query storage buffer information by ID
     *
     * @param[in] po				Pointer to shader program object
     * @param[in] id				Storage buffer ID
     * @param[in] nameBufSize		Size of nameBuf
     * @param[out] paramCount		Count of parameter contained in the buffer
     * @param[out] bufferSize		Data size of the constant buffer
     * @param[out] length			Actual characters written
     * @param[out] nameBuf			Receives the name of the buffer
     *
     * @return Return current storage buffer binding or -1 if not found
     */
    virtual qint32 getStorageBufferInfoByID(QSSGRenderBackendShaderProgramObject po,
                                            quint32 id,
                                            quint32 nameBufSize,
                                            qint32 *paramCount,
                                            qint32 *bufferSize,
                                            qint32 *length,
                                            char *nameBuf) = 0;

    /**
     * @brief Bind a storage buffer for usage in the current active shader program
     *
     * @param[in] index				Constant ID
     * @param[in] bo				Pointer to storage buffer object
     *
     * @return No return
     */
    virtual void programSetStorageBuffer(quint32 index, QSSGRenderBackendBufferObject bo) = 0;

    /**
     * @brief Set constant value
     *
     * @param[in] po				Pointer program object
     * @param[in] id				Constant ID
     * @param[in] type				Constant data type (QVector4D, QVector3D,...)
     * @param[in] count				Element count
     * @param[in] value				Pointer to constant value
     * @param[in] transpose			Transpose a matrix
     *
     * @return No return
     */
    virtual void setConstantValue(QSSGRenderBackendShaderProgramObject po,
                                  quint32 id,
                                  QSSGRenderShaderDataType type,
                                  qint32 count,
                                  const void *value,
                                  bool transpose = false) = 0;

    /**
     * @brief Draw the current active vertex buffer
     *
     * @param[in] drawMode	Draw mode (Triangles, ....)
     * @param[in] start		Start vertex
     * @param[in] count		Vertex count
     *
     * @return no return.
     */
    virtual void draw(QSSGRenderDrawMode drawMode, quint32 start, quint32 count) = 0;

    /**
     * @brief Draw the current active index buffer
     *
     * @param[in] drawMode	Draw mode (Triangles, ....)
     * @param[in] count		Index count
     * @param[in] type		Index type (quint16, quint8)
     * @param[in] indices	Pointer to index buffer. In the case of buffer objects
     *						this is an offset into the active index buffer
     *object.
     *
     * @return no return.
     */
    virtual void drawIndexed(QSSGRenderDrawMode drawMode, quint32 count, QSSGRenderComponentType type, const void *indices) = 0;

    /**
     * @brief Read a pixel rectangle from render target (from bottom left)
     *
     * @param[in]  rto		Pointer to render target object
     * @param[in]  x		Windows X start coord
     * @param[in]  y		Windows Y start coord
     * @param[in]  width	Read width dim
     * @param[in]  height	Read height dim
     * @param[out] pixels	Returned pixel data
     *
     * @return No return
     */
    virtual void readPixel(QSSGRenderBackendRenderTargetObject rto,
                           qint32 x,
                           qint32 y,
                           qint32 width,
                           qint32 height,
                           QSSGRenderReadPixelFormat inFormat,
                           QSSGByteRef pixels) = 0;

    virtual QSurfaceFormat format() const = 0;

protected:
    /// struct for what the backend supports
    typedef struct QSSGRenderBackendSupport
    {
        union {
            struct
            {
                bool bDXTImagesSupported : 1; ///< compressed images supported
                bool bAnistropySupported : 1; ///< anistropic filtering supported
                bool bTextureSwizzleSupported : 1; ///< texture swizzle supported
                bool bDepthStencilSupported : 1; ///< depth stencil textures are supported
                bool bFPRenderTargetsSupported : 1; ///< floating point render targets are
                /// supported
                bool bConstantBufferSupported : 1; ///< Constant (uniform) buffers are supported
                bool bMsTextureSupported : 1; ///< Multisample textures are esupported
                bool bFastBlitsSupported : 1; ///< The hardware supports fast memor blits
                bool bTessellationSupported : 1; ///< Hardware supports tessellation
                bool bComputeSupported : 1; ///< Hardware supports compute shader
                bool bGeometrySupported : 1; ///< Hardware supports geometry shader
                bool bTimerQuerySupported : 1; ///< Hardware supports timer queries
                bool bProgramInterfaceSupported : 1; ///< API supports program interface queries
                bool bStorageBufferSupported : 1; ///< Shader storage buffers are supported
                /// supported
                bool bShaderImageLoadStoreSupported : 1; ///< Shader image load / store
                /// operations are supported
                bool bProgramPipelineSupported : 1; ///< Driver supports separate programs
                bool bNVAdvancedBlendSupported : 1; ///< Advanced blend modes supported
                bool bNVBlendCoherenceSupported : 1; ///< Advanced blend done coherently
                /// supported
                bool bGPUShader5ExtensionSupported : 1;
                bool bKHRAdvancedBlendSupported : 1; ///< Advanced blend modes supported
                bool bKHRBlendCoherenceSupported : 1; ///< Advanced blend done coherently
                bool bVertexArrayObjectSupported : 1;
                bool bStandardDerivativesSupported : 1;
                bool bTextureLodSupported : 1;
            } bits;

            quint32 u32Values;
        } caps;
    } QSSGRenderBackendSupportBits;

    QSSGRenderBackendSupportBits m_backendSupport; ///< holds the backend support bits
};

QT_END_NAMESPACE

#endif
