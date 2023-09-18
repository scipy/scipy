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

#ifndef QSSG_RENDER_EFFECT_SYSTEM_COMMANDS_H
#define QSSG_RENDER_EFFECT_SYSTEM_COMMANDS_H

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

QT_BEGIN_NAMESPACE
namespace dynamic {

enum class CommandType
{
    Unknown = 0,
    AllocateBuffer,
    BindTarget,
    BindBuffer,
    BindShader,
    ApplyInstanceValue,
    ApplyBufferValue,
    // Apply the depth buffer as an input texture.
    ApplyDepthValue,
    Render, // Render to current FBO
    ApplyBlending,
    ApplyRenderState, // apply a render state
    ApplyBlitFramebuffer,
    ApplyValue,
    DepthStencil,
    AllocateImage,
    ApplyImageValue,
    AllocateDataBuffer,
    ApplyDataBufferValue,
    ApplyCullMode,
};

// All commands need at least two constructors.  One for when they are created that should
// setup all their member variables and one for when we are copying commands from an outside
// entity into the effect system.  We have to re-register strings in that case because we
// can't assume the outside entity was using the same string table we are...
struct QSSGCommand
{
    CommandType m_type;
    QSSGCommand(CommandType inType) : m_type(inType) {}
    QSSGCommand() : m_type(CommandType::Unknown) {}
    // Implemented in UICRenderEffectSystem.cpp
//    static quint32 getSizeofCommand(const QSSGCommand &inCommand);
//    static void copyConstructCommand(quint8 *inDataBuffer, const QSSGCommand &inCommand);
};

enum class AllocateBufferFlagValues
{
    None = 0,
    SceneLifetime = 1,
};

struct QSSGAllocateBufferFlags : public QFlags<AllocateBufferFlagValues>
{
    QSSGAllocateBufferFlags(quint32 inValues) : QFlags(inValues) {}
    QSSGAllocateBufferFlags() {}
    void setSceneLifetime(bool inValue) { setFlag(AllocateBufferFlagValues::SceneLifetime, inValue); }
    // If isSceneLifetime is unset the buffer is assumed to be frame lifetime and will be
    // released after this render operation.
    bool isSceneLifetime() const { return this->operator&(AllocateBufferFlagValues::SceneLifetime); }
};

struct QSSGAllocateBuffer : public QSSGCommand
{
    QByteArray m_name;
    QSSGRenderTextureFormat m_format = QSSGRenderTextureFormat::RGBA8;
    QSSGRenderTextureMagnifyingOp m_filterOp = QSSGRenderTextureMagnifyingOp::Linear;
    QSSGRenderTextureCoordOp m_texCoordOp = QSSGRenderTextureCoordOp::ClampToEdge;
    float m_sizeMultiplier = 1.0f;
    QSSGAllocateBufferFlags m_bufferFlags;
    QSSGAllocateBuffer() : QSSGCommand(CommandType::AllocateBuffer) {}
    QSSGAllocateBuffer(const QByteArray &inName,
                         QSSGRenderTextureFormat inFormat,
                         QSSGRenderTextureMagnifyingOp inFilterOp,
                         QSSGRenderTextureCoordOp inCoordOp,
                         float inMultiplier,
                         QSSGAllocateBufferFlags inFlags)
        : QSSGCommand(CommandType::AllocateBuffer)
        , m_name(inName)
        , m_format(inFormat)
        , m_filterOp(inFilterOp)
        , m_texCoordOp(inCoordOp)
        , m_sizeMultiplier(inMultiplier)
        , m_bufferFlags(inFlags)
    {
    }
    QSSGAllocateBuffer(const QSSGAllocateBuffer &inOther)
        : QSSGCommand(CommandType::AllocateBuffer)
        , m_name(inOther.m_name)
        , m_format(inOther.m_format)
        , m_filterOp(inOther.m_filterOp)
        , m_texCoordOp(inOther.m_texCoordOp)
        , m_sizeMultiplier(inOther.m_sizeMultiplier)
        , m_bufferFlags(inOther.m_bufferFlags)
    {
    }
};

struct QSSGAllocateImage : public QSSGAllocateBuffer
{
    QSSGRenderImageAccessType m_access = QSSGRenderImageAccessType::ReadWrite;

    QSSGAllocateImage() : QSSGAllocateBuffer() { m_type = CommandType::AllocateImage; }
    QSSGAllocateImage(QByteArray &inName,
                        QSSGRenderTextureFormat inFormat,
                        QSSGRenderTextureMagnifyingOp inFilterOp,
                        QSSGRenderTextureCoordOp inCoordOp,
                        float inMultiplier,
                        QSSGAllocateBufferFlags inFlags,
                        QSSGRenderImageAccessType inAccess)
        : QSSGAllocateBuffer(inName, inFormat, inFilterOp, inCoordOp, inMultiplier, inFlags), m_access(inAccess)
    {
        m_type = CommandType::AllocateImage;
    }

    QSSGAllocateImage(const QSSGAllocateImage &inOther)
        : QSSGAllocateBuffer(inOther.m_name, inOther.m_format, inOther.m_filterOp, inOther.m_texCoordOp, inOther.m_sizeMultiplier, inOther.m_bufferFlags)
        , m_access(inOther.m_access)
    {
        m_type = CommandType::AllocateImage;
    }
};

struct QSSGAllocateDataBuffer : public QSSGCommand
{
    QByteArray m_name;
    QSSGRenderBufferType m_dataBufferType;
    QByteArray m_wrapName;
    QSSGRenderBufferType m_dataBufferWrapType;
    float m_size;
    QSSGAllocateBufferFlags m_bufferFlags;

    QSSGAllocateDataBuffer() : QSSGCommand(CommandType::AllocateDataBuffer) {}

    QSSGAllocateDataBuffer(const QByteArray &inName,
                             QSSGRenderBufferType inBufferType,
                             const QByteArray &inWrapName,
                             QSSGRenderBufferType inBufferWrapType,
                             float inSize,
                             QSSGAllocateBufferFlags inFlags)
        : QSSGCommand(CommandType::AllocateDataBuffer)
        , m_name(inName)
        , m_dataBufferType(inBufferType)
        , m_wrapName(inWrapName)
        , m_dataBufferWrapType(inBufferWrapType)
        , m_size(inSize)
        , m_bufferFlags(inFlags)
    {
    }

    QSSGAllocateDataBuffer(const QSSGAllocateDataBuffer &inOther)
        : QSSGCommand(CommandType::AllocateDataBuffer)
        , m_name(inOther.m_name)
        , m_dataBufferType(inOther.m_dataBufferType)
        , m_wrapName(inOther.m_wrapName)
        , m_dataBufferWrapType(inOther.m_dataBufferWrapType)
        , m_size(inOther.m_size)
        , m_bufferFlags(inOther.m_bufferFlags)
    {
    }
};

struct QSSGBindTarget : public QSSGCommand
{
    QSSGRenderTextureFormat m_outputFormat;

    explicit QSSGBindTarget(QSSGRenderTextureFormat inFormat = QSSGRenderTextureFormat::RGBA8)
        : QSSGCommand(CommandType::BindTarget), m_outputFormat(inFormat)
    {
    }
    QSSGBindTarget(const QSSGBindTarget &inOther)
        : QSSGCommand(CommandType::BindTarget), m_outputFormat(inOther.m_outputFormat)
    {
    }
};

struct QSSGBindBuffer : public QSSGCommand
{
    QByteArray m_bufferName;
    bool m_needsClear;
    QSSGBindBuffer(const QByteArray &inBufName, bool inNeedsClear)
        : QSSGCommand(CommandType::BindBuffer), m_bufferName(inBufName), m_needsClear(inNeedsClear)
    {
    }
    QSSGBindBuffer(const QSSGBindBuffer &inOther)
        : QSSGCommand(CommandType::BindBuffer), m_bufferName(inOther.m_bufferName), m_needsClear(inOther.m_needsClear)
    {
    }
};

struct QSSGBindShader : public QSSGCommand
{
    QByteArray m_shaderPath;
    // One GLSL file can hold multiple shaders in the case of multipass effects.
    // This makes it significantly easier for authors to reason about the shader
    // but it means we need to #define a preprocessor token to indicate which
    // effect we intend to compile at this point.
    QByteArray m_shaderDefine;
    QSSGBindShader(const QByteArray &inShaderPath, const QByteArray &inShaderDefine = QByteArray())
        : QSSGCommand(CommandType::BindShader), m_shaderPath(inShaderPath), m_shaderDefine(inShaderDefine)
    {
    }
    QSSGBindShader() : QSSGCommand(CommandType::BindShader) {}
    QSSGBindShader(const QSSGBindShader &inOther)
        : QSSGCommand(CommandType::BindShader), m_shaderPath(inOther.m_shaderPath), m_shaderDefine(inOther.m_shaderDefine)
    {
    }
};

// The value sits immediately after the 'this' object
// in memory.
// If propertyName is not valid then we attempt to apply all of the effect property values
// to the shader, ignoring ones that don't match up.
struct QSSGApplyInstanceValue : public QSSGCommand
{
    // Name of value to apply in shader
    QByteArray m_propertyName;
    // type of value
    QSSGRenderShaderDataType m_valueType;
    // offset in the effect data section of value.
    quint32 m_valueOffset;
    QSSGApplyInstanceValue(const QByteArray &inName, QSSGRenderShaderDataType inValueType, quint32 inValueOffset)
        : QSSGCommand(CommandType::ApplyInstanceValue), m_propertyName(inName), m_valueType(inValueType), m_valueOffset(inValueOffset)
    {
    }
    // Default will attempt to apply all effect values to the currently bound shader
    QSSGApplyInstanceValue()
        : QSSGCommand(CommandType::ApplyInstanceValue), m_valueType(QSSGRenderShaderDataType::Unknown), m_valueOffset(0)
    {
    }
    QSSGApplyInstanceValue(const QSSGApplyInstanceValue &inOther)
        : QSSGCommand(CommandType::ApplyInstanceValue)
        , m_propertyName(inOther.m_propertyName)
        , m_valueType(inOther.m_valueType)
        , m_valueOffset(inOther.m_valueOffset)
    {
    }
};

struct QSSGApplyValue : public QSSGCommand
{
    QByteArray m_propertyName;
    QVariant m_value;
    explicit QSSGApplyValue(const QByteArray &inName)
        : QSSGCommand(CommandType::ApplyValue), m_propertyName(inName)
    {
    }
    // Default will attempt to apply all effect values to the currently bound shader
    QSSGApplyValue() : QSSGCommand(CommandType::ApplyValue) {}

    QSSGApplyValue(const QSSGApplyValue &inOther)
        : QSSGCommand(CommandType::ApplyValue)
        , m_propertyName(inOther.m_propertyName)
        , m_value(inOther.m_value)
    {
    }
};

// bind a buffer to a given shader parameter.
struct QSSGApplyBufferValue : public QSSGCommand
{
    // If no buffer name is given then the special buffer [source]
    // is assumed.
    QByteArray m_bufferName;
    // If no param name is given, the buffer is bound to the
    // input texture parameter (texture0).
    QByteArray m_paramName;

    QSSGApplyBufferValue(const QByteArray &bufferName, const QByteArray &shaderParam)
        : QSSGCommand(CommandType::ApplyBufferValue), m_bufferName(bufferName), m_paramName(shaderParam)
    {
    }
    QSSGApplyBufferValue(const QSSGApplyBufferValue &inOther)
        : QSSGCommand(CommandType::ApplyBufferValue), m_bufferName(inOther.m_bufferName), m_paramName(inOther.m_paramName)
    {
    }
};

// bind a buffer to a given shader parameter.
struct QSSGApplyImageValue : public QSSGCommand
{
    QByteArray m_imageName; ///< name which the image was allocated
    QByteArray m_paramName; ///< must match the name in the shader
    bool m_bindAsTexture; ///< bind image as texture
    bool m_needSync; ///< if true we add a memory barrier before usage

    QSSGApplyImageValue(const QByteArray &bufferName, const QByteArray &shaderParam, bool inBindAsTexture, bool inNeedSync)
        : QSSGCommand(CommandType::ApplyImageValue)
        , m_imageName(bufferName)
        , m_paramName(shaderParam)
        , m_bindAsTexture(inBindAsTexture)
        , m_needSync(inNeedSync)
    {
    }
    QSSGApplyImageValue(const QSSGApplyImageValue &inOther)
        : QSSGCommand(CommandType::ApplyImageValue)
        , m_imageName(inOther.m_imageName)
        , m_paramName(inOther.m_paramName)
        , m_bindAsTexture(inOther.m_bindAsTexture)
        , m_needSync(inOther.m_needSync)
    {
    }
};

// bind a buffer to a given shader parameter.
struct QSSGApplyDataBufferValue : public QSSGCommand
{
    QByteArray m_paramName; ///< must match the name in the shader
    QSSGRenderBufferType m_bindAs; ///< to which target we bind this buffer

    QSSGApplyDataBufferValue(const QByteArray &inShaderParam, QSSGRenderBufferType inBufferType)
        : QSSGCommand(CommandType::ApplyDataBufferValue), m_paramName(inShaderParam), m_bindAs(inBufferType)
    {
    }
    QSSGApplyDataBufferValue(const QSSGApplyDataBufferValue &inOther)
        : QSSGCommand(CommandType::ApplyDataBufferValue), m_paramName(inOther.m_paramName), m_bindAs(inOther.m_bindAs)
    {
    }
};

struct QSSGApplyDepthValue : public QSSGCommand
{
    // If no param name is given, the buffer is bound to the
    // input texture parameter (texture0).
    QByteArray m_paramName;
    QSSGApplyDepthValue(const QByteArray &param) : QSSGCommand(CommandType::ApplyDepthValue), m_paramName(param) {}
    QSSGApplyDepthValue(const QSSGApplyDepthValue &inOther)
        : QSSGCommand(CommandType::ApplyDepthValue), m_paramName(inOther.m_paramName)
    {
    }
};

struct QSSGRender : public QSSGCommand
{
    explicit QSSGRender() : QSSGCommand(CommandType::Render) { }

    QSSGRender(const QSSGRender &)
        : QSSGCommand(CommandType::Render)
    {
    }
};

struct QSSGApplyBlending : public QSSGCommand
{
    QSSGRenderSrcBlendFunc m_srcBlendFunc;
    QSSGRenderDstBlendFunc m_dstBlendFunc;

    QSSGApplyBlending(QSSGRenderSrcBlendFunc inSrcBlendFunc, QSSGRenderDstBlendFunc inDstBlendFunc)
        : QSSGCommand(CommandType::ApplyBlending), m_srcBlendFunc(inSrcBlendFunc), m_dstBlendFunc(inDstBlendFunc)
    {
    }

    QSSGApplyBlending(const QSSGApplyBlending &inOther)
        : QSSGCommand(CommandType::ApplyBlending), m_srcBlendFunc(inOther.m_srcBlendFunc), m_dstBlendFunc(inOther.m_dstBlendFunc)
    {
    }
};

struct QSSGApplyRenderState : public QSSGCommand
{
    QSSGRenderState m_renderState;
    bool m_enabled;

    QSSGApplyRenderState(QSSGRenderState inRenderStateValue, bool inEnabled)
        : QSSGCommand(CommandType::ApplyRenderState), m_renderState(inRenderStateValue), m_enabled(inEnabled)
    {
    }

    QSSGApplyRenderState(const QSSGApplyRenderState &inOther)
        : QSSGCommand(CommandType::ApplyRenderState), m_renderState(inOther.m_renderState), m_enabled(inOther.m_enabled)
    {
    }
};

struct QSSGApplyCullMode : public QSSGCommand
{
    QSSGCullFaceMode m_cullMode;

    QSSGApplyCullMode(QSSGCullFaceMode cullMode)
        : QSSGCommand(CommandType::ApplyCullMode), m_cullMode(cullMode)
    {
    }

    QSSGApplyCullMode(const QSSGApplyCullMode &inOther)
        : QSSGCommand(CommandType::ApplyCullMode), m_cullMode(inOther.m_cullMode)
    {
    }
};

struct QSSGApplyBlitFramebuffer : public QSSGCommand
{
    // If no buffer name is given then the special buffer [source]
    // is assumed. Which is the default render target
    QByteArray m_sourceBufferName;
    // If no buffer name is given then the special buffer [dest]
    // is assumed. Which is the default render target
    QByteArray m_destBufferName;

    QSSGApplyBlitFramebuffer(const QByteArray &inSourceBufferName, const QByteArray &inDestBufferName)
        : QSSGCommand(CommandType::ApplyBlitFramebuffer), m_sourceBufferName(inSourceBufferName), m_destBufferName(inDestBufferName)
    {
    }

    QSSGApplyBlitFramebuffer(const QSSGApplyBlitFramebuffer &inOther)
        : QSSGCommand(CommandType::ApplyBlitFramebuffer)
        , m_sourceBufferName(inOther.m_sourceBufferName)
        , m_destBufferName(inOther.m_destBufferName)
    {
    }
};

enum class QSSGDepthStencilFlagValue
{
    NoFlagValue = 0,
    ClearStencil = 1 << 0,
    ClearDepth = 1 << 1,
};

struct QSSGDepthStencilFlags : public QFlags<QSSGDepthStencilFlagValue>
{
    bool hasClearStencil() const { return operator&(QSSGDepthStencilFlagValue::ClearStencil); }
    void setClearStencil(bool value) { setFlag(QSSGDepthStencilFlagValue::ClearStencil, value); }

    bool hasClearDepth() const { return operator&(QSSGDepthStencilFlagValue::ClearDepth); }
    void setClearDepth(bool value) { setFlag(QSSGDepthStencilFlagValue::ClearDepth, value); }
};

struct QSSGDepthStencil : public QSSGCommand
{
    QByteArray m_bufferName;
    QSSGDepthStencilFlags m_glags;
    QSSGRenderStencilOp m_stencilFailOperation = QSSGRenderStencilOp::Keep;
    QSSGRenderStencilOp m_depthPassOperation = QSSGRenderStencilOp::Keep;
    QSSGRenderStencilOp m_depthFailOperation = QSSGRenderStencilOp::Keep;
    QSSGRenderBoolOp m_stencilFunction = QSSGRenderBoolOp::Equal;
    quint32 m_reference = 0;
    quint32 m_mask = std::numeric_limits<quint32>::max();

    QSSGDepthStencil() : QSSGCommand(CommandType::DepthStencil) {}

    QSSGDepthStencil(const QByteArray &bufName,
                       QSSGDepthStencilFlags flags,
                       QSSGRenderStencilOp inStencilOp,
                       QSSGRenderStencilOp inDepthPassOp,
                       QSSGRenderStencilOp inDepthFailOp,
                       QSSGRenderBoolOp inStencilFunc,
                       quint32 value,
                       quint32 mask)
        : QSSGCommand(CommandType::DepthStencil)
        , m_bufferName(bufName)
        , m_glags(flags)
        , m_stencilFailOperation(inStencilOp)
        , m_depthPassOperation(inDepthPassOp)
        , m_depthFailOperation(inDepthFailOp)
        , m_stencilFunction(inStencilFunc)
        , m_reference(value)
        , m_mask(mask)
    {
    }

    QSSGDepthStencil(const QSSGDepthStencil &inOther)
        : QSSGCommand(CommandType::DepthStencil)
        , m_bufferName(inOther.m_bufferName)
        , m_glags(inOther.m_glags)
        , m_stencilFailOperation(inOther.m_stencilFailOperation)
        , m_depthPassOperation(inOther.m_depthPassOperation)
        , m_depthFailOperation(inOther.m_depthFailOperation)
        , m_stencilFunction(inOther.m_stencilFunction)
        , m_reference(inOther.m_reference)
        , m_mask(inOther.m_mask)
    {
    }
    QSSGDepthStencil& operator=(const QSSGDepthStencil&) = default;
};
}
QT_END_NAMESPACE

#endif
