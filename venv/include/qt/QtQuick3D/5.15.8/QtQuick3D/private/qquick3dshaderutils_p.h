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

#ifndef QQUICK3DSHADERUTILS_H
#define QQUICK3DSHADERUTILS_H

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

#include <QtQuick3D/qtquick3dglobal.h>
#include <QtQuick3D/private/qquick3dobject_p.h>
#include <QtQuick3D/private/qquick3dtexture_p.h>
#include <QtQuick3D/private/qquick3dmaterial_p.h>

#include <QtQuick3DRender/private/qssgrenderbasetypes_p.h>

#include <QtQuick3DRuntimeRender/private/qssgrenderdynamicobjectsystemcommands_p.h>

QT_BEGIN_NAMESPACE


namespace QSSGShaderUtils
{
void addSnapperSampler(const QByteArray &texName, QByteArray &shaderPrefix);
QByteArray resolveShader(const QByteArray &shader, QByteArray &shaderPath, const QObject *qmlObj);
QByteArray mergeShaderCode(const QByteArray &shared,
                           const QByteArray &uniforms,
                           const QByteArray &textures,
                           const QByteArray &vertex,
                           const QByteArray &geometry,
                           const QByteArray &fragment);
}

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsTextureInput : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQuick3DTexture *texture READ texture WRITE setTexture)
    Q_PROPERTY(bool enabled MEMBER enabled)

public:
    QQuick3DShaderUtilsTextureInput() = default;
    virtual ~QQuick3DShaderUtilsTextureInput() = default;
    QQuick3DTexture *m_texture = nullptr;
    bool enabled = true;
    QByteArray name;
    QQuick3DTexture *texture() const
    {
        return m_texture;
    }

public Q_SLOTS:
    void setTexture(QQuick3DTexture *texture);

Q_SIGNALS:
    void textureDirty(QQuick3DShaderUtilsTextureInput *texture);
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsBuffer : public QObject
{
    Q_OBJECT
    Q_PROPERTY(TextureFormat format READ format WRITE setFormat)
    Q_PROPERTY(TextureFilterOperation textureFilterOperation READ textureFilterOperation WRITE setTextureFilterOperation)
    Q_PROPERTY(TextureCoordOperation textureCoordOperation READ textureCoordOperation WRITE setTextureCoordOperation)
    Q_PROPERTY(float sizeMultiplier MEMBER sizeMultiplier)
    Q_PROPERTY(AllocateBufferFlagValues bufferFlags READ bufferFlags WRITE setBufferFlags)
    Q_PROPERTY(QByteArray name MEMBER name)
public:
    QQuick3DShaderUtilsBuffer() = default;
    ~QQuick3DShaderUtilsBuffer() override = default;

    enum class TextureFilterOperation
    {
        Unknown = 0,
        Nearest,
        Linear
    };
    Q_ENUM(TextureFilterOperation)

    enum class TextureCoordOperation
    {
        Unknown = 0,
        ClampToEdge,
        MirroredRepeat,
        Repeat
    };
    Q_ENUM(TextureCoordOperation)

    enum class AllocateBufferFlagValues
    {
        None = 0,
        SceneLifetime = 1
    };
    Q_ENUM(AllocateBufferFlagValues)

    enum class TextureFormat {
        Unknown = 0,
        R8,
        R16,
        R16F,
        R32I,
        R32UI,
        R32F,
        RG8,
        RGBA8,
        RGB8,
        SRGB8,
        SRGB8A8,
        RGB565,
        RGBA16F,
        RG16F,
        RG32F,
        RGB32F,
        RGBA32F,
        R11G11B10,
        RGB9E5,
        Depth16,
        Depth24,
        Depth32,
        Depth24Stencil8
    };
    Q_ENUM(TextureFormat)

    dynamic::QSSGAllocateBuffer command {};
    TextureFilterOperation textureFilterOperation() const { return TextureFilterOperation(command.m_filterOp); }
    void setTextureFilterOperation(TextureFilterOperation op) { command.m_filterOp = QSSGRenderTextureMagnifyingOp(op); }

    TextureCoordOperation textureCoordOperation() const { return TextureCoordOperation(command.m_texCoordOp); }
    void setTextureCoordOperation(TextureCoordOperation texCoordOp) { command.m_texCoordOp = QSSGRenderTextureCoordOp(texCoordOp); }
    float &sizeMultiplier = command.m_sizeMultiplier;
    dynamic::QSSGCommand *getCommand() { return &command; }

    TextureFormat format() const;
    void setFormat(TextureFormat format);

    AllocateBufferFlagValues bufferFlags() const { return AllocateBufferFlagValues(int(command.m_bufferFlags)); }
    void setBufferFlags(AllocateBufferFlagValues flag) { command.m_bufferFlags = quint32(flag);}

    QByteArray &name = command.m_name;

    static QSSGRenderTextureFormat::Format mapTextureFormat(QQuick3DShaderUtilsBuffer::TextureFormat fmt);
    static QQuick3DShaderUtilsBuffer::TextureFormat mapRenderTextureFormat(QSSGRenderTextureFormat::Format fmt);
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsRenderCommand : public QObject
{
    Q_OBJECT
public:
    QQuick3DShaderUtilsRenderCommand() = default;
    ~QQuick3DShaderUtilsRenderCommand() override = default;
    virtual dynamic::QSSGCommand *getCommand() { Q_ASSERT(0); return nullptr; }
    virtual int bufferCount() const { return 0; }
    virtual QQuick3DShaderUtilsBuffer *bufferAt(int idx) const { Q_UNUSED(idx) return nullptr; }
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsBufferInput : public QQuick3DShaderUtilsRenderCommand
{
    Q_OBJECT
    Q_PROPERTY(QQuick3DShaderUtilsBuffer *buffer READ buffer WRITE setBuffer)
    Q_PROPERTY(QByteArray param MEMBER param)
public:
    QQuick3DShaderUtilsBufferInput() = default;
    ~QQuick3DShaderUtilsBufferInput() override = default;
    dynamic::QSSGApplyBufferValue command { QByteArray(), QByteArray() };
    QByteArray &param = command.m_paramName;
    dynamic::QSSGCommand *getCommand() override { return &command; }

    int bufferCount() const override { return (m_buffer != nullptr) ? 1 : 0; }
    QQuick3DShaderUtilsBuffer *bufferAt(int idx) const override
    {
        Q_ASSERT(idx < 1 && idx >= 0);
        return (m_buffer && idx == 0) ? m_buffer : nullptr;
    }

    QQuick3DShaderUtilsBuffer *buffer() const { return m_buffer; }
    void setBuffer(QQuick3DShaderUtilsBuffer *buffer) {
        if (m_buffer == buffer)
            return;

        if (buffer) {
            Q_ASSERT(!buffer->name.isEmpty());
            command.m_bufferName = buffer->name;
        }
        m_buffer = buffer;
    }

    QQuick3DShaderUtilsBuffer *m_buffer = nullptr;

};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsBufferBlit : public QQuick3DShaderUtilsRenderCommand
{
    Q_OBJECT
    Q_PROPERTY(QQuick3DShaderUtilsBuffer *source READ source WRITE setSource)
    Q_PROPERTY(QQuick3DShaderUtilsBuffer *destination READ destination WRITE setDestination)
public:
    QQuick3DShaderUtilsBufferBlit() = default;
    ~QQuick3DShaderUtilsBufferBlit() override = default;
    dynamic::QSSGApplyBlitFramebuffer command { QByteArray(), QByteArray() };
    dynamic::QSSGCommand *getCommand() override { return &command; }

    int bufferCount() const override {
        if (m_source != nullptr && m_destination != nullptr)
            return 2;
        if (m_source || m_destination)
            return 1;

        return 0;
    }

    QQuick3DShaderUtilsBuffer *bufferAt(int idx) const override
    {
        Q_ASSERT(idx < 2 && idx >= 0);
        if (idx == 0)
            return m_source ? m_source : m_destination;
        if (idx == 1 && m_destination)
            return m_destination;

        return nullptr;
    }

    QQuick3DShaderUtilsBuffer *source() const { return m_source; }
    void setSource(QQuick3DShaderUtilsBuffer *src)
    {
        if (src == m_source)
            return;

        if (src) {
            Q_ASSERT(!src->name.isEmpty());
            command.m_sourceBufferName = src->name;
        }
        m_source = src;
    }

    QQuick3DShaderUtilsBuffer *destination() const { return m_destination; }
    void setDestination(QQuick3DShaderUtilsBuffer *dest)
    {
        if (dest == m_destination)
            return;

        if (dest) {
            Q_ASSERT(!dest->name.isEmpty());
            command.m_destBufferName = dest->name;
        }
        m_destination = dest;
    }

    QQuick3DShaderUtilsBuffer *m_source = nullptr;
    QQuick3DShaderUtilsBuffer *m_destination = nullptr;
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsBlending : public QQuick3DShaderUtilsRenderCommand
{
    Q_OBJECT
    Q_PROPERTY(SrcBlending srcBlending READ srcBlending WRITE setSrcBlending)
    Q_PROPERTY(DestBlending destBlending READ destBlending WRITE setDestBlending)

public:
    enum class SrcBlending
    {
        Unknown = 0,
        Zero,
        One,
        SrcColor,
        OneMinusSrcColor,
        DstColor,
        OneMinusDstColor,
        SrcAlpha,
        OneMinusSrcAlpha,
        DstAlpha,
        OneMinusDstAlpha,
        ConstantColor,
        OneMinusConstantColor,
        ConstantAlpha,
        OneMinusConstantAlpha,
        SrcAlphaSaturate
    };
    Q_ENUM(SrcBlending)

    enum class DestBlending
    {
        Unknown = 0,
        Zero,
        One,
        SrcColor,
        OneMinusSrcColor,
        DstColor,
        OneMinusDstColor,
        SrcAlpha,
        OneMinusSrcAlpha,
        DstAlpha,
        OneMinusDstAlpha,
        ConstantColor,
        OneMinusConstantColor,
        ConstantAlpha,
        OneMinusConstantAlpha
    };
    Q_ENUM(DestBlending)

    QQuick3DShaderUtilsBlending() = default;
    ~QQuick3DShaderUtilsBlending() override = default;
    dynamic::QSSGApplyBlending command { QSSGRenderSrcBlendFunc::Unknown, QSSGRenderDstBlendFunc::Unknown };
    DestBlending destBlending() const
    {
        return DestBlending(command.m_dstBlendFunc);
    }
    SrcBlending srcBlending() const
    {
        return SrcBlending(command.m_srcBlendFunc);
    }

    dynamic::QSSGCommand *getCommand() override { return &command; }

public Q_SLOTS:
    void setDestBlending(DestBlending destBlending)
    {
        command.m_dstBlendFunc = QSSGRenderDstBlendFunc(destBlending);
    }
    void setSrcBlending(SrcBlending srcBlending)
    {
        command.m_srcBlendFunc= QSSGRenderSrcBlendFunc(srcBlending);
    }
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsRenderState : public QQuick3DShaderUtilsRenderCommand
{
    Q_OBJECT
    Q_PROPERTY(RenderState renderState READ renderState WRITE setRenderState)
    Q_PROPERTY(bool enabled MEMBER enabled)

public:
    enum class RenderState
    {
        Unknown = 0,
        Blend,
        CullFace,
        DepthTest,
        StencilTest,
        ScissorTest,
        DepthWrite,
        Multisample
    };
    Q_ENUM(RenderState)

    QQuick3DShaderUtilsRenderState() = default;
    ~QQuick3DShaderUtilsRenderState() override = default;
    dynamic::QSSGApplyRenderState command { QSSGRenderState::Unknown, false };
    bool &enabled = command.m_enabled;
    RenderState renderState() const
    {
        return RenderState(command.m_renderState);
    }

    dynamic::QSSGCommand *getCommand() override { return &command; }
public Q_SLOTS:
    void setRenderState(RenderState renderState)
    {
        command.m_renderState = QSSGRenderState(renderState);
    }
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsCullMode : public QQuick3DShaderUtilsRenderCommand
{
    Q_OBJECT
    Q_PROPERTY(QQuick3DMaterial::CullMode cullMode READ cullMode WRITE setCullMode)

public:
    QQuick3DShaderUtilsCullMode() = default;
    ~QQuick3DShaderUtilsCullMode() override = default;

    dynamic::QSSGApplyCullMode command { QSSGCullFaceMode::Back };

    QQuick3DMaterial::CullMode cullMode() const
    {
        return QQuick3DMaterial::CullMode(command.m_cullMode);
    }
    dynamic::QSSGCommand *getCommand() override { return &command; }
public Q_SLOTS:
    void setCullMode(QQuick3DMaterial::CullMode cullMode)
    {
        command.m_cullMode = QSSGCullFaceMode(cullMode);
    }
};

class Q_QUICK3D_EXPORT QQuick3DShaderApplyDepthValue : public QQuick3DShaderUtilsRenderCommand
{
    Q_OBJECT
    Q_PROPERTY(QByteArray param MEMBER param)

public:
    QQuick3DShaderApplyDepthValue() = default;
    ~QQuick3DShaderApplyDepthValue() override = default;

    dynamic::QSSGApplyDepthValue command { QByteArray() };

    dynamic::QSSGCommand *getCommand() override { return &command; }
    QByteArray &param = command.m_paramName;
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsApplyValue : public QQuick3DShaderUtilsRenderCommand
{
    Q_OBJECT
    Q_PROPERTY(QByteArray target MEMBER target)
    Q_PROPERTY(QVariant value MEMBER value)

public:
    QQuick3DShaderUtilsApplyValue() = default;
    ~QQuick3DShaderUtilsApplyValue() override = default;
    dynamic::QSSGCommand *getCommand() override { return &command; }
    dynamic::QSSGApplyValue command { };
    QVariant &value = command.m_value;
    QByteArray &target = command.m_propertyName;
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsShader : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QByteArray shader MEMBER shader)
    Q_PROPERTY(Stage stage MEMBER stage)
public:
    QQuick3DShaderUtilsShader() = default;
    virtual ~QQuick3DShaderUtilsShader() = default;
    enum class Stage : quint8
    {
        Shared,
        Vertex,
        Fragment,
        Geometry,
        Compute
    };
    Q_ENUM(Stage)

    QByteArray shader;
    Stage stage = Stage::Shared;
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsRenderPass : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QQmlListProperty<QQuick3DShaderUtilsRenderCommand> commands READ commands)
    Q_PROPERTY(QQuick3DShaderUtilsBuffer *output MEMBER outputBuffer)
    Q_PROPERTY(QQmlListProperty<QQuick3DShaderUtilsShader> shaders READ shaders)
public:
    QQuick3DShaderUtilsRenderPass() = default;
    ~QQuick3DShaderUtilsRenderPass() override = default;

    static void qmlAppendCommand(QQmlListProperty<QQuick3DShaderUtilsRenderCommand> *list, QQuick3DShaderUtilsRenderCommand *command);
    static QQuick3DShaderUtilsRenderCommand *qmlCommandAt(QQmlListProperty<QQuick3DShaderUtilsRenderCommand> *list, int index);
    static int qmlCommandCount(QQmlListProperty<QQuick3DShaderUtilsRenderCommand> *list);
    static void qmlCommandClear(QQmlListProperty<QQuick3DShaderUtilsRenderCommand> *list);

    static void qmlAppendShader(QQmlListProperty<QQuick3DShaderUtilsShader> *list, QQuick3DShaderUtilsShader *shader);
    static QQuick3DShaderUtilsShader *qmlShaderAt(QQmlListProperty<QQuick3DShaderUtilsShader> *list, int index);
    static int qmlShaderCount(QQmlListProperty<QQuick3DShaderUtilsShader> *list);
    static void qmlShaderClear(QQmlListProperty<QQuick3DShaderUtilsShader> *list);

    QQmlListProperty<QQuick3DShaderUtilsRenderCommand> commands();
    QVector<QQuick3DShaderUtilsRenderCommand *> m_commands;
    QQuick3DShaderUtilsBuffer *outputBuffer = nullptr;
    QQmlListProperty<QQuick3DShaderUtilsShader> shaders();
    QQuick3DShaderUtilsShader *shader(QQuick3DShaderUtilsShader::Stage stage) const;

private:
    QVarLengthArray<QQuick3DShaderUtilsShader *, 5> m_shaders;
};

class Q_QUICK3D_EXPORT QQuick3DShaderUtilsShaderInfo : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QByteArray version MEMBER version)
    Q_PROPERTY(QByteArray type MEMBER type)
    Q_PROPERTY(qint32 shaderKey MEMBER shaderKey)
public:
    QQuick3DShaderUtilsShaderInfo() = default;
    ~QQuick3DShaderUtilsShaderInfo() override = default;
    QByteArray version;
    QByteArray type; // I.e., GLSL

    enum class MaterialShaderKeyValues
    {
        Diffuse = 1 << 0,
        Specular = 1 << 1,
        Cutout = 1 << 2,
        Refraction = 1 << 3,
        Transparent = 1 << 4,
        Displace = 1 << 5,
        Transmissive = 1 << 6,
        Glossy = Diffuse | Specular
    };
    Q_ENUM(MaterialShaderKeyValues)
    Q_DECLARE_FLAGS(MaterialShaderKeyFlags, MaterialShaderKeyValues)

    qint32 shaderKey {0};
    bool isValid() const { return !(version.isEmpty() && type.isEmpty()); }
};

QT_END_NAMESPACE

#endif // QQUICK3DSHADERUTILS_H
