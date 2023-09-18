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

#ifndef QSSG_RENDER_SHADER_CONSTANT_H
#define QSSG_RENDER_SHADER_CONSTANT_H

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

#include <QtQuick3DRender/private/qssgrendercontext_p.h>

#include <QtCore/QByteArray>

#include <QtGui/QMatrix4x4>
#include <QtGui/QMatrix3x3>

#include <limits>
#include <type_traits>

QT_BEGIN_NAMESPACE

///< forward declarations
class QSSGRenderContext;
class QSSGRenderConstantBuffer;

///< A shader constant belongs to a program
class Q_QUICK3DRENDER_EXPORT QSSGRenderShaderConstantBase
{
public:
    QAtomicInt ref;
    QByteArray m_name; ///< register constant name
    qint32 m_location; ///< constant index
    qint32 m_elementCount; ///< constant element count for arrays
    QSSGRenderShaderDataType m_type; ///< constant type
    qint32 m_binding; ///< sampler/imnage binding point

public:
    QSSGRenderShaderConstantBase(const QByteArray &name,
                                   qint32 location,
                                   qint32 elementCount,
                                   QSSGRenderShaderDataType type,
                                   qint32 binding)
        : m_name(name), m_location(location), m_elementCount(elementCount), m_type(type), m_binding(binding)
    {
    }

    virtual ~QSSGRenderShaderConstantBase() {}

    QSSGRenderShaderDataType getShaderConstantType() const { return m_type; }

    bool isCompatibleType(QSSGRenderShaderDataType type) const
    {
        if (m_type == type)
            return true;

        if (m_type == QSSGRenderShaderDataType::Vec4
            && type == QSSGRenderShaderDataType::Rgba) {
            return true;
        }

        return false;
    }

    virtual void release() = 0;
};

///< A general class for shader types
template<typename TDataType>
class QSSGRenderShaderConstant : public QSSGRenderShaderConstantBase
{
public:
    TDataType m_value; ///< constant value

public:
    QSSGRenderShaderConstant(const QByteArray &name,
                               qint32 location,
                               qint32 elementCount,
                               QSSGRenderShaderDataType type,
                               qint32 binding)
        : QSSGRenderShaderConstantBase(name, location, elementCount, type, binding)
    {
        if (QTypeInfo<TDataType>::isComplex)
            new (&m_value) TDataType();
        else
            memset(static_cast<void *>(&m_value), 0, sizeof(TDataType));
    }

    void release() override {}
};

///< A specialized class for textures
template<>
class QSSGRenderShaderConstant<QSSGRenderTexture2D *> : public QSSGRenderShaderConstantBase
{
public:
    quint32 m_value; ///< constant value

public:
    QSSGRenderShaderConstant(const QByteArray &name,
                               qint32 location,
                               qint32 elementCount,
                               QSSGRenderShaderDataType type,
                               qint32 binding)
        : QSSGRenderShaderConstantBase(name, location, elementCount, type, binding)
    {
        m_value = std::numeric_limits<quint32>::max();
    }

    void release() override {}
};

///< A specialized class for textures
template<>
class QSSGRenderShaderConstant<QSSGRenderTexture2D **> : public QSSGRenderShaderConstantBase
{
public:
    QVector<quint32> m_value; ///< constant value

public:
    QSSGRenderShaderConstant(const QByteArray &name,
                               qint32 location,
                               qint32 elementCount,
                               QSSGRenderShaderDataType type,
                               qint32 binding)
        : QSSGRenderShaderConstantBase(name, location, elementCount, type, binding)
    {
        m_value.resize(elementCount);
        m_value.fill(std::numeric_limits<quint32>::max());
    }

    void release() override {}
};

///< A specialized class for cubemap textures
template<>
class QSSGRenderShaderConstant<QSSGRenderTextureCube *> : public QSSGRenderShaderConstantBase
{
public:
    quint32 m_value; ///< constant value

public:
    QSSGRenderShaderConstant(const QByteArray &name,
                               qint32 location,
                               qint32 elementCount,
                               QSSGRenderShaderDataType type,
                               qint32 binding)
        : QSSGRenderShaderConstantBase(name, location, elementCount, type, binding)
    {
        m_value = std::numeric_limits<quint32>::max();
    }

    void release() override {}
};

///< A specialized class for cubemap textures
template<>
class QSSGRenderShaderConstant<QSSGRenderTextureCube **> : public QSSGRenderShaderConstantBase
{
public:
    QVector<quint32> m_value; ///< constant value

public:
    QSSGRenderShaderConstant(const QByteArray &name,
                               qint32 location,
                               qint32 elementCount,
                               QSSGRenderShaderDataType type,
                               qint32 binding)
        : QSSGRenderShaderConstantBase(name, location, elementCount, type, binding)
    {
        m_value.resize(elementCount);
        m_value.fill(std::numeric_limits<quint32>::max());
    }

    void release() override {}
};

///< A specialized class for texture image buffer
template<>
class QSSGRenderShaderConstant<QSSGRenderImage2D *> : public QSSGRenderShaderConstantBase
{
public:
    quint32 m_value; ///< constant value

public:
    QSSGRenderShaderConstant(const QByteArray &name,
                               qint32 location,
                               qint32 elementCount,
                               QSSGRenderShaderDataType type,
                               qint32 binding)
        : QSSGRenderShaderConstantBase(name, location, elementCount, type, binding)
    {
        m_value = std::numeric_limits<quint32>::max();
    }

    void release() override {}
};

///< Base for any buffer ( constant, texture, ... ) which is used by this program
class QSSGRenderShaderBufferBase
{
public:
    QAtomicInt ref;
    QSSGRef<QSSGRenderContext> m_context; ///< pointer to context
    QByteArray m_name; ///< buffer name
    quint32 m_location; ///< program buffer block location
    quint32 m_binding; ///< program buffer binding
    qint32 m_size; ///< buffer size

public:
    QSSGRenderShaderBufferBase(const QSSGRef<QSSGRenderContext> &context, const QByteArray &name, qint32 location, qint32 binding, qint32 size)
        : m_context(context), m_name(name), m_location(location), m_binding(binding), m_size(size)
    {
    }

    virtual ~QSSGRenderShaderBufferBase() {}

    virtual void release() = 0;

    virtual void validate(const QSSGRef<QSSGRenderShaderProgram> &inShader) = 0;
    virtual void bindToProgram(const QSSGRef<QSSGRenderShaderProgram> &inShader) = 0;
    virtual void update() = 0;
};

class QSSGRenderShaderConstantBuffer : public QSSGRenderShaderBufferBase
{
public:
    qint32 m_paramCount; ///< count of parameters contained in the constant buffer
    QSSGRef<QSSGRenderConstantBuffer> m_constBuffer; ///< pointer to constant buffer

public:
    QSSGRenderShaderConstantBuffer(const QSSGRef<QSSGRenderContext> &context,
                                     const QByteArray &name,
                                     quint32 location,
                                     qint32 binding,
                                     qint32 size,
                                     qint32 count,
                                     const QSSGRef<QSSGRenderConstantBuffer> &pCB)
        : QSSGRenderShaderBufferBase(context, name, location, binding, size), m_paramCount(count), m_constBuffer(pCB)
    {
    }

    void release() override {}

    void validate(const QSSGRef<QSSGRenderShaderProgram> &inShader) override;

    void update() override
    {
        if (m_constBuffer)
            m_constBuffer->update();
    }

    void bindToProgram(const QSSGRef<QSSGRenderShaderProgram> &inShader) override;
};

class QSSGRenderShaderStorageBuffer : public QSSGRenderShaderBufferBase
{
public:
    qint32 m_paramCount; ///< count of parameters contained in the constant buffer
    QSSGRef<QSSGRenderStorageBuffer> m_storageBuffer; ///< pointer to storage buffer

public:
    QSSGRenderShaderStorageBuffer(const QSSGRef<QSSGRenderContext> &context,
                                    const QByteArray &name,
                                    quint32 location,
                                    qint32 binding,
                                    qint32 size,
                                    qint32 count,
                                    const QSSGRef<QSSGRenderStorageBuffer> &pSB)
        : QSSGRenderShaderBufferBase(context, name, location, binding, size), m_paramCount(count), m_storageBuffer(pSB)
    {
    }

    void release() override {}

    void validate(const QSSGRef<QSSGRenderShaderProgram> &/*inShader*/) override;

    void update() override
    {
        if (m_storageBuffer)
            m_storageBuffer->update();
    }

    void bindToProgram(const QSSGRef<QSSGRenderShaderProgram> &/*inShader*/) override;
};

QT_END_NAMESPACE

#endif
