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

#ifndef QSSG_RENDER_CONSTANT_BUFFER_H
#define QSSG_RENDER_CONSTANT_BUFFER_H

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
#include <QtQuick3DRender/private/qssgrenderdatabuffer_p.h>
#include <QtQuick3DRender/private/qssgrenderbackend_p.h>

#include <QtCore/QString>
#include <QtCore/QHash>

QT_BEGIN_NAMESPACE

// forward declaration
class QSSGRenderContext;
class ConstantBufferParamEntry;
class QSSGRenderShaderProgram;

///< Constant (uniform) buffer representation
class Q_QUICK3DRENDER_EXPORT QSSGRenderConstantBuffer : public QSSGRenderDataBuffer
{
public:
    struct ParamHandle
    {
        QByteArray name;
        uint key = 0;
        inline static ParamHandle create(const QByteArray &name)
        {
            ParamHandle h;
            h.name = name;
            h.key = qHash(name);
            return h;
        }
    };

    enum class Param
    {
        AoProperties,
        AoProperties2,
        AoScreenConst,
        ShadowProperties,
        UvToEyeConst
    };

    template <Param>
    struct ParamData {
        static QByteArray name();
        static ParamHandle handle();
    };


    /**
     * @brief constructor
     *
     * @param[in] context		Pointer to context
     * @param[in] bufferName	Name of the buffer. Must match the name used in programs
     * @param[in] size			Size of the buffer
     * @param[in] usage			Usage of the buffer (e.g. static, dynamic...)
     * @param[in] data			A pointer to the buffer data that is allocated by the
     * application.
     *
     * @return No return.
     */
    QSSGRenderConstantBuffer(const QSSGRef<QSSGRenderContext> &context,
                               const QByteArray &bufferName,
                               QSSGRenderBufferUsageType usageType,
                               QSSGByteView data);

    ///< destructor
    virtual ~QSSGRenderConstantBuffer() override;

    /**
     * @brief bind the buffer bypasses the context state
     *
     * @return no return.
     */
    void bind() override;

    /**
     * @brief bind the buffer to a shader program
     *
     * @param[in] inShader		Pointer to active program
     * @param[in] blockIndex	Index of the constant block within the program
     * @param[in] binding		Binding point of constant buffer
     *
     * @return no return.
     */
    void bindToShaderProgram(const QSSGRef<QSSGRenderShaderProgram> &inShader, quint32 blockIndex, quint32 binding);

    /**
     * @brief update the buffer to hardware
     *
     * @return no return.
     */
    void update();

    /**
     * @brief setup constant buffer
     *
     * @param[in] pProgram		Pointer to the shader program
     * @param[in] index			Index of the constant buffer within the program
     * @param[in] bufSize		Size of the constant buffer
     * @param[in] paramCount	Parameter entry count of the constant buffer
     *
     * @return return if successful
     */
    bool setupBuffer(const QSSGRenderShaderProgram *program, qint32 index, qint32 bufSize, qint32 paramCount);

    /**
     * @brief add a parameter to the constant buffer
     *
     * @param[in] name		Name of the parameter (must match the name in the shader
     * program)
     * @param[in] type		Type of the parameter like Mat44
     * @param[in] count		One or size of array
     *
     * @return no return
     */
    void addParam(const ParamHandle &handle, QSSGRenderShaderDataType type, qint32 count);

    /**
     * @brief update a parameter in the constant buffer
     *
     * @param[in] name		Name of the parameter (must match the name in the shader
     * program)
     * @param[in] value		New value
     *
     * @return no return
     */
    void updateParam(const ParamHandle &handle, QSSGByteView value);

    /**
     * @brief update a piece of memory directly within the constant buffer
     *
     * Note: When you use this function you should know what you are doing.
     *		 The memory layout within C++ must exactly match the memory layout in the
     *shader.
     *		 We use std140 layout which guarantees a specific layout behavior across all
     *HW vendors.
     *		 How the memory layout is computed can be found in the GL spec.
     *
     * @param[in] offset	offset into constant buffer
     * @param[in] data		pointer to new data
     *
     * @return no return
     */
    void updateRaw(quint32 offset, QSSGByteView data);

    /**
     * @brief get the buffer name
     *
     * @return the buffer name
     */
    QByteArray name() const { return m_name; }

private:
    /**
     * @brief Create a parameter entry
     *
     * @param[in] name		Name of the parameter (must match the name in the shader
     * program)
     * @param[in] type		Type of the parameter like Mat44
     * @param[in] count		One or size of array
     * @param[in] offset	Offset of the parameter in the memory buffer
     *
     * @return return new Entry
     */
    ConstantBufferParamEntry *createParamEntry(const QByteArray &name, QSSGRenderShaderDataType type, qint32 count, qint32 offset);

    /**
     * @brief get size of a uniform type
     *
     * @param[in] type		type of uniform
     *
     * @return return uniform size
     */
    qint32 uniformTypeSize(QSSGRenderShaderDataType type);

    inline void setDirty(quint32 start, quint32 size)
    {
        m_rangeStart = qMin(m_rangeStart, start);
        m_rangeEnd = qMax(m_rangeEnd, start + size);
    }

    using RenderConstantBufferEntryMap = QHash<ParamHandle, ConstantBufferParamEntry *>;

    QByteArray m_name; ///< buffer name
    RenderConstantBufferEntryMap m_constantBufferEntryMap; ///< holds the entries of a constant buffer
    quint32 m_currentOffset; ///< holds the current offset
    quint32 m_currentSize; ///< holds the current size
    bool m_hwBufferInitialized; ///< true if the hardware version of the buffer is initialized
    quint32 m_rangeStart = 0; ///< start offset of the range to update
    quint32 m_rangeEnd = std::numeric_limits<quint32>::max(); ///< end of the range to update
    qint32 m_maxBlockSize; ///< maximum size for a single constant buffer
    QByteArray m_shadowCopy; ///< host copy of the data in the GPU
};

inline bool operator==(const QSSGRenderConstantBuffer::ParamHandle &h1, const QSSGRenderConstantBuffer::ParamHandle &h2)
{
    return (h1.name == h2.name);
}

template<>
struct QSSGRenderConstantBuffer::ParamData<QSSGRenderConstantBuffer::Param::AoProperties>
{
    static QByteArray name() { return QByteArrayLiteral("aoProperties"); }
    static ParamHandle handle() { return ParamHandle::create(name()); }
};

template<>
struct QSSGRenderConstantBuffer::ParamData<QSSGRenderConstantBuffer::Param::AoProperties2>
{
    static QByteArray name() { return QByteArrayLiteral("aoProperties2"); }
    static ParamHandle handle() { return ParamHandle::create(name()); }
};
template<>
struct QSSGRenderConstantBuffer::ParamData<QSSGRenderConstantBuffer::Param::AoScreenConst>
{
    static QByteArray name() { return QByteArrayLiteral("aoScreenConst"); }
    static ParamHandle handle() { return ParamHandle::create(name()); }
};
template<>
struct QSSGRenderConstantBuffer::ParamData<QSSGRenderConstantBuffer::Param::ShadowProperties>
{
    static QByteArray name() { return QByteArrayLiteral("shadowProperties"); }
    static ParamHandle handle() { return ParamHandle::create(name()); }
};
template<>
struct QSSGRenderConstantBuffer::ParamData<QSSGRenderConstantBuffer::Param::UvToEyeConst>
{
    static QByteArray name() { return QByteArrayLiteral("uvToEyeConst"); }
    static ParamHandle handle() { return ParamHandle::create(name()); }
};


QT_END_NAMESPACE

#endif
