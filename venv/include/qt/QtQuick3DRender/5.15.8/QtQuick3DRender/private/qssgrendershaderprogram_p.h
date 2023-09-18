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

#ifndef QSSG_RENDER_SHADER_PROGRAM_H
#define QSSG_RENDER_SHADER_PROGRAM_H

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

#include <QtQuick3DRender/private/qssgrendershaderconstant_p.h>

#include <QtCore/QString>
#include <QtCore/QHash>

QT_BEGIN_NAMESPACE

///< forward declarations
class QSSGRenderContext;
class QSSGRenderVertexShader;
class QSSGRenderFragmentShader;
class QSSGRenderTessControlShader;
class QSSGRenderTessEvaluationShader;
class QSSGRenderGeometryShader;
class QSSGRenderShaderConstantBase;
class QSSGRenderShaderBufferBase;
class QSSGRenderComputeShader;
class QColor;

typedef QHash<QByteArray, QSSGRef<QSSGRenderShaderConstantBase>> TShaderConstantMap;
typedef QHash<QByteArray, QSSGRef<QSSGRenderShaderBufferBase>> TShaderBufferMap;

///< A shader program is an object composed of a multiple shaders (vertex, fragment,
/// geometry,....)
class Q_QUICK3DRENDER_EXPORT QSSGRenderShaderProgram
{
    Q_DISABLE_COPY(QSSGRenderShaderProgram)
public:
    QAtomicInt ref;

    enum class ProgramType
    {
        Graphics,
        Compute
    };

private:
    const QSSGRef<QSSGRenderContext> m_context; ///< pointer to context
    const QSSGRef<QSSGRenderBackend> m_backend; ///< pointer to backend
    const char *m_programName; /// Name of the program
    QSSGRenderBackend::QSSGRenderBackendShaderProgramObject m_handle; ///< opaque backend handle
    TShaderConstantMap m_constants; ///< map of shader constants
    TShaderBufferMap m_shaderBuffers; ///< map of shader buffers
    ProgramType m_programType; ///< shader type
    QByteArray m_errorMessage; ///< contains the error message if linking fails

    template<typename TShaderObject>
    void attach(TShaderObject *pShader);
    template<typename TShaderObject>
    void detach(TShaderObject *pShader);

    void getShaderParameters();

    QSSGRenderShaderProgram(const QSSGRef<QSSGRenderContext> &context, const char *programName, bool separableProgram);
public:
    ~QSSGRenderShaderProgram();

    /**
     * @brief link a program
     *
     * @return true if successfully linked.
     */
    bool link();
    bool link(quint32 format, const QByteArray &binary);

    ProgramType programType() const { return m_programType; }

    /**
     * @brief Get Error Message
     *
     *
     * @return error message.
     */
    QByteArray errorMessage();

    /**
     * @brief Query constant class
     *
     * @param[in] constantName	Pointer to constant name
     *
     * @return return a pointer to a constant class.
     */
    QSSGRef<QSSGRenderShaderConstantBase> shaderConstant(const QByteArray &constantName) const;

    /**
     * @brief Query a shader buffer (constant, ... )
     *
     * @param[in] bufferName	Pointer to constant name
     *
     * @return return a pointer to a constant class.
     */
    QSSGRef<QSSGRenderShaderBufferBase> shaderBuffer(const QByteArray &bufferName) const;

    const QSSGRef<QSSGRenderBackend> &backend() const { return m_backend; }

    // concrete set functions
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, qint32 inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const qint32_2 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const qint32_3 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const qint32_4 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, bool inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const bool_2 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const bool_3 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const bool_4 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const float &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const QVector2D &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const QVector3D &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const QVector4D &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const QColor &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const quint32 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const quint32_2 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const quint32_3 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const quint32_4 &inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const QMatrix3x3 &inValue, const qint32 inCount, bool inTranspose = false);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const QMatrix4x4 &inValue, const qint32 inCount, bool inTranspose = false);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, const QSSGDataView<QMatrix4x4> inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, QSSGRenderTexture2D *inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, QSSGRenderTexture2D **inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, QSSGRenderTextureCube *inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, QSSGRenderTextureCube **inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, QSSGRenderImage2D *inValue, const qint32 inCount);
    void setConstantValue(QSSGRenderShaderConstantBase *inConstant, QSSGRenderDataBuffer *inValue, const qint32 inCount);

    /**
     * @brief Template to set constant value via name
     *
     * @param[in] inConstantName	Pointer to constant name
     * @param[in] inValue			Pointer to data
     * @param[in] inCount			Number of elements (array count)
     *
     * @return return a pointer to a constant class.
     */
    template<typename TDataType>
    void setPropertyValue(const char *inConstantName, const TDataType &inValue, const qint32 inCount = 1)
    {
        const QSSGRef<QSSGRenderShaderConstantBase> &theConstant = shaderConstant(inConstantName);

        if (theConstant) {
            if (theConstant->getShaderConstantType() == QSSGDataTypeToShaderDataTypeMap<TDataType>::getType()) {
                setConstantValue(theConstant.data(), inValue, inCount);
            } else {
                // Types don't match or property not found
                Q_ASSERT(false);
            }
        }
    }

    /**
     * @brief Template to set constant value shader constant object
     *
     * @param[in] inConstant	Pointer shader constant object
     * @param[in] inValue		Pointer to data
     * @param[in] inCount		Number of elements (array count)
     *
     * @return return a pointer to a constant class.
     */
    template<typename TDataType>
    void setPropertyValue(QSSGRenderShaderConstantBase *inConstant, const TDataType &inValue, const qint32 inCount = 1)
    {
        if (inConstant) {
            if (inConstant->isCompatibleType(QSSGDataTypeToShaderDataTypeMap<TDataType>::getType())) {
                setConstantValue(inConstant, inValue, inCount);
            } else {
                // Types don't match or property not found
                Q_ASSERT(false);
            }
        }
    }

    void bindComputeInput(QSSGRenderDataBuffer *inBuffer, quint32 inIndex);

    /**
     * @brief get the backend object handle
     *
     * @return the backend object handle.
     */
    QSSGRenderBackend::QSSGRenderBackendShaderProgramObject handle() const
    {
        return m_handle;
    }

    void getProgramBinary(quint32 &outFormat, QByteArray &outBinary) const
    {
        m_backend->getProgramBinary(m_handle, outFormat, outBinary);
    }

    /**
     * @brief get the context object
     *
     * @return context which this shader belongs to.
     */
    const QSSGRef<QSSGRenderContext> &renderContext();

    /**
     * @brief Create a shader program
     *
     * @param[in] context						Pointer to context
     * @param[in] programName					Name of the program
     * @param[in] vertShaderSource				Vertex shader source code
     * @param[in] fragShaderSource				Fragment shader source code
     * @param[in] tessControlShaderSource		tessellation control shader source code
     * @param[in] tessEvaluationShaderSource	tessellation evaluation shader source code
     * @param[in] separateProgram				True if this will we a separate
     * program
     * @param[in] type							Binary program type
     * @param[in] binaryProgram					True if program is binary
     *
     * @return a render result
     */
    static QSSGRenderVertFragCompilationResult create(
            const QSSGRef<QSSGRenderContext> &context,
            const char *programName,
            QSSGByteView vertShaderSource,
            QSSGByteView fragShaderSource,
            QSSGByteView tessControlShaderSource = QSSGByteView(),
            QSSGByteView tessEvaluationShaderSource = QSSGByteView(),
            QSSGByteView geometryShaderSource = QSSGByteView(),
            bool separateProgram = false,
            QSSGRenderShaderProgramBinaryType type = QSSGRenderShaderProgramBinaryType::Unknown,
            bool binaryProgram = false);

    /**
     * @brief Create a shader program
     *
     * @param[in] context               Pointer to context
     * @param[in] programName           Name of the program
     * @param[in] format                Format of the program
     * @param[in] binary                Program binary
     *
     * @return a render result
     */
    static QSSGRenderVertFragCompilationResult create(
            const QSSGRef<QSSGRenderContext> &context, const char *programName,
            quint32 format, const QByteArray &binary);

    /**
     * @brief Create a compute shader program
     *
     * @param[in] context						Pointer to context
     * @param[in] programName					Name of the program
     * @param[in] computeShaderSource			Compute shader source code
     *
     * @return a render result
     */
    static QSSGRenderVertFragCompilationResult createCompute(const QSSGRef<QSSGRenderContext> &context,
                                                               const char *programName,
                                                               QSSGByteView computeShaderSource);
};

// Helper class to cache the lookup of shader properties and apply them quickly in a typesafe
// way.
template<typename TDataType>
struct QSSGRenderCachedShaderProperty
{
    QSSGRef<QSSGRenderShaderProgram> shader; ///< pointer to shader program
    QSSGRef<QSSGRenderShaderConstantBase> constant; ///< poiner to shader constant object

    QSSGRenderCachedShaderProperty(const QByteArray &inConstantName, const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : shader(inShader)
    {
        const QSSGRef<QSSGRenderShaderConstantBase> &theConstant = inShader->shaderConstant(inConstantName);
        if (theConstant) {
            if (Q_LIKELY(theConstant->getShaderConstantType() == QSSGDataTypeToShaderDataTypeMap<TDataType>::getType())) {
                constant = theConstant;
            } else {
                // Property types do not match, this probably indicates that the shader changed
                // while the
                // code creating this object did not change.
                Q_ASSERT(false);
            }
        }
    }

    QSSGRenderCachedShaderProperty() = default;

    void set(const TDataType &inValue)
    {
        // TODO: Make sure the raw pointer her is ok
        if (constant)
            shader->setPropertyValue(constant.data(), inValue);
    }

    bool isValid() const { return constant != nullptr; }
};

template<typename TDataType, int size>
struct QSSGRenderCachedShaderPropertyArray
{
    QSSGRef<QSSGRenderShaderProgram> shader; ///< pointer to shader program
    QSSGRef<QSSGRenderShaderConstantBase> constant; ///< poiner to shader constant object
    TDataType m_array[size];

    QSSGRenderCachedShaderPropertyArray(const QByteArray &inConstantName, const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : shader(inShader)
    {
        memset(m_array, 0, sizeof(m_array));
        QSSGRef<QSSGRenderShaderConstantBase> theConstant = inShader->shaderConstant(inConstantName);
        if (Q_LIKELY(theConstant)) {
            if (theConstant->m_elementCount > 1 && theConstant->m_elementCount <= size
                && theConstant->getShaderConstantType() == QSSGDataTypeToShaderDataTypeMap<TDataType *>::getType()) {
                constant = theConstant;
            } else {
                // Property types do not match, this probably indicates that the shader changed
                // while the code creating this object did not change.
                Q_ASSERT(false);
            }
        }
    }

    QSSGRenderCachedShaderPropertyArray()
    {
        memset(m_array, 0, sizeof(m_array));
    }

    void set(int count)
    {
        if (constant)
            shader->setPropertyValue(constant.data(), static_cast<TDataType *>(m_array), qMin(size, count));
    }

    bool isValid() const { return constant != nullptr; }
};

// Helper class to cache the lookup of shader properties and apply them quickly in a typesafe
// way.
template<typename TDataType>
struct QSSGRenderCachedShaderBuffer
{
    QSSGRef<QSSGRenderShaderProgram> shader; ///< pointer to shader program
    QSSGRef<TDataType> shaderBuffer; ///< poiner to shader buffer object

    QSSGRenderCachedShaderBuffer(const QByteArray &inShaderBufferName, const QSSGRef<QSSGRenderShaderProgram> &inShader)
        : shader(inShader)
    {
        TDataType *theShaderBuffer = static_cast<TDataType *>(inShader->shaderBuffer(inShaderBufferName).get());
        if (theShaderBuffer)
            shaderBuffer = theShaderBuffer;
    }
    QSSGRenderCachedShaderBuffer() = default;

    void set()
    {
        if (shaderBuffer) {
            shaderBuffer->validate(shader);
            shaderBuffer->update();
            shaderBuffer->bindToProgram(shader);
        }
    }

    bool isValid() const { return shaderBuffer != 0; }
};

QT_END_NAMESPACE

#endif
