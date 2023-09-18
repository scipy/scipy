/****************************************************************************
**
** Copyright (C) 2016 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DRENDER_RENDER_UNIFORM_P_H
#define QT3DRENDER_RENDER_UNIFORM_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <qt3drender_global.h>
#include <Qt3DCore/qnodeid.h>
#include <Qt3DCore/private/matrix4x4_p.h>
#include <Qt3DCore/private/vector3d_p.h>
#include <Qt3DCore/private/vector4d_p.h>
#include <Qt3DRender/private/qt3drender_global_p.h>
#include <QMatrix4x4>
#include <QVector2D>
#include <QVector3D>
#include <QColor>

#include <QDebug>
#include <string.h>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
namespace Render {

enum UniformType {
    Float = 0,
    Vec2,
    Vec3,
    Vec4,
    Double,
    DVec2,
    DVec3,
    DVec4,
    Int,
    IVec2,
    IVec3,
    IVec4,
    UInt,
    UIVec2,
    UIVec3,
    UIVec4,
    Bool,
    BVec2,
    BVec3,
    BVec4,
    Mat2,
    Mat3,
    Mat4,
    Mat2x3,
    Mat3x2,
    Mat2x4,
    Mat4x2,
    Mat3x4,
    Mat4x3,
    Sampler,
    Image,
    Unknown
};

class Q_3DRENDERSHARED_PRIVATE_EXPORT UniformValue
{
public:
    enum ValueType {
        ScalarValue,
        NodeId,
        TextureValue,
        BufferValue,
        ShaderImageValue
    };

    // UniformValue implicitely converts doubles to floats to ensure
    // correct rendering behavior for the cases where Qt3D parameters created from
    // a double or QVariant(double) are used to fill uniform values that in reality
    // should be floats. This occur especially with QML where qreal might be double
    // and not float. Otherwise, at when filling the uniforms, calling constData<float>
    // on something that contains a double will result in wrong values

    UniformValue()
        : m_data(4)
    {
        memset(m_data.data(), 0, m_data.size() * sizeof(float));
    }

    UniformValue(int i) : UniformValue() { data<int>()[0] = i; }
    UniformValue(uint i) : UniformValue() { data<uint>()[0] = i; }
    UniformValue(float f) : UniformValue() { data<float>()[0] = f; }
    UniformValue(double d) : UniformValue() { data<float>()[0] = d; } // Double to float conversion
    UniformValue(bool b) : UniformValue() { data<bool>()[0] = b; }
    UniformValue(const QVector2D &vec2) : UniformValue() { memcpy(m_data.data(), &vec2, sizeof(QVector2D)); }
    UniformValue(const Vector3D &vec3) : UniformValue() { memcpy(m_data.data(), &vec3, sizeof(Vector3D)); }
    UniformValue(const Vector4D &vec4) : m_data(sizeof(Vector4D) / sizeof(float)) { memcpy(m_data.data(), &vec4, sizeof(Vector4D)); }

    UniformValue(const QMatrix3x3 &mat33)
        : m_data(9)
    {
        // Use constData because we want column-major layout
        memcpy(m_data.data(), mat33.constData(), 9 * sizeof(float));
    }

    // We don t want the QMatrix4x4 builder to use sizeof since QMatrix4x4 contains a type flag
#if defined(__SSE2__) || defined(__AVX2__)
    UniformValue(const Matrix4x4 &mat44)
        : m_data(sizeof(Matrix4x4) / sizeof(float))
    {
        // Use constData because we want column-major layout
        memcpy(m_data.data(), &mat44, sizeof(Matrix4x4));
    }
#else
    UniformValue(const QMatrix4x4 &mat44)
        : m_data(16)
    {
        // Use constData because we want column-major layout
        memcpy(m_data.data(), mat44.constData(), 16 * sizeof(float));
    }
#endif

    UniformValue(const QVector<QMatrix4x4> &v)
        : m_data(16 * v.size())
    {
        int offset = 0;
        const int byteSize = 16 * sizeof(float);
        float *data = m_data.data();
        for (const auto &m : v) {
            memcpy(data + offset, m.constData(), byteSize);
            offset += 16;
        }
    }

    // Reserve data to be filled in later
    UniformValue(int byteSize, ValueType valueType)
        : m_data(byteSize / sizeof(float))
        , m_valueType(valueType)
    {
    }

    // For nodes which will later be replaced by a Texture or Buffer
    UniformValue(Qt3DCore::QNodeId id)
        : m_data(sizeof(Qt3DCore::QNodeId) / sizeof(float))
    {
        m_valueType = NodeId;
        memcpy(m_data.data(), &id, sizeof(Qt3DCore::QNodeId));
    }

    ValueType valueType() const { return m_valueType; }
    UniformType storedType() const { return m_storedType; }

    template<typename T>
    void setData(const QVector<T> &v)
    {
        m_data.resize(v.size() * sizeof(T) / sizeof(float));
        m_valueType = ScalarValue;
        float *data = m_data.data();
        memcpy(data, v.constData(), v.size() * sizeof(T));
    }

    static UniformValue fromVariant(const QVariant &variant);

    int byteSize() const { return m_data.size() * sizeof(float); }

    template<typename T>
    const T *constData() const
    {
        return reinterpret_cast<const T *>(m_data.constData());
    }

    template<typename T>
    T *data()
    {
        return reinterpret_cast<T *>(m_data.data());
    }

    bool operator==(const UniformValue &other) const
    {
        return other.m_data == m_data;
    }

    bool operator!=(const UniformValue &other) const
    {
        return !(*this == other);
    }
private:
    // Allocate 16 floats on stack
    // For larger elements, heap allocation will be used
    QVarLengthArray<float, 16> m_data;

    ValueType m_valueType = ScalarValue;

    // TODO: Replace this hack see QTBUG-57510
    UniformType m_storedType = Unknown;
};

template<>
Q_3DRENDERSHARED_PRIVATE_EXPORT void UniformValue::setData<QMatrix4x4>(const QVector<QMatrix4x4> &v);

} // namespace Render
} // namespace Qt3DRender

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DRender::Render::UniformType) // LCOV_EXCL_LINE

#endif // QT3DRENDER_RENDER_UNIFORM_P_H
