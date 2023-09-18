/****************************************************************************
**
** Copyright (C) 2016 Paul Lemire <paul.lemire350@gmail.com>
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

#ifndef QT3DCORE_MATRIX4X4_SSE_P_H
#define QT3DCORE_MATRIX4X4_SSE_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt3D API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <Qt3DCore/private/vector4d_p.h>
#include <Qt3DCore/private/vector3d_p.h>
#include <private/qsimd_p.h>
#include <QMatrix4x4>

#ifdef QT_COMPILER_SUPPORTS_SSE2

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class Matrix4x4_SSE
{
public:

    Q_ALWAYS_INLINE Matrix4x4_SSE() { setToIdentity(); }
    explicit Q_ALWAYS_INLINE Matrix4x4_SSE(Qt::Initialization) {}

    // QMatrix4x4::constData returns in column major order
    explicit Q_ALWAYS_INLINE Matrix4x4_SSE(const QMatrix4x4 &mat)
    {
        // data may not be properly aligned, using unaligned loads
        const float *data = mat.constData();
        m_col1 = _mm_loadu_ps(data);
        m_col2 = _mm_loadu_ps(data + 4);
        m_col3 = _mm_loadu_ps(data + 8);
        m_col4 = _mm_loadu_ps(data + 12);
    }

    // Assumes data is 16 bytes aligned (and in column major order)
    explicit Q_ALWAYS_INLINE Matrix4x4_SSE(float *data)
    {
        m_col1 = _mm_load_ps(data);
        m_col2 = _mm_load_ps(data + 4);
        m_col3 = _mm_load_ps(data + 8);
        m_col4 = _mm_load_ps(data + 12);
    }

    // In (row major) but we store in column major order
    explicit Q_ALWAYS_INLINE Matrix4x4_SSE(float m11, float m12, float m13, float m14,
                                           float m21, float m22, float m23, float m24,
                                           float m31, float m32, float m33, float m34,
                                           float m41, float m42, float m43, float m44)
    {
        m_col1 = _mm_set_ps(m41, m31, m21, m11);
        m_col2 = _mm_set_ps(m42, m32, m22, m12);
        m_col3 = _mm_set_ps(m43, m33, m23, m13);
        m_col4 = _mm_set_ps(m44, m34, m24, m14);
    }

    Q_ALWAYS_INLINE void setToIdentity()
    {
        m_col1 = _mm_set_ss(1.0f);
        m_col2 = _mm_set_ps(0.0f, 0.0f, 1.0f, 0.0f);
        m_col3 = _mm_set_ps(0.0f, 1.0f, 0.0f, 0.0f);
        m_col4 = _mm_set_ps(1.0f, 0.0f, 0.0f, 0.0f);
    }

    Q_ALWAYS_INLINE Matrix4x4_SSE operator*(const Matrix4x4_SSE &other) const
    {
        Matrix4x4_SSE c(Qt::Uninitialized);

        const __m128 c1 = m_col1;
        const __m128 c2 = m_col2;
        const __m128 c3 = m_col3;
        const __m128 c4 = m_col4;

        // c11, c21, c31, c41
        // 1) (m11 x n11), (m11 x n21), (m11 x n31), (m11 x n41)
        // 2) (m11 x n11) + (m21 x n12), (m11 x n21) + (m21 x n22), (m11 x n31) + (m21 x n32), (m11 x n41) + (m21 x n42)
        // 3) (m11 x n11) + (m21 x n21) + (m31 x n13), (m11 x n21) + (m21 x n22) + (m31 x n 23), (m11 x n31) + (m21 x n32) + (m31 x n33), (m11 x n41) + (m21 x n42) (m31 x n43)
        // 4) (m11 x n11) + (m21 x n21) + (m31 x n13) + (m41 x n14), (m11 x n21) + (m21 x n22) + (m31 x n 23) + (m41 x n24), (m11 x n31) + (m21 x n32) + (m31 x n33) + (m41 x n34), (m11 x n41) + (m21 x n42) (m31 x n43) + (m41 x n44)
        __m128 tmp = _mm_mul_ps(_mm_set1_ps(other.m11()), c1);
        tmp = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m21()), c2), tmp);
        tmp = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m31()), c3), tmp);
        c.m_col1 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m41()), c4), tmp);

        // c21, c22, c23, c24
        tmp = _mm_mul_ps(_mm_set1_ps(other.m12()), c1);
        tmp = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m22()), c2), tmp);
        tmp = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m32()), c3), tmp);
        c.m_col2 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m42()), c4), tmp);

        // c31, c32, c33, c34
        tmp = _mm_mul_ps(_mm_set1_ps(other.m13()), c1);
        tmp = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m23()), c2), tmp);
        tmp = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m33()), c3), tmp);
        c.m_col3 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m43()), c4), tmp);

        // c41, c42, c43, c44
        tmp = _mm_mul_ps(_mm_set1_ps(other.m14()), c1);
        tmp = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m24()), c2), tmp);
        tmp = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m34()), c3), tmp);
        c.m_col4 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(other.m44()), c4), tmp);

        return c;
    }

    Q_ALWAYS_INLINE Matrix4x4_SSE operator-(const Matrix4x4_SSE &other) const
    {
        Matrix4x4_SSE c(Qt::Uninitialized);

        c.m_col1 = _mm_sub_ps(m_col1, other.m_col1);
        c.m_col2 = _mm_sub_ps(m_col2, other.m_col2);
        c.m_col3 = _mm_sub_ps(m_col3, other.m_col3);
        c.m_col4 = _mm_sub_ps(m_col4, other.m_col4);

        return c;
    }

    Q_ALWAYS_INLINE Matrix4x4_SSE operator+(const Matrix4x4_SSE &other) const
    {
        Matrix4x4_SSE c(Qt::Uninitialized);

        c.m_col1 = _mm_add_ps(m_col1, other.m_col1);
        c.m_col2 = _mm_add_ps(m_col2, other.m_col2);
        c.m_col3 = _mm_add_ps(m_col3, other.m_col3);
        c.m_col4 = _mm_add_ps(m_col4, other.m_col4);

        return c;
    }

    Q_ALWAYS_INLINE Matrix4x4_SSE &operator*=(const Matrix4x4_SSE &other)
    {
        *this = *this * other;
        return *this;
    }

    Q_ALWAYS_INLINE Matrix4x4_SSE &operator-=(const Matrix4x4_SSE &other)
    {
        *this = *this - other;
        return *this;
    }

    Q_ALWAYS_INLINE Matrix4x4_SSE &operator+=(const Matrix4x4_SSE &other)
    {
        *this = *this + other;
        return *this;
    }

    Q_ALWAYS_INLINE Matrix4x4_SSE transposed() const
    {
        Matrix4x4_SSE c(Qt::Uninitialized);

        // ~113 instructions
        // 0b11011101 == 0xdd
        // 0b10001000 == 0x88
        const __m128 tmp1 = _mm_shuffle_ps(m_col1, m_col2, 0xdd);
        const __m128 tmp2 = _mm_shuffle_ps(m_col1, m_col2, 0x88);
        const __m128 tmp3 = _mm_shuffle_ps(m_col3, m_col4, 0xdd);
        const __m128 tmp4 = _mm_shuffle_ps(m_col3, m_col4, 0x88);
        c.m_col1 = _mm_shuffle_ps(tmp2, tmp4, 0x88);
        c.m_col2 = _mm_shuffle_ps(tmp1, tmp3, 0x88);
        c.m_col3 = _mm_shuffle_ps(tmp2, tmp4, 0xdd);
        c.m_col4 = _mm_shuffle_ps(tmp1, tmp3, 0xdd);

        return c;
    }

    Q_ALWAYS_INLINE Matrix4x4_SSE inverted() const
    {
        // TO DO: Optimize
        const QMatrix4x4 mat = toQMatrix4x4();
        return Matrix4x4_SSE(mat.inverted());
    }

    Q_ALWAYS_INLINE bool operator==(const Matrix4x4_SSE &other) const
    {
        // 0b1111 == 0xf
        return (_mm_movemask_ps(_mm_cmpeq_ps(m_col1, other.m_col1)) == 0xf &&
                _mm_movemask_ps(_mm_cmpeq_ps(m_col2, other.m_col2)) == 0xf &&
                _mm_movemask_ps(_mm_cmpeq_ps(m_col3, other.m_col3)) == 0xf &&
                _mm_movemask_ps(_mm_cmpeq_ps(m_col4, other.m_col4)) == 0xf);
    }

    Q_ALWAYS_INLINE bool operator!=(const Matrix4x4_SSE &other) const
    {
        return !(*this == other);
    }

    Q_ALWAYS_INLINE float m11() const { return _mm_cvtss_f32(m_col1); }
    Q_ALWAYS_INLINE float m12() const { return _mm_cvtss_f32(m_col2); }
    Q_ALWAYS_INLINE float m13() const { return _mm_cvtss_f32(m_col3); }
    Q_ALWAYS_INLINE float m14() const { return _mm_cvtss_f32(m_col4); }

    Q_ALWAYS_INLINE float m21() const
    {
        // 0b01010101 = 0x55
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col1, m_col1, 0x55));
    }
    Q_ALWAYS_INLINE float m22() const
    {
        // 0b01010101 = 0x55
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col2, m_col2, 0x55));
    }
    Q_ALWAYS_INLINE float m23() const
    {
        // 0b01010101 = 0x55
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col3, m_col3, 0x55));
    }
    Q_ALWAYS_INLINE float m24() const
    {
        // 0b01010101 = 0x55
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col4, m_col4, 0x55));
    }

    Q_ALWAYS_INLINE float m31() const
    {
        // 0b10101010 = 0xaa
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col1, m_col1, 0xaa));
    }
    Q_ALWAYS_INLINE float m32() const
    {
        // 0b10101010 = 0xaa
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col2, m_col2, 0xaa));
    }
    Q_ALWAYS_INLINE float m33() const
    {
        // 0b10101010 = 0xaa
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col3, m_col3, 0xaa));
    }
    Q_ALWAYS_INLINE float m34() const
    {
        // 0b10101010 = 0xaa
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col4, m_col4, 0xaa));
    }

    Q_ALWAYS_INLINE float m41() const
    {
        // 0b11111111 = 0xff
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col1, m_col1, 0xff));
    }
    Q_ALWAYS_INLINE float m42() const
    {
        // 0b11111111 = 0xff
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col2, m_col2, 0xff));
    }
    Q_ALWAYS_INLINE float m43() const
    {
        // 0b11111111 = 0xff
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col3, m_col3, 0xff));
    }
    Q_ALWAYS_INLINE float m44() const
    {
        // 0b11111111 = 0xff
        return _mm_cvtss_f32(_mm_shuffle_ps(m_col4, m_col4, 0xff));
    }

    Q_ALWAYS_INLINE Vector4D row(int index) const
    {
        switch (index) {
        case 0:
            return Vector4D(m11(), m12(), m13(), m14());
        case 1:
            return Vector4D(m21(), m22(), m23(), m24());
        case 2:
            return Vector4D(m31(), m32(), m33(), m34());
        case 3:
            return Vector4D(m41(), m42(), m43(), m44());
        default:
            Q_UNREACHABLE();
            return Vector4D();
        }
    }

    Q_ALWAYS_INLINE Vector4D column(int index) const
    {
        Vector4D c(Qt::Uninitialized);
        switch (index) {
        case 0:
            c.m_xyzw = m_col1;
            break;
        case 1:
            c.m_xyzw = m_col2;
            break;
        case 2:
            c.m_xyzw = m_col3;
            break;
        case 3:
            c.m_xyzw = m_col4;
            break;
        default:
            Q_UNREACHABLE();
            return Vector4D();
        }
        return c;
    }

    Q_ALWAYS_INLINE QMatrix4x4 toQMatrix4x4() const { return QMatrix4x4(m11(), m12(), m13(), m14(),
                                                                        m21(), m22(), m23(), m24(),
                                                                        m31(), m32(), m33(), m34(),
                                                                        m41(), m42(), m43(), m44()); }

    Q_ALWAYS_INLINE Vector3D_SSE map(const Vector3D_SSE &point) const
    {
        return *this * point;
    }

    Q_ALWAYS_INLINE Vector4D_SSE map(const Vector4D_SSE &point) const
    {
        return *this * point;
    }

    Q_ALWAYS_INLINE Vector3D_SSE mapVector(const Vector3D_SSE &vector) const
    {
        const Vector3D_SSE row1(m11(), m12(), m13());
        const Vector3D_SSE row2(m21(), m22(), m23());
        const Vector3D_SSE row3(m31(), m32(), m33());

        return Vector3D(Vector3D_SSE::dotProduct(row1, vector),
                        Vector3D_SSE::dotProduct(row2, vector),
                        Vector3D_SSE::dotProduct(row3, vector));
    }

    friend Q_ALWAYS_INLINE Vector4D operator*(const Vector4D &vector, const Matrix4x4_SSE &matrix);
    friend Q_ALWAYS_INLINE Vector4D operator*(const Matrix4x4_SSE &matrix, const Vector4D &vector);

    friend Q_ALWAYS_INLINE Vector3D operator*(const Vector3D &vector, const Matrix4x4_SSE &matrix);
    friend Q_ALWAYS_INLINE Vector3D operator*(const Matrix4x4_SSE &matrix, const Vector3D &vector);

    friend Q_3DCORE_PRIVATE_EXPORT QDebug operator<<(QDebug dbg, const Matrix4x4_SSE &m);

private:
    // Internally we will store the matrix as indicated below
    //  Q_DECL_ALIGN(16)  // aligned on 16 bytes boundary for SSE (column major)
    //        struct
    //        {
    //            float m_m11, m_m21, m_m31, m_m41;
    //            float m_m12, m_m22, m_m32, m_m42;
    //            float m_m13, m_m23, m_m33, m_m43;
    //            float m_m14, m_m24, m_m34, m_m44;
    //        };
    //        struct
    //        {
    //            float m[16];
    //        };
    __m128 m_col1;
    __m128 m_col2;
    __m128 m_col3;
    __m128 m_col4;
};

Q_ALWAYS_INLINE Vector4D operator*(const Vector4D &vector, const Matrix4x4_SSE &matrix)
{
    const __m128 vCol1 = _mm_mul_ps(matrix.m_col1, vector.m_xyzw);
    const __m128 vCol2 = _mm_mul_ps(matrix.m_col2, vector.m_xyzw);
    const __m128 vCol3 = _mm_mul_ps(matrix.m_col3, vector.m_xyzw);
    const __m128 vCol4 = _mm_mul_ps(matrix.m_col4, vector.m_xyzw);


    // 0b01000100 == 0x44
    // 0b11101110 == 0xee

    // vCol1.x, vCol1.y, vCol2.x, vCol2.y
    __m128 tmp1 = _mm_shuffle_ps(vCol1, vCol2, 0x44);
    // vCol1.z, vCol1.w, vCol2.z, vCol2.w
    __m128 tmp2 = _mm_shuffle_ps(vCol1, vCol2, 0xee);

    // vCol1.x + vCol1.z, vCol1.y + vCol1.w, vCol2.x + vCol2.z, vCol2.y + vCol2.w,
    const __m128 tmpSum01 = _mm_add_ps(tmp1, tmp2);

    // vCol3.x, vCol3.y, vCol4.x, vCol4.y
    tmp1 = _mm_shuffle_ps(vCol3, vCol4, 0x44);
    // vCol3.z, vCol3.w, vCol4.z, vCol4.w
    tmp2 = _mm_shuffle_ps(vCol3, vCol4, 0xee);

    // vCol3.x + vCol3.z, vCol3.y + vCol3.w, vCol4.x + vCol4.z, vCol4.y + vCol4.w,
    const __m128 tmpSum02 = _mm_add_ps(tmp1, tmp2);

    // 0b10001000 == 0x88
    // 0b11011101 == 0xdd

    // vCol1.x + vCol1.z, vCol2.x + vCol2.z, vCol3.x + vCol3.z, vCol4.x + vCol4.z,
    tmp1 = _mm_shuffle_ps(tmpSum01, tmpSum02, 0x88);
    // vCol1.y + vCol1.w, vCol2.y + vCol2.w, vCol3.y + vCol3.w, vCol4.y + vCol4.w,
    tmp2 = _mm_shuffle_ps(tmpSum01, tmpSum02, 0xdd);

    Vector4D v(Qt::Uninitialized);
    v.m_xyzw = _mm_add_ps(tmp1, tmp2);
    return v;
}

Q_ALWAYS_INLINE Vector4D operator*(const Matrix4x4_SSE &matrix, const Vector4D &vector)
{
    const Matrix4x4_SSE transposed = matrix.transposed();
    return vector * transposed;
}

Q_ALWAYS_INLINE Vector3D operator*(const Vector3D &vector, const Matrix4x4_SSE &matrix)
{
    const __m128 vec4 = _mm_set_ps(1.0f, vector.z(), vector.y(), vector.x());

    const __m128 vCol1 = _mm_mul_ps(matrix.m_col1, vec4);
    const __m128 vCol2 = _mm_mul_ps(matrix.m_col2, vec4);
    const __m128 vCol3 = _mm_mul_ps(matrix.m_col3, vec4);
    const __m128 vCol4 = _mm_mul_ps(matrix.m_col4, vec4);

    // 0b01000100 == 0x44
    // 0b11101110 == 0xee

    // vCol1.x, vCol1.y, vCol2.x, vCol2.y
    __m128 tmp1 = _mm_shuffle_ps(vCol1, vCol2, 0x44);
    // vCol1.z, vCol1.w, vCol2.z, vCol2.w
    __m128 tmp2 = _mm_shuffle_ps(vCol1, vCol2, 0xee);

    // vCol1.x + vCol1.z, vCol1.y + vCol1.w, vCol2.x + vCol2.z, vCol2.y + vCol2.w,
    const __m128 tmpSum01 = _mm_add_ps(tmp1, tmp2);

    // vCol3.x, vCol3.y, vCol4.x, vCol4.y
    tmp1 = _mm_shuffle_ps(vCol3, vCol4, 0x44);
    // vCol3.z, vCol3.w, vCol4.z, vCol4.w
    tmp2 = _mm_shuffle_ps(vCol3, vCol4, 0xee);

    // vCol3.x + vCol3.z, vCol3.y + vCol3.w, vCol4.x + vCol4.z, vCol4.y + vCol4.w,
    const __m128 tmpSum02 = _mm_add_ps(tmp1, tmp2);

    // 0b10001000 == 0x88
    // 0b11011101 == 0xdd

    // vCol1.x + vCol1.z, vCol2.x + vCol2.z, vCol3.x + vCol3.z, vCol4.x + vCol4.z,
    tmp1 = _mm_shuffle_ps(tmpSum01, tmpSum02, 0x88);
    // vCol1.y + vCol1.w, vCol2.y + vCol2.w, vCol3.y + vCol3.w, vCol4.y + vCol4.w,
    tmp2 = _mm_shuffle_ps(tmpSum01, tmpSum02, 0xdd);

    const __m128 result = _mm_add_ps(tmp1, tmp2);
    // 0b11111111 = 0xff
    const __m128 divisor = _mm_shuffle_ps(result, result, 0xff);
    Vector3D v(Qt::Uninitialized);
    v.m_xyzw = _mm_div_ps(result, divisor);
    return v;
}

Q_ALWAYS_INLINE Vector3D operator*(const Matrix4x4_SSE &matrix, const Vector3D &vector)
{
    const Matrix4x4_SSE transposed = matrix.transposed();
    return vector * transposed;
}

} // Qt3DCore


Q_DECLARE_TYPEINFO(Qt3DCore::Matrix4x4_SSE, Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DCore::Matrix4x4_SSE)

#endif // QT_COMPILER_SUPPORTS_SSE2

#endif // QT3DCORE_MATRIX4X4_SSE_P_H
