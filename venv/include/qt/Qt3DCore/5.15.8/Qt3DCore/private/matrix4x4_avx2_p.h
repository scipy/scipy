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

#ifndef QT3DCORE_MATRIX4X4_AVX2_P_H
#define QT3DCORE_MATRIX4X4_AVX2_P_H

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

#ifdef QT_COMPILER_SUPPORTS_AVX2

// Some GCC versions don't have _mm256_set_m128 available
// Work around that
#ifndef _mm256_set_m128
#define _mm256_set_m128(va, vb) \
        _mm256_insertf128_ps(_mm256_castps128_ps256(vb), va, 1)
#endif

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class Matrix4x4_AVX2
{
public:

    Q_ALWAYS_INLINE Matrix4x4_AVX2() { setToIdentity(); }
    explicit Q_ALWAYS_INLINE Matrix4x4_AVX2(Qt::Initialization) {}

    // Assumes data is 32 bytes aligned (and in column major order)
    explicit Q_ALWAYS_INLINE Matrix4x4_AVX2(float *data)
    {
        m_col12 = _mm256_load_ps(data);
        m_col34 = _mm256_load_ps(data + 8);
    }

    // QMatrix4x4::constData returns in column major order
    explicit Q_ALWAYS_INLINE Matrix4x4_AVX2(const QMatrix4x4 &mat)
    {
        // data may not be properly aligned, using unaligned loads
        const float *data = mat.constData();
        m_col12 = _mm256_loadu_ps(data);
        m_col34 = _mm256_loadu_ps(data + 8);
    }

    // In (row major) but we store in column major order
    explicit Q_ALWAYS_INLINE Matrix4x4_AVX2(float m11, float m12, float m13, float m14,
                                            float m21, float m22, float m23, float m24,
                                            float m31, float m32, float m33, float m34,
                                            float m41, float m42, float m43, float m44)
    {
        m_col12 = _mm256_set_ps(m42, m32, m22, m12, m41, m31, m21, m11);
        m_col34 = _mm256_set_ps(m44, m34, m24, m14, m43, m33, m23, m13);
    }

    Q_ALWAYS_INLINE void setToIdentity()
    {
        // 23 instructions
        m_col12 = _mm256_set_ps(0.0f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f);
        m_col34 = _mm256_set_ps(1.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f, 0.0f);

        // 23 instructions
        // 1, 0, 0, 0
        //        __m128 vec = _mm_set_ss(1.0f);
        //        // 0, 1, 0, 0
        //        // 0b01010001 == 0x51
        //        __m128 tmp = _mm_permute_ps(vec, 0x51);

        //        // 1, 0, 0, 0, 0, 1, 0, 0
        //        m_col12 = _mm256_set_m128(tmp, vec);

        //        // 0, 0, 1, 0
        //        // 0b01000101 == 0x45
        //        tmp = _mm_permute_ps(vec, 0x45);

        //        // 0, 0, 0, 1
        //        // 0b00010101 == 0x15
        //        vec = _mm_permute_ps(vec, 0x15);

        //        // 0, 0, 1, 0, 0, 0, 0, 1
        //        m_col34 = _mm256_set_m128(vec, tmp);

        // Using a static identity matrix and assigning it is 27 instructions
    }

    Q_ALWAYS_INLINE Matrix4x4_AVX2 operator*(const Matrix4x4_AVX2 &other) const
    {
        // Shuffling: (Latency 1)
        // (8 bits -> first two pairs used to select from the first vector, second pairs from second vector)

        //       00  01  10  11, 00  01  10  11
        // v1 = m11 m12 m13 m14 m21 m22 m23 m24
        // v2 = m11 m12 m13 m14 m21 m22 m23 m24

        // shuffled with 00 00 00 00
        // v1[0] v1[0] v2[0] v2[0] v1[0] v1[0] v2[0] v2[0]
        // -> m11 m11 m11 m11 m21 m21 m21 m21

        // Broadcasting: (Latency 1)
        // -> n11 n12 n13 n14 broadcasted
        // n11 n12 n13 n14 n11 n12 n13 n14

        // Multiplying (Latency 5):
        // m11 m11 m11 m11 m21 m21 m21 m21 *  n11 n12 n13 n14 n11 n12 n13 n14

        // -> m11 n11, m11 n12, m11 n13, m11 n14
        //    m21 n11, m21 n12, m21 n13, m21 n14

        //       00  01  10  11, 00  01  10  11
        // v1 = m11 m12 m13 m14 m21 m22 m23 m24
        // v2 = m11 m12 m13 m14 m21 m22 m23 m24
        const __m256 otherCol12 = other.m_col12;
        const __m256 otherCol34 = other.m_col34;
        const __m128 col1 = _mm256_extractf128_ps(m_col12, 0);
        const __m128 col2 = _mm256_extractf128_ps(m_col12, 1);
        const __m128 col3 = _mm256_extractf128_ps(m_col34, 0);
        const __m128 col4 = _mm256_extractf128_ps(m_col34, 1);

        //    const __m256 col12 = _mm256_load_ps(m);
        //    const __m256 col34 = _mm256_load_ps(m + 8);
        //    const __m128 otherCol1 = _mm_load_ps(other.m);
        //    const __m128 otherCol2 = _mm_load_ps(other.m + 4);
        //    const __m128 otherCol3 = _mm_load_ps(other.m + 8);
        //    const __m128 otherCol4 = _mm_load_ps(other.m + 12);

        __m256 tmp = _mm256_mul_ps(_mm256_shuffle_ps(otherCol12, otherCol12, 0x00), _mm256_broadcast_ps(&col1));

        // shuffled with 01 01 01 01
        // v1[1] v1[1] v2[1] v2[1] v1[1] v1[1] v2[1] v2[1]
        // -> m12 m12 m12 m12 m22 m22 m22 m22


        //  00  01  10  11,  00  01  10  11
        // m11 m12 m13 m14, m21 m22 m23 m24 shuffled with 01 01 01 01
        // -> m12 m12 m12 m12 m22 m22 m22 m22 x n21 n22 n23 n24 n21 n22 n23 n24

        // -> m12 n21, m12 n22, m12 n23, m12 n24,
        // -> m22 n21, m22 n22, m22 n23, m22 n24
        tmp = _mm256_add_ps(_mm256_mul_ps(_mm256_shuffle_ps(otherCol12, otherCol12, 0x55), _mm256_broadcast_ps(&col2)), tmp);

        // m11 m12 m13 m14 m11 m12 m13 m14 shuffled with 10 10 10 10
        // m13 m13 m13 m13, m23 m23 m23 m23

        // Multiplying with other.col3
        // -> m13 n31, m13 n32, m13 n33, m13 n34
        // -> m23 n31, m23 n32, m23 n33, m23 n34
        tmp = _mm256_add_ps(_mm256_mul_ps(_mm256_shuffle_ps(otherCol12, otherCol12, 0xaa), _mm256_broadcast_ps(&col3)), tmp);

        // m11 m12 m13 m14 m11 m12 m13 m14 shuffled with 11 11 11 11
        // m14 m14 m14 m14 m24 m24 m24 m24

        // -> m14 n41, m14 n42, m14 n43, m14 n44
        // -> m24 n41, m24 n42, m24 n43, m24 n44
        tmp = _mm256_add_ps(_mm256_mul_ps(_mm256_shuffle_ps(otherCol12, otherCol12, 0xff), _mm256_broadcast_ps(&col4)), tmp);

        // Which finally gives
        // c11 -> m11 n11 + m12 n21 + m13 n31 + m14 n41,
        // c12 -> m11 n12 + m12 n22 + m13 n32 + m14 n42
        // c13 -> m11 n13 + m12 n23 + m13 n33 + m14 n43
        // c14 -> m11 n14 + m12 n24 + m13 n34 + m14 n44

        // c21 -> m21 n11 + m22 n21 + m23 n31 + m24 n41,
        // c12 -> m21 n12 + m22 n22 + m23 n32 + m24 n42
        // c13 -> m21 n13 + m22 n23 + m23 n33 + m24 n43
        // c14 -> m21 n14 + m22 n24 + m23 n34 + m24 n44

        __m256 tmp2 = _mm256_mul_ps(_mm256_shuffle_ps(otherCol34, otherCol34, 0x00), _mm256_broadcast_ps(&col1));
        tmp2 = _mm256_add_ps(_mm256_mul_ps(_mm256_shuffle_ps(otherCol34, otherCol34, 0x55), _mm256_broadcast_ps(&col2)), tmp2);
        tmp2 = _mm256_add_ps(_mm256_mul_ps(_mm256_shuffle_ps(otherCol34, otherCol34, 0xaa), _mm256_broadcast_ps(&col3)), tmp2);
        tmp2 = _mm256_add_ps(_mm256_mul_ps(_mm256_shuffle_ps(otherCol34, otherCol34, 0xff), _mm256_broadcast_ps(&col4)), tmp2);

        Matrix4x4_AVX2 c(Qt::Uninitialized);
        c.m_col12 = tmp;
        c.m_col34 = tmp2;
        return c;
    }

    Q_ALWAYS_INLINE Matrix4x4_AVX2 operator-(const Matrix4x4_AVX2 &other) const
    {
        Matrix4x4_AVX2 c(Qt::Uninitialized);

        c.m_col12 = _mm256_sub_ps(m_col12, other.m_col12);
        c.m_col34 = _mm256_sub_ps(m_col34, other.m_col34);
        return c;
    }

    Q_ALWAYS_INLINE Matrix4x4_AVX2 operator+(const Matrix4x4_AVX2 &other) const
    {
        Matrix4x4_AVX2 c(Qt::Uninitialized);

        c.m_col12 = _mm256_add_ps(m_col12, other.m_col12);
        c.m_col34 = _mm256_add_ps(m_col34, other.m_col34);
        return c;
    }

    Q_ALWAYS_INLINE Matrix4x4_AVX2 &operator*=(const Matrix4x4_AVX2 &other)
    {
        *this = *this * other;
        return *this;
    }

    Q_ALWAYS_INLINE Matrix4x4_AVX2 &operator-=(const Matrix4x4_AVX2 &other)
    {
        *this = *this - other;
        return *this;
    }

    Q_ALWAYS_INLINE Matrix4x4_AVX2 &operator+=(const Matrix4x4_AVX2 &other)
    {
        *this = *this + other;
        return *this;
    }

    Q_ALWAYS_INLINE Matrix4x4_AVX2 transposed() const
    {
        Matrix4x4_AVX2 c(Qt::Uninitialized);
        const __m128 col1 = _mm256_extractf128_ps(m_col12, 0);
        const __m128 col2 = _mm256_extractf128_ps(m_col12, 1);
        const __m128 col3 = _mm256_extractf128_ps(m_col34, 0);
        const __m128 col4 = _mm256_extractf128_ps(m_col34, 1);

        // ~117 instructions
        // Matrix4x4_AVX2 c = *this;
        // _MM_TRANSPOSE4_PS(c.m_col1, c.m_col2, c.m_col3, c.m_col4);

        // ~131 instructions - AVX2
        // const __m256i indexes = _mm256_set_epi32(7, 3, 6, 2, 5, 1, 4, 0);
        // c.m_col12 = _mm256_permutevar8x32_ps(_mm256_unpacklo_ps(m_col12, m_col34), indexes);
        // c.m_col34 = _mm256_permutevar8x32_ps(_mm256_unpackhi_ps(m_col12, m_col34), indexes);

        // ~193 instructions
        // c.m_col12 = _mm256_setr_ps(m_m11, m_m21, m_m31, m_m41, m_m12, m_m22, m_m32, m_m42);
        // c.m_col34 = _mm256_setr_ps(m_m13, m_m23, m_m33, m_m43, m_m14, m_m24, m_m34, m_m44);

        // ~113 instructions
        //    union {
        //        struct
        //        {
        //            __m256 twin;
        //        };
        //        struct
        //        {
        //            __m128 col1;
        //            __m128 col2;
        //        };
        //    } u;

        //    u.twin = _mm256_shuffle_ps(m_col12, m_col34, 0b01000100);
        //    c.m_col1 = _mm_permute_ps(_mm_shuffle_ps(u.col1, u.col2, 0b10001000), 0b11011000);
        //    c.m_col2 = _mm_permute_ps(_mm_shuffle_ps(u.col1, u.col2, 0b11011101), 0b11011000);

        //    u.twin = _mm256_shuffle_ps(m_col12, m_col34, 0b11101110);
        //    c.m_col3 = _mm_permute_ps(_mm_shuffle_ps(u.col1, u.col2, 0b10001000), 0b11011000);
        //    c.m_col4 = _mm_permute_ps(_mm_shuffle_ps(u.col1, u.col2, 0b11011101), 0b11011000);

        // ~113 instructions
        // 0b11011101 == 0xdd
        // 0b10001000 == 0x88
        const __m128 tmp1 = _mm_shuffle_ps(col1, col2, 0xdd);
        const __m128 tmp2 = _mm_shuffle_ps(col1, col2, 0x88);
        const __m128 tmp3 = _mm_shuffle_ps(col3, col4, 0xdd);
        const __m128 tmp4 = _mm_shuffle_ps(col3, col4, 0x88);
        c.m_col12 = _mm256_set_m128(_mm_shuffle_ps(tmp1, tmp3, 0x88), _mm_shuffle_ps(tmp2, tmp4, 0x88));
        c.m_col34 = _mm256_set_m128(_mm_shuffle_ps(tmp1, tmp3, 0xdd), _mm_shuffle_ps(tmp2, tmp4, 0xdd));

        return c;
    }

    Q_ALWAYS_INLINE Matrix4x4_AVX2 inverted() const
    {
        // TO DO: Optimize
        const QMatrix4x4 mat = toQMatrix4x4();
        return Matrix4x4_AVX2(mat.inverted());
    }

    Q_ALWAYS_INLINE bool operator==(const Matrix4x4_AVX2 &other) const
    {
        // cmp returns (-1, -1, -1, -1, -1, -1, -1, -1) if the two m256 are equals
        // movemask takes the most significant bits (8x 1 in this case) which equals 0xff
        return (_mm256_movemask_ps(_mm256_cmp_ps(m_col12, other.m_col12, _CMP_EQ_OQ)) == 0xff &&
                _mm256_movemask_ps(_mm256_cmp_ps(m_col34, other.m_col34, _CMP_EQ_OQ)) == 0xff);

    }

    Q_ALWAYS_INLINE bool operator!=(const Matrix4x4_AVX2 &other) const
    {
        return !(*this == other);
    }

    // For some reason _mm256_cvtss_f32 doesn't seem to be defined
    Q_ALWAYS_INLINE float m11() const { return _mm_cvtss_f32(_mm256_extractf128_ps(m_col12, 0)); }
    Q_ALWAYS_INLINE float m12() const { return _mm_cvtss_f32(_mm256_extractf128_ps(m_col12, 1)); }
    Q_ALWAYS_INLINE float m13() const { return _mm_cvtss_f32(_mm256_extractf128_ps(m_col34, 0)); }
    Q_ALWAYS_INLINE float m14() const { return _mm_cvtss_f32(_mm256_extractf128_ps(m_col34, 1)); }

    Q_ALWAYS_INLINE float m21() const
    {
        // 0b01010101 = 0x55
        const __m128 v = _mm256_extractf128_ps(m_col12, 0);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0x55));
    }
    Q_ALWAYS_INLINE float m22() const
    {
        // 0b01010101 = 0x55
        const __m128 v = _mm256_extractf128_ps(m_col12, 1);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0x55));
    }
    Q_ALWAYS_INLINE float m23() const
    {
        // 0b01010101 = 0x55
        const __m128 v = _mm256_extractf128_ps(m_col34, 0);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0x55));
    }
    Q_ALWAYS_INLINE float m24() const
    {
        // 0b01010101 = 0x55
        const __m128 v = _mm256_extractf128_ps(m_col34, 1);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0x55));
    }

    Q_ALWAYS_INLINE float m31() const
    {
        // 0b10101010 = 0xaa
        const __m128 v = _mm256_extractf128_ps(m_col12, 0);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0xaa));
    }
    Q_ALWAYS_INLINE float m32() const
    {
        // 0b10101010 = 0xaa
        const __m128 v = _mm256_extractf128_ps(m_col12, 1);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0xaa));
    }
    Q_ALWAYS_INLINE float m33() const
    {
        // 0b10101010 = 0xaa
        const __m128 v = _mm256_extractf128_ps(m_col34, 0);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0xaa));
    }
    Q_ALWAYS_INLINE float m34() const
    {
        // 0b10101010 = 0xaa
        const __m128 v = _mm256_extractf128_ps(m_col34, 1);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0xaa));
    }

    Q_ALWAYS_INLINE float m41() const
    {
        // 0b11111111 = 0xff
        const __m128 v = _mm256_extractf128_ps(m_col12, 0);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0xff));
    }
    Q_ALWAYS_INLINE float m42() const
    {
        // 0b11111111 = 0xff
        const __m128 v = _mm256_extractf128_ps(m_col12, 1);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0xff));
    }
    Q_ALWAYS_INLINE float m43() const
    {
        // 0b11111111 = 0xff
        const __m128 v = _mm256_extractf128_ps(m_col34, 0);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0xff));
    }
    Q_ALWAYS_INLINE float m44() const
    {
        // 0b11111111 = 0xff
        const __m128 v = _mm256_extractf128_ps(m_col34, 1);
        return _mm_cvtss_f32(_mm_shuffle_ps(v, v, 0xff));
    }
    Q_ALWAYS_INLINE QMatrix4x4 toQMatrix4x4() const { return QMatrix4x4(m11(), m12(), m13(), m14(),
                                                                        m21(), m22(), m23(), m24(),
                                                                        m31(), m32(), m33(), m34(),
                                                                        m41(), m42(), m43(), m44()); }

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
            c.m_xyzw = _mm256_extractf128_ps(m_col12, 0);
            break;
        case 1:
            c.m_xyzw = _mm256_extractf128_ps(m_col12, 1);
            break;
        case 2:
            c.m_xyzw = _mm256_extractf128_ps(m_col34, 0);
            break;
        case 3:
            c.m_xyzw = _mm256_extractf128_ps(m_col34, 1);
            break;
        default:
            Q_UNREACHABLE();
            return Vector4D();
        }
        return c;
    }

    Q_ALWAYS_INLINE Vector3D_SSE map(const Vector3D_SSE &point) const
    {
        return *this * point;
    }

    Q_ALWAYS_INLINE Vector4D_SSE map(const Vector4D_SSE &point) const
    {
        return *this * point;
    }

    Vector3D_SSE mapVector(const Vector3D_SSE &vector) const
    {
        const Vector3D_SSE row1(m11(), m12(), m13());
        const Vector3D_SSE row2(m21(), m22(), m23());
        const Vector3D_SSE row3(m31(), m32(), m33());

        return Vector3D(Vector3D_SSE::dotProduct(row1, vector),
                        Vector3D_SSE::dotProduct(row2, vector),
                        Vector3D_SSE::dotProduct(row3, vector));
    }

    friend Vector4D operator*(const Vector4D &vector, const Matrix4x4_AVX2 &matrix);
    friend Vector4D operator*(const Matrix4x4_AVX2 &matrix, const Vector4D &vector);

    friend Vector3D operator*(const Vector3D &vector, const Matrix4x4_AVX2 &matrix);
    friend Vector3D operator*(const Matrix4x4_AVX2 &matrix, const Vector3D &vector);

    friend Q_3DCORE_PRIVATE_EXPORT QDebug operator<<(QDebug dbg, const Matrix4x4_AVX2 &m);

private:
    // column major order
    // aligned on 32 bytes boundaries for AVX, compatible with 16 bytes boundary for SSE
    //    union Q_DECL_ALIGN(32)
    //    {
    //        float m[16];
    //        struct
    //        {
    //            float m_m11, m_m21, m_m31, m_m41;
    //            float m_m12, m_m22, m_m32, m_m42;
    //            float m_m13, m_m23, m_m33, m_m43;
    //            float m_m14, m_m24, m_m34, m_m44;
    //        };
    //    };
    __m256 m_col12;
    __m256 m_col34;
};

Q_ALWAYS_INLINE Vector4D operator*(const Vector4D &vector, const Matrix4x4_AVX2 &matrix)
{
    const __m256 vecMultiplier = _mm256_broadcast_ps(&vector.m_xyzw);
    // a1 a2 a3 a4 b1 b2 b3 b4, c1 c2 c3 c4 d1 d2 d3 d4
    // a1 + a2, a3 + a4, c1 + c2, c3 + c4
    // b1 + b2, b3 + b3, d1 + d2, d3 + d4

    const __m256 partialSum = _mm256_hadd_ps(_mm256_mul_ps(matrix.m_col12, vecMultiplier),
                                             _mm256_mul_ps(matrix.m_col34, vecMultiplier));

    Vector4D v(Qt::Uninitialized);
    // a12 + a34, b12 + b34, c12 + c34, d12 + d34
    // _mm256_permute4x64_pd is AVX2
    // 0b11011000 == 0xd8
    const __m256 shuffledSum = _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(partialSum), 0xd8));
    v.m_xyzw = _mm_hadd_ps(_mm256_extractf128_ps(shuffledSum, 0), _mm256_extractf128_ps(shuffledSum, 1));
    return v;
}

Q_ALWAYS_INLINE Vector4D operator*(const Matrix4x4_AVX2 &matrix, const Vector4D &vector)
{
    const Matrix4x4_AVX2 transposed = matrix.transposed();
    return vector * transposed;
}

Q_ALWAYS_INLINE Vector3D operator*(const Vector3D &vector, const Matrix4x4_AVX2 &matrix)
{
    const __m128 vec4 = _mm_set_ps(1.0f, vector.z(), vector.y(), vector.x());
    const __m256 vecMultiplier = _mm256_broadcast_ps(&vec4);
    // a1 a2 a3 a4 b1 b2 b3 b4, c1 c2 c3 c4 d1 d2 d3 d4
    // a1 + a2, a3 + a4, c1 + c2, c3 + c4
    // b1 + b2, b3 + b3, d1 + d2, d3 + d4
    const __m256 partialSum = _mm256_hadd_ps(_mm256_mul_ps(matrix.m_col12, vecMultiplier),
                                             _mm256_mul_ps(matrix.m_col34, vecMultiplier));

    // _mm256_permute4x64_pd is AVX2
    // 0b11011000 == 0xd8
    const __m256 shuffledSum = _mm256_castpd_ps(_mm256_permute4x64_pd(_mm256_castps_pd(partialSum), 0xd8));
    // a12 + a34, b12 + b34, c12 + c34, d12 + d34
    const __m128 result = _mm_hadd_ps(_mm256_extractf128_ps(shuffledSum, 0), _mm256_extractf128_ps(shuffledSum, 1));
    // 0b11111111 = 0xff
    const __m128 divisor = _mm_shuffle_ps(result, result, 0xff);

    Vector3D v(Qt::Uninitialized);
    v.m_xyzw = _mm_div_ps(result, divisor);;
    return v;
}

Q_ALWAYS_INLINE Vector3D operator*(const Matrix4x4_AVX2 &matrix, const Vector3D &vector)
{
    const Matrix4x4_AVX2 transposed = matrix.transposed();
    return vector * transposed;
}

} // Qt3DCore

Q_DECLARE_TYPEINFO(Qt3DCore::Matrix4x4_AVX2, Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DCore::Matrix4x4_AVX2)

#endif // QT_COMPILER_SUPPORTS_AVX

#endif // QT3DCORE_MATRIX4X4_AVX2_P_H
