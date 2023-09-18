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

#ifndef QT3DCORE_VECTOR4D_SSE_P_H
#define QT3DCORE_VECTOR4D_SSE_P_H

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

#include <Qt3DCore/private/vector3d_p.h>
#include <QtGui/qvector4d.h>

#ifdef QT_COMPILER_SUPPORTS_SSE2

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class Matrix4x4_SSE;
class Matrix4x4_AVX2;

class Vector4D_SSE
{
public:
    Q_ALWAYS_INLINE Vector4D_SSE()
        : m_xyzw(_mm_setzero_ps())
    {
    }

    explicit Q_ALWAYS_INLINE Vector4D_SSE(Qt::Initialization) {}

    explicit Q_ALWAYS_INLINE Vector4D_SSE(float x, float y, float z, float w)
        : m_xyzw(_mm_set_ps(w, z, y, x))
    {
    }

    explicit Q_ALWAYS_INLINE Vector4D_SSE(QVector4D v)
        : m_xyzw(_mm_set_ps(v.w(), v.z(), v.y(), v.x()))
    {
    }

    explicit Q_ALWAYS_INLINE Vector4D_SSE(const Vector3D_SSE &vec3, float w = 0.0f)
        : m_xyzw(vec3.m_xyzw)
    {
        setW(w);
    }

    explicit Q_ALWAYS_INLINE Vector4D_SSE(QVector3D v, float w = 0.0f)
        : m_xyzw(_mm_set_ps(w, v.z(), v.y(), v.x()))
    {
    }

    Q_ALWAYS_INLINE Vector4D_SSE &operator+=(Vector4D_SSE vector)
    {
        m_xyzw = _mm_add_ps(m_xyzw, vector.m_xyzw);
        return *this;
    }

    Q_ALWAYS_INLINE Vector4D_SSE &operator-=(Vector4D_SSE vector)
    {
        m_xyzw = _mm_sub_ps(m_xyzw, vector.m_xyzw);
        return *this;
    }

    Q_ALWAYS_INLINE Vector4D_SSE &operator*=(Vector4D_SSE vector)
    {
        m_xyzw = _mm_mul_ps(m_xyzw, vector.m_xyzw);
        return *this;
    }

    Q_ALWAYS_INLINE Vector4D_SSE &operator/=(Vector4D_SSE vector)
    {
        m_xyzw = _mm_div_ps(m_xyzw, vector.m_xyzw);
        return *this;
    }

    Q_ALWAYS_INLINE Vector4D_SSE &operator*=(float factor)
    {
        m_xyzw = _mm_mul_ps(m_xyzw, _mm_set1_ps(factor));
        return *this;
    }

    Q_ALWAYS_INLINE Vector4D_SSE &operator/=(float factor)
    {
        m_xyzw = _mm_div_ps(m_xyzw, _mm_set1_ps(factor));
        return *this;
    }

    Q_ALWAYS_INLINE bool operator==(Vector4D_SSE other) const
    {
        // 0b1111 == 0xf
        return (_mm_movemask_ps(_mm_cmpeq_ps(m_xyzw, other.m_xyzw)) == 0xf);
    }

    Q_ALWAYS_INLINE bool operator!=(Vector4D_SSE other) const
    {
        return !(*this == other);
    }

    Q_ALWAYS_INLINE QVector4D toQVector4D() const
    {
        return QVector4D(x(), y(), z(), w());
    }

    // TODO: Uncomment when we introduce Vector3D_SSE
    //Q_ALWAYS_INLINE Vector3D_SSE toVector3D() const { return Vector3D_SSE(*this); }

    Q_ALWAYS_INLINE float lengthSquared() const
    {
        return dotProduct(*this, *this);
    }

    Q_ALWAYS_INLINE float length() const
    {
        return sqrt(dotProduct(*this, *this));
    }

    Q_ALWAYS_INLINE void normalize()
    {
        const float len = length();
        m_xyzw = _mm_div_ps(m_xyzw, _mm_set_ps1(len));
    }

    Q_ALWAYS_INLINE Vector4D_SSE normalized() const
    {
        Vector4D_SSE v = *this;
        v.normalize();
        return v;
    }

    Q_ALWAYS_INLINE bool isNull() const
    {
        // 0b1111 == 0xf
        return _mm_movemask_ps(_mm_cmpeq_ps(m_xyzw, _mm_setzero_ps())) == 0xf;
    }

    Q_ALWAYS_INLINE float x() const { return _mm_cvtss_f32(m_xyzw); }

    Q_ALWAYS_INLINE float y() const
    {
        // 0b01010101 = 0x55
        return _mm_cvtss_f32(_mm_shuffle_ps(m_xyzw, m_xyzw, 0x55));
    }

    Q_ALWAYS_INLINE float z() const
    {
        // 0b10101010 = 0xaa
        return _mm_cvtss_f32(_mm_unpackhi_ps(m_xyzw, m_xyzw));
    }

    Q_ALWAYS_INLINE float w() const
    {
        // 0b11111111 = 0xff
        return _mm_cvtss_f32(_mm_shuffle_ps(m_xyzw, m_xyzw, 0xff));
    }

    Q_ALWAYS_INLINE void setX(float x)
    {
        m_xyzw = _mm_move_ss(m_xyzw, _mm_set_ss(x));
    }

    Q_ALWAYS_INLINE void setY(float y)
    {
        // m_xyzw = a, b, c, d

        // y, y, y, y
        const __m128 yVec = _mm_set_ps1(y);

        // y, y, a, a
        // 0b00000000 == 0x0
        const __m128 yaVec = _mm_shuffle_ps(yVec, m_xyzw, 0x0);

        // a, y, c, d
        // 0b11100010 == 0xe2
        m_xyzw = _mm_shuffle_ps(yaVec, m_xyzw, 0xe2);
    }

    Q_ALWAYS_INLINE void setZ(float z)
    {
        // m_xyzw = a, b, c, d

        // z, z, z, z
        const __m128 zVec = _mm_set_ps1(z);

        // z, z, d, d
        // 0b11110000 == 0xf0
        const __m128 zdVec = _mm_shuffle_ps(zVec, m_xyzw, 0xf0);

        // a, b, z, d
        // 0b10000100 == 0x84
        m_xyzw = _mm_shuffle_ps(m_xyzw, zdVec, 0x84);
    }

    Q_ALWAYS_INLINE void setW(float w)
    {
#ifdef __SSE4_1__
        const __m128 wVec = _mm_set_ss(w);
        // insert element 0 of wVec into position 3 in vec3, don't zero anything
        m_xyzw = _mm_insert_ps(m_xyzw, wVec, 0x30);
#else
        // m_xyzw = a, b, c, d

        // w, w, w, w
        const __m128 wVec = _mm_set_ps1(w);

        // c, c, w, w
        const __m128 cwVec = _mm_shuffle_ps(m_xyzw, wVec, _MM_SHUFFLE(0, 0, 2, 2));

        // a, b, c, w
        m_xyzw = _mm_shuffle_ps(m_xyzw, cwVec, _MM_SHUFFLE(2, 0, 1, 0));
#endif
    }

    Q_ALWAYS_INLINE float operator[](int idx) const
    {
        Q_DECL_ALIGN(16) float vec[4];
        _mm_store_ps(vec, m_xyzw);
        return vec[idx];
    }

    struct DigitWrapper
    {
        explicit DigitWrapper(int idx, Vector4D_SSE *vec)
            : m_vec(vec)
            , m_idx(idx)
        {}

        operator float() const
        {
            switch (m_idx) {
            case 0:
                return m_vec->x();
            case 1:
                return m_vec->y();
            case 2:
                return m_vec->z();
            case 3:
                return m_vec->w();
            default:
                Q_UNREACHABLE();
                return 0.0f;
            }
        }
        void operator =(float value)
        {
            switch (m_idx) {
            case 0:
                m_vec->setX(value);
                break;
            case 1:
                m_vec->setY(value);
                break;
            case 2:
                m_vec->setZ(value);
                break;
            case 3:
                m_vec->setW(value);
                break;
            default:
                Q_UNREACHABLE();
            }
        }

    private:
        Vector4D_SSE *m_vec;
        const int m_idx;
    };

    Q_ALWAYS_INLINE DigitWrapper operator[](int idx)
    {
        return DigitWrapper(idx, this);
    }

    static Q_ALWAYS_INLINE float dotProduct(Vector4D_SSE a, Vector4D_SSE b)
    {
#if defined(__SSE4_1__)
        // 0b11111111 = 0xff
        return _mm_cvtss_f32(_mm_dp_ps(a.m_xyzw, b.m_xyzw, 0xff));
#elif defined(__SSE3__)
        const __m128 mult = _mm_mul_ps(a.m_xyzw, b.m_xyzw);
        // a + b, c + d, a + d, c + d
        const __m128 partialSum = _mm_hadd_ps(mult, mult);
        // c + d, ......
        // 0x00000001 =
        const __m128 partialSumShuffle = _mm_shuffle_ps(partialSum, partialSum, 0x1);
        return _mm_cvtss_f32(_mm_hadd_ps(partialSum, partialSumShuffle));
#else
        const __m128 mult = _mm_mul_ps(a.m_xyzw, b.m_xyzw);
        // (multX, multY, 0, 0) + (multZ, multW, 0, 0) -> (multX + multZ, multY + multW, 0, 0)
        // 0b00001110 == 0xe
        const __m128 shuffled = _mm_shuffle_ps(mult, mult, 0xe);
        __m128 result = _mm_add_ps(shuffled, mult);
        // (multX + multZ, 0, 0, 0) + (multY + multW, 0, 0, 0);
        // 0b00000001 == 0x1
        const __m128 shuffled2 = _mm_shuffle_ps(result, result, 0x1);
        result = _mm_add_ps(result, shuffled2);
        return _mm_cvtss_f32(result);
#endif
    }

    friend class Matrix4x4_SSE;

#ifdef __AVX2__
    friend class Matrix4x4_AVX2;
    friend Vector4D_SSE operator*(const Vector4D_SSE &vector, const Matrix4x4_AVX2 &matrix);
    friend Vector4D_SSE operator*(const Matrix4x4_AVX2 &matrix, const Vector4D_SSE &vector);
#endif

    friend class Vector3D_SSE;
    friend Vector4D_SSE operator*(const Vector4D_SSE &vector, const Matrix4x4_SSE &matrix);
    friend Vector4D_SSE operator*(const Matrix4x4_SSE  &matrix, const Vector4D_SSE &vector);

    friend Q_ALWAYS_INLINE const Vector4D_SSE operator+(Vector4D_SSE v1, Vector4D_SSE v2) { return v1 += v2; }
    friend Q_ALWAYS_INLINE const Vector4D_SSE operator-(Vector4D_SSE v1, Vector4D_SSE v2) { return v1 -= v2; }
    friend Q_ALWAYS_INLINE const Vector4D_SSE operator*(float factor, Vector4D_SSE vector) { return vector *= factor; }
    friend Q_ALWAYS_INLINE const Vector4D_SSE operator*(Vector4D_SSE vector, float factor) { return vector *= factor; }
    friend Q_ALWAYS_INLINE const Vector4D_SSE operator*(Vector4D_SSE v1, Vector4D_SSE v2) { return v1 *= v2; }
    friend Q_ALWAYS_INLINE const Vector4D_SSE operator-(Vector4D_SSE vector)
    {
        Vector4D_SSE c(Qt::Uninitialized);

        c.m_xyzw = _mm_xor_ps(vector.m_xyzw, _mm_set1_ps(-0.0f));

        return c;
    }

    friend Q_ALWAYS_INLINE const Vector4D_SSE operator/(Vector4D_SSE vector, float divisor) { return vector /= divisor; }
    friend Q_ALWAYS_INLINE const Vector4D_SSE operator/(Vector4D_SSE vector, Vector4D_SSE divisor) { return vector /= divisor; }

    friend Q_3DCORE_PRIVATE_EXPORT QDebug operator<<(QDebug dbg, const Vector4D_SSE &v);
    friend Q_ALWAYS_INLINE bool qFuzzyCompare(const Vector4D_SSE& v1, const Vector4D_SSE& v2)
    {
        return ::qFuzzyCompare(v1.x(), v2.x()) &&
               ::qFuzzyCompare(v1.y(), v2.y()) &&
               ::qFuzzyCompare(v1.z(), v2.z()) &&
               ::qFuzzyCompare(v1.w(), v2.w());
    }

private:
    // Q_DECL_ALIGN(16) float m[4];// for SSE support
    __m128 m_xyzw;
};

} // Qt3DCore

Q_DECLARE_TYPEINFO(Qt3DCore::Vector4D_SSE, Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DCore::Vector4D_SSE)

#endif // QT_COMPILER_SUPPORTS_SSE2

#endif // QT3DCORE_VECTOR4D_SSE_P_H
