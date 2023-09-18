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

#ifndef QT3DCORE_VECTOR3D_SSE_P_H
#define QT3DCORE_VECTOR3D_SSE_P_H

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

#include <Qt3DCore/private/qt3dcore_global_p.h>
#include <QtCore/private/qsimd_p.h>
#include <QtCore/QtGlobal>
#include <QtGui/qvector3d.h>
#include <QDebug>
#include <math.h>

#ifdef QT_COMPILER_SUPPORTS_SSE2

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

class Matrix4x4_SSE;
class Matrix4x4_AVX2;
class Vector4D_SSE;

class Vector3D_SSE
{
public:

    Q_ALWAYS_INLINE Vector3D_SSE()
        : m_xyzw(_mm_setzero_ps())
    {
    }

    explicit Q_ALWAYS_INLINE Vector3D_SSE(Qt::Initialization) {}

    explicit Q_ALWAYS_INLINE Vector3D_SSE(float x, float y, float z)
        : m_xyzw(_mm_set_ps(0.0f, z, y, x))
    {
    }

    explicit Q_ALWAYS_INLINE Vector3D_SSE(QVector3D v)
        : m_xyzw(_mm_set_ps(0.0f, v.z(), v.y(), v.x()))
    {
    }

    explicit Q_3DCORE_PRIVATE_EXPORT Vector3D_SSE(const Vector4D_SSE &v);

    Q_ALWAYS_INLINE Vector3D_SSE &operator+=(Vector3D_SSE vector)
    {
        m_xyzw = _mm_add_ps(m_xyzw, vector.m_xyzw);
        return *this;
    }

    Q_ALWAYS_INLINE Vector3D_SSE &operator-=(Vector3D_SSE vector)
    {
        m_xyzw = _mm_sub_ps(m_xyzw, vector.m_xyzw);
        return *this;
    }

    Q_ALWAYS_INLINE Vector3D_SSE &operator*=(Vector3D_SSE vector)
    {
        m_xyzw = _mm_mul_ps(m_xyzw, vector.m_xyzw);
        return *this;
    }

    Q_ALWAYS_INLINE Vector3D_SSE &operator/=(Vector3D_SSE vector)
    {
        m_xyzw = _mm_div_ps(m_xyzw, vector.m_xyzw);
        return *this;
    }

    Q_ALWAYS_INLINE Vector3D_SSE &operator*=(float factor)
    {
        m_xyzw = _mm_mul_ps(m_xyzw, _mm_set1_ps(factor));
        return *this;
    }

    Q_ALWAYS_INLINE Vector3D_SSE &operator/=(float factor)
    {
        m_xyzw = _mm_div_ps(m_xyzw, _mm_set1_ps(factor));
        return *this;
    }

    Q_ALWAYS_INLINE bool operator==(Vector3D_SSE other) const
    {
        // 0b111 == 0x7
        return ((_mm_movemask_ps(_mm_cmpeq_ps(m_xyzw, other.m_xyzw)) & 0x7) == 0x7);
    }

    Q_ALWAYS_INLINE bool operator!=(Vector3D_SSE other) const
    {
        return !(*this == other);
    }

    Q_ALWAYS_INLINE QVector3D toQVector3D() const
    {
        return QVector3D(x(), y(), z());
    }

    Q_ALWAYS_INLINE float lengthSquared() const
    {
        return Qt3DCore::Vector3D_SSE::dotProduct(*this, *this);
    }

    Q_ALWAYS_INLINE float length() const
    {
        return sqrt(Qt3DCore::Vector3D_SSE::dotProduct(*this, *this));
    }

    Q_ALWAYS_INLINE float distanceToPoint(const Vector3D_SSE &point) const
    {
        return (*this - point).length();
    }

    Q_ALWAYS_INLINE void normalize()
    {
        const float len = length();
        m_xyzw = _mm_div_ps(m_xyzw, _mm_set_ps1(len));
    }

    Q_ALWAYS_INLINE Vector3D_SSE normalized() const
    {
        Vector3D_SSE v = *this;
        v.normalize();
        return v;
    }

    Q_ALWAYS_INLINE bool isNull() const
    {
        // Ignore last bit
        // 0b111 = 0x7
        return ((_mm_movemask_ps(_mm_cmpeq_ps(m_xyzw, _mm_set_ps1(0.0f))) & 0x7) == 0x7);
    }

#if defined(__AVX2__) && defined(QT_COMPILER_SUPPORTS_AVX2)
    Q_3DCORE_PRIVATE_EXPORT Vector3D_SSE unproject(const Matrix4x4_AVX2 &modelView, const Matrix4x4_AVX2 &projection, const QRect &viewport) const;
    Q_3DCORE_PRIVATE_EXPORT Vector3D_SSE project(const Matrix4x4_AVX2 &modelView, const Matrix4x4_AVX2 &projection, const QRect &viewport) const;
#else
    Q_3DCORE_PRIVATE_EXPORT Vector3D_SSE unproject(const Matrix4x4_SSE &modelView, const Matrix4x4_SSE &projection, const QRect &viewport) const;
    Q_3DCORE_PRIVATE_EXPORT Vector3D_SSE project(const Matrix4x4_SSE &modelView, const Matrix4x4_SSE &projection, const QRect &viewport) const;
#endif

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

    Q_ALWAYS_INLINE float operator[](int idx) const
    {
        switch (idx) {
        case 0:
            return x();
        case 1:
            return y();
        case 2:
            return z();
        default:
            Q_UNREACHABLE();
            return 0.0f;
        }
    }

    struct DigitWrapper
    {
        explicit DigitWrapper(int idx, Vector3D_SSE *vec)
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
            default:
                Q_UNREACHABLE();
            }
        }

    private:
        Vector3D_SSE *m_vec;
        const int m_idx;
    };

    Q_ALWAYS_INLINE DigitWrapper operator[](int idx)
    {
        return DigitWrapper(idx, this);
    }

    static Q_ALWAYS_INLINE float dotProduct(Vector3D_SSE a, Vector3D_SSE b)
    {
#if defined(__SSE4_1__)
        // 0b01111111 = 0x7f
        return _mm_cvtss_f32(_mm_dp_ps(a.m_xyzw, b.m_xyzw, 0x7f));
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

        // (multX, 0, 0, 0) + (multY, 0, 0, 0) -> (multX + multY, 0, 0, 0)
        // 0b11111101 == 0xfd
        const __m128 shuffled = _mm_shuffle_ps(mult, mult, 0xfd);
        // (multX + multY, 0, 0, 0) + (multZ, 0, 0, 0);
        // 0b11111110 == 0xfe
        const __m128 shuffled2 = _mm_shuffle_ps(mult, mult, 0xfe);
        const __m128 result = _mm_add_ps(_mm_add_ps(shuffled, mult), shuffled2);
        return _mm_cvtss_f32(result);
#endif
    }

    static Q_ALWAYS_INLINE Vector3D_SSE crossProduct(Vector3D_SSE a, Vector3D_SSE b)
    {
        // a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x
        // (a.y, a.z, a.z, a.x, a.x, a.y) (b.z, b.y, b.x, b.z, b.y, b.x)
        // (a.y, a.z, a.x) * (b.z, b.x, b.y) - (a.z, a.x, a.y) (b.y, b.z, b.x)

        // 0b11001001 == 0xc9
        const __m128 a1 = _mm_shuffle_ps(a.m_xyzw, a.m_xyzw, 0xc9);
        const __m128 b2 = _mm_shuffle_ps(b.m_xyzw, b.m_xyzw, 0xc9);
        // 0b11010010 == 0xd2
        const __m128 a2 = _mm_shuffle_ps(a.m_xyzw, a.m_xyzw, 0xd2);
        const __m128 b1 = _mm_shuffle_ps(b.m_xyzw, b.m_xyzw, 0xd2);

        Vector3D_SSE v(Qt::Uninitialized);
        v.m_xyzw = _mm_sub_ps(_mm_mul_ps(a1, b1), _mm_mul_ps(a2, b2));
        return v;
    }

    friend class Vector4D_SSE;

#if defined(__AVX2__) && defined(QT_COMPILER_SUPPORTS_AVX2)
    friend class Matrix4x4_AVX2;
    friend Vector3D_SSE operator*(const Vector3D_SSE &vector, const Matrix4x4_AVX2 &matrix);
    friend Vector3D_SSE operator*(const Matrix4x4_AVX2 &matrix, const Vector3D_SSE &vector);
#endif

    friend class Matrix4x4_SSE;
    friend Vector3D_SSE operator*(const Vector3D_SSE &vector, const Matrix4x4_SSE &matrix);
    friend Vector3D_SSE operator*(const Matrix4x4_SSE &matrix, const Vector3D_SSE &vector);

    friend Q_ALWAYS_INLINE const Vector3D_SSE operator+(Vector3D_SSE v1, Vector3D_SSE v2) { return v1 += v2; }
    friend Q_ALWAYS_INLINE const Vector3D_SSE operator-(Vector3D_SSE v1, Vector3D_SSE v2) { return v1 -= v2; }
    friend Q_ALWAYS_INLINE const Vector3D_SSE operator*(float factor, Vector3D_SSE vector) { return vector *= factor; }
    friend Q_ALWAYS_INLINE const Vector3D_SSE operator*(Vector3D_SSE vector, float factor) { return vector *= factor; }
    friend Q_ALWAYS_INLINE const Vector3D_SSE operator*(Vector3D_SSE v1, Vector3D_SSE v2) { return v1 *= v2; }
    friend Q_ALWAYS_INLINE const Vector3D_SSE operator-(Vector3D_SSE vector)
    {
        Vector3D_SSE c(Qt::Uninitialized);

        c.m_xyzw = _mm_xor_ps(vector.m_xyzw, _mm_set1_ps(-0.0f));

        return c;
    }

    friend Q_ALWAYS_INLINE const Vector3D_SSE operator/(Vector3D_SSE vector, float divisor) { return vector /= divisor; }
    friend Q_ALWAYS_INLINE const Vector3D_SSE operator/(Vector3D_SSE vector, Vector3D_SSE divisor) { return vector /= divisor; }

    friend Q_3DCORE_PRIVATE_EXPORT QDebug operator<<(QDebug dbg, const Vector3D_SSE &v);
    friend Q_ALWAYS_INLINE bool qFuzzyCompare(const Vector3D_SSE& v1, const Vector3D_SSE& v2)
    {
        return ::qFuzzyCompare(v1.x(), v2.x()) &&
               ::qFuzzyCompare(v1.y(), v2.y()) &&
               ::qFuzzyCompare(v1.z(), v2.z());
    }

private:
    // Q_DECL_ALIGN(16) float m[4];// for SSE support
    __m128 m_xyzw;
};

} // Qt3DCore

Q_DECLARE_TYPEINFO(Qt3DCore::Vector3D_SSE, Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DCore::Vector3D_SSE)

#endif // QT_COMPILER_SUPPORTS_SSE2

#endif // QT3DCORE_VECTOR3D_SSE_P_H
