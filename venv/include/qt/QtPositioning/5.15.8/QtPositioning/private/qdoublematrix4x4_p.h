/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtGui module of the Qt Toolkit.
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

#ifndef QDOUBLEMATRIX4X4_H
#define QDOUBLEMATRIX4X4_H

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

#include <QtPositioning/private/qpositioningglobal_p.h>
#include <QtPositioning/private/qdoublevector3d_p.h>
#include <QtCore/QDebug>
#include <QtCore/qmetatype.h>
#include <QtCore/QRectF>

QT_BEGIN_NAMESPACE

/*
 * This class is a copy/paste/replace of QMatrix4x4
 * No algorithm has been changed.
 * Some methods have been removed.
 */

class Q_POSITIONING_PRIVATE_EXPORT QDoubleMatrix4x4
{
public:
    inline QDoubleMatrix4x4() { setToIdentity(); }
    explicit QDoubleMatrix4x4(Qt::Initialization) : flagBits(General) {}
    explicit QDoubleMatrix4x4(const double *values);
    inline QDoubleMatrix4x4(double m11, double m12, double m13, double m14,
                      double m21, double m22, double m23, double m24,
                      double m31, double m32, double m33, double m34,
                      double m41, double m42, double m43, double m44);

    QDoubleMatrix4x4(const double *values, int cols, int rows);

    inline const double& operator()(int row, int column) const;
    inline double& operator()(int row, int column);

    inline bool isAffine() const;

    inline bool isIdentity() const;
    inline void setToIdentity();

    inline void fill(double value);

    double determinant() const;
    QDoubleMatrix4x4 inverted(bool *invertible = nullptr) const;
    QDoubleMatrix4x4 transposed() const;

    inline QDoubleMatrix4x4& operator+=(const QDoubleMatrix4x4& other);
    inline QDoubleMatrix4x4& operator-=(const QDoubleMatrix4x4& other);
    inline QDoubleMatrix4x4& operator*=(const QDoubleMatrix4x4& other);
    inline QDoubleMatrix4x4& operator*=(double factor);
    QDoubleMatrix4x4& operator/=(double divisor);
    inline bool operator==(const QDoubleMatrix4x4& other) const;
    inline bool operator!=(const QDoubleMatrix4x4& other) const;

    friend QDoubleMatrix4x4 operator+(const QDoubleMatrix4x4& m1, const QDoubleMatrix4x4& m2);
    friend QDoubleMatrix4x4 operator-(const QDoubleMatrix4x4& m1, const QDoubleMatrix4x4& m2);
    friend QDoubleMatrix4x4 operator*(const QDoubleMatrix4x4& m1, const QDoubleMatrix4x4& m2);

    friend QDoubleVector3D operator*(const QDoubleMatrix4x4& matrix, const QDoubleVector3D& vector);
    friend QDoubleVector3D operator*(const QDoubleVector3D& vector, const QDoubleMatrix4x4& matrix);

    friend QPoint operator*(const QPoint& point, const QDoubleMatrix4x4& matrix);
    friend QPointF operator*(const QPointF& point, const QDoubleMatrix4x4& matrix);
    friend QDoubleMatrix4x4 operator-(const QDoubleMatrix4x4& matrix);
    friend QPoint operator*(const QDoubleMatrix4x4& matrix, const QPoint& point);
    friend QPointF operator*(const QDoubleMatrix4x4& matrix, const QPointF& point);
    friend QDoubleMatrix4x4 operator*(double factor, const QDoubleMatrix4x4& matrix);
    friend QDoubleMatrix4x4 operator*(const QDoubleMatrix4x4& matrix, double factor);
    friend Q_POSITIONING_PRIVATE_EXPORT QDoubleMatrix4x4 operator/(const QDoubleMatrix4x4& matrix, double divisor);

    friend inline bool qFuzzyCompare(const QDoubleMatrix4x4& m1, const QDoubleMatrix4x4& m2);


    void scale(const QDoubleVector3D& vector);
    void translate(const QDoubleVector3D& vector);
    void rotate(double angle, const QDoubleVector3D& vector);

    void scale(double x, double y);
    void scale(double x, double y, double z);
    void scale(double factor);
    void translate(double x, double y);
    void translate(double x, double y, double z);
    void rotate(double angle, double x, double y, double z = 0.0f);

    void ortho(const QRect& rect);
    void ortho(const QRectF& rect);
    void ortho(double left, double right, double bottom, double top, double nearPlane, double farPlane);
    void frustum(double left, double right, double bottom, double top, double nearPlane, double farPlane);
    void perspective(double verticalAngle, double aspectRatio, double nearPlane, double farPlane);

    void lookAt(const QDoubleVector3D& eye, const QDoubleVector3D& center, const QDoubleVector3D& up);

    void viewport(const QRectF &rect);
    void viewport(double left, double bottom, double width, double height, double nearPlane = 0.0f, double farPlane = 1.0f);
    void flipCoordinates();

    void copyDataTo(double *values) const;

    QPoint map(const QPoint& point) const;
    QPointF map(const QPointF& point) const;

    QDoubleVector3D map(const QDoubleVector3D& point) const;
    QDoubleVector3D mapVector(const QDoubleVector3D& vector) const;

    QRect mapRect(const QRect& rect) const;
    QRectF mapRect(const QRectF& rect) const;

    inline double *data();
    inline const double *data() const { return *m; }
    inline const double *constData() const { return *m; }

    void optimize();

#ifndef QT_NO_DEBUG_STREAM
    friend Q_POSITIONING_PRIVATE_EXPORT QDebug operator<<(QDebug dbg, const QDoubleMatrix4x4 &m);
#endif

private:
    double m[4][4];          // Column-major order to match OpenGL.
    int flagBits;           // Flag bits from the enum below.

    // When matrices are multiplied, the flag bits are or-ed together.
    enum {
        Identity        = 0x0000, // Identity matrix
        Translation     = 0x0001, // Contains a translation
        Scale           = 0x0002, // Contains a scale
        Rotation2D      = 0x0004, // Contains a rotation about the Z axis
        Rotation        = 0x0008, // Contains an arbitrary rotation
        Perspective     = 0x0010, // Last row is different from (0, 0, 0, 1)
        General         = 0x001f  // General matrix, unknown contents
    };

    // Construct without initializing identity matrix.
    explicit QDoubleMatrix4x4(int) { }

    QDoubleMatrix4x4 orthonormalInverse() const;

    void projectedRotate(double angle, double x, double y, double z);
};

Q_DECLARE_TYPEINFO(QDoubleMatrix4x4, Q_MOVABLE_TYPE);

inline QDoubleMatrix4x4::QDoubleMatrix4x4
        (double m11, double m12, double m13, double m14,
         double m21, double m22, double m23, double m24,
         double m31, double m32, double m33, double m34,
         double m41, double m42, double m43, double m44)
{
    m[0][0] = m11; m[0][1] = m21; m[0][2] = m31; m[0][3] = m41;
    m[1][0] = m12; m[1][1] = m22; m[1][2] = m32; m[1][3] = m42;
    m[2][0] = m13; m[2][1] = m23; m[2][2] = m33; m[2][3] = m43;
    m[3][0] = m14; m[3][1] = m24; m[3][2] = m34; m[3][3] = m44;
    flagBits = General;
}

inline const double& QDoubleMatrix4x4::operator()(int aRow, int aColumn) const
{
    Q_ASSERT(aRow >= 0 && aRow < 4 && aColumn >= 0 && aColumn < 4);
    return m[aColumn][aRow];
}

inline double& QDoubleMatrix4x4::operator()(int aRow, int aColumn)
{
    Q_ASSERT(aRow >= 0 && aRow < 4 && aColumn >= 0 && aColumn < 4);
    flagBits = General;
    return m[aColumn][aRow];
}

Q_POSITIONING_PRIVATE_EXPORT QDoubleMatrix4x4 operator/(const QDoubleMatrix4x4& matrix, double divisor);

inline bool QDoubleMatrix4x4::isAffine() const
{
    return m[0][3] == 0.0f && m[1][3] == 0.0f && m[2][3] == 0.0f && m[3][3] == 1.0f;
}

inline bool QDoubleMatrix4x4::isIdentity() const
{
    if (flagBits == Identity)
        return true;
    if (m[0][0] != 1.0f || m[0][1] != 0.0f || m[0][2] != 0.0f)
        return false;
    if (m[0][3] != 0.0f || m[1][0] != 0.0f || m[1][1] != 1.0f)
        return false;
    if (m[1][2] != 0.0f || m[1][3] != 0.0f || m[2][0] != 0.0f)
        return false;
    if (m[2][1] != 0.0f || m[2][2] != 1.0f || m[2][3] != 0.0f)
        return false;
    if (m[3][0] != 0.0f || m[3][1] != 0.0f || m[3][2] != 0.0f)
        return false;
    return (m[3][3] == 1.0f);
}

inline void QDoubleMatrix4x4::setToIdentity()
{
    m[0][0] = 1.0f;
    m[0][1] = 0.0f;
    m[0][2] = 0.0f;
    m[0][3] = 0.0f;
    m[1][0] = 0.0f;
    m[1][1] = 1.0f;
    m[1][2] = 0.0f;
    m[1][3] = 0.0f;
    m[2][0] = 0.0f;
    m[2][1] = 0.0f;
    m[2][2] = 1.0f;
    m[2][3] = 0.0f;
    m[3][0] = 0.0f;
    m[3][1] = 0.0f;
    m[3][2] = 0.0f;
    m[3][3] = 1.0f;
    flagBits = Identity;
}

inline void QDoubleMatrix4x4::fill(double value)
{
    m[0][0] = value;
    m[0][1] = value;
    m[0][2] = value;
    m[0][3] = value;
    m[1][0] = value;
    m[1][1] = value;
    m[1][2] = value;
    m[1][3] = value;
    m[2][0] = value;
    m[2][1] = value;
    m[2][2] = value;
    m[2][3] = value;
    m[3][0] = value;
    m[3][1] = value;
    m[3][2] = value;
    m[3][3] = value;
    flagBits = General;
}

inline QDoubleMatrix4x4& QDoubleMatrix4x4::operator+=(const QDoubleMatrix4x4& other)
{
    m[0][0] += other.m[0][0];
    m[0][1] += other.m[0][1];
    m[0][2] += other.m[0][2];
    m[0][3] += other.m[0][3];
    m[1][0] += other.m[1][0];
    m[1][1] += other.m[1][1];
    m[1][2] += other.m[1][2];
    m[1][3] += other.m[1][3];
    m[2][0] += other.m[2][0];
    m[2][1] += other.m[2][1];
    m[2][2] += other.m[2][2];
    m[2][3] += other.m[2][3];
    m[3][0] += other.m[3][0];
    m[3][1] += other.m[3][1];
    m[3][2] += other.m[3][2];
    m[3][3] += other.m[3][3];
    flagBits = General;
    return *this;
}

inline QDoubleMatrix4x4& QDoubleMatrix4x4::operator-=(const QDoubleMatrix4x4& other)
{
    m[0][0] -= other.m[0][0];
    m[0][1] -= other.m[0][1];
    m[0][2] -= other.m[0][2];
    m[0][3] -= other.m[0][3];
    m[1][0] -= other.m[1][0];
    m[1][1] -= other.m[1][1];
    m[1][2] -= other.m[1][2];
    m[1][3] -= other.m[1][3];
    m[2][0] -= other.m[2][0];
    m[2][1] -= other.m[2][1];
    m[2][2] -= other.m[2][2];
    m[2][3] -= other.m[2][3];
    m[3][0] -= other.m[3][0];
    m[3][1] -= other.m[3][1];
    m[3][2] -= other.m[3][2];
    m[3][3] -= other.m[3][3];
    flagBits = General;
    return *this;
}

inline QDoubleMatrix4x4& QDoubleMatrix4x4::operator*=(const QDoubleMatrix4x4& other)
{
    flagBits |= other.flagBits;

    if (flagBits < Rotation2D) {
        m[3][0] += m[0][0] * other.m[3][0];
        m[3][1] += m[1][1] * other.m[3][1];
        m[3][2] += m[2][2] * other.m[3][2];

        m[0][0] *= other.m[0][0];
        m[1][1] *= other.m[1][1];
        m[2][2] *= other.m[2][2];
        return *this;
    }

    double m0, m1, m2;
    m0 = m[0][0] * other.m[0][0]
            + m[1][0] * other.m[0][1]
            + m[2][0] * other.m[0][2]
            + m[3][0] * other.m[0][3];
    m1 = m[0][0] * other.m[1][0]
            + m[1][0] * other.m[1][1]
            + m[2][0] * other.m[1][2]
            + m[3][0] * other.m[1][3];
    m2 = m[0][0] * other.m[2][0]
            + m[1][0] * other.m[2][1]
            + m[2][0] * other.m[2][2]
            + m[3][0] * other.m[2][3];
    m[3][0] = m[0][0] * other.m[3][0]
            + m[1][0] * other.m[3][1]
            + m[2][0] * other.m[3][2]
            + m[3][0] * other.m[3][3];
    m[0][0] = m0;
    m[1][0] = m1;
    m[2][0] = m2;

    m0 = m[0][1] * other.m[0][0]
            + m[1][1] * other.m[0][1]
            + m[2][1] * other.m[0][2]
            + m[3][1] * other.m[0][3];
    m1 = m[0][1] * other.m[1][0]
            + m[1][1] * other.m[1][1]
            + m[2][1] * other.m[1][2]
            + m[3][1] * other.m[1][3];
    m2 = m[0][1] * other.m[2][0]
            + m[1][1] * other.m[2][1]
            + m[2][1] * other.m[2][2]
            + m[3][1] * other.m[2][3];
    m[3][1] = m[0][1] * other.m[3][0]
            + m[1][1] * other.m[3][1]
            + m[2][1] * other.m[3][2]
            + m[3][1] * other.m[3][3];
    m[0][1] = m0;
    m[1][1] = m1;
    m[2][1] = m2;

    m0 = m[0][2] * other.m[0][0]
            + m[1][2] * other.m[0][1]
            + m[2][2] * other.m[0][2]
            + m[3][2] * other.m[0][3];
    m1 = m[0][2] * other.m[1][0]
            + m[1][2] * other.m[1][1]
            + m[2][2] * other.m[1][2]
            + m[3][2] * other.m[1][3];
    m2 = m[0][2] * other.m[2][0]
            + m[1][2] * other.m[2][1]
            + m[2][2] * other.m[2][2]
            + m[3][2] * other.m[2][3];
    m[3][2] = m[0][2] * other.m[3][0]
            + m[1][2] * other.m[3][1]
            + m[2][2] * other.m[3][2]
            + m[3][2] * other.m[3][3];
    m[0][2] = m0;
    m[1][2] = m1;
    m[2][2] = m2;

    m0 = m[0][3] * other.m[0][0]
            + m[1][3] * other.m[0][1]
            + m[2][3] * other.m[0][2]
            + m[3][3] * other.m[0][3];
    m1 = m[0][3] * other.m[1][0]
            + m[1][3] * other.m[1][1]
            + m[2][3] * other.m[1][2]
            + m[3][3] * other.m[1][3];
    m2 = m[0][3] * other.m[2][0]
            + m[1][3] * other.m[2][1]
            + m[2][3] * other.m[2][2]
            + m[3][3] * other.m[2][3];
    m[3][3] = m[0][3] * other.m[3][0]
            + m[1][3] * other.m[3][1]
            + m[2][3] * other.m[3][2]
            + m[3][3] * other.m[3][3];
    m[0][3] = m0;
    m[1][3] = m1;
    m[2][3] = m2;
    return *this;
}

inline QDoubleMatrix4x4& QDoubleMatrix4x4::operator*=(double factor)
{
    m[0][0] *= factor;
    m[0][1] *= factor;
    m[0][2] *= factor;
    m[0][3] *= factor;
    m[1][0] *= factor;
    m[1][1] *= factor;
    m[1][2] *= factor;
    m[1][3] *= factor;
    m[2][0] *= factor;
    m[2][1] *= factor;
    m[2][2] *= factor;
    m[2][3] *= factor;
    m[3][0] *= factor;
    m[3][1] *= factor;
    m[3][2] *= factor;
    m[3][3] *= factor;
    flagBits = General;
    return *this;
}

inline bool QDoubleMatrix4x4::operator==(const QDoubleMatrix4x4& other) const
{
    return m[0][0] == other.m[0][0] &&
           m[0][1] == other.m[0][1] &&
           m[0][2] == other.m[0][2] &&
           m[0][3] == other.m[0][3] &&
           m[1][0] == other.m[1][0] &&
           m[1][1] == other.m[1][1] &&
           m[1][2] == other.m[1][2] &&
           m[1][3] == other.m[1][3] &&
           m[2][0] == other.m[2][0] &&
           m[2][1] == other.m[2][1] &&
           m[2][2] == other.m[2][2] &&
           m[2][3] == other.m[2][3] &&
           m[3][0] == other.m[3][0] &&
           m[3][1] == other.m[3][1] &&
           m[3][2] == other.m[3][2] &&
           m[3][3] == other.m[3][3];
}

inline bool QDoubleMatrix4x4::operator!=(const QDoubleMatrix4x4& other) const
{
    return m[0][0] != other.m[0][0] ||
           m[0][1] != other.m[0][1] ||
           m[0][2] != other.m[0][2] ||
           m[0][3] != other.m[0][3] ||
           m[1][0] != other.m[1][0] ||
           m[1][1] != other.m[1][1] ||
           m[1][2] != other.m[1][2] ||
           m[1][3] != other.m[1][3] ||
           m[2][0] != other.m[2][0] ||
           m[2][1] != other.m[2][1] ||
           m[2][2] != other.m[2][2] ||
           m[2][3] != other.m[2][3] ||
           m[3][0] != other.m[3][0] ||
           m[3][1] != other.m[3][1] ||
           m[3][2] != other.m[3][2] ||
           m[3][3] != other.m[3][3];
}

inline QDoubleMatrix4x4 operator+(const QDoubleMatrix4x4& m1, const QDoubleMatrix4x4& m2)
{
    QDoubleMatrix4x4 m(1);
    m.m[0][0] = m1.m[0][0] + m2.m[0][0];
    m.m[0][1] = m1.m[0][1] + m2.m[0][1];
    m.m[0][2] = m1.m[0][2] + m2.m[0][2];
    m.m[0][3] = m1.m[0][3] + m2.m[0][3];
    m.m[1][0] = m1.m[1][0] + m2.m[1][0];
    m.m[1][1] = m1.m[1][1] + m2.m[1][1];
    m.m[1][2] = m1.m[1][2] + m2.m[1][2];
    m.m[1][3] = m1.m[1][3] + m2.m[1][3];
    m.m[2][0] = m1.m[2][0] + m2.m[2][0];
    m.m[2][1] = m1.m[2][1] + m2.m[2][1];
    m.m[2][2] = m1.m[2][2] + m2.m[2][2];
    m.m[2][3] = m1.m[2][3] + m2.m[2][3];
    m.m[3][0] = m1.m[3][0] + m2.m[3][0];
    m.m[3][1] = m1.m[3][1] + m2.m[3][1];
    m.m[3][2] = m1.m[3][2] + m2.m[3][2];
    m.m[3][3] = m1.m[3][3] + m2.m[3][3];
    m.flagBits = QDoubleMatrix4x4::General;
    return m;
}

inline QDoubleMatrix4x4 operator-(const QDoubleMatrix4x4& m1, const QDoubleMatrix4x4& m2)
{
    QDoubleMatrix4x4 m(1);
    m.m[0][0] = m1.m[0][0] - m2.m[0][0];
    m.m[0][1] = m1.m[0][1] - m2.m[0][1];
    m.m[0][2] = m1.m[0][2] - m2.m[0][2];
    m.m[0][3] = m1.m[0][3] - m2.m[0][3];
    m.m[1][0] = m1.m[1][0] - m2.m[1][0];
    m.m[1][1] = m1.m[1][1] - m2.m[1][1];
    m.m[1][2] = m1.m[1][2] - m2.m[1][2];
    m.m[1][3] = m1.m[1][3] - m2.m[1][3];
    m.m[2][0] = m1.m[2][0] - m2.m[2][0];
    m.m[2][1] = m1.m[2][1] - m2.m[2][1];
    m.m[2][2] = m1.m[2][2] - m2.m[2][2];
    m.m[2][3] = m1.m[2][3] - m2.m[2][3];
    m.m[3][0] = m1.m[3][0] - m2.m[3][0];
    m.m[3][1] = m1.m[3][1] - m2.m[3][1];
    m.m[3][2] = m1.m[3][2] - m2.m[3][2];
    m.m[3][3] = m1.m[3][3] - m2.m[3][3];
    m.flagBits = QDoubleMatrix4x4::General;
    return m;
}

inline QDoubleMatrix4x4 operator*(const QDoubleMatrix4x4& m1, const QDoubleMatrix4x4& m2)
{
    int flagBits = m1.flagBits | m2.flagBits;
    if (flagBits < QDoubleMatrix4x4::Rotation2D) {
        QDoubleMatrix4x4 m = m1;
        m.m[3][0] += m.m[0][0] * m2.m[3][0];
        m.m[3][1] += m.m[1][1] * m2.m[3][1];
        m.m[3][2] += m.m[2][2] * m2.m[3][2];

        m.m[0][0] *= m2.m[0][0];
        m.m[1][1] *= m2.m[1][1];
        m.m[2][2] *= m2.m[2][2];
        m.flagBits = flagBits;
        return m;
    }

    QDoubleMatrix4x4 m(1);
    m.m[0][0] = m1.m[0][0] * m2.m[0][0]
              + m1.m[1][0] * m2.m[0][1]
              + m1.m[2][0] * m2.m[0][2]
              + m1.m[3][0] * m2.m[0][3];
    m.m[0][1] = m1.m[0][1] * m2.m[0][0]
              + m1.m[1][1] * m2.m[0][1]
              + m1.m[2][1] * m2.m[0][2]
              + m1.m[3][1] * m2.m[0][3];
    m.m[0][2] = m1.m[0][2] * m2.m[0][0]
              + m1.m[1][2] * m2.m[0][1]
              + m1.m[2][2] * m2.m[0][2]
              + m1.m[3][2] * m2.m[0][3];
    m.m[0][3] = m1.m[0][3] * m2.m[0][0]
              + m1.m[1][3] * m2.m[0][1]
              + m1.m[2][3] * m2.m[0][2]
              + m1.m[3][3] * m2.m[0][3];

    m.m[1][0] = m1.m[0][0] * m2.m[1][0]
              + m1.m[1][0] * m2.m[1][1]
              + m1.m[2][0] * m2.m[1][2]
              + m1.m[3][0] * m2.m[1][3];
    m.m[1][1] = m1.m[0][1] * m2.m[1][0]
              + m1.m[1][1] * m2.m[1][1]
              + m1.m[2][1] * m2.m[1][2]
              + m1.m[3][1] * m2.m[1][3];
    m.m[1][2] = m1.m[0][2] * m2.m[1][0]
              + m1.m[1][2] * m2.m[1][1]
              + m1.m[2][2] * m2.m[1][2]
              + m1.m[3][2] * m2.m[1][3];
    m.m[1][3] = m1.m[0][3] * m2.m[1][0]
              + m1.m[1][3] * m2.m[1][1]
              + m1.m[2][3] * m2.m[1][2]
              + m1.m[3][3] * m2.m[1][3];

    m.m[2][0] = m1.m[0][0] * m2.m[2][0]
              + m1.m[1][0] * m2.m[2][1]
              + m1.m[2][0] * m2.m[2][2]
              + m1.m[3][0] * m2.m[2][3];
    m.m[2][1] = m1.m[0][1] * m2.m[2][0]
              + m1.m[1][1] * m2.m[2][1]
              + m1.m[2][1] * m2.m[2][2]
              + m1.m[3][1] * m2.m[2][3];
    m.m[2][2] = m1.m[0][2] * m2.m[2][0]
              + m1.m[1][2] * m2.m[2][1]
              + m1.m[2][2] * m2.m[2][2]
              + m1.m[3][2] * m2.m[2][3];
    m.m[2][3] = m1.m[0][3] * m2.m[2][0]
              + m1.m[1][3] * m2.m[2][1]
              + m1.m[2][3] * m2.m[2][2]
              + m1.m[3][3] * m2.m[2][3];

    m.m[3][0] = m1.m[0][0] * m2.m[3][0]
              + m1.m[1][0] * m2.m[3][1]
              + m1.m[2][0] * m2.m[3][2]
              + m1.m[3][0] * m2.m[3][3];
    m.m[3][1] = m1.m[0][1] * m2.m[3][0]
              + m1.m[1][1] * m2.m[3][1]
              + m1.m[2][1] * m2.m[3][2]
              + m1.m[3][1] * m2.m[3][3];
    m.m[3][2] = m1.m[0][2] * m2.m[3][0]
              + m1.m[1][2] * m2.m[3][1]
              + m1.m[2][2] * m2.m[3][2]
              + m1.m[3][2] * m2.m[3][3];
    m.m[3][3] = m1.m[0][3] * m2.m[3][0]
              + m1.m[1][3] * m2.m[3][1]
              + m1.m[2][3] * m2.m[3][2]
              + m1.m[3][3] * m2.m[3][3];
    m.flagBits = flagBits;
    return m;
}

inline QDoubleVector3D operator*(const QDoubleVector3D& vector, const QDoubleMatrix4x4& matrix)
{
    double x, y, z, w;
    x = vector.x() * matrix.m[0][0] +
        vector.y() * matrix.m[0][1] +
        vector.z() * matrix.m[0][2] +
        matrix.m[0][3];
    y = vector.x() * matrix.m[1][0] +
        vector.y() * matrix.m[1][1] +
        vector.z() * matrix.m[1][2] +
        matrix.m[1][3];
    z = vector.x() * matrix.m[2][0] +
        vector.y() * matrix.m[2][1] +
        vector.z() * matrix.m[2][2] +
        matrix.m[2][3];
    w = vector.x() * matrix.m[3][0] +
        vector.y() * matrix.m[3][1] +
        vector.z() * matrix.m[3][2] +
        matrix.m[3][3];
    if (w == 1.0f)
        return QDoubleVector3D(x, y, z);
    else
        return QDoubleVector3D(x / w, y / w, z / w);
}

inline QDoubleVector3D operator*(const QDoubleMatrix4x4& matrix, const QDoubleVector3D& vector)
{
    double x, y, z, w;
    if (matrix.flagBits == QDoubleMatrix4x4::Identity) {
        return vector;
    } else if (matrix.flagBits < QDoubleMatrix4x4::Rotation2D) {
        // Translation | Scale
        return QDoubleVector3D(vector.x() * matrix.m[0][0] + matrix.m[3][0],
                         vector.y() * matrix.m[1][1] + matrix.m[3][1],
                         vector.z() * matrix.m[2][2] + matrix.m[3][2]);
    } else if (matrix.flagBits < QDoubleMatrix4x4::Rotation) {
        // Translation | Scale | Rotation2D
        return QDoubleVector3D(vector.x() * matrix.m[0][0] + vector.y() * matrix.m[1][0] + matrix.m[3][0],
                         vector.x() * matrix.m[0][1] + vector.y() * matrix.m[1][1] + matrix.m[3][1],
                         vector.z() * matrix.m[2][2] + matrix.m[3][2]);
    } else {
        x = vector.x() * matrix.m[0][0] +
            vector.y() * matrix.m[1][0] +
            vector.z() * matrix.m[2][0] +
            matrix.m[3][0];
        y = vector.x() * matrix.m[0][1] +
            vector.y() * matrix.m[1][1] +
            vector.z() * matrix.m[2][1] +
            matrix.m[3][1];
        z = vector.x() * matrix.m[0][2] +
            vector.y() * matrix.m[1][2] +
            vector.z() * matrix.m[2][2] +
            matrix.m[3][2];
        w = vector.x() * matrix.m[0][3] +
            vector.y() * matrix.m[1][3] +
            vector.z() * matrix.m[2][3] +
            matrix.m[3][3];
        if (w == 1.0f)
            return QDoubleVector3D(x, y, z);
        else
            return QDoubleVector3D(x / w, y / w, z / w);
    }
}

inline QPoint operator*(const QPoint& point, const QDoubleMatrix4x4& matrix)
{
    double xin, yin;
    double x, y, w;
    xin = point.x();
    yin = point.y();
    x = xin * matrix.m[0][0] +
        yin * matrix.m[0][1] +
        matrix.m[0][3];
    y = xin * matrix.m[1][0] +
        yin * matrix.m[1][1] +
        matrix.m[1][3];
    w = xin * matrix.m[3][0] +
        yin * matrix.m[3][1] +
        matrix.m[3][3];
    if (w == 1.0f)
        return QPoint(qRound(x), qRound(y));
    else
        return QPoint(qRound(x / w), qRound(y / w));
}

inline QPointF operator*(const QPointF& point, const QDoubleMatrix4x4& matrix)
{
    double xin, yin;
    double x, y, w;
    xin = point.x();
    yin = point.y();
    x = xin * matrix.m[0][0] +
        yin * matrix.m[0][1] +
        matrix.m[0][3];
    y = xin * matrix.m[1][0] +
        yin * matrix.m[1][1] +
        matrix.m[1][3];
    w = xin * matrix.m[3][0] +
        yin * matrix.m[3][1] +
        matrix.m[3][3];
    if (w == 1.0f) {
        return QPointF(double(x), double(y));
    } else {
        return QPointF(double(x / w), double(y / w));
    }
}

inline QPoint operator*(const QDoubleMatrix4x4& matrix, const QPoint& point)
{
    double xin, yin;
    double x, y, w;
    xin = point.x();
    yin = point.y();
    if (matrix.flagBits == QDoubleMatrix4x4::Identity) {
        return point;
    } else if (matrix.flagBits < QDoubleMatrix4x4::Rotation2D) {
        // Translation | Scale
        return QPoint(qRound(xin * matrix.m[0][0] + matrix.m[3][0]),
                      qRound(yin * matrix.m[1][1] + matrix.m[3][1]));
    } else if (matrix.flagBits < QDoubleMatrix4x4::Perspective) {
        return QPoint(qRound(xin * matrix.m[0][0] + yin * matrix.m[1][0] + matrix.m[3][0]),
                      qRound(xin * matrix.m[0][1] + yin * matrix.m[1][1] + matrix.m[3][1]));
    } else {
        x = xin * matrix.m[0][0] +
            yin * matrix.m[1][0] +
            matrix.m[3][0];
        y = xin * matrix.m[0][1] +
            yin * matrix.m[1][1] +
            matrix.m[3][1];
        w = xin * matrix.m[0][3] +
            yin * matrix.m[1][3] +
            matrix.m[3][3];
        if (w == 1.0f)
            return QPoint(qRound(x), qRound(y));
        else
            return QPoint(qRound(x / w), qRound(y / w));
    }
}

inline QPointF operator*(const QDoubleMatrix4x4& matrix, const QPointF& point)
{
    double xin, yin;
    double x, y, w;
    xin = point.x();
    yin = point.y();
    if (matrix.flagBits == QDoubleMatrix4x4::Identity) {
        return point;
    } else if (matrix.flagBits < QDoubleMatrix4x4::Rotation2D) {
        // Translation | Scale
        return QPointF(xin * matrix.m[0][0] + matrix.m[3][0],
                       yin * matrix.m[1][1] + matrix.m[3][1]);
    } else if (matrix.flagBits < QDoubleMatrix4x4::Perspective) {
        return QPointF(xin * matrix.m[0][0] + yin * matrix.m[1][0] + matrix.m[3][0],
                       xin * matrix.m[0][1] + yin * matrix.m[1][1] + matrix.m[3][1]);
    } else {
        x = xin * matrix.m[0][0] +
            yin * matrix.m[1][0] +
            matrix.m[3][0];
        y = xin * matrix.m[0][1] +
            yin * matrix.m[1][1] +
            matrix.m[3][1];
        w = xin * matrix.m[0][3] +
            yin * matrix.m[1][3] +
            matrix.m[3][3];
        if (w == 1.0f) {
            return QPointF(double(x), double(y));
        } else {
            return QPointF(double(x / w), double(y / w));
        }
    }
}

inline QDoubleMatrix4x4 operator-(const QDoubleMatrix4x4& matrix)
{
    QDoubleMatrix4x4 m(1);
    m.m[0][0] = -matrix.m[0][0];
    m.m[0][1] = -matrix.m[0][1];
    m.m[0][2] = -matrix.m[0][2];
    m.m[0][3] = -matrix.m[0][3];
    m.m[1][0] = -matrix.m[1][0];
    m.m[1][1] = -matrix.m[1][1];
    m.m[1][2] = -matrix.m[1][2];
    m.m[1][3] = -matrix.m[1][3];
    m.m[2][0] = -matrix.m[2][0];
    m.m[2][1] = -matrix.m[2][1];
    m.m[2][2] = -matrix.m[2][2];
    m.m[2][3] = -matrix.m[2][3];
    m.m[3][0] = -matrix.m[3][0];
    m.m[3][1] = -matrix.m[3][1];
    m.m[3][2] = -matrix.m[3][2];
    m.m[3][3] = -matrix.m[3][3];
    m.flagBits = QDoubleMatrix4x4::General;
    return m;
}

inline QDoubleMatrix4x4 operator*(double factor, const QDoubleMatrix4x4& matrix)
{
    QDoubleMatrix4x4 m(1);
    m.m[0][0] = matrix.m[0][0] * factor;
    m.m[0][1] = matrix.m[0][1] * factor;
    m.m[0][2] = matrix.m[0][2] * factor;
    m.m[0][3] = matrix.m[0][3] * factor;
    m.m[1][0] = matrix.m[1][0] * factor;
    m.m[1][1] = matrix.m[1][1] * factor;
    m.m[1][2] = matrix.m[1][2] * factor;
    m.m[1][3] = matrix.m[1][3] * factor;
    m.m[2][0] = matrix.m[2][0] * factor;
    m.m[2][1] = matrix.m[2][1] * factor;
    m.m[2][2] = matrix.m[2][2] * factor;
    m.m[2][3] = matrix.m[2][3] * factor;
    m.m[3][0] = matrix.m[3][0] * factor;
    m.m[3][1] = matrix.m[3][1] * factor;
    m.m[3][2] = matrix.m[3][2] * factor;
    m.m[3][3] = matrix.m[3][3] * factor;
    m.flagBits = QDoubleMatrix4x4::General;
    return m;
}

inline QDoubleMatrix4x4 operator*(const QDoubleMatrix4x4& matrix, double factor)
{
    QDoubleMatrix4x4 m(1);
    m.m[0][0] = matrix.m[0][0] * factor;
    m.m[0][1] = matrix.m[0][1] * factor;
    m.m[0][2] = matrix.m[0][2] * factor;
    m.m[0][3] = matrix.m[0][3] * factor;
    m.m[1][0] = matrix.m[1][0] * factor;
    m.m[1][1] = matrix.m[1][1] * factor;
    m.m[1][2] = matrix.m[1][2] * factor;
    m.m[1][3] = matrix.m[1][3] * factor;
    m.m[2][0] = matrix.m[2][0] * factor;
    m.m[2][1] = matrix.m[2][1] * factor;
    m.m[2][2] = matrix.m[2][2] * factor;
    m.m[2][3] = matrix.m[2][3] * factor;
    m.m[3][0] = matrix.m[3][0] * factor;
    m.m[3][1] = matrix.m[3][1] * factor;
    m.m[3][2] = matrix.m[3][2] * factor;
    m.m[3][3] = matrix.m[3][3] * factor;
    m.flagBits = QDoubleMatrix4x4::General;
    return m;
}

inline bool qFuzzyCompare(const QDoubleMatrix4x4& m1, const QDoubleMatrix4x4& m2)
{
    return qFuzzyCompare(m1.m[0][0], m2.m[0][0]) &&
           qFuzzyCompare(m1.m[0][1], m2.m[0][1]) &&
           qFuzzyCompare(m1.m[0][2], m2.m[0][2]) &&
           qFuzzyCompare(m1.m[0][3], m2.m[0][3]) &&
           qFuzzyCompare(m1.m[1][0], m2.m[1][0]) &&
           qFuzzyCompare(m1.m[1][1], m2.m[1][1]) &&
           qFuzzyCompare(m1.m[1][2], m2.m[1][2]) &&
           qFuzzyCompare(m1.m[1][3], m2.m[1][3]) &&
           qFuzzyCompare(m1.m[2][0], m2.m[2][0]) &&
           qFuzzyCompare(m1.m[2][1], m2.m[2][1]) &&
           qFuzzyCompare(m1.m[2][2], m2.m[2][2]) &&
           qFuzzyCompare(m1.m[2][3], m2.m[2][3]) &&
           qFuzzyCompare(m1.m[3][0], m2.m[3][0]) &&
           qFuzzyCompare(m1.m[3][1], m2.m[3][1]) &&
           qFuzzyCompare(m1.m[3][2], m2.m[3][2]) &&
           qFuzzyCompare(m1.m[3][3], m2.m[3][3]);
}

inline QPoint QDoubleMatrix4x4::map(const QPoint& point) const
{
    return *this * point;
}

inline QPointF QDoubleMatrix4x4::map(const QPointF& point) const
{
    return *this * point;
}

inline QDoubleVector3D QDoubleMatrix4x4::map(const QDoubleVector3D& point) const
{
    return *this * point;
}

inline QDoubleVector3D QDoubleMatrix4x4::mapVector(const QDoubleVector3D& vector) const
{
    if (flagBits < Scale) {
        // Translation
        return vector;
    } else if (flagBits < Rotation2D) {
        // Translation | Scale
        return QDoubleVector3D(vector.x() * m[0][0],
                         vector.y() * m[1][1],
                         vector.z() * m[2][2]);
    } else {
        return QDoubleVector3D(vector.x() * m[0][0] +
                         vector.y() * m[1][0] +
                         vector.z() * m[2][0],
                         vector.x() * m[0][1] +
                         vector.y() * m[1][1] +
                         vector.z() * m[2][1],
                         vector.x() * m[0][2] +
                         vector.y() * m[1][2] +
                         vector.z() * m[2][2]);
    }
}

inline double *QDoubleMatrix4x4::data()
{
    // We have to assume that the caller will modify the matrix elements,
    // so we flip it over to "General" mode.
    flagBits = General;
    return *m;
}

inline void QDoubleMatrix4x4::viewport(const QRectF &rect)
{
    viewport(rect.x(), rect.y(), rect.width(), rect.height());
}

#ifndef QT_NO_DEBUG_STREAM
Q_POSITIONING_PRIVATE_EXPORT QDebug operator<<(QDebug dbg, const QDoubleMatrix4x4 &m);
#endif

#ifndef QT_NO_DATASTREAM
Q_POSITIONING_PRIVATE_EXPORT QDataStream &operator<<(QDataStream &, const QDoubleMatrix4x4 &);
Q_POSITIONING_PRIVATE_EXPORT QDataStream &operator>>(QDataStream &, QDoubleMatrix4x4 &);
#endif


QT_END_NAMESPACE


#endif // QDOUBLEMATRIX4X4_H
