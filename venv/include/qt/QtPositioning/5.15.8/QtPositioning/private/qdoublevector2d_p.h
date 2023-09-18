/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtPositioning module of the Qt Toolkit.
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

#ifndef QDOUBLEVECTOR2D_P_H
#define QDOUBLEVECTOR2D_P_H

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

#ifdef QT_BUILD_LOCATION_LIB
#include <QVector2D>
#endif

#include "qpositioningglobal_p.h"
#include <QtCore/qmetatype.h>
#include <QPointF>

QT_BEGIN_NAMESPACE

class QDoubleVector3D;

class Q_POSITIONING_PRIVATE_EXPORT QDoubleVector2D
{
public:
    Q_DECL_CONSTEXPR inline QDoubleVector2D();
    Q_DECL_CONSTEXPR inline QDoubleVector2D(double xpos, double ypos);
    Q_DECL_CONSTEXPR explicit inline QDoubleVector2D(const QPointF &p);
    explicit QDoubleVector2D(const QDoubleVector3D &vector);

    Q_DECL_CONSTEXPR inline double manhattanLength() const;
    inline bool isNull() const;
    inline bool isFinite() const;

    Q_DECL_CONSTEXPR inline double x() const;
    Q_DECL_CONSTEXPR inline double y() const;

    inline void setX(double x);
    inline void setY(double y);

    double length() const;
    Q_DECL_CONSTEXPR inline double lengthSquared() const;

    QDoubleVector2D normalized() const;
    void normalize();

    inline QDoubleVector2D &operator+=(const QDoubleVector2D &vector);
    inline QDoubleVector2D &operator-=(const QDoubleVector2D &vector);
    inline QDoubleVector2D &operator*=(double factor);
    inline QDoubleVector2D &operator*=(const QDoubleVector2D &vector);
    inline QDoubleVector2D &operator/=(double divisor);
    inline QDoubleVector2D &operator/=(const QDoubleVector2D &vector);

    Q_DECL_CONSTEXPR static inline double dotProduct(const QDoubleVector2D &v1, const QDoubleVector2D &v2)
    { return v1.xp * v2.xp + v1.yp * v2.yp; }


    friend Q_DECL_CONSTEXPR inline bool operator==(const QDoubleVector2D &v1, const QDoubleVector2D &v2);
    friend Q_DECL_CONSTEXPR inline bool operator!=(const QDoubleVector2D &v1, const QDoubleVector2D &v2);
    friend Q_DECL_CONSTEXPR inline const QDoubleVector2D operator+(const QDoubleVector2D &v1, const QDoubleVector2D &v2);
    friend Q_DECL_CONSTEXPR inline const QDoubleVector2D operator-(const QDoubleVector2D &v1, const QDoubleVector2D &v2);
    friend Q_DECL_CONSTEXPR inline const QDoubleVector2D operator*(double factor, const QDoubleVector2D &vector);
    friend Q_DECL_CONSTEXPR inline const QDoubleVector2D operator*(const QDoubleVector2D &vector, double factor);
    friend Q_DECL_CONSTEXPR inline const QDoubleVector2D operator*(const QDoubleVector2D &v1, const QDoubleVector2D &v2);
    friend Q_DECL_CONSTEXPR inline const QDoubleVector2D operator-(const QDoubleVector2D &vector);
    friend Q_DECL_CONSTEXPR inline const QDoubleVector2D operator/(const QDoubleVector2D &vector, double divisor);

    friend Q_DECL_CONSTEXPR inline bool qFuzzyCompare(const QDoubleVector2D &v1, const QDoubleVector2D &v2);

    QDoubleVector3D toVector3D() const;
    Q_DECL_CONSTEXPR inline QPointF toPointF() const;

private:
    double xp, yp;

    friend class QDoubleVector3D;
};

Q_DECLARE_TYPEINFO(QDoubleVector2D, Q_MOVABLE_TYPE);

Q_DECL_CONSTEXPR inline QDoubleVector2D::QDoubleVector2D() : xp(0.0), yp(0.0) {}

Q_DECL_CONSTEXPR inline QDoubleVector2D::QDoubleVector2D(double xpos, double ypos) : xp(xpos), yp(ypos) {}

Q_DECL_CONSTEXPR inline QDoubleVector2D::QDoubleVector2D(const QPointF &p) : xp(p.x()), yp(p.y()) { }

Q_DECL_CONSTEXPR inline double QDoubleVector2D::manhattanLength() const
{
    return qAbs(x())+qAbs(y());
}

inline bool QDoubleVector2D::isNull() const
{
    return qIsNull(xp) && qIsNull(yp);
}

inline bool QDoubleVector2D::isFinite() const
{
    return qIsFinite(xp) && qIsFinite(yp);
}

Q_DECL_CONSTEXPR inline double QDoubleVector2D::x() const { return xp; }
Q_DECL_CONSTEXPR inline double QDoubleVector2D::y() const { return yp; }

inline void QDoubleVector2D::setX(double aX) { xp = aX; }
inline void QDoubleVector2D::setY(double aY) { yp = aY; }

Q_DECL_CONSTEXPR inline double QDoubleVector2D::lengthSquared() const
{ return xp * xp + yp * yp; }

inline QDoubleVector2D &QDoubleVector2D::operator+=(const QDoubleVector2D &vector)
{
    xp += vector.xp;
    yp += vector.yp;
    return *this;
}

inline QDoubleVector2D &QDoubleVector2D::operator-=(const QDoubleVector2D &vector)
{
    xp -= vector.xp;
    yp -= vector.yp;
    return *this;
}

inline QDoubleVector2D &QDoubleVector2D::operator*=(double factor)
{
    xp *= factor;
    yp *= factor;
    return *this;
}

inline QDoubleVector2D &QDoubleVector2D::operator*=(const QDoubleVector2D &vector)
{
    xp *= vector.xp;
    yp *= vector.yp;
    return *this;
}

inline QDoubleVector2D &QDoubleVector2D::operator/=(double divisor)
{
    xp /= divisor;
    yp /= divisor;
    return *this;
}

inline QDoubleVector2D &QDoubleVector2D::operator/=(const QDoubleVector2D &vector)
{
    xp /= vector.xp;
    yp /= vector.yp;
    return *this;
}

Q_DECL_CONSTEXPR inline bool operator==(const QDoubleVector2D &v1, const QDoubleVector2D &v2)
{
    return v1.xp == v2.xp && v1.yp == v2.yp;
}

Q_DECL_CONSTEXPR inline bool operator!=(const QDoubleVector2D &v1, const QDoubleVector2D &v2)
{
    return v1.xp != v2.xp || v1.yp != v2.yp;
}

Q_DECL_CONSTEXPR inline const QDoubleVector2D operator+(const QDoubleVector2D &v1, const QDoubleVector2D &v2)
{
    return QDoubleVector2D(v1.xp + v2.xp, v1.yp + v2.yp);
}

Q_DECL_CONSTEXPR inline const QDoubleVector2D operator-(const QDoubleVector2D &v1, const QDoubleVector2D &v2)
{
    return QDoubleVector2D(v1.xp - v2.xp, v1.yp - v2.yp);
}

Q_DECL_CONSTEXPR inline const QDoubleVector2D operator*(double factor, const QDoubleVector2D &vector)
{
    return QDoubleVector2D(vector.xp * factor, vector.yp * factor);
}

Q_DECL_CONSTEXPR inline const QDoubleVector2D operator*(const QDoubleVector2D &vector, double factor)
{
    return QDoubleVector2D(vector.xp * factor, vector.yp * factor);
}

Q_DECL_CONSTEXPR inline const QDoubleVector2D operator*(const QDoubleVector2D &v1, const QDoubleVector2D &v2)
{
    return QDoubleVector2D(v1.xp * v2.xp, v1.yp * v2.yp);
}

Q_DECL_CONSTEXPR inline const QDoubleVector2D operator-(const QDoubleVector2D &vector)
{
    return QDoubleVector2D(-vector.xp, -vector.yp);
}

Q_DECL_CONSTEXPR inline const QDoubleVector2D operator/(const QDoubleVector2D &vector, double divisor)
{
    return QDoubleVector2D(vector.xp / divisor, vector.yp / divisor);
}

Q_DECL_CONSTEXPR inline bool qFuzzyCompare(const QDoubleVector2D &v1, const QDoubleVector2D &v2)
{
    return qFuzzyCompare(v1.xp, v2.xp) && qFuzzyCompare(v1.yp, v2.yp);
}

Q_DECL_CONSTEXPR inline QPointF QDoubleVector2D::toPointF() const
{
    return QPointF(qreal(xp), qreal(yp));
}

#ifndef QT_NO_DEBUG_STREAM
Q_POSITIONING_EXPORT QDebug operator<<(QDebug dbg, const QDoubleVector2D &vector);
#endif

#ifndef QT_NO_DATASTREAM
Q_POSITIONING_EXPORT QDataStream &operator<<(QDataStream &, const QDoubleVector2D &);
Q_POSITIONING_EXPORT QDataStream &operator>>(QDataStream &, QDoubleVector2D &);
#endif

QT_END_NAMESPACE

#endif
