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

#ifndef QSSGPLANE_H
#define QSSGPLANE_H

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

#include <QtQuick3DUtils/private/qtquick3dutilsglobal_p.h>
#include <QtGui/QVector3D>
#include <QtCore/qmath.h>

QT_BEGIN_NAMESPACE

/**
\brief Representation of a plane.

 Plane equation used: n.dot(v) + d = 0
*/
class Q_QUICK3DUTILS_EXPORT QSSGPlane
{
public:
    QSSGPlane() = default;

    /**
    \brief Constructor from a normal and a distance
    */
    Q_ALWAYS_INLINE QSSGPlane(float nx, float ny, float nz, float distance) : n(nx, ny, nz), d(distance) {}

    /**
    \brief Constructor from a normal and a distance
    */
    Q_ALWAYS_INLINE QSSGPlane(const QVector3D &normal, float distance) : n(normal), d(distance) {}

    /**
    \brief Constructor from a point on the plane and a normal
    */
    Q_ALWAYS_INLINE QSSGPlane(const QVector3D &point, const QVector3D &normal)
        : n(normal), d(-QVector3D::dotProduct(point, n)) // p satisfies normal.dot(p) + d = 0
    {
    }

    /**
    \brief Constructor from three points
    */
    Q_ALWAYS_INLINE QSSGPlane(const QVector3D &p0, const QVector3D &p1, const QVector3D &p2)
    {
        n = QVector3D::crossProduct(p1 - p0, p2 - p0).normalized();
        d = QVector3D::dotProduct(-p0, n);
    }

    Q_ALWAYS_INLINE float distance(const QVector3D &p) const { return QVector3D::dotProduct(p, n) + d; }

    Q_ALWAYS_INLINE bool contains(const QVector3D &p) const { return qAbs(distance(p)) < (1.0e-7f); }

    /**
    \brief projects p into the plane
    */
    Q_ALWAYS_INLINE QVector3D project(const QVector3D &p) const { return p - n * distance(p); }

    /**
    \brief find an arbitrary point in the plane
    */
    Q_ALWAYS_INLINE QVector3D pointInPlane() const { return -n * d; }

    /**
    \brief equivalent plane with unit normal
    */

    Q_ALWAYS_INLINE void normalize();

    QVector3D n; //!< The normal to the plane
    float d = 0.0f; //!< The distance from the origin
};

QT_END_NAMESPACE

#endif // QSSGPLANE_H
