/****************************************************************************
**
** Copyright (C) 2015 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QCOLLISIONQUERYRESULT_P_H
#define QT3DRENDER_QCOLLISIONQUERYRESULT_P_H

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

#include <Qt3DRender/qt3drender_global.h>
#include <Qt3DCore/qnodeid.h>
#include <Qt3DCore/private/vector3d_p.h>
#include <QVector>
#include <QSharedData>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
namespace RayCasting {

typedef int QQueryHandle;
class QCollisionQueryResultPrivate;

class Q_3DRENDERSHARED_EXPORT QCollisionQueryResult
{
public:
    struct Hit {
        enum HitType {
            Entity,
            Point,
            Edge,
            Triangle
        };

        Hit()
            : m_type(Entity)
            , m_distance(-1.f)
            , m_primitiveIndex(0)
        {
            m_vertexIndex[0] = m_vertexIndex[1] = m_vertexIndex[2] = 0;
        }

        Hit(Qt3DCore::QNodeId entity, const Vector3D &intersection, float distance, const Vector3D &uvw)
            : m_entityId(entity)
            , m_type(Entity)
            , m_intersection(intersection)
            , m_distance(distance)
            , m_primitiveIndex(0U)
            , m_uvw(uvw)
        {
        }

        Qt3DCore::QNodeId m_entityId;
        HitType m_type;
        Vector3D m_intersection;
        float m_distance;
        uint m_primitiveIndex;
        uint m_vertexIndex[3];
        Vector3D m_uvw;
    };

    QCollisionQueryResult();
    QCollisionQueryResult(const QCollisionQueryResult &);
    ~QCollisionQueryResult();

    QCollisionQueryResult &operator=(const QCollisionQueryResult &);
#ifdef Q_COMPILER_RVALUE_REFS
    QCollisionQueryResult &operator=(QCollisionQueryResult &&other) Q_DECL_NOTHROW
    {
        swap(other);
        return *this;
    }
#endif

    void swap(QCollisionQueryResult &other) Q_DECL_NOTHROW
    {
        qSwap(d_ptr, other.d_ptr);
    }

    QQueryHandle handle() const;
    QVector<Hit> hits() const;
    QVector<Qt3DCore::QNodeId> entitiesHit() const;

private:
    friend class QAbstractCollisionQueryService;

    explicit QCollisionQueryResult(QCollisionQueryResultPrivate &p);

    QSharedDataPointer<QCollisionQueryResultPrivate> d_ptr;
    // Q_DECLARE_PRIVATE equivalent for shared data pointers
    QCollisionQueryResultPrivate *d_func();
    inline const QCollisionQueryResultPrivate *d_func() const
    {
        return d_ptr.constData();
    }
};
QT3D_DECLARE_TYPEINFO_2(Qt3DRender, RayCasting, QCollisionQueryResult::Hit, Q_PRIMITIVE_TYPE)
QT3D_DECLARE_SHARED_2(Qt3DRender, RayCasting, QCollisionQueryResult)

class QCollisionQueryResultPrivate : public QSharedData
{
public:
    explicit QCollisionQueryResultPrivate();
    explicit QCollisionQueryResultPrivate(const QCollisionQueryResultPrivate &copy);

    void setHandle(const QQueryHandle &handle);
    void addEntityHit(Qt3DCore::QNodeId entity, const Vector3D& intersection, float distance,
                      const Vector3D& uvw);

    QQueryHandle m_handle;
    QVector<QCollisionQueryResult::Hit> m_hits;
};

inline bool operator==(const QCollisionQueryResult::Hit& left, const QCollisionQueryResult::Hit& right)
{
    return left.m_entityId == right.m_entityId;
}

} // RayCasting
} // Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QCOLLISIONQUERYRESULT_P_H
