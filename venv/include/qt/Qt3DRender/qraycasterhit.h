/****************************************************************************
**
** Copyright (C) 2018 Klaralvdalens Datakonsult AB (KDAB).
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

#ifndef QT3DRENDER_QRAYCASTERHIT_H
#define QT3DRENDER_QRAYCASTERHIT_H

#include <Qt3DCore/qcomponent.h>
#include <Qt3DRender/qt3drender_global.h>

#include <QtCore/QSharedData>
#include <QtGui/QVector3D>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QRayCasterHitData;
class QAbstractRayCasterPrivate;

class Q_3DRENDERSHARED_EXPORT QRayCasterHit
{
    Q_GADGET
public:
    enum HitType {
        TriangleHit,
        LineHit,
        PointHit,
        EntityHit
    };
    Q_ENUM(HitType)

    QRayCasterHit();
    explicit QRayCasterHit(QRayCasterHit::HitType type, Qt3DCore::QNodeId id, float distance,
                  const QVector3D &localIntersect, const QVector3D &worldIntersect,
                  uint primitiveIndex, uint v1, uint v2, uint v3);
    QRayCasterHit(const QRayCasterHit &other);
    ~QRayCasterHit();

    QRayCasterHit &operator =(const QRayCasterHit &other);

    HitType type() const;
    Qt3DCore::QNodeId entityId() const;
    Qt3DCore::QEntity *entity() const;
    float distance() const;
    QVector3D localIntersection() const;
    QVector3D worldIntersection() const;
    uint primitiveIndex() const;
    uint vertex1Index() const;
    uint vertex2Index() const;
    uint vertex3Index() const;

private:
    friend class QAbstractRayCasterPrivate;
    void setEntity(Qt3DCore::QEntity *entity) const;

    QSharedDataPointer<QRayCasterHitData> d;
};

} // Qt3D

QT_END_NAMESPACE

#endif // QT3DRENDER_QRAYCASTERHIT_H
