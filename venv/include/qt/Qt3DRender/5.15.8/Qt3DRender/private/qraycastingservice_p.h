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

#ifndef QT3DRENDER_QRAYCASTINGSERVICE_P_H
#define QT3DRENDER_QRAYCASTINGSERVICE_P_H

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

#include <Qt3DCore/qt3dcore_global.h>
#include <Qt3DRender/private/qabstractcollisionqueryservice_p.h>
#include <Qt3DRender/private/qray3d_p.h>

#include <QHash>
#include <QtConcurrent>
#include <QAtomicInt>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {
namespace RayCasting {

class QBoundingVolumeProvider;
class QRayCastingServicePrivate;

typedef QFuture<QCollisionQueryResult> FutureQueryResult;

class Q_3DRENDERSHARED_EXPORT QRayCastingService : public QAbstractCollisionQueryService
{
public:
    QRayCastingService();

    QQueryHandle query(const QRay3D &ray, QueryMode mode, QBoundingVolumeProvider *provider) override;
    QCollisionQueryResult::Hit query(const QRay3D &ray, const QBoundingVolume *volume) override;

    QCollisionQueryResult fetchResult(const QQueryHandle &handle) override;
    QVector<QCollisionQueryResult> fetchAllResults() const override;

protected:
    QRayCastingService(QRayCastingServicePrivate &dd);

private:
    Q_DISABLE_COPY(QRayCastingService)
    Q_DECLARE_PRIVATE(QRayCastingService)
};

class QRayCastingServicePrivate : public QAbstractCollisionQueryServicePrivate
{
public:
    QRayCastingServicePrivate(const QString &description);

    QCollisionQueryResult collides(const QRay3D &ray,
                                   QBoundingVolumeProvider *provider,
                                   QAbstractCollisionQueryService::QueryMode mode,
                                   const QQueryHandle &handle);
    QCollisionQueryResult::Hit collides(const QRay3D &ray,
                                        const QBoundingVolume *volume);

    Q_DECLARE_PUBLIC(QRayCastingService)

    struct Query
    {
        QQueryHandle handle;
        QRay3D ray;
        QRayCastingService::QueryMode mode;
    };

#if QT_CONFIG(concurrent)
    QHash<QQueryHandle, FutureQueryResult> m_results;
#else
    QHash<QQueryHandle, QCollisionQueryResult> m_results;
#endif
    QAtomicInt m_handlesCount;
};

} // namespace RayCasting
} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QRAYCASTINGSERVICE_P_H
