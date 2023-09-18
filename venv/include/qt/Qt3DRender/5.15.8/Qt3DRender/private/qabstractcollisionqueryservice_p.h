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

#ifndef QT3DRENDER_QABSTRACTCOLLISIONQUERYSERVICE_P_H
#define QT3DRENDER_QABSTRACTCOLLISIONQUERYSERVICE_P_H

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

#include <QVector>

#include <Qt3DRender/qt3drender_global.h>
#include <Qt3DCore/private/qservicelocator_p.h>
#include <Qt3DCore/qnodeid.h>
#include <Qt3DCore/private/qabstractserviceprovider_p.h>
#include <Qt3DRender/private/qcollisionqueryresult_p.h>

QT_BEGIN_NAMESPACE

class QAbstractCollisionQueryServicePrivate;
namespace Qt3DRender {
namespace RayCasting {

class QRay3D;
class QBoundingVolume;
class QBoundingVolumeProvider;

class QAbstractCollisionQueryServicePrivate : public Qt3DCore::QAbstractServiceProviderPrivate
{
public:
    QAbstractCollisionQueryServicePrivate(const QString &description)
        : QAbstractServiceProviderPrivate(Qt3DCore::QServiceLocator::CollisionService, description)
    {}
};

class Q_3DRENDERSHARED_EXPORT QAbstractCollisionQueryService : public Qt3DCore::QAbstractServiceProvider
{
    Q_OBJECT
public:
    enum QueryMode {
        FirstHit,
        AllHits
    };

    virtual QQueryHandle query(const QRay3D &ray, QueryMode mode, QBoundingVolumeProvider *provider) = 0;
    virtual QCollisionQueryResult::Hit query(const QRay3D &ray, const QBoundingVolume* volume) = 0;

    virtual QCollisionQueryResult fetchResult(const QQueryHandle &handle) = 0;
    virtual QVector<QCollisionQueryResult> fetchAllResults() const = 0;

protected:
    QAbstractCollisionQueryService(const QString &description = QString());
    QAbstractCollisionQueryService(QAbstractCollisionQueryServicePrivate &dd);

    void setResultHandle(QCollisionQueryResult &result, const QQueryHandle &handle);
    void addEntityHit(QCollisionQueryResult &result, Qt3DCore::QNodeId entity, const Vector3D &intersection,
                      float distance, const Vector3D &uvw);

private:
    Q_DECLARE_PRIVATE(QAbstractCollisionQueryService)
};

} // RayCasting
} // Qt3DRender

QT_END_NAMESPACE

Q_DECLARE_METATYPE(Qt3DRender::RayCasting::QAbstractCollisionQueryService::QueryMode) // LCOV_EXCL_LINE

#endif // QT3DRENDER_QABSTRACTCOLLISIONQUERYSERVICE_P_H
