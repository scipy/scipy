/****************************************************************************
**
** Copyright (C) 2018 Klaralvdalens Datakonsult AB (KDAB).
** Contact: http://www.qt-project.org/legal
**
** This file is part of the Qt3D module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QT3DRENDER_QABSTRACTRAYCASTER_P_H
#define QT3DRENDER_QABSTRACTRAYCASTER_P_H

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

#include <Qt3DCore/private/qcomponent_p.h>
#include <Qt3DCore/qnodeid.h>
#include <Qt3DRender/QAbstractRayCaster>
#include <Qt3DRender/private/qt3drender_global_p.h>

#include <QtGui/QVector3D>

QT_BEGIN_NAMESPACE

namespace Qt3DRender {

class QLayer;

/*!
    \internal
*/
class Q_3DRENDERSHARED_PRIVATE_EXPORT QAbstractRayCasterPrivate : public Qt3DCore::QComponentPrivate
{
public:
    QAbstractRayCasterPrivate();

    enum RayCasterType {
        WorldSpaceRayCaster,
        ScreenScapeRayCaster
    };

    static QAbstractRayCasterPrivate *get(QAbstractRayCaster *obj);
    static const QAbstractRayCasterPrivate *get(const QAbstractRayCaster *obj);
    static void updateHitEntites(QAbstractRayCaster::Hits &hits, Qt3DCore::QScene *scene);

    RayCasterType m_rayCasterType = WorldSpaceRayCaster;
    QAbstractRayCaster::RunMode m_runMode = QAbstractRayCaster::SingleShot;
    QAbstractRayCaster::Hits m_hits;
    QPoint m_position;
    QVector3D m_origin;
    QVector3D m_direction = QVector3D(0., 0., 1.f);
    float m_length = 1.f;
    QAbstractRayCaster::FilterMode m_filterMode = QAbstractRayCaster::AcceptAnyMatchingLayers;
    QVector<QLayer*> m_layers;

    virtual void dispatchHits(const QAbstractRayCaster::Hits &hits);

    Q_DECLARE_PUBLIC(QAbstractRayCaster)
};

struct QAbstractRayCasterData
{
    QAbstractRayCaster::RunMode runMode;
    QAbstractRayCasterPrivate::RayCasterType casterType;
    QPoint position;
    QVector3D origin;
    QVector3D direction;
    float length = 1.f;
    QAbstractRayCaster::FilterMode filterMode;
    Qt3DCore::QNodeIdVector layerIds;
};

} // namespace Qt3DRender

QT_END_NAMESPACE

#endif // QT3DRENDER_QABSTRACTRAYCASTER_P_H
