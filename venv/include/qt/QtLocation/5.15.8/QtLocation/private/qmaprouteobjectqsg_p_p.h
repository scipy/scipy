/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtLocation module of the Qt Toolkit.
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

#ifndef QMAPROUTEOBJECTQSG_P_P_H
#define QMAPROUTEOBJECTQSG_P_P_H

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

#include <QtLocation/private/qlocationglobal_p.h>
#include <QtLocation/private/qmappolylineobjectqsg_p_p.h>
#include <QtLocation/private/qmaprouteobject_p_p.h>
#include <QtLocation/private/qdeclarativegeoroute_p.h>
#include <QtLocation/private/qmaprouteobject_p.h>
#include <QtLocation/private/qqsgmapobject_p.h>
#include <QtCore/qscopedvaluerollback.h>
#include <QtCore/qscopedpointer.h>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QMapRouteObjectPrivateQSG : public QMapRouteObjectPrivate, public QQSGMapObject
{
public:
    QMapRouteObjectPrivateQSG(QGeoMapObject *q);
    QMapRouteObjectPrivateQSG(const QMapRouteObjectPrivate &other);
    ~QMapRouteObjectPrivateQSG() override;

    // QQSGMapObject
    void updateGeometry() override;
    QSGNode *updateMapObjectNode(QSGNode *oldNode,
                                 VisibleNode **visibleNode,
                                 QSGNode *root,
                                 QQuickWindow *window) override;

    // QMapRouteObjectPrivate interface
    void setRoute(const QDeclarativeGeoRoute *route) override;

    // QGeoMapObjectPrivate interface
    QGeoMapObjectPrivate *clone() override;
    void setMap(QGeoMap *map) override;
    void setVisible(bool visible) override;
    virtual QGeoShape geoShape() const override;

    // Data Members
    QScopedPointer<QMapPolylineObjectPrivateQSG> m_polyline;
};

QT_END_NAMESPACE

#endif // QMAPROUTEOBJECTQSG_P_P_H
