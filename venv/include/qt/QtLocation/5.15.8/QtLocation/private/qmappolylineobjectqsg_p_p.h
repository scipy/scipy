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

#ifndef QMAPPOLYLINEOBJECTOBJECTSOVERLAY_H
#define QMAPPOLYLINEOBJECTOBJECTSOVERLAY_H

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
#include <QtLocation/private/qmappolylineobject_p_p.h>
#include <QtLocation/private/qdeclarativepolylinemapitem_p.h>
#include <QtLocation/private/qmappolylineobject_p.h>
#include <QtLocation/private/qqsgmapobject_p.h>
#include <QtCore/qscopedvaluerollback.h>
#include <QtPositioning/private/qgeopath_p.h>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QMapPolylineObjectPrivateQSG : public QMapPolylineObjectPrivateDefault, public QQSGMapObject
{
public:
    QMapPolylineObjectPrivateQSG(QGeoMapObject *q);
    QMapPolylineObjectPrivateQSG(const QMapPolylineObjectPrivate &other);
    ~QMapPolylineObjectPrivateQSG() override;

    QList<QDoubleVector2D> projectPath();

    // QQSGMapObject
    void markSourceDirty();
    void updateGeometry() override;
    QSGNode *updateMapObjectNode(QSGNode *oldNode,
                                 VisibleNode **visibleNode,
                                 QSGNode *root,
                                 QQuickWindow *window) override;

    // QGeoMapPolylinePrivate interface
    QList<QGeoCoordinate> path() const override;
    void setPath(const QList<QGeoCoordinate> &path) override;
    void setColor(const QColor &color) override;
    void setWidth(qreal width) override;

    // QGeoMapObjectPrivate
    QGeoMapObjectPrivate *clone() override;
    virtual QGeoShape geoShape() const override;

    unsigned int zoomForLOD(int zoom) const;

    // Data Members
    QDoubleVector2D m_leftBoundMercator;
    QGeoMapPolylineGeometryOpenGL m_borderGeometry;
    MapPolylineNodeOpenGLExtruded *m_polylinenode = nullptr;
};

QT_END_NAMESPACE

#endif // QMAPPOLYLINEOBJECTOBJECTSOVERLAY_H
