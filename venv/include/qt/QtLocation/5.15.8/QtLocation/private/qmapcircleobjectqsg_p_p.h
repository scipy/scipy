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

#ifndef QMAPCIRCLEOBJECTQSG_P_H
#define QMAPCIRCLEOBJECTQSG_P_H

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
#include <QtLocation/private/qgeomapobject_p_p.h>
#include <QtLocation/private/qdeclarativecirclemapitem_p_p.h>
#include <QtLocation/private/qdeclarativepolygonmapitem_p.h>
#include <QtLocation/private/qmapcircleobject_p.h>
#include <QtLocation/private/qmapcircleobject_p_p.h>
#include <QtLocation/private/qqsgmapobject_p.h>
#include <QtCore/qscopedvaluerollback.h>
#include <QGeoCoordinate>
#include <QColor>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QMapCircleObjectPrivateQSG : public QMapCircleObjectPrivateDefault, public QQSGMapObject
{
public:
    QMapCircleObjectPrivateQSG(QGeoMapObject *q);
    QMapCircleObjectPrivateQSG(const QMapCircleObjectPrivate &other);
    ~QMapCircleObjectPrivateQSG() override;

    // QQSGMapObject
    void updateGeometry() override;
    void updateGeometryCPU();
    void updateGeometryGL();
    QSGNode *updateMapObjectNode(QSGNode *oldNode,
                                 VisibleNode **visibleNode,
                                 QSGNode *root,
                                 QQuickWindow *window) override;
    QSGNode *updateMapObjectNodeCPU(QSGNode *oldNode,
                                 VisibleNode **visibleNode,
                                 QSGNode *root,
                                 QQuickWindow *window);
    QSGNode *updateMapObjectNodeGL(QSGNode *oldNode,
                                 VisibleNode **visibleNode,
                                 QSGNode *root,
                                 QQuickWindow *window);

    // QGeoMapCirclePrivate interface
    void setCenter(const QGeoCoordinate &center) override;
    void setRadius(qreal radius) override;
    void setColor(const QColor &color) override;
    void setBorderColor(const QColor &color) override;
    void setBorderWidth(qreal width) override;

    // QGeoMapObjectPrivate
    QGeoMapObjectPrivate *clone() override;

    void switchToGL();
    void switchToCPU();

public:
    // Data Members
struct CircleDataCPU {
    MapPolygonNode *m_node = nullptr;
    QList<QDoubleVector2D> m_circlePath;
    QGeoCoordinate m_leftBound;
    QGeoMapCircleGeometry m_geometry;
    QGeoMapPolylineGeometry m_borderGeometry;
    bool m_updatingGeometry = false;

    void updateCirclePath(const QGeoCoordinate &center, qreal radius, const QGeoProjectionWebMercator &p);
};
struct CircleDataGL {
    QList<QGeoCoordinate> m_circlePath;
    QGeoCoordinate m_leftBound;
    QDoubleVector2D m_leftBoundMercator;
    QGeoMapPolygonGeometryOpenGL m_geometry;
    QGeoMapPolylineGeometryOpenGL m_borderGeometry;
    QDeclarativePolygonMapItemPrivateOpenGL::RootNode *m_rootNode = nullptr;
    MapPolygonNodeGL *m_node = nullptr;
    MapPolylineNodeOpenGLExtruded *m_polylinenode = nullptr;

    void updateCirclePath(const QGeoCoordinate &center, qreal radius, const QGeoProjectionWebMercator &p);
    void markSourceDirty()
    {
        m_geometry.markSourceDirty();
        m_borderGeometry.markSourceDirty();
    }
};
    QScopedPointer<CircleDataCPU> m_dataCPU;
    QScopedPointer<CircleDataGL>  m_dataGL;
};

QT_END_NAMESPACE

#endif // QMAPCIRCLEOBJECT_P_P_H

