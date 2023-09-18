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

#ifndef QMAPPOLYGONOBJECTQSG_P_P_H
#define QMAPPOLYGONOBJECTQSG_P_P_H

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
#include <QtLocation/private/qmappolygonobject_p.h>
#include <QtLocation/private/qmappolygonobject_p_p.h>
#include <QtLocation/private/qqsgmapobject_p.h>
#include <QtLocation/private/qdeclarativepolygonmapitem_p_p.h>
#include <QtCore/qscopedvaluerollback.h>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QMapPolygonObjectPrivateQSG : public QMapPolygonObjectPrivateDefault, public QQSGMapObject
{
public:
    QMapPolygonObjectPrivateQSG(QGeoMapObject *q);
    QMapPolygonObjectPrivateQSG(const QMapPolygonObjectPrivate &other);
    ~QMapPolygonObjectPrivateQSG() override;

    QList<QDoubleVector2D> projectPath();

    // QQSGMapObject
    void markSourceDirty();
    void updateGeometry() override;
    QSGNode *updateMapObjectNode(QSGNode *oldNode,
                                 VisibleNode **visibleNode,
                                 QSGNode *root,
                                 QQuickWindow *window) override;

    // QGeoMapPolygonPrivate interface
    void setPath(const QList<QGeoCoordinate> &path) override;
    void setFillColor(const QColor &color) override;
    void setBorderColor(const QColor &color) override;
    void setBorderWidth(qreal width) override;

    // QGeoMapObjectPrivate
    QGeoMapObjectPrivate *clone() override;
    virtual void setGeoShape(const QGeoShape &shape) override;

    // Data Members
    QDoubleVector2D m_leftBoundMercator;
    QGeoMapPolygonGeometryOpenGL m_geometry;
    QGeoMapPolylineGeometryOpenGL m_borderGeometry;
    QDeclarativePolygonMapItemPrivateOpenGL::RootNode *m_rootNode = nullptr;
    MapPolygonNodeGL *m_node = nullptr;
    MapPolylineNodeOpenGLExtruded *m_polylinenode = nullptr;
};

QT_END_NAMESPACE

#endif // QMAPPOLYGONOBJECTQSG_P_P_H
