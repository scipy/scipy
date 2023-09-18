/****************************************************************************
**
** Copyright (C) 2020 Paolo Angelelli <paolo.angelelli@gmail.com>
** Copyright (C) 2020 The Qt Company Ltd.
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

#ifndef QDECLARATIVERECTANGLEMAPITEM_P_P_H
#define QDECLARATIVERECTANGLEMAPITEM_P_P_H

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
#include <QtLocation/private/qdeclarativepolygonmapitem_p_p.h>
#include <QtLocation/private/qdeclarativerectanglemapitem_p.h>
#include <QtPositioning/private/qwebmercator_p.h>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeRectangleMapItemPrivate
{
public:
    QDeclarativeRectangleMapItemPrivate(QDeclarativeRectangleMapItem &rect) : m_rect(rect)
    {

    }
    QDeclarativeRectangleMapItemPrivate(QDeclarativeRectangleMapItemPrivate &other) : m_rect(other.m_rect)
    {
    }

    virtual ~QDeclarativeRectangleMapItemPrivate();
    virtual void onLinePropertiesChanged() = 0;
    virtual void markSourceDirtyAndUpdate() = 0;
    virtual void onMapSet() = 0;
    virtual void onGeoGeometryChanged() = 0;
    virtual void onItemGeometryChanged() = 0;
    virtual void updatePolish() = 0;
    virtual void afterViewportChanged() = 0;
    virtual QSGNode * updateMapItemPaintNode(QSGNode *oldNode, QQuickItem::UpdatePaintNodeData *data) = 0;
    virtual bool contains(const QPointF &point) const = 0;

    QDeclarativeRectangleMapItem &m_rect;
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeRectangleMapItemPrivateCPU: public QDeclarativeRectangleMapItemPrivate
{
public:
    QDeclarativeRectangleMapItemPrivateCPU(QDeclarativeRectangleMapItem &rect) : QDeclarativeRectangleMapItemPrivate(rect)
    {
    }

    QDeclarativeRectangleMapItemPrivateCPU(QDeclarativeRectangleMapItemPrivate &other)
    : QDeclarativeRectangleMapItemPrivate(other)
    {
    }

    ~QDeclarativeRectangleMapItemPrivateCPU() override;

    void onLinePropertiesChanged() override
    {
        // mark dirty just in case we're a width change
        markSourceDirtyAndUpdate();
    }
    virtual void markSourceDirtyAndUpdate() override
    {
        m_geometry.markSourceDirty();
        m_borderGeometry.markSourceDirty();
        m_rect.polishAndUpdate();
    }
    virtual void onMapSet() override
    {
        markSourceDirtyAndUpdate();
    }
    virtual void onGeoGeometryChanged() override
    {
        markSourceDirtyAndUpdate();
    }
    virtual void onItemGeometryChanged() override
    {
        m_geometry.setPreserveGeometry(true, m_rect.m_rectangle.topLeft());
        m_borderGeometry.setPreserveGeometry(true, m_rect.m_rectangle.topLeft());
        markSourceDirtyAndUpdate();
    }
    virtual void afterViewportChanged() override
    {
        m_geometry.setPreserveGeometry(true, m_rect.m_rectangle.topLeft());
        m_borderGeometry.setPreserveGeometry(true, m_rect.m_rectangle.topLeft());
        markSourceDirtyAndUpdate();
    }
    virtual void updatePolish() override
    {
        if (!m_rect.topLeft().isValid() || !m_rect.bottomRight().isValid()) {
            m_geometry.clear();
            m_borderGeometry.clear();
            m_rect.setWidth(0);
            m_rect.setHeight(0);
            return;
        }

        const QGeoProjectionWebMercator &p = static_cast<const QGeoProjectionWebMercator&>(m_rect.map()->geoProjection());

        QScopedValueRollback<bool> rollback(m_rect.m_updatingGeometry);
        m_rect.m_updatingGeometry = true;

        const QList<QGeoCoordinate> perimeter = path(m_rect.m_rectangle);
        const QList<QDoubleVector2D> pathMercator_ = pathMercator(perimeter);
        m_geometry.setPreserveGeometry(true, m_rect.m_rectangle.topLeft());
        m_geometry.updateSourcePoints(*m_rect.map(), pathMercator_);
        m_geometry.updateScreenPoints(*m_rect.map(), m_rect.m_border.width());

        QList<QGeoMapItemGeometry *> geoms;
        geoms << &m_geometry;
        m_borderGeometry.clear();

        if (m_rect.m_border.color().alpha() != 0 && m_rect.m_border.width() > 0) {
            QList<QDoubleVector2D> closedPath = pathMercator_;
            closedPath << closedPath.first();

            m_borderGeometry.setPreserveGeometry(true, m_rect.m_rectangle.topLeft());
            const QGeoCoordinate &geometryOrigin = m_geometry.origin();

            m_borderGeometry.srcPoints_.clear();
            m_borderGeometry.srcPointTypes_.clear();

            QDoubleVector2D borderLeftBoundWrapped;
            QList<QList<QDoubleVector2D > > clippedPaths = m_borderGeometry.clipPath(*m_rect.map(), closedPath, borderLeftBoundWrapped);
            if (clippedPaths.size()) {
                borderLeftBoundWrapped = p.geoToWrappedMapProjection(geometryOrigin);
                m_borderGeometry.pathToScreen(*m_rect.map(), clippedPaths, borderLeftBoundWrapped);
                m_borderGeometry.updateScreenPoints(*m_rect.map(), m_rect.m_border.width());

                geoms << &m_borderGeometry;
            } else {
                m_borderGeometry.clear();
            }
        }

        QRectF combined = QGeoMapItemGeometry::translateToCommonOrigin(geoms);
        m_rect.setWidth(combined.width()  + 2 * m_rect.m_border.width()); // ToDo: fix this! 2 is incorrect
        m_rect.setHeight(combined.height()  + 2 * m_rect.m_border.width());

        m_rect.setPositionOnMap(m_geometry.origin(), m_geometry.firstPointOffset());
    }

    virtual QSGNode * updateMapItemPaintNode(QSGNode *oldNode, QQuickItem::UpdatePaintNodeData *data) override
    {
        Q_UNUSED(data);
        if (!m_node || !oldNode) {
            m_node = new MapPolygonNode();
            if (oldNode) {
                delete oldNode;
                oldNode = nullptr;
            }
        } else {
            m_node = static_cast<MapPolygonNode *>(oldNode);
        }

        //TODO: update only material
        if (m_geometry.isScreenDirty() || m_borderGeometry.isScreenDirty() || m_rect.m_dirtyMaterial) {
            m_node->update(m_rect.m_color, m_rect.m_border.color(), &m_geometry, &m_borderGeometry);
            m_geometry.setPreserveGeometry(false);
            m_borderGeometry.setPreserveGeometry(false);
            m_geometry.markClean();
            m_borderGeometry.markClean();
            m_rect.m_dirtyMaterial = false;
        }
        return m_node;
    }
    virtual bool contains(const QPointF &point) const override
    {
        return (m_geometry.contains(point) || m_borderGeometry.contains(point));
    }

    static QList<QGeoCoordinate> path(const QGeoRectangle &rect)
    {
        QList<QGeoCoordinate> res;
        res << rect.topLeft();
        res << QGeoCoordinate(rect.topLeft().latitude(), rect.bottomRight().longitude());
        res << rect.bottomRight();
        res << QGeoCoordinate(rect.bottomRight().latitude(), rect.topLeft().longitude());
        return res;
    }

    static QList<QGeoCoordinate> perimeter(const QGeoRectangle &rect)
    {
        QList<QGeoCoordinate> res;
        res << rect.topLeft();
        res << QGeoCoordinate(rect.topLeft().latitude(), rect.bottomRight().longitude());
        res << rect.bottomRight();
        res << QGeoCoordinate(rect.bottomRight().latitude(), rect.topLeft().longitude());
        res << res.first();
        return res;
    }

    static QList<QDoubleVector2D> pathMercator(const QList<QGeoCoordinate> &p)
    {
        QList<QDoubleVector2D> res;
        for (const auto &c: p)
            res << QWebMercator::coordToMercator(c);
        return res;
    }

    QGeoMapPolygonGeometry m_geometry;
    QGeoMapPolylineGeometry m_borderGeometry;
    MapPolygonNode *m_node = nullptr;
};

#if QT_CONFIG(opengl)
class Q_LOCATION_PRIVATE_EXPORT QDeclarativeRectangleMapItemPrivateOpenGL: public QDeclarativeRectangleMapItemPrivate
{
public:
    QDeclarativeRectangleMapItemPrivateOpenGL(QDeclarativeRectangleMapItem &rect) : QDeclarativeRectangleMapItemPrivate(rect)
    {
    }

    QDeclarativeRectangleMapItemPrivateOpenGL(QDeclarativeRectangleMapItemPrivate &other)
    : QDeclarativeRectangleMapItemPrivate(other)
    {
    }

    ~QDeclarativeRectangleMapItemPrivateOpenGL() override;

    void markScreenDirtyAndUpdate()
    {
        // preserveGeometry is cleared in updateMapItemPaintNode
        m_geometry.markScreenDirty();
        m_borderGeometry.markScreenDirty();
        m_rect.polishAndUpdate();
    }
    void onLinePropertiesChanged() override
    {
        m_rect.m_dirtyMaterial = true;
        afterViewportChanged();
    }
    virtual void markSourceDirtyAndUpdate() override
    {
        m_geometry.markSourceDirty();
        m_borderGeometry.markSourceDirty();
        m_rect.polishAndUpdate();
    }
    void preserveGeometry()
    {
        m_geometry.setPreserveGeometry(true, m_rect.m_rectangle.topLeft());
        m_borderGeometry.setPreserveGeometry(true, m_rect.m_rectangle.topLeft());
    }
    virtual void onMapSet() override
    {
        markSourceDirtyAndUpdate();
    }
    virtual void onGeoGeometryChanged() override
    {
        preserveGeometry();
        markSourceDirtyAndUpdate();
    }
    virtual void onItemGeometryChanged() override
    {
        onGeoGeometryChanged();
    }
    virtual void afterViewportChanged() override
    {
        preserveGeometry();
        markScreenDirtyAndUpdate();
    }
    virtual void updatePolish() override
    {
        if (!m_rect.topLeft().isValid() || !m_rect.bottomRight().isValid()) {
            m_geometry.clear();
            m_borderGeometry.clear();
            m_rect.setWidth(0);
            m_rect.setHeight(0);
            return;
        }

        QScopedValueRollback<bool> rollback(m_rect.m_updatingGeometry);
        m_rect.m_updatingGeometry = true;
        const qreal lineWidth = m_rect.m_border.width();
        const QColor &lineColor = m_rect.m_border.color();
        const QColor &fillColor = m_rect.color();
        if (fillColor.alpha() != 0) {
            m_geometry.updateSourcePoints(*m_rect.map(), m_rect.m_rectangle);
            m_geometry.markScreenDirty();
            m_geometry.updateScreenPoints(*m_rect.map(), lineWidth, lineColor);
        } else {
            m_geometry.clearBounds();
        }

        QGeoMapItemGeometry * geom = &m_geometry;
        m_borderGeometry.clearScreen();
        if (lineColor.alpha() != 0 && lineWidth > 0) {
            m_borderGeometry.updateSourcePoints(*m_rect.map(), m_rect.m_rectangle);
            m_borderGeometry.markScreenDirty();
            m_borderGeometry.updateScreenPoints(*m_rect.map(), lineWidth);
            geom = &m_borderGeometry;
        }
        m_rect.setWidth(geom->sourceBoundingBox().width());
        m_rect.setHeight(geom->sourceBoundingBox().height());
        m_rect.setPosition(1.0 * geom->firstPointOffset() - QPointF(lineWidth * 0.5,lineWidth * 0.5));
    }

    virtual QSGNode * updateMapItemPaintNode(QSGNode *oldNode, QQuickItem::UpdatePaintNodeData *data) override
    {
        Q_UNUSED(data);

        if (!m_rootNode || !oldNode) {
            m_rootNode = new QDeclarativePolygonMapItemPrivateOpenGL::RootNode();
            m_node = new MapPolygonNodeGL();
            m_rootNode->appendChildNode(m_node);
            m_polylinenode = new MapPolylineNodeOpenGLExtruded();
            m_rootNode->appendChildNode(m_polylinenode);
            m_rootNode->markDirty(QSGNode::DirtyNodeAdded);
            if (oldNode)
                delete oldNode;
        } else {
            m_rootNode = static_cast<QDeclarativePolygonMapItemPrivateOpenGL::RootNode *>(oldNode);
        }

        const QGeoMap *map = m_rect.map();
        const QMatrix4x4 &combinedMatrix = map->geoProjection().qsgTransform();
        const QDoubleVector3D &cameraCenter = map->geoProjection().centerMercator();

        if (m_borderGeometry.isScreenDirty()) {
            /* Do the border update first */
            m_polylinenode->update(m_rect.m_border.color(),
                                   float(m_rect.m_border.width()),
                                   &m_borderGeometry,
                                   combinedMatrix,
                                   cameraCenter,
                                   Qt::SquareCap,
                                   true,
                                   30); // No LOD for rectangles
            m_borderGeometry.setPreserveGeometry(false);
            m_borderGeometry.markClean();
        } else {
            m_polylinenode->setSubtreeBlocked(true);
        }
        if (m_geometry.isScreenDirty()) {
            m_node->update(m_rect.m_color,
                         &m_geometry,
                         combinedMatrix,
                         cameraCenter);
            m_geometry.setPreserveGeometry(false);
            m_geometry.markClean();
        } else {
            m_node->setSubtreeBlocked(true);
        }

        m_rootNode->setSubtreeBlocked(false);
        return m_rootNode;
    }
    virtual bool contains(const QPointF &point) const override
    {
        const qreal lineWidth = m_rect.m_border.width();
        const QColor &lineColor = m_rect.m_border.color();
        const QRectF &bounds = (lineColor.alpha() != 0 && lineWidth > 0) ? m_borderGeometry.sourceBoundingBox() : m_geometry.sourceBoundingBox();
        if (bounds.contains(point)) {
            QDeclarativeGeoMap *m = m_rect.quickMap();
            if (m) {
                const QGeoCoordinate crd = m->toCoordinate(m->mapFromItem(&m_rect, point));
                return  m_rect.m_rectangle.contains(crd) || m_borderGeometry.contains(m_rect.mapToItem(m_rect.quickMap(), point),
                                                                                      m_rect.border()->width(),
                                                                                      static_cast<const QGeoProjectionWebMercator&>(m_rect.map()->geoProjection()));
            } else {
                return  true;
            }
        }
        return false;
    }

    QGeoMapPolygonGeometryOpenGL m_geometry;
    QGeoMapPolylineGeometryOpenGL m_borderGeometry;
    QDeclarativePolygonMapItemPrivateOpenGL::RootNode *m_rootNode = nullptr;
    MapPolygonNodeGL *m_node = nullptr;
    MapPolylineNodeOpenGLExtruded *m_polylinenode = nullptr;
};
#endif // QT_CONFIG(opengl)

QT_END_NAMESPACE

#endif // QDECLARATIVERECTANGLEMAPITEM_P_P_H

