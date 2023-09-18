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

#ifndef QDECLARATIVECIRCLEMAPITEM_P_P_H
#define QDECLARATIVECIRCLEMAPITEM_P_P_H

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
#include <QtLocation/private/qdeclarativecirclemapitem_p.h>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QGeoMapCircleGeometry : public QGeoMapPolygonGeometry
{
public:
    QGeoMapCircleGeometry();

    void updateScreenPointsInvert(const QList<QDoubleVector2D> &circlePath, const QGeoMap &map);
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeCircleMapItemPrivate
{
public:
    static const int CircleSamples = 128; // ToDo: make this radius && ZL dependent?

    QDeclarativeCircleMapItemPrivate(QDeclarativeCircleMapItem &circle) : m_circle(circle)
    {

    }
    QDeclarativeCircleMapItemPrivate(QDeclarativeCircleMapItemPrivate &other) : m_circle(other.m_circle)
    {
    }

    virtual ~QDeclarativeCircleMapItemPrivate();
    virtual void onLinePropertiesChanged() = 0;
    virtual void markSourceDirtyAndUpdate() = 0;
    virtual void onMapSet() = 0;
    virtual void onGeoGeometryChanged() = 0;
    virtual void onItemGeometryChanged() = 0;
    virtual void updatePolish() = 0;
    virtual void afterViewportChanged() = 0;
    virtual QSGNode * updateMapItemPaintNode(QSGNode *oldNode, QQuickItem::UpdatePaintNodeData *data) = 0;
    virtual bool contains(const QPointF &point) const = 0;

    void updateCirclePath()
    {
        if (!m_circle.map() || m_circle.map()->geoProjection().projectionType() != QGeoProjection::ProjectionWebMercator)
            return;

        const QGeoProjectionWebMercator &p = static_cast<const QGeoProjectionWebMercator&>(m_circle.map()->geoProjection());
        QList<QGeoCoordinate> path;
        calculatePeripheralPoints(path, m_circle.center(), m_circle.radius(), CircleSamples, m_leftBound);
        m_circlePath.clear();
        for (const QGeoCoordinate &c : path)
            m_circlePath << p.geoToMapProjection(c);
    }

    static bool crossEarthPole(const QGeoCoordinate &center, qreal distance);

    static bool preserveCircleGeometry(QList<QDoubleVector2D> &path, const QGeoCoordinate &center,
                                qreal distance, const QGeoProjectionWebMercator &p);
    static void updateCirclePathForRendering(QList<QDoubleVector2D> &path, const QGeoCoordinate &center,
                                      qreal distance, const QGeoProjectionWebMercator &p);

    static void calculatePeripheralPoints(QList<QGeoCoordinate> &path, const QGeoCoordinate &center,
                                   qreal distance, int steps, QGeoCoordinate &leftBound);

    QDeclarativeCircleMapItem &m_circle;
    QList<QDoubleVector2D> m_circlePath;
    QGeoCoordinate m_leftBound;
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeCircleMapItemPrivateCPU: public QDeclarativeCircleMapItemPrivate
{
public:

    QDeclarativeCircleMapItemPrivateCPU(QDeclarativeCircleMapItem &circle) : QDeclarativeCircleMapItemPrivate(circle)
    {
    }

    QDeclarativeCircleMapItemPrivateCPU(QDeclarativeCircleMapItemPrivate &other)
    : QDeclarativeCircleMapItemPrivate(other)
    {
    }

    ~QDeclarativeCircleMapItemPrivateCPU() override;

    void onLinePropertiesChanged() override
    {
        // mark dirty just in case we're a width change
        markSourceDirtyAndUpdate();
    }
    void markSourceDirtyAndUpdate() override
    {
        // preserveGeometry is cleared in updateMapItemPaintNode
        m_geometry.markSourceDirty();
        m_borderGeometry.markSourceDirty();
        m_circle.polishAndUpdate();
    }
    void onMapSet() override
    {
        updateCirclePath();
        markSourceDirtyAndUpdate();
    }
    void onGeoGeometryChanged() override
    {
        updateCirclePath();
        markSourceDirtyAndUpdate();
    }
    void onItemGeometryChanged() override
    {
        onGeoGeometryChanged();
    }
    void afterViewportChanged() override
    {
        markSourceDirtyAndUpdate();
    }
    void updatePolish() override
    {
        if (!m_circle.m_circle.isValid()) {
            m_geometry.clear();
            m_borderGeometry.clear();
            m_circle.setWidth(0);
            m_circle.setHeight(0);
            return;
        }

        const QGeoProjectionWebMercator &p = static_cast<const QGeoProjectionWebMercator&>(m_circle.map()->geoProjection());
        QScopedValueRollback<bool> rollback(m_circle.m_updatingGeometry);
        m_circle.m_updatingGeometry = true;

        QList<QDoubleVector2D> circlePath = m_circlePath;

        int pathCount = circlePath.size();
        bool preserve = preserveCircleGeometry(circlePath, m_circle.m_circle.center(), m_circle.m_circle.radius(), p);
        // using leftBound_ instead of the analytically calculated circle_.boundingGeoRectangle().topLeft());
        // to fix QTBUG-62154
        m_geometry.setPreserveGeometry(true, m_leftBound); // to set the geoLeftBound_
        m_geometry.setPreserveGeometry(preserve, m_leftBound);

        bool invertedCircle = false;
        if (crossEarthPole(m_circle.m_circle.center(), m_circle.m_circle.radius()) && circlePath.size() == pathCount) {
            m_geometry.updateScreenPointsInvert(circlePath, *m_circle.map()); // invert fill area for really huge circles
            invertedCircle = true;
        } else {
            m_geometry.updateSourcePoints(*m_circle.map(), circlePath);
            m_geometry.updateScreenPoints(*m_circle.map(), m_circle.m_border.width());
        }

        m_borderGeometry.clear();
        QList<QGeoMapItemGeometry *> geoms;
        geoms << &m_geometry;

        if (m_circle.m_border.color() != Qt::transparent && m_circle.m_border.width() > 0) {
            QList<QDoubleVector2D> closedPath = circlePath;
            closedPath << closedPath.first();

            if (invertedCircle) {
                closedPath = m_circlePath;
                closedPath << closedPath.first();
                std::reverse(closedPath.begin(), closedPath.end());
            }

            m_borderGeometry.setPreserveGeometry(true, m_leftBound);
            m_borderGeometry.setPreserveGeometry(preserve, m_leftBound);

            // Use srcOrigin_ from fill geometry after clipping to ensure that translateToCommonOrigin won't fail.
            const QGeoCoordinate &geometryOrigin = m_geometry.origin();

            m_borderGeometry.srcPoints_.clear();
            m_borderGeometry.srcPointTypes_.clear();

            QDoubleVector2D borderLeftBoundWrapped;
            QList<QList<QDoubleVector2D > > clippedPaths = m_borderGeometry.clipPath(*m_circle.map(), closedPath, borderLeftBoundWrapped);
            if (clippedPaths.size()) {
                borderLeftBoundWrapped = p.geoToWrappedMapProjection(geometryOrigin);
                m_borderGeometry.pathToScreen(*m_circle.map(), clippedPaths, borderLeftBoundWrapped);
                m_borderGeometry.updateScreenPoints(*m_circle.map(), m_circle.m_border.width());
                geoms << &m_borderGeometry;
            } else {
                m_borderGeometry.clear();
            }
        }

        QRectF combined = QGeoMapItemGeometry::translateToCommonOrigin(geoms);

        if (invertedCircle || !preserve) {
            m_circle.setWidth(combined.width());
            m_circle.setHeight(combined.height());
        } else {
            m_circle.setWidth(combined.width() + 2 * m_circle.m_border.width()); // ToDo: Fix this!
            m_circle.setHeight(combined.height() + 2 * m_circle.m_border.width());
        }

        // No offsetting here, even in normal case, because first point offset is already translated
        m_circle.setPositionOnMap(m_geometry.origin(), m_geometry.firstPointOffset());
    }

    QSGNode * updateMapItemPaintNode(QSGNode *oldNode, QQuickItem::UpdatePaintNodeData *data) override
    {
        Q_UNUSED(data);
        if (!m_node || !oldNode) { // Apparently the QSG might delete the nodes if they become invisible
            m_node = new MapPolygonNode();
            if (oldNode) {
                delete oldNode;
                oldNode = nullptr;
            }
        } else {
            m_node = static_cast<MapPolygonNode *>(oldNode);
        }

        //TODO: update only material
        if (m_geometry.isScreenDirty() || m_borderGeometry.isScreenDirty() || m_circle.m_dirtyMaterial) {
            m_node->update(m_circle.m_color, m_circle.m_border.color(), &m_geometry, &m_borderGeometry);
            m_geometry.setPreserveGeometry(false);
            m_borderGeometry.setPreserveGeometry(false);
            m_geometry.markClean();
            m_borderGeometry.markClean();
            m_circle.m_dirtyMaterial = false;
        }
        return m_node;
    }
    bool contains(const QPointF &point) const override
    {
        return (m_geometry.contains(point) || m_borderGeometry.contains(point));
    }

    QGeoMapCircleGeometry m_geometry;
    QGeoMapPolylineGeometry m_borderGeometry;
    MapPolygonNode *m_node = nullptr;
};

#if QT_CONFIG(opengl)
class Q_LOCATION_PRIVATE_EXPORT QDeclarativeCircleMapItemPrivateOpenGL: public QDeclarativeCircleMapItemPrivate
{
public:
    QDeclarativeCircleMapItemPrivateOpenGL(QDeclarativeCircleMapItem &circle) : QDeclarativeCircleMapItemPrivate(circle)
    {
    }

    QDeclarativeCircleMapItemPrivateOpenGL(QDeclarativeCircleMapItemPrivate &other)
    : QDeclarativeCircleMapItemPrivate(other)
    {
    }

    ~QDeclarativeCircleMapItemPrivateOpenGL() override;

    void onLinePropertiesChanged() override
    {
        m_circle.m_dirtyMaterial = true;
        afterViewportChanged();
    }
    void markScreenDirtyAndUpdate()
    {
        // preserveGeometry is cleared in updateMapItemPaintNode
        m_geometry.markScreenDirty();
        m_borderGeometry.markScreenDirty();
        m_circle.polishAndUpdate();
    }
    virtual void markSourceDirtyAndUpdate() override
    {
        updateCirclePath();
        preserveGeometry();
        m_geometry.markSourceDirty();
        m_borderGeometry.markSourceDirty();
        m_circle.polishAndUpdate();
    }
    void preserveGeometry()
    {
        m_geometry.setPreserveGeometry(true, m_leftBound);
        m_borderGeometry.setPreserveGeometry(true, m_leftBound);
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
        onGeoGeometryChanged();
    }
    virtual void afterViewportChanged() override
    {
        preserveGeometry();
        markScreenDirtyAndUpdate();
    }
    virtual void updatePolish() override
    {
        if (m_circle.m_circle.isEmpty()) {
            m_geometry.clear();
            m_borderGeometry.clear();
            m_circle.setWidth(0);
            m_circle.setHeight(0);
            return;
        }

        QScopedValueRollback<bool> rollback(m_circle.m_updatingGeometry);
        m_circle.m_updatingGeometry = true;
        const qreal lineWidth = m_circle.m_border.width();
        const QColor &lineColor = m_circle.m_border.color();
        const QColor &fillColor = m_circle.color();
        if (fillColor.alpha() != 0) {
            m_geometry.updateSourcePoints(*m_circle.map(), m_circlePath);
            m_geometry.markScreenDirty();
            m_geometry.updateScreenPoints(*m_circle.map(), lineWidth, lineColor);
        } else {
            m_geometry.clearBounds();
        }

        QGeoMapItemGeometry * geom = &m_geometry;
        m_borderGeometry.clearScreen();
        if (lineColor.alpha() != 0 && lineWidth > 0) {
            m_borderGeometry.updateSourcePoints(*m_circle.map(), m_circle.m_circle);
            m_borderGeometry.markScreenDirty();
            m_borderGeometry.updateScreenPoints(*m_circle.map(), lineWidth);
            geom = &m_borderGeometry;
        }
        m_circle.setWidth(geom->sourceBoundingBox().width());
        m_circle.setHeight(geom->sourceBoundingBox().height());
        m_circle.setPosition(1.0 * geom->firstPointOffset() - QPointF(lineWidth * 0.5,lineWidth * 0.5));
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

        const QGeoMap *map = m_circle.map();
        const QMatrix4x4 &combinedMatrix = map->geoProjection().qsgTransform();
        const QDoubleVector3D &cameraCenter = map->geoProjection().centerMercator();

        if (m_borderGeometry.isScreenDirty()) {
            /* Do the border update first */
            m_polylinenode->update(m_circle.m_border.color(),
                                   float(m_circle.m_border.width()),
                                   &m_borderGeometry,
                                   combinedMatrix,
                                   cameraCenter,
                                   Qt::SquareCap,
                                   true,
                                   30); // No LOD for circles
            m_borderGeometry.setPreserveGeometry(false);
            m_borderGeometry.markClean();
        } else {
            m_polylinenode->setSubtreeBlocked(true);
        }
        if (m_geometry.isScreenDirty()) {
            m_node->update(m_circle.m_color,
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
        const qreal lineWidth = m_circle.m_border.width();
        const QColor &lineColor = m_circle.m_border.color();
        const QRectF &bounds = (lineColor.alpha() != 0 && lineWidth > 0) ? m_borderGeometry.sourceBoundingBox() : m_geometry.sourceBoundingBox();
        if (bounds.contains(point)) {
            QDeclarativeGeoMap *m = m_circle.quickMap();
            if (m) {
                const QGeoCoordinate crd = m->toCoordinate(m->mapFromItem(&m_circle, point));
                return  m_circle.m_circle.contains(crd) || m_borderGeometry.contains(m_circle.mapToItem(m_circle.quickMap(), point),
                                                                                     m_circle.border()->width(),
                                                                                     static_cast<const QGeoProjectionWebMercator&>(m_circle.map()->geoProjection()));
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

#endif // QDECLARATIVECIRCLEMAPITEM_P_P_H
