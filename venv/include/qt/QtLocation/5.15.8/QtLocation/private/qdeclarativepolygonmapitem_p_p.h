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

#ifndef QDECLARATIVEPOLYGONMAPITEM_P_P_H
#define QDECLARATIVEPOLYGONMAPITEM_P_P_H

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
#include <QtLocation/private/qgeomapitemgeometry_p.h>
#include <QtLocation/private/qdeclarativegeomapitembase_p.h>
#include <QtLocation/private/qdeclarativepolylinemapitem_p.h>
#include <QtLocation/private/qdeclarativegeomapitemutils_p.h>
#include <QtLocation/private/qdeclarativepolygonmapitem_p.h>
#include <QtLocation/private/qdeclarativepolylinemapitem_p_p.h>
#include <QSGGeometryNode>
#include <QSGFlatColorMaterial>
#include <QtPositioning/QGeoPath>
#include <QtPositioning/QGeoRectangle>
#include <QtPositioning/QGeoPolygon>
#include <QtPositioning/private/qdoublevector2d_p.h>
#include <QSGFlatColorMaterial>
#include <QSGSimpleMaterial>
#include <QtGui/QMatrix4x4>
#include <QColor>
#include <QList>
#include <QVector>
#include <QtCore/QScopedValueRollback>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QGeoMapPolygonGeometry : public QGeoMapItemGeometry
{
public:
    QGeoMapPolygonGeometry();

    inline void setAssumeSimple(bool value) { assumeSimple_ = value; }

    void updateSourcePoints(const QGeoMap &map,
                            const QList<QDoubleVector2D> &path);

    void updateScreenPoints(const QGeoMap &map, qreal strokeWidth = 0.0);

protected:
    QPainterPath srcPath_;
    bool assumeSimple_;
};

#if QT_CONFIG(opengl)
class Q_LOCATION_PRIVATE_EXPORT QGeoMapPolygonGeometryOpenGL : public QGeoMapItemGeometry
{
public:
    typedef struct {
        QList<QDoubleVector2D> wrappedBboxes;
    } WrappedPolygon;
    QGeoMapPolygonGeometryOpenGL();
    ~QGeoMapPolygonGeometryOpenGL() override {}

    // Temporary method for compatibility in MapCircleObject. Remove when MapObjects are ported.
    void updateSourcePoints(const QGeoMap &map,
                            const QList<QDoubleVector2D> &path);

    void updateSourcePoints(const QGeoMap &map,
                            const QList<QGeoCoordinate> &perimeter);

    void updateSourcePoints(const QGeoMap &map,
                            const QGeoPolygon &poly);

    void updateSourcePoints(const QGeoMap &map,
                            const QGeoRectangle &rect);

    void updateScreenPoints(const QGeoMap &map, qreal strokeWidth = 0.0, const QColor &strokeColor = Qt::transparent);
    void updateQuickGeometry(const QGeoProjectionWebMercator &p, qreal strokeWidth = 0.0);

    void allocateAndFillPolygon(QSGGeometry *geom) const
    {


        const QVector<QDeclarativeGeoMapItemUtils::vec2> &vx = m_screenVertices;
        const QVector<quint32> &ix = m_screenIndices;

        geom->allocate(vx.size(), ix.size());
        if (geom->indexType() == QSGGeometry::UnsignedShortType) {
            quint16 *its = geom->indexDataAsUShort();
            for (int i = 0; i < ix.size(); ++i)
                its[i] = ix[i];
        } else if (geom->indexType() == QSGGeometry::UnsignedIntType) {
            quint32 *its = geom->indexDataAsUInt();
            for (int i = 0; i < ix.size(); ++i)
                its[i] = ix[i];
        }

        QSGGeometry::Point2D *pts = geom->vertexDataAsPoint2D();
        for (int i = 0; i < vx.size(); ++i)
            pts[i].set(vx[i].x, vx[i].y);
    }

    QVector<QDeclarativeGeoMapItemUtils::vec2> m_screenVertices;
    QVector<quint32> m_screenIndices;
    QDoubleVector2D m_bboxLeftBoundWrapped;
    QVector<WrappedPolygon> m_wrappedPolygons;
    int m_wrapOffset;
};

class Q_LOCATION_PRIVATE_EXPORT MapPolygonShader : public QSGMaterialShader
{
public:
    MapPolygonShader();

    const char *vertexShader() const override {
        return
        "attribute highp vec4 vertex;               \n"
        "uniform highp mat4 qt_Matrix;              \n"
        "uniform highp mat4 mapProjection;          \n"
        "uniform highp vec3 center;                 \n"
        "uniform highp vec3 center_lowpart;         \n"
        "uniform lowp float wrapOffset;             \n"
        "vec4 wrapped(in vec4 v) { return vec4(v.x + wrapOffset, v.y, 0.0, 1.0); }\n"
        "void main() {                              \n"
        "    vec4 vtx = wrapped(vertex) - vec4(center, 0.0);   \n"
        "    vtx = vtx - vec4(center_lowpart, 0.0);   \n"
        "    gl_Position = qt_Matrix * mapProjection * vtx;      \n"
        "}";
    }

    const char *fragmentShader() const override {
        return
        "uniform lowp vec4 color;                   \n"
        "void main() {                              \n"
        "    gl_FragColor = color;                  \n"
        "}";
    }

    void updateState(const RenderState &state, QSGMaterial *newEffect, QSGMaterial *oldEffect) override;
    char const *const *attributeNames() const override
    {
        static char const *const attr[] = { "vertex", nullptr };
        return attr;
    }

private:
    void initialize() override
    {
        m_matrix_id = program()->uniformLocation("qt_Matrix");
        m_color_id = program()->uniformLocation("color");
        m_mapProjection_id = program()->uniformLocation("mapProjection");
        m_center_id = program()->uniformLocation("center");
        m_center_lowpart_id = program()->uniformLocation("center_lowpart");
        m_wrapOffset_id = program()->uniformLocation("wrapOffset");
    }
    int m_center_id;
    int m_center_lowpart_id;
    int m_mapProjection_id;
    int m_matrix_id;
    int m_color_id;
    int m_wrapOffset_id;
};
#endif // QT_CONFIG(opengl)

class Q_LOCATION_PRIVATE_EXPORT MapPolygonMaterial : public QSGFlatColorMaterial
{
public:
    MapPolygonMaterial()
        : QSGFlatColorMaterial()
    {
        // Passing RequiresFullMatrix is essential in order to prevent the
        // batch renderer from baking in simple, translate-only transforms into
        // the vertex data. The shader will rely on the fact that
        // vertexCoord.xy is the Shape-space coordinate and so no modifications
        // are welcome.
        setFlag(Blending | RequiresFullMatrix | CustomCompileStep);
    }

    QSGMaterialShader *createShader() const override;

    void setGeoProjection(const QMatrix4x4 &p)
    {
        m_geoProjection = p;
    }

    QMatrix4x4 geoProjection() const
    {
        return m_geoProjection;
    }

    void setCenter(const QDoubleVector3D &c)
    {
        m_center = c;
    }

    QDoubleVector3D center() const
    {
        return m_center;
    }

    int wrapOffset() const
    {
        return m_wrapOffset;
    }

    void setWrapOffset(int wrapOffset)
    {
        m_wrapOffset = wrapOffset;
    }

    int compare(const QSGMaterial *other) const override;
    QSGMaterialType *type() const override;

protected:
    QMatrix4x4 m_geoProjection;
    QDoubleVector3D m_center;
    int m_wrapOffset = 0;
};

class Q_LOCATION_PRIVATE_EXPORT MapPolygonNode : public MapItemGeometryNode
{

public:
    MapPolygonNode();
    ~MapPolygonNode() override;

    void update(const QColor &fillColor, const QColor &borderColor,
                const QGeoMapItemGeometry *fillShape,
                const QGeoMapItemGeometry *borderShape);
private:
    QSGFlatColorMaterial fill_material_;
    MapPolylineNode *border_;
    QSGGeometry geometry_;
};

#if QT_CONFIG(opengl)
class Q_LOCATION_PRIVATE_EXPORT MapPolygonNodeGL : public MapItemGeometryNode
{

public:
    MapPolygonNodeGL();
    ~MapPolygonNodeGL() override;

    void update(const QColor &fillColor,
                const QGeoMapPolygonGeometryOpenGL *fillShape,
                const QMatrix4x4 &geoProjection,
                const QDoubleVector3D &center);

    MapPolygonMaterial fill_material_;
    QSGGeometry geometry_;
};
#endif // QT_CONFIG(opengl)

class Q_LOCATION_PRIVATE_EXPORT QDeclarativePolygonMapItemPrivate
{
public:
    QDeclarativePolygonMapItemPrivate(QDeclarativePolygonMapItem &polygon) : m_poly(polygon)
    {

    }
    QDeclarativePolygonMapItemPrivate(QDeclarativePolygonMapItemPrivate &other) : m_poly(other.m_poly)
    {
    }

    virtual ~QDeclarativePolygonMapItemPrivate();
    virtual void onLinePropertiesChanged() = 0;
    virtual void markSourceDirtyAndUpdate() = 0;
    virtual void onMapSet() = 0;
    virtual void onGeoGeometryChanged() = 0;
    virtual void onGeoGeometryUpdated() = 0;
    virtual void onItemGeometryChanged() = 0;
    virtual void updatePolish() = 0;
    virtual void afterViewportChanged() = 0;
    virtual QSGNode * updateMapItemPaintNode(QSGNode *oldNode, QQuickItem::UpdatePaintNodeData *data) = 0;
    virtual bool contains(const QPointF &point) const = 0;

    QDeclarativePolygonMapItem &m_poly;
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativePolygonMapItemPrivateCPU: public QDeclarativePolygonMapItemPrivate
{
public:
    QDeclarativePolygonMapItemPrivateCPU(QDeclarativePolygonMapItem &polygon) : QDeclarativePolygonMapItemPrivate(polygon)
    {
    }

    QDeclarativePolygonMapItemPrivateCPU(QDeclarativePolygonMapItemPrivate &other)
    : QDeclarativePolygonMapItemPrivate(other)
    {
    }

    ~QDeclarativePolygonMapItemPrivateCPU() override;
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
        m_poly.polishAndUpdate();
    }
    void regenerateCache()
    {
        if (!m_poly.map() || m_poly.map()->geoProjection().projectionType() != QGeoProjection::ProjectionWebMercator)
            return;
        const QGeoProjectionWebMercator &p = static_cast<const QGeoProjectionWebMercator&>(m_poly.map()->geoProjection());
        m_geopathProjected.clear();
        m_geopathProjected.reserve(m_poly.m_geopoly.size());
        for (const QGeoCoordinate &c : m_poly.m_geopoly.path())
            m_geopathProjected << p.geoToMapProjection(c);
    }
    void updateCache()
    {
        if (!m_poly.map() || m_poly.map()->geoProjection().projectionType() != QGeoProjection::ProjectionWebMercator)
            return;
        const QGeoProjectionWebMercator &p = static_cast<const QGeoProjectionWebMercator&>(m_poly.map()->geoProjection());
        m_geopathProjected << p.geoToMapProjection(m_poly.m_geopoly.path().last());
    }
    void preserveGeometry()
    {
        m_geometry.setPreserveGeometry(true, m_poly.m_geopoly.boundingGeoRectangle().topLeft());
        m_borderGeometry.setPreserveGeometry(true, m_poly.m_geopoly.boundingGeoRectangle().topLeft());
    }
    void afterViewportChanged() override
    {
        // preserveGeometry is cleared in updateMapItemPaintNode
        preserveGeometry();
        markSourceDirtyAndUpdate();
    }
    void onMapSet() override
    {
        regenerateCache();
        markSourceDirtyAndUpdate();
    }
    void onGeoGeometryChanged() override
    {
        regenerateCache();
        preserveGeometry();
        markSourceDirtyAndUpdate();
    }
    void onGeoGeometryUpdated() override
    {
        updateCache();
        preserveGeometry();
        markSourceDirtyAndUpdate();
    }
    void onItemGeometryChanged() override
    {
        onGeoGeometryChanged();
    }
    void updatePolish() override
    {
        if (m_poly.m_geopoly.path().length() == 0) { // Possibly cleared
            m_geometry.clear();
            m_borderGeometry.clear();
            m_poly.setWidth(0);
            m_poly.setHeight(0);
            return;
        }
        const QGeoMap *map = m_poly.map();
        const qreal borderWidth = m_poly.m_border.width();
        const QGeoProjectionWebMercator &p = static_cast<const QGeoProjectionWebMercator&>(map->geoProjection());
        QScopedValueRollback<bool> rollback(m_poly.m_updatingGeometry);
        m_poly.m_updatingGeometry = true;

        m_geometry.updateSourcePoints(*map, m_geopathProjected);
        m_geometry.updateScreenPoints(*map, borderWidth);

        QList<QGeoMapItemGeometry *> geoms;
        geoms << &m_geometry;
        m_borderGeometry.clear();

        if (m_poly.m_border.color().alpha() != 0 && borderWidth > 0) {
            QList<QDoubleVector2D> closedPath = m_geopathProjected;
            closedPath << closedPath.first();

            m_borderGeometry.setPreserveGeometry(true, m_poly.m_geopoly.boundingGeoRectangle().topLeft());

            const QGeoCoordinate &geometryOrigin = m_geometry.origin();

            m_borderGeometry.srcPoints_.clear();
            m_borderGeometry.srcPointTypes_.clear();

            QDoubleVector2D borderLeftBoundWrapped;
            QList<QList<QDoubleVector2D > > clippedPaths = m_borderGeometry.clipPath(*map, closedPath, borderLeftBoundWrapped);
            if (clippedPaths.size()) {
                borderLeftBoundWrapped = p.geoToWrappedMapProjection(geometryOrigin);
                m_borderGeometry.pathToScreen(*map, clippedPaths, borderLeftBoundWrapped);
                m_borderGeometry.updateScreenPoints(*map, borderWidth);

                geoms <<  &m_borderGeometry;
            } else {
                m_borderGeometry.clear();
            }
        }

        QRectF combined = QGeoMapItemGeometry::translateToCommonOrigin(geoms);
        m_poly.setWidth(combined.width() + 2 * borderWidth);
        m_poly.setHeight(combined.height() + 2 * borderWidth);

        m_poly.setPositionOnMap(m_geometry.origin(), -1 * m_geometry.sourceBoundingBox().topLeft()
                                                + QPointF(borderWidth, borderWidth));
    }
    QSGNode *updateMapItemPaintNode(QSGNode *oldNode, QQuickItem::UpdatePaintNodeData *data) override
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
        if (m_geometry.isScreenDirty()
                || m_borderGeometry.isScreenDirty()
                || m_poly.m_dirtyMaterial
                || !oldNode) {
            m_node->update(m_poly.m_color,
                         m_poly.m_border.color(),
                         &m_geometry,
                         &m_borderGeometry);
            m_geometry.setPreserveGeometry(false);
            m_borderGeometry.setPreserveGeometry(false);
            m_geometry.markClean();
            m_borderGeometry.markClean();
            m_poly.m_dirtyMaterial = false;
        }
        return m_node;
    }
    bool contains(const QPointF &point) const override
    {
        return (m_geometry.contains(point) || m_borderGeometry.contains(point));
    }

    QList<QDoubleVector2D> m_geopathProjected;
    QGeoMapPolygonGeometry m_geometry;
    QGeoMapPolylineGeometry m_borderGeometry;
    MapPolygonNode *m_node = nullptr;
};

#if QT_CONFIG(opengl)
class Q_LOCATION_PRIVATE_EXPORT QDeclarativePolygonMapItemPrivateOpenGL: public QDeclarativePolygonMapItemPrivate
{
public:
    struct RootNode : public QSGNode /*QSGTransformNode*/, public VisibleNode
    {
        RootNode() { }

        bool isSubtreeBlocked() const override
        {
            return subtreeBlocked();
        }
    };

    QDeclarativePolygonMapItemPrivateOpenGL(QDeclarativePolygonMapItem &polygon) : QDeclarativePolygonMapItemPrivate(polygon)
    {
    }

    QDeclarativePolygonMapItemPrivateOpenGL(QDeclarativePolygonMapItemPrivate &other)
    : QDeclarativePolygonMapItemPrivate(other)
    {
    }

    ~QDeclarativePolygonMapItemPrivateOpenGL() override;

    void markScreenDirtyAndUpdate()
    {
        // preserveGeometry is cleared in updateMapItemPaintNode
        m_geometry.markScreenDirty();
        m_borderGeometry.markScreenDirty();
        m_poly.polishAndUpdate();
    }
    void onLinePropertiesChanged() override
    {
        m_poly.m_dirtyMaterial = true;
        afterViewportChanged();
    }
    void markSourceDirtyAndUpdate() override
    {
        // preserveGeometry is cleared in updateMapItemPaintNode
        m_geometry.markSourceDirty();
        m_borderGeometry.markSourceDirty();
        m_poly.polishAndUpdate();
    }
    void preserveGeometry()
    {
        m_geometry.setPreserveGeometry(true, m_poly.m_geopoly.boundingGeoRectangle().topLeft());
        m_borderGeometry.setPreserveGeometry(true, m_poly.m_geopoly.boundingGeoRectangle().topLeft());
    }
    void afterViewportChanged() override // This is called when the camera changes, or visibleArea changes.
    {
        // preserveGeometry is cleared in updateMapItemPaintNode
        preserveGeometry();
        markScreenDirtyAndUpdate();
    }
    void onMapSet() override
    {
        markSourceDirtyAndUpdate();
    }
    void onGeoGeometryChanged() override
    {
        preserveGeometry();
        markSourceDirtyAndUpdate();
    }
    void onGeoGeometryUpdated() override
    {
        preserveGeometry();
        markSourceDirtyAndUpdate();
    }
    void onItemGeometryChanged() override
    {
        onGeoGeometryChanged();
    }
    void updatePolish() override
    {
        if (m_poly.m_geopoly.path().length() == 0) { // Possibly cleared
            m_geometry.clear();
            m_borderGeometry.clear();
            m_poly.setWidth(0);
            m_poly.setHeight(0);
            return;
        }

        QScopedValueRollback<bool> rollback(m_poly.m_updatingGeometry);
        m_poly.m_updatingGeometry = true;
        const qreal lineWidth = m_poly.m_border.width();
        const QColor &lineColor = m_poly.m_border.color();
        const QColor &fillColor = m_poly.color();
        if (fillColor.alpha() != 0) {
            m_geometry.updateSourcePoints(*m_poly.map(), m_poly.m_geopoly);
            m_geometry.markScreenDirty();
            m_geometry.updateScreenPoints(*m_poly.map(), lineWidth, lineColor);
        } else {
            m_geometry.clearBounds();
        }

        QGeoMapItemGeometry * geom = &m_geometry;
        m_borderGeometry.clearScreen();
        if (lineColor.alpha() != 0 && lineWidth > 0) {
            m_borderGeometry.updateSourcePoints(*m_poly.map(), m_poly.m_geopoly);
            m_borderGeometry.markScreenDirty();
            m_borderGeometry.updateScreenPoints(*m_poly.map(), lineWidth);
            geom = &m_borderGeometry;
        }
        m_poly.setWidth(geom->sourceBoundingBox().width());
        m_poly.setHeight(geom->sourceBoundingBox().height());
        m_poly.setPosition(1.0 * geom->firstPointOffset() - QPointF(lineWidth * 0.5,lineWidth * 0.5));
    }
    QSGNode * updateMapItemPaintNode(QSGNode *oldNode, QQuickItem::UpdatePaintNodeData *data) override
    {
        Q_UNUSED(data);

        if (!m_rootNode || !oldNode) {
            m_rootNode = new RootNode();
            m_node = new MapPolygonNodeGL();
            m_rootNode->appendChildNode(m_node);
            m_polylinenode = new MapPolylineNodeOpenGLExtruded();
            m_rootNode->appendChildNode(m_polylinenode);
            m_rootNode->markDirty(QSGNode::DirtyNodeAdded);
            if (oldNode)
                delete oldNode;
        } else {
            m_rootNode = static_cast<RootNode *>(oldNode);
        }

        const QGeoMap *map = m_poly.map();
        const QMatrix4x4 &combinedMatrix = map->geoProjection().qsgTransform();
        const QDoubleVector3D &cameraCenter = map->geoProjection().centerMercator();

        if (m_borderGeometry.isScreenDirty()) {
            /* Do the border update first */
            m_polylinenode->update(m_poly.m_border.color(),
                                   float(m_poly.m_border.width()),
                                   &m_borderGeometry,
                                   combinedMatrix,
                                   cameraCenter,
                                   Qt::SquareCap,
                                   true,
                                   30); // No LOD for polygons just yet.
                                        // First figure out what to do with holes.
            m_borderGeometry.setPreserveGeometry(false);
            m_borderGeometry.markClean();
        } else {
            m_polylinenode->setSubtreeBlocked(true);
        }
        if (m_geometry.isScreenDirty()) {
            m_node->update(m_poly.m_color,
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
    bool contains(const QPointF &point) const override
    {
        const qreal lineWidth = m_poly.m_border.width();
        const QColor &lineColor = m_poly.m_border.color();
        const QRectF &bounds = (lineColor.alpha() != 0 && lineWidth > 0) ? m_borderGeometry.sourceBoundingBox() : m_geometry.sourceBoundingBox();
        if (bounds.contains(point)) {
            QDeclarativeGeoMap *m = m_poly.quickMap();
            if (m) {
                const QGeoCoordinate crd = m->toCoordinate(m->mapFromItem(&m_poly, point));
                return  m_poly.m_geopoly.contains(crd) || m_borderGeometry.contains(m_poly.mapToItem(m_poly.quickMap(), point),
                                                                              m_poly.border()->width(),
                                                                              static_cast<const QGeoProjectionWebMercator&>(m_poly.map()->geoProjection()));
            } else {
                return  true;
            }
        }
        return false;
    }

    QGeoMapPolygonGeometryOpenGL m_geometry;
    QGeoMapPolylineGeometryOpenGL m_borderGeometry;
    RootNode *m_rootNode = nullptr;
    MapPolygonNodeGL *m_node = nullptr;
    MapPolylineNodeOpenGLExtruded *m_polylinenode = nullptr;
};
#endif // QT_CONFIG(opengl)

QT_END_NAMESPACE

#endif // QDECLARATIVEPOLYGONMAPITEM_P_P_H
