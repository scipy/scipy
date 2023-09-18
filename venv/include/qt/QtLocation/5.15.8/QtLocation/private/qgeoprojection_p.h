/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
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

#ifndef QGEOPROJECTION_H
#define QGEOPROJECTION_H

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
#include <QtLocation/private/qgeocameradata_p.h>
#include <QtPositioning/private/qdoublematrix4x4_p.h>
#include <QtPositioning/QGeoShape>
#include <QMatrix4x4>
#include <QTransform>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QGeoProjection
{
public:
    enum ProjectionGroup {
        ProjectionOther,
        ProjectionCylindrical,
        ProjectionPseudocylindrical,
        ProjectionAzimuthal,
        ProjectionPseudoazimuthal,
        ProjectionConic,
        ProjectionPseudoconic
        //Polyhedral
        //Retroazimuthal
    };

    enum Datum {
        DatumUnknown,
        DatumWGS84,
        DatumSphere
    };

    enum ProjectionType {
        ProjectionUnknown,
        ProjectionGeneralPerspective,
        ProjectionWebMercator
    };

    QGeoProjection();
    virtual ~QGeoProjection();

    virtual void setVisibleArea(const QRectF &visibleArea) = 0;
    virtual void setViewportSize(const QSize &size) = 0;
    virtual void setCameraData(const QGeoCameraData &cameraData, bool force = true) = 0;
    virtual QGeoCameraData cameraData() const = 0;

    // returns the minimum zoom at the current viewport size
    virtual double minimumZoom() const = 0;
    virtual double maximumCenterLatitudeAtZoom(const QGeoCameraData &cameraData) const = 0;
    virtual double minimumCenterLatitudeAtZoom(const QGeoCameraData &cameraData) const = 0;

    virtual QGeoCoordinate itemPositionToCoordinate(const QDoubleVector2D &pos, bool clipToViewport = true) const = 0;
    virtual QDoubleVector2D coordinateToItemPosition(const QGeoCoordinate &coordinate, bool clipToViewport = true) const = 0;

    virtual ProjectionGroup projectionGroup() const = 0;
    virtual Datum datum() const = 0;
    virtual ProjectionType projectionType() const = 0;

    // Returns the new map center after anchoring coordinate to anchorPoint on the screen
    virtual QGeoCoordinate anchorCoordinateToPoint(const QGeoCoordinate &coordinate, const QPointF &anchorPoint) const;

    virtual QGeoShape visibleRegion() const;
    virtual bool setBearing(qreal bearing, const QGeoCoordinate &coordinate);
    virtual QMatrix4x4 projectionTransformation() const = 0; // This brings a mercator coord into the correct viewport coordinate.
    virtual QMatrix4x4 projectionTransformation_centered() const = 0; // Same as projectionTransformation, but the center of the camera is around 0,0.
                                                                      // Requires subsequent shifting of the geometry to fit such camera.
    virtual const QMatrix4x4 &qsgTransform() const = 0;
    virtual QDoubleVector3D centerMercator() const = 0;

    void setItemToWindowTransform(const QTransform &itemToWindowTransform);
    virtual QTransform itemToWindowTransform() const;

    QTransform       m_itemToWindowTransform;
    mutable bool     m_qsgTransformDirty = true;
};

class Q_LOCATION_PRIVATE_EXPORT QGeoProjectionWebMercator : public QGeoProjection
{
public:
    QGeoProjectionWebMercator();
    ~QGeoProjectionWebMercator();

    // From QGeoProjection
    double minimumZoom() const override;
    QMatrix4x4 projectionTransformation() const override;
    QMatrix4x4 projectionTransformation_centered() const override;
    const QMatrix4x4 &qsgTransform() const override;
    QDoubleVector3D centerMercator() const override;

    double maximumCenterLatitudeAtZoom(const QGeoCameraData &cameraData) const override;
    double minimumCenterLatitudeAtZoom(const QGeoCameraData &cameraData) const override;

    void setVisibleArea(const QRectF &visibleArea) override;
    void setViewportSize(const QSize &size) override;
    void setCameraData(const QGeoCameraData &cameraData, bool force = true) override;
    QGeoCameraData cameraData() const override;

    QGeoCoordinate itemPositionToCoordinate(const QDoubleVector2D &pos, bool clipToViewport = true) const override;
    QDoubleVector2D coordinateToItemPosition(const QGeoCoordinate &coordinate, bool clipToViewport = true) const override;

    QGeoProjection::ProjectionGroup projectionGroup() const override;
    QGeoProjection::Datum datum() const override;
    QGeoProjection::ProjectionType projectionType() const override;

    QGeoCoordinate anchorCoordinateToPoint(const QGeoCoordinate &coordinate, const QPointF &anchorPoint) const override;
    bool setBearing(qreal bearing, const QGeoCoordinate &coordinate) override;

    QGeoShape visibleRegion() const override;

    // Specific to QGeoProjectionWebMercator
    double mapWidth() const; // The size of the underlying map, at the current zoom level.
    double mapHeight() const;

    QDoubleVector2D geoToMapProjection(const QGeoCoordinate &coordinate) const;
    QGeoCoordinate mapProjectionToGeo(const QDoubleVector2D &projection) const;

    int projectionWrapFactor(const QDoubleVector2D &projection) const;
    QDoubleVector2D wrapMapProjection(const QDoubleVector2D &projection) const;
    QDoubleVector2D unwrapMapProjection(const QDoubleVector2D &wrappedProjection) const;

    QDoubleVector2D wrappedMapProjectionToItemPosition(const QDoubleVector2D &wrappedProjection) const;
    QDoubleVector2D itemPositionToWrappedMapProjection(const QDoubleVector2D &itemPosition) const;

    QDoubleVector2D geoToWrappedMapProjection(const QGeoCoordinate &coordinate) const;
    QGeoCoordinate wrappedMapProjectionToGeo(const QDoubleVector2D &wrappedProjection) const;
    QMatrix4x4 quickItemTransformation(const QGeoCoordinate &coordinate, const QPointF &anchorPoint, qreal zoomLevel) const;

    bool isProjectable(const QDoubleVector2D &wrappedProjection) const;
    QList<QDoubleVector2D> visibleGeometry() const;
    QList<QDoubleVector2D> visibleGeometryExpanded() const;
    QList<QDoubleVector2D> projectableGeometry() const;

    inline QDoubleVector2D viewportToWrappedMapProjection(const QDoubleVector2D &itemPosition) const;
    inline QDoubleVector2D viewportToWrappedMapProjection(const QDoubleVector2D &itemPosition, double &s) const;

    QPair<QGeoCoordinate, qreal> fitViewportToGeoRectangle(const QGeoRectangle &rectangle,
                                                           const QMargins &margins) const;

private:
    void setupCamera();
    void updateVisibleRegion();

public:
    struct Line2D
    {
        Line2D();
        Line2D(const QDoubleVector2D &linePoint, const QDoubleVector2D &lineDirection);

        bool isValid() const;

        QDoubleVector2D m_point;
        QDoubleVector2D m_direction;
    };

    struct Plane
    {
        Plane();
        Plane(const QDoubleVector3D &planePoint, const QDoubleVector3D &planeNormal);

        QDoubleVector3D lineIntersection(const QDoubleVector3D &linePoint, const QDoubleVector3D &lineDirection) const;
        inline QDoubleVector3D lineIntersection(const QDoubleVector3D &linePoint, const QDoubleVector3D &lineDirection, double &s) const;
        Line2D planeXYIntersection() const;
        bool isValid() const;

        QDoubleVector3D m_point;
        QDoubleVector3D m_normal;
    };

#ifdef QT_LOCATION_DEBUG
public:
#else
protected:
#endif
    QGeoCameraData m_cameraData;
    double m_mapEdgeSize;
    double m_minimumZoom;
    // mercator to camera transform for coordinates (not tiles!)
    double m_cameraCenterXMercator;
    double m_cameraCenterYMercator;

    // cameraToScreen transform
    double m_viewportWidth; // in pixels
    double m_viewportHeight; // in pixels
    double m_1_viewportWidth;
    double m_1_viewportHeight;

    QDoubleMatrix4x4 m_cameraMatrix;
    QDoubleMatrix4x4 m_cameraMatrix0;
    QDoubleMatrix4x4 m_transformation;
    QDoubleMatrix4x4 m_transformation0;
    QDoubleMatrix4x4 m_quickItemTransformation;
    QDoubleVector3D  m_eye;
    QDoubleVector3D  m_up;
    QDoubleVector3D  m_center;
    QDoubleVector3D  m_view;
    QDoubleVector3D  m_viewNormalized;
    QDoubleVector3D  m_side;
    QDoubleVector3D  m_centerNearPlane;
    double           m_sideLengthPixels; // map edge size at integer zoom level
    double           m_aperture;
    double           m_nearPlane;
    double           m_farPlane;
    double           m_halfWidth;
    double           m_halfHeight;
    double           m_minimumUnprojectableY;
    double           m_verticalEstateToSkip;

    // For the clipping region
    QDoubleVector3D  m_centerMercator;
    QDoubleVector3D  m_eyeMercator;
    QDoubleVector3D  m_eyeMercator0;
    QDoubleVector3D  m_viewMercator;
    QDoubleVector3D  m_upMercator;
    QDoubleVector3D  m_sideMercator;
    QDoubleVector3D  m_centerNearPlaneMercator;
    double           m_nearPlaneMercator;
    Line2D           m_nearPlaneMapIntersection;

    QList<QDoubleVector2D> m_visibleRegion;
    QList<QDoubleVector2D> m_visibleRegionExpanded;
    QList<QDoubleVector2D> m_projectableRegion;
    bool             m_visibleRegionDirty;

    mutable QMatrix4x4 m_qsgTransform;
    QRectF           m_visibleArea;

    Q_DISABLE_COPY(QGeoProjectionWebMercator)
};

QT_END_NAMESPACE

#endif // QGEOPROJECTION_H
