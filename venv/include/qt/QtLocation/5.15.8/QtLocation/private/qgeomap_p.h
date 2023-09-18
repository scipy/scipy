/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
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
#ifndef QGEOMAP_P_H
#define QGEOMAP_P_H

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
#include <QtLocation/private/qgeomaptype_p.h>
#include <QtLocation/private/qgeocameracapabilities_p.h>
#include <QtCore/QObject>
#include <QtPositioning/private/qdoublevector2d_p.h>
#include <QtLocation/private/qgeoprojection_p.h>
#include <QtLocation/qgeoroute.h>
#include <QTransform>

QT_BEGIN_NAMESPACE

class QGeoMappingManagerEngine;
class QGeoMapPrivate;
class QGeoMapController;
class QGeoCoordinate;
class QSGNode;
class QQuickWindow;
class QGeoMapParameter;
class QDeclarativeGeoMapItemBase;
class QGeoMapObject;
class QDeclarativeGeoMap;

class Q_LOCATION_PRIVATE_EXPORT QGeoMap : public QObject
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QGeoMap)

    Q_ENUMS(Capability)
    Q_FLAGS(Capabilities)
public:
    enum ItemType {
        NoItem = 0x0000,
        MapRectangle = 0x0001,
        MapCircle = 0x0002,
        MapPolyline = 0x0004,
        MapPolygon = 0x0008,
        MapQuickItem = 0x0010,
        CustomMapItem = 0x8000
    };

    Q_DECLARE_FLAGS(ItemTypes, ItemType)

    enum Capability {
        SupportsNothing = 0x0000,
        SupportsVisibleRegion = 0x0001,
        SupportsSetBearing = 0x0002,
        SupportsAnchoringCoordinate = 0x0004,
        SupportsFittingViewportToGeoRectangle = 0x0008,
        SupportsVisibleArea = 0x0010,
    };

    Q_DECLARE_FLAGS(Capabilities, Capability)

    virtual ~QGeoMap();

    // Sets the display size
    void setViewportSize(const QSize& viewportSize);
    QSize viewportSize() const;
    int viewportWidth() const;
    int viewportHeight() const;


    const QGeoCameraData &cameraData() const;
    QGeoCameraCapabilities cameraCapabilities() const;
    virtual Capabilities capabilities() const;

    void setActiveMapType(const QGeoMapType mapType);
    const QGeoMapType activeMapType() const;

    // returns the minimum zoom at the current viewport size
    double minimumZoom() const;
    double maximumCenterLatitudeAtZoom(const QGeoCameraData &cameraData) const;
    double minimumCenterLatitudeAtZoom(const QGeoCameraData &cameraData) const;

    // returns the size of the underlying map, at the current zoom level. Unrelated to width()/height()/size().
    double mapWidth() const;
    double mapHeight() const;

    const QGeoProjection &geoProjection() const;

    virtual void prefetchData();
    virtual void clearData();

    void addParameter(QGeoMapParameter *param);
    void removeParameter(QGeoMapParameter *param);
    void clearParameters();

    ItemTypes supportedMapItemTypes() const;

    void addMapItem(QDeclarativeGeoMapItemBase *item);
    void removeMapItem(QDeclarativeGeoMapItemBase *item);
    void clearMapItems();

    virtual bool createMapObjectImplementation(QGeoMapObject *obj);
    QList<QGeoMapObject *> mapObjects() const;


    virtual QString copyrightsStyleSheet() const;
    virtual void setAcceptedGestures(bool pan, bool flick, bool pinch, bool rotate, bool tilt);
    virtual bool handleEvent(QEvent *event);

    virtual bool setBearing(qreal bearing, const QGeoCoordinate &coordinate);
    virtual QGeoShape visibleRegion() const;
    virtual bool anchorCoordinateToPoint(const QGeoCoordinate &coordinate, const QPointF &anchorPoint);
    virtual bool fitViewportToGeoRectangle(const QGeoRectangle &rectangle, const QMargins &borders);

    virtual void setCopyrightVisible(bool visible);
    virtual void removeMapObject(QGeoMapObject *obj);
    virtual QList<QObject *> mapObjectsAt(const QGeoCoordinate &coordinate) const;
    virtual void setItemToWindowTransform(const QTransform &itemToWindowTransform);

    void setVisibleArea(const QRectF &visibleArea);
    QRectF visibleArea() const;

protected:
    QGeoMap(QGeoMapPrivate &dd, QObject *parent = 0);
    void setCameraData(const QGeoCameraData &cameraData);
    void setCameraCapabilities(const QGeoCameraCapabilities &cameraCapabilities);
    virtual QSGNode *updateSceneGraph(QSGNode *node, QQuickWindow *window) = 0;

Q_SIGNALS:
    void cameraDataChanged(const QGeoCameraData &cameraData);
    void sgNodeChanged();
    void activeMapTypeChanged();
    void cameraCapabilitiesChanged(const QGeoCameraCapabilities &oldCameraCapabilities);
    void copyrightsChanged(const QImage &copyrightsImage);
    void copyrightsChanged(const QString &copyrightsHtml);
    void copyrightsStyleSheetChanged(const QString &styleSheet);
    void visibleAreaChanged();

private:
    Q_DISABLE_COPY(QGeoMap)
    friend class QDeclarativeGeoMap; //updateSceneGraph
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QGeoMap::ItemTypes)

QT_END_NAMESPACE

#endif // QGEOMAP_P_H
