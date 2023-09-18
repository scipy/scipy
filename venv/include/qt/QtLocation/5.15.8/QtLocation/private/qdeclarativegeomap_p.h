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

#ifndef QDECLARATIVEGEOMAP_H
#define QDECLARATIVEGEOMAP_H

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
#include <QtLocation/private/qdeclarativegeomapitemview_p.h>
#include <QtLocation/private/qquickgeomapgesturearea_p.h>
#include <QtLocation/private/qdeclarativegeomapitemgroup_p.h>
#include <QtLocation/qgeoserviceprovider.h>
#include <QtLocation/private/qgeocameradata_p.h>
#include <QtLocation/private/qgeocameracapabilities_p.h>
#include <QtQuick/QQuickItem>
#include <QtCore/QList>
#include <QtCore/QPointer>
#include <QtGui/QColor>
#include <QtPositioning/qgeorectangle.h>
#include <QtLocation/private/qgeomap_p.h>
#include <QtQuick/private/qquickitemchangelistener_p.h>

QT_BEGIN_NAMESPACE

class QDeclarativeGeoServiceProvider;
class QDeclarativeGeoMapType;
class QDeclarativeGeoMapCopyrightNotice;
class QDeclarativeGeoMapParameter;

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeGeoMap : public QQuickItem
{
    Q_OBJECT
    Q_ENUMS(QGeoServiceProvider::Error)
    Q_PROPERTY(QQuickGeoMapGestureArea *gesture READ gesture CONSTANT)
    Q_PROPERTY(QDeclarativeGeoServiceProvider *plugin READ plugin WRITE setPlugin NOTIFY pluginChanged)
    Q_PROPERTY(qreal minimumZoomLevel READ minimumZoomLevel WRITE setMinimumZoomLevel NOTIFY minimumZoomLevelChanged)
    Q_PROPERTY(qreal maximumZoomLevel READ maximumZoomLevel WRITE setMaximumZoomLevel NOTIFY maximumZoomLevelChanged)
    Q_PROPERTY(qreal zoomLevel READ zoomLevel WRITE setZoomLevel NOTIFY zoomLevelChanged)

    Q_PROPERTY(qreal tilt READ tilt WRITE setTilt NOTIFY tiltChanged)
    Q_PROPERTY(qreal minimumTilt READ minimumTilt WRITE setMinimumTilt NOTIFY minimumTiltChanged)
    Q_PROPERTY(qreal maximumTilt READ maximumTilt WRITE setMaximumTilt NOTIFY maximumTiltChanged)

    Q_PROPERTY(qreal bearing READ bearing WRITE setBearing NOTIFY bearingChanged)

    Q_PROPERTY(qreal fieldOfView READ fieldOfView WRITE setFieldOfView NOTIFY fieldOfViewChanged)
    Q_PROPERTY(qreal minimumFieldOfView READ minimumFieldOfView WRITE setMinimumFieldOfView NOTIFY minimumFieldOfViewChanged)
    Q_PROPERTY(qreal maximumFieldOfView READ maximumFieldOfView WRITE setMaximumFieldOfView NOTIFY minimumFieldOfViewChanged)

    Q_PROPERTY(QDeclarativeGeoMapType *activeMapType READ activeMapType WRITE setActiveMapType NOTIFY activeMapTypeChanged)
    Q_PROPERTY(QQmlListProperty<QDeclarativeGeoMapType> supportedMapTypes READ supportedMapTypes NOTIFY supportedMapTypesChanged)
    Q_PROPERTY(QGeoCoordinate center READ center WRITE setCenter NOTIFY centerChanged)
    Q_PROPERTY(QList<QObject *> mapItems READ mapItems NOTIFY mapItemsChanged)
    Q_PROPERTY(QList<QObject *> mapParameters READ mapParameters)
    Q_PROPERTY(QGeoServiceProvider::Error error READ error NOTIFY errorChanged)
    Q_PROPERTY(QString errorString READ errorString NOTIFY errorChanged)
    Q_PROPERTY(QGeoShape visibleRegion READ visibleRegion WRITE setVisibleRegion NOTIFY visibleRegionChanged)
    Q_PROPERTY(bool copyrightsVisible READ copyrightsVisible WRITE setCopyrightsVisible NOTIFY copyrightsVisibleChanged)
    Q_PROPERTY(QColor color READ color WRITE setColor NOTIFY colorChanged)
    Q_PROPERTY(bool mapReady READ mapReady NOTIFY mapReadyChanged)
    Q_PROPERTY(QRectF visibleArea READ visibleArea WRITE setVisibleArea NOTIFY visibleAreaChanged  REVISION 12)
    Q_INTERFACES(QQmlParserStatus)

public:

    explicit QDeclarativeGeoMap(QQuickItem *parent = 0);
    ~QDeclarativeGeoMap();

    void setPlugin(QDeclarativeGeoServiceProvider *plugin);
    QDeclarativeGeoServiceProvider *plugin() const;

    void setActiveMapType(QDeclarativeGeoMapType *mapType);
    QDeclarativeGeoMapType *activeMapType() const;

    void setMinimumZoomLevel(qreal minimumZoomLevel, bool userSet = true);
    qreal minimumZoomLevel() const;
    qreal implicitMinimumZoomLevel() const;
    qreal effectiveMinimumZoomLevel() const;

    void setMaximumZoomLevel(qreal maximumZoomLevel, bool userSet = true);
    qreal maximumZoomLevel() const;

    void setZoomLevel(qreal zoomLevel);
    qreal zoomLevel() const;

    void setBearing(qreal bearing);
    qreal bearing() const;

    void setTilt(qreal tilt);
    qreal tilt() const;
    void setMinimumTilt(qreal minimumTilt, bool userSet = true);
    qreal minimumTilt() const;
    void setMaximumTilt(qreal maximumTilt, bool userSet = true);
    qreal maximumTilt() const;

    void setFieldOfView(qreal fieldOfView);
    qreal fieldOfView() const;
    void setMinimumFieldOfView(qreal minimumFieldOfView, bool userSet = true);
    qreal minimumFieldOfView() const;
    void setMaximumFieldOfView(qreal maximumFieldOfView, bool userSet = true);
    qreal maximumFieldOfView() const;

    void setCenter(const QGeoCoordinate &center);
    QGeoCoordinate center() const;

    void setVisibleRegion(const QGeoShape &shape);
    QGeoShape visibleRegion() const;

    void setCopyrightsVisible(bool visible);
    bool copyrightsVisible() const;

    void setColor(const QColor &color);
    QColor color() const;

    QRectF visibleArea() const;
    void setVisibleArea(const QRectF &visibleArea);

    bool mapReady() const;

    QQmlListProperty<QDeclarativeGeoMapType> supportedMapTypes();

    Q_INVOKABLE void setBearing(qreal bearing, const QGeoCoordinate &coordinate);
    Q_INVOKABLE void alignCoordinateToPoint(const QGeoCoordinate &coordinate, const QPointF &point);

    Q_INVOKABLE void removeMapItem(QDeclarativeGeoMapItemBase *item);
    Q_INVOKABLE void addMapItem(QDeclarativeGeoMapItemBase *item);

    Q_INVOKABLE void addMapItemGroup(QDeclarativeGeoMapItemGroup *itemGroup);
    Q_INVOKABLE void removeMapItemGroup(QDeclarativeGeoMapItemGroup *itemGroup);

    Q_INVOKABLE void removeMapItemView(QDeclarativeGeoMapItemView *itemView);
    Q_INVOKABLE void addMapItemView(QDeclarativeGeoMapItemView *itemView);

    Q_INVOKABLE void clearMapItems();
    QList<QObject *> mapItems();

    Q_INVOKABLE void addMapParameter(QDeclarativeGeoMapParameter *parameter);
    Q_INVOKABLE void removeMapParameter(QDeclarativeGeoMapParameter *parameter);
    Q_INVOKABLE void clearMapParameters();
    QList<QObject *> mapParameters();

    void addMapObject(QGeoMapObject *object); // Not invokable as currently meant to be used through a main MapObjectView
    void removeMapObject(QGeoMapObject *object);
    void clearMapObjects();
    QList<QGeoMapObject *> mapObjects();


    Q_INVOKABLE QGeoCoordinate toCoordinate(const QPointF &position, bool clipToViewPort = true) const;
    Q_INVOKABLE QPointF fromCoordinate(const QGeoCoordinate &coordinate, bool clipToViewPort = true) const;

    QQuickGeoMapGestureArea *gesture();

    Q_INVOKABLE void fitViewportToMapItems(const QVariantList &items = {});
    Q_INVOKABLE void fitViewportToVisibleMapItems();
    Q_INVOKABLE void pan(int dx, int dy);
    Q_INVOKABLE void prefetchData(); // optional hint for prefetch
    Q_INVOKABLE void clearData();
    Q_REVISION(13) Q_INVOKABLE void fitViewportToGeoShape(const QGeoShape &shape, QVariant margins);
    void fitViewportToGeoShape(const QGeoShape &shape, const QMargins &borders = QMargins(10, 10, 10, 10));

    QString errorString() const;
    QGeoServiceProvider::Error error() const;
    QGeoMap* map() const;

    // From QQuickItem
    void itemChange(ItemChange, const ItemChangeData &) override;

Q_SIGNALS:
    void pluginChanged(QDeclarativeGeoServiceProvider *plugin);
    void zoomLevelChanged(qreal zoomLevel);
    void centerChanged(const QGeoCoordinate &coordinate);
    void activeMapTypeChanged();
    void supportedMapTypesChanged();
    void minimumZoomLevelChanged();
    void maximumZoomLevelChanged();
    void mapItemsChanged();
    void errorChanged();
    void copyrightLinkActivated(const QString &link);
    void copyrightsVisibleChanged(bool visible);
    void colorChanged(const QColor &color);
    void bearingChanged(qreal bearing);
    void tiltChanged(qreal tilt);
    void fieldOfViewChanged(qreal fieldOfView);
    void minimumTiltChanged(qreal minimumTilt);
    void maximumTiltChanged(qreal maximumTilt);
    void minimumFieldOfViewChanged(qreal minimumFieldOfView);
    void maximumFieldOfViewChanged(qreal maximumFieldOfView);
    void copyrightsChanged(const QImage &copyrightsImage);
    void copyrightsChanged(const QString &copyrightsHtml);
    void mapReadyChanged(bool ready);
    Q_REVISION(11) void mapObjectsChanged();
    void visibleAreaChanged();
    Q_REVISION(14) void visibleRegionChanged();

protected:
    void mousePressEvent(QMouseEvent *event) override ;
    void mouseMoveEvent(QMouseEvent *event) override ;
    void mouseReleaseEvent(QMouseEvent *event) override ;
    void mouseUngrabEvent() override ;
    void touchUngrabEvent() override;
    void touchEvent(QTouchEvent *event) override ;
#if QT_CONFIG(wheelevent)
    void wheelEvent(QWheelEvent *event) override ;
#endif

    bool childMouseEventFilter(QQuickItem *item, QEvent *event) override;
    bool sendMouseEvent(QMouseEvent *event);
    bool sendTouchEvent(QTouchEvent *event);

    void componentComplete() override;
    QSGNode *updatePaintNode(QSGNode *, UpdatePaintNodeData *) override;
    void geometryChanged(const QRectF &newGeometry, const QRectF &oldGeometry) override;

    void setError(QGeoServiceProvider::Error error, const QString &errorString);
    void initialize();
    void setZoomLevel(qreal zoomLevel, bool overzoom);
    bool addMapChild(QObject *child);
    bool removeMapChild(QObject *child);
    bool isGroupNested(QDeclarativeGeoMapItemGroup *group);

    bool addMapItem_real(QDeclarativeGeoMapItemBase *item);
    bool removeMapItem_real(QDeclarativeGeoMapItemBase *item);
    bool addMapItemGroup_real(QDeclarativeGeoMapItemGroup *itemGroup);
    bool removeMapItemGroup_real(QDeclarativeGeoMapItemGroup *itemGroup);
    bool addMapItemView_real(QDeclarativeGeoMapItemView *itemView);
    bool removeMapItemView_real(QDeclarativeGeoMapItemView *itemView);
    void updateItemToWindowTransform();
    void onSGNodeChanged();

private Q_SLOTS:
    void mappingManagerInitialized();
    void pluginReady();
    void onSupportedMapTypesChanged();
    void onCameraCapabilitiesChanged(const QGeoCameraCapabilities &oldCameraCapabilities);
    void onAttachedCopyrightNoticeVisibilityChanged();
    void onCameraDataChanged(const QGeoCameraData &cameraData);

private:
    void setupMapView(QDeclarativeGeoMapItemView *view);
    void populateMap();
    void populateParameters();
    void fitViewportToMapItemsRefine(const QList<QPointer<QDeclarativeGeoMapItemBase> > &mapItems, bool refine, bool onlyVisible);
    bool isInteractive();
    void attachCopyrightNotice(bool initialVisibility);
    void detachCopyrightNotice(bool currentVisibility);
    QMargins mapMargins() const;

private:
    QDeclarativeGeoServiceProvider *m_plugin;
    QGeoMappingManager *m_mappingManager;
    QDeclarativeGeoMapType *m_activeMapType;
    QList<QDeclarativeGeoMapType *> m_supportedMapTypes;
    QList<QDeclarativeGeoMapItemView *> m_mapViews;
    QQuickGeoMapGestureArea *m_gestureArea;
    QPointer<QGeoMap> m_map;
    QPointer<QDeclarativeGeoMapCopyrightNotice> m_copyrights;
    QList<QPointer<QDeclarativeGeoMapItemBase> > m_mapItems;
    QList<QPointer<QDeclarativeGeoMapItemGroup> > m_mapItemGroups;
    QString m_errorString;
    QGeoServiceProvider::Error m_error;
    QGeoRectangle m_visibleRegion;
    QColor m_color;
    QGeoCameraData m_cameraData;
    bool m_componentCompleted;
    bool m_pendingFitViewport;
    bool m_copyrightsVisible;
    double m_maximumViewportLatitude;
    double m_minimumViewportLatitude = 0.0;
    bool m_initialized;
    bool m_sgNodeHasChanged = false;
    QList<QDeclarativeGeoMapParameter *> m_mapParameters;
    QList<QGeoMapObject*> m_pendingMapObjects; // Used only in the initialization phase
    QGeoCameraCapabilities m_cameraCapabilities;
    qreal m_userMinimumZoomLevel;
    qreal m_userMaximumZoomLevel;

    qreal m_minimumTilt;
    qreal m_maximumTilt;
    qreal m_userMinimumTilt;
    qreal m_userMaximumTilt;

    qreal m_minimumFieldOfView;
    qreal m_maximumFieldOfView;
    qreal m_userMinimumFieldOfView;
    qreal m_userMaximumFieldOfView;

    int m_copyNoticesVisible = 0;
    qreal m_maxChildZ = 0;
    QRectF m_visibleArea;


    friend class QDeclarativeGeoMapItem;
    friend class QDeclarativeGeoMapItemView;
    friend class QQuickGeoMapGestureArea;
    friend class QDeclarativeGeoMapCopyrightNotice;
    Q_DISABLE_COPY(QDeclarativeGeoMap)
};


QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativeGeoMap)

#endif
