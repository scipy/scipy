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

#ifndef QGEOSERVICEPROVIDER_H
#define QGEOSERVICEPROVIDER_H

#include <QtCore/QVariant>
#include <QtCore/QString>
#include <QtCore/QObject>
#include <QtLocation/qlocationglobal.h>

QT_BEGIN_NAMESPACE

class QLocale;
class QStringList;
class QGeoCodingManager;
class QGeoMappingManager;
class QGeoRoutingManager;
class QPlaceManager;
class QNavigationManager;
class QGeoCodingManagerEngine;
class QGeoMappingManagerEngine;
class QGeoRoutingManagerEngine;
class QPlaceManagerEngine;
class QNavigationManagerEngine;
class QGeoServiceProviderPrivate;
class QQmlEngine;

class Q_LOCATION_EXPORT QGeoServiceProvider : public QObject
{
    Q_OBJECT
    Q_ENUMS(Error)
public:
    enum Error {
        NoError,
        NotSupportedError,
        UnknownParameterError,
        MissingRequiredParameterError,
        ConnectionError,
        LoaderError
    };

    enum RoutingFeature {
        NoRoutingFeatures               = 0,
        OnlineRoutingFeature            = (1<<0),
        OfflineRoutingFeature           = (1<<1),
        LocalizedRoutingFeature         = (1<<2),
        RouteUpdatesFeature             = (1<<3),
        AlternativeRoutesFeature        = (1<<4),
        ExcludeAreasRoutingFeature      = (1<<5),
        AnyRoutingFeatures              = ~(0)
    };

    enum GeocodingFeature {
        NoGeocodingFeatures             = 0,
        OnlineGeocodingFeature          = (1<<0),
        OfflineGeocodingFeature         = (1<<1),
        ReverseGeocodingFeature         = (1<<2),
        LocalizedGeocodingFeature       = (1<<3),
        AnyGeocodingFeatures            = ~(0)
    };

    enum MappingFeature {
        NoMappingFeatures               = 0,
        OnlineMappingFeature            = (1<<0),
        OfflineMappingFeature           = (1<<1),
        LocalizedMappingFeature         = (1<<2),
        AnyMappingFeatures              = ~(0)
    };

    enum PlacesFeature {
        NoPlacesFeatures                = 0,
        OnlinePlacesFeature             = (1<<0),
        OfflinePlacesFeature            = (1<<1),
        SavePlaceFeature                = (1<<2),
        RemovePlaceFeature              = (1<<3),
        SaveCategoryFeature             = (1<<4),
        RemoveCategoryFeature           = (1<<5),
        PlaceRecommendationsFeature     = (1<<6),
        SearchSuggestionsFeature        = (1<<7),
        LocalizedPlacesFeature          = (1<<8),
        NotificationsFeature            = (1<<9),
        PlaceMatchingFeature            = (1<<10),
        AnyPlacesFeatures               = ~(0)
    };

    enum NavigationFeature {
        NoNavigationFeatures            = 0,
        OnlineNavigationFeature         = (1<<0),
        OfflineNavigationFeature        = (1<<1),
        AnyNavigationFeatures           = ~(0)
    };

    Q_DECLARE_FLAGS(RoutingFeatures, RoutingFeature)
    Q_FLAGS(RoutingFeatures)

    Q_DECLARE_FLAGS(GeocodingFeatures, GeocodingFeature)
    Q_FLAGS(GeocodingFeatures)

    Q_DECLARE_FLAGS(MappingFeatures, MappingFeature)
    Q_FLAGS(MappingFeatures)

    Q_DECLARE_FLAGS(PlacesFeatures, PlacesFeature)
    Q_FLAGS(PlacesFeatures)

    Q_DECLARE_FLAGS(NavigationFeatures, NavigationFeature)
    Q_FLAGS(NavigationFeatures)

    static QStringList availableServiceProviders();
    QGeoServiceProvider(const QString &providerName,
                        const QVariantMap &parameters = QVariantMap(),
                        bool allowExperimental = false);

    ~QGeoServiceProvider();

    RoutingFeatures routingFeatures() const;
    GeocodingFeatures geocodingFeatures() const;
    MappingFeatures mappingFeatures() const;
    PlacesFeatures placesFeatures() const;
    NavigationFeatures navigationFeatures() const;

    QGeoCodingManager *geocodingManager() const;
    QGeoMappingManager *mappingManager() const;
    QGeoRoutingManager *routingManager() const;
    QPlaceManager *placeManager() const;
    QNavigationManager *navigationManager() const;

    Error error() const;
    QString errorString() const;

    Error mappingError() const;
    QString mappingErrorString() const;
    Error geocodingError() const;
    QString geocodingErrorString() const;
    Error routingError() const;
    QString routingErrorString() const;
    Error placesError() const;
    QString placesErrorString() const;
    Error navigationError() const;
    QString navigationErrorString() const;

    void setParameters(const QVariantMap &parameters);
    void setLocale(const QLocale &locale);
    void setAllowExperimental(bool allow);
    void setQmlEngine(QQmlEngine *engine);

private:
    QGeoServiceProviderPrivate *d_ptr;
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QGeoServiceProvider::RoutingFeatures)
Q_DECLARE_OPERATORS_FOR_FLAGS(QGeoServiceProvider::GeocodingFeatures)
Q_DECLARE_OPERATORS_FOR_FLAGS(QGeoServiceProvider::MappingFeatures)
Q_DECLARE_OPERATORS_FOR_FLAGS(QGeoServiceProvider::PlacesFeatures)
Q_DECLARE_OPERATORS_FOR_FLAGS(QGeoServiceProvider::NavigationFeatures)

QT_END_NAMESPACE

#endif
