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

#ifndef QDECLARATIVEQGEOSERVICEPROVIDER_H
#define QDECLARATIVEQGEOSERVICEPROVIDER_H

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

#include <QtCore/QMap>
#include <QtCore/QString>
#include <QtCore/QStringList>
#include <QtCore/QVariant>
#include <QtQml/qqml.h>
#include <QtQml/QQmlParserStatus>
#include <QtQml/QQmlListProperty>
#include <QtLocation/QGeoServiceProvider>
#include <QtPositioningQuick/private/qdeclarativepluginparameter_p.h>

QT_BEGIN_NAMESPACE

class QDeclarativeGeoServiceProviderRequirements;

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeGeoServiceProvider : public QObject, public QQmlParserStatus
{
    Q_OBJECT
    Q_ENUMS(RoutingFeature)
    Q_ENUMS(GeocodingFeature)
    Q_ENUMS(MappingFeature)
    Q_ENUMS(PlacesFeature)

    Q_PROPERTY(QString name READ name WRITE setName NOTIFY nameChanged)
    Q_PROPERTY(QStringList availableServiceProviders READ availableServiceProviders CONSTANT)
    Q_PROPERTY(QQmlListProperty<QDeclarativePluginParameter> parameters READ parameters)
    Q_PROPERTY(QDeclarativeGeoServiceProviderRequirements *required READ requirements WRITE setRequirements)
    Q_PROPERTY(QStringList locales READ locales WRITE setLocales NOTIFY localesChanged)
    Q_PROPERTY(QStringList preferred READ preferred WRITE setPreferred NOTIFY preferredChanged)
    Q_PROPERTY(bool allowExperimental READ allowExperimental WRITE setAllowExperimental NOTIFY allowExperimentalChanged)
    Q_PROPERTY(bool isAttached READ isAttached NOTIFY attached)

    Q_CLASSINFO("DefaultProperty", "parameters")
    Q_INTERFACES(QQmlParserStatus)

public:
    explicit QDeclarativeGeoServiceProvider(QObject *parent = nullptr);
    ~QDeclarativeGeoServiceProvider();

    enum RoutingFeature {
        NoRoutingFeatures               = QGeoServiceProvider::NoRoutingFeatures,
        OnlineRoutingFeature            = QGeoServiceProvider::OnlineRoutingFeature,
        OfflineRoutingFeature           = QGeoServiceProvider::OfflineRoutingFeature,
        LocalizedRoutingFeature         = QGeoServiceProvider::LocalizedRoutingFeature,
        RouteUpdatesFeature             = QGeoServiceProvider::RouteUpdatesFeature,
        AlternativeRoutesFeature        = QGeoServiceProvider::AlternativeRoutesFeature,
        ExcludeAreasRoutingFeature      = QGeoServiceProvider::ExcludeAreasRoutingFeature,
        AnyRoutingFeatures              = QGeoServiceProvider::AnyRoutingFeatures
    };

    enum GeocodingFeature {
        NoGeocodingFeatures             = QGeoServiceProvider::NoGeocodingFeatures,
        OnlineGeocodingFeature          = QGeoServiceProvider::OnlineGeocodingFeature,
        OfflineGeocodingFeature         = QGeoServiceProvider::OfflineGeocodingFeature,
        ReverseGeocodingFeature         = QGeoServiceProvider::ReverseGeocodingFeature,
        LocalizedGeocodingFeature       = QGeoServiceProvider::LocalizedGeocodingFeature,
        AnyGeocodingFeatures            = QGeoServiceProvider::AnyGeocodingFeatures
    };

    enum MappingFeature {
        NoMappingFeatures               = QGeoServiceProvider::NoMappingFeatures,
        OnlineMappingFeature            = QGeoServiceProvider::OnlineMappingFeature,
        OfflineMappingFeature           = QGeoServiceProvider::OfflineMappingFeature,
        LocalizedMappingFeature         = QGeoServiceProvider::LocalizedMappingFeature,
        AnyMappingFeatures              = QGeoServiceProvider::AnyMappingFeatures
    };

    enum PlacesFeature {
        NoPlacesFeatures                = QGeoServiceProvider::NoPlacesFeatures,
        OnlinePlacesFeature             = QGeoServiceProvider::OnlinePlacesFeature,
        OfflinePlacesFeature            = QGeoServiceProvider::OfflinePlacesFeature,
        SavePlaceFeature                = QGeoServiceProvider::SavePlaceFeature,
        RemovePlaceFeature              = QGeoServiceProvider::RemovePlaceFeature,
        SaveCategoryFeature             = QGeoServiceProvider::SaveCategoryFeature,
        RemoveCategoryFeature           = QGeoServiceProvider::RemoveCategoryFeature,
        PlaceRecommendationsFeature     = QGeoServiceProvider::PlaceRecommendationsFeature,
        SearchSuggestionsFeature        = QGeoServiceProvider::SearchSuggestionsFeature,
        LocalizedPlacesFeature          = QGeoServiceProvider::LocalizedPlacesFeature,
        NotificationsFeature            = QGeoServiceProvider::NotificationsFeature,
        PlaceMatchingFeature            = QGeoServiceProvider::PlaceMatchingFeature,
        AnyPlacesFeatures               = QGeoServiceProvider::AnyPlacesFeatures
    };

    enum NavigationFeature {
        NoNavigationFeatures            = QGeoServiceProvider::NoNavigationFeatures,
        OnlineNavigationFeature         = QGeoServiceProvider::OnlineNavigationFeature,
        OfflineNavigationFeature        = QGeoServiceProvider::OfflineNavigationFeature,
        AnyNavigationFeatures           = QGeoServiceProvider::AnyNavigationFeatures
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

    // From QQmlParserStatus
    virtual void classBegin() {}
    virtual void componentComplete();

    void setName(const QString &name);
    QString name() const;

    QQmlListProperty<QDeclarativePluginParameter> parameters();
    QVariantMap parameterMap() const;

    QStringList availableServiceProviders();

    QDeclarativeGeoServiceProviderRequirements *requirements() const;
    void setRequirements(QDeclarativeGeoServiceProviderRequirements *req);

    QStringList preferred() const;
    void setPreferred(const QStringList &val);

    QGeoServiceProvider *sharedGeoServiceProvider() const;

    Q_INVOKABLE bool supportsRouting(const RoutingFeatures &feature = AnyRoutingFeatures) const;
    Q_INVOKABLE bool supportsGeocoding(const GeocodingFeatures &feature = AnyGeocodingFeatures) const;
    Q_INVOKABLE bool supportsMapping(const MappingFeatures &feature = AnyMappingFeatures) const;
    Q_INVOKABLE bool supportsPlaces(const PlacesFeatures &feature = AnyPlacesFeatures) const;
    Q_REVISION(11) Q_INVOKABLE bool supportsNavigation(const NavigationFeature &feature = AnyNavigationFeatures) const;

    QStringList locales() const;
    void setLocales(const QStringList &locales);

    bool isAttached() const;

    void setAllowExperimental(bool allow);
    bool allowExperimental() const;

Q_SIGNALS:
    void nameChanged(const QString &name);
    void localesChanged();
    void attached();
    void preferredChanged(const QStringList &preferences);
    void allowExperimentalChanged(bool allow);

private:
    bool parametersReady();
    void tryAttach();
    static void parameter_append(QQmlListProperty<QDeclarativePluginParameter> *prop, QDeclarativePluginParameter *mapObject);
    static int parameter_count(QQmlListProperty<QDeclarativePluginParameter> *prop);
    static QDeclarativePluginParameter *parameter_at(QQmlListProperty<QDeclarativePluginParameter> *prop, int index);
    static void parameter_clear(QQmlListProperty<QDeclarativePluginParameter> *prop);

    QGeoServiceProvider *sharedProvider_;
    QString name_;
    QList<QDeclarativePluginParameter *> parameters_;
    QDeclarativeGeoServiceProviderRequirements *required_;
    bool complete_;
    bool experimental_;
    QStringList locales_;
    QStringList prefer_;
    Q_DISABLE_COPY(QDeclarativeGeoServiceProvider)
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeGeoServiceProviderRequirements : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QDeclarativeGeoServiceProvider::MappingFeatures mapping
               READ mappingRequirements WRITE setMappingRequirements
               NOTIFY mappingRequirementsChanged)
    Q_PROPERTY(QDeclarativeGeoServiceProvider::RoutingFeatures routing
               READ routingRequirements WRITE setRoutingRequirements
               NOTIFY routingRequirementsChanged)
    Q_PROPERTY(QDeclarativeGeoServiceProvider::GeocodingFeatures geocoding
               READ geocodingRequirements WRITE setGeocodingRequirements
               NOTIFY geocodingRequirementsChanged)
    Q_PROPERTY(QDeclarativeGeoServiceProvider::PlacesFeatures places
               READ placesRequirements WRITE setPlacesRequirements
               NOTIFY placesRequirementsChanged)
    Q_PROPERTY(QDeclarativeGeoServiceProvider::NavigationFeatures navigation
               READ navigationRequirements WRITE setNavigationRequirements
               NOTIFY navigationRequirementsChanged)

public:
    explicit QDeclarativeGeoServiceProviderRequirements(QObject *parent = 0);
    ~QDeclarativeGeoServiceProviderRequirements();

    QDeclarativeGeoServiceProvider::MappingFeatures mappingRequirements() const;
    void setMappingRequirements(const QDeclarativeGeoServiceProvider::MappingFeatures &features);

    QDeclarativeGeoServiceProvider::RoutingFeatures routingRequirements() const;
    void setRoutingRequirements(const QDeclarativeGeoServiceProvider::RoutingFeatures &features);

    QDeclarativeGeoServiceProvider::GeocodingFeatures geocodingRequirements() const;
    void setGeocodingRequirements(const QDeclarativeGeoServiceProvider::GeocodingFeatures &features);

    QDeclarativeGeoServiceProvider::PlacesFeatures placesRequirements() const;
    void setPlacesRequirements(const QDeclarativeGeoServiceProvider::PlacesFeatures &features);

    QDeclarativeGeoServiceProvider::NavigationFeatures navigationRequirements() const;
    void setNavigationRequirements(const QDeclarativeGeoServiceProvider::NavigationFeatures &features);

    Q_INVOKABLE bool matches(const QGeoServiceProvider *provider) const;

    bool operator == (const QDeclarativeGeoServiceProviderRequirements &rhs) const;

Q_SIGNALS:
    void mappingRequirementsChanged(const QDeclarativeGeoServiceProvider::MappingFeatures &features);
    void routingRequirementsChanged(const QDeclarativeGeoServiceProvider::RoutingFeatures &features);
    void geocodingRequirementsChanged(const QDeclarativeGeoServiceProvider::GeocodingFeatures &features);
    void placesRequirementsChanged(const QDeclarativeGeoServiceProvider::PlacesFeatures &features);
    void navigationRequirementsChanged(const QDeclarativeGeoServiceProvider::NavigationFeatures &features);

    void requirementsChanged();

private:
    QDeclarativeGeoServiceProvider::MappingFeatures mapping_;
    QDeclarativeGeoServiceProvider::RoutingFeatures routing_;
    QDeclarativeGeoServiceProvider::GeocodingFeatures geocoding_;
    QDeclarativeGeoServiceProvider::PlacesFeatures places_;
    QDeclarativeGeoServiceProvider::NavigationFeatures navigation_;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativeGeoServiceProviderRequirements)
QML_DECLARE_TYPE(QDeclarativeGeoServiceProvider)

#endif
