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

#ifndef QGEOSERVICEPROVIDER_P_H
#define QGEOSERVICEPROVIDER_P_H

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

#include "qgeoserviceprovider.h"

#include <QHash>
#include <QJsonObject>
#include <QJsonArray>
#include <QLocale>

QT_BEGIN_NAMESPACE

class QGeoCodingManager;
class QGeoRoutingManager;
class QGeoMappingManager;

class QGeoServiceProviderFactory;
class QGeoServiceProviderFactoryV2;
class QGeoServiceProviderFactoryV3;
class QQmlEngine;

class QGeoServiceProviderPrivate
{
public:
    QGeoServiceProviderPrivate();
    ~QGeoServiceProviderPrivate();

    void loadMeta();
    void loadPlugin(const QVariantMap &parameters);
    void unload();
    void filterParameterMap();

    /* helper templates for generating the feature and manager accessors */
    template <class Manager, class Engine>
    Manager *manager(QGeoServiceProvider::Error *error,
                     QString *errorString, Manager **manager);
    template <class Flags>
    Flags features(const char *enumName);

    QGeoServiceProviderFactory *factory;
    QGeoServiceProviderFactoryV2 *factoryV2 = nullptr;
    QGeoServiceProviderFactoryV3 *factoryV3 = nullptr;
    QJsonObject metaData;

    QVariantMap parameterMap;
    QVariantMap cleanedParameterMap;

    bool experimental;

    QGeoCodingManager *geocodingManager;
    QGeoRoutingManager *routingManager;
    QGeoMappingManager *mappingManager;
    QPlaceManager *placeManager;
    QNavigationManager *navigationManager = nullptr;
    QQmlEngine *qmlEngine = nullptr;

    QGeoServiceProvider::Error geocodeError;
    QGeoServiceProvider::Error routingError;
    QGeoServiceProvider::Error mappingError;
    QGeoServiceProvider::Error placeError;
    QGeoServiceProvider::Error navigationError = QGeoServiceProvider::NoError;

    QString geocodeErrorString;
    QString routingErrorString;
    QString mappingErrorString;
    QString placeErrorString;
    QString navigationErrorString;

    QGeoServiceProvider::Error error;
    QString errorString;

    QString providerName;

    QLocale locale;
    bool localeSet;

    static QMultiHash<QString, QJsonObject> plugins(bool reload = false);
    static void loadPluginMetadata(QMultiHash<QString, QJsonObject> &list);
};

QT_END_NAMESPACE

#endif
