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

#ifndef QGEOSERVICEPROVIDERFACTORY_H
#define QGEOSERVICEPROVIDERFACTORY_H

#include <QtLocation/QGeoServiceProvider>

#include <QtCore/QtPlugin>
#include <QtCore/QMap>
#include <QtCore/QString>

QT_BEGIN_NAMESPACE
class QQmlEngine;

class Q_LOCATION_EXPORT QGeoServiceProviderFactory
{
public:
    virtual ~QGeoServiceProviderFactory() {}

    virtual QGeoCodingManagerEngine *createGeocodingManagerEngine(const QVariantMap &parameters,
            QGeoServiceProvider::Error *error,
            QString *errorString) const;
    virtual QGeoMappingManagerEngine *createMappingManagerEngine(const QVariantMap &parameters,
            QGeoServiceProvider::Error *error,
            QString *errorString) const;
    virtual QGeoRoutingManagerEngine *createRoutingManagerEngine(const QVariantMap &parameters,
            QGeoServiceProvider::Error *error,
            QString *errorString) const;
    virtual QPlaceManagerEngine *createPlaceManagerEngine(const QVariantMap &parameters,
            QGeoServiceProvider::Error *error,
            QString *errorString) const;
};

Q_DECLARE_INTERFACE(QGeoServiceProviderFactory,
                    "org.qt-project.qt.geoservice.serviceproviderfactory/5.0")

class Q_LOCATION_EXPORT QGeoServiceProviderFactoryV2 : public QGeoServiceProviderFactory
{
public:
    virtual ~QGeoServiceProviderFactoryV2() {}

    virtual QNavigationManagerEngine *createNavigationManagerEngine(const QVariantMap &parameters,
            QGeoServiceProvider::Error *error,
            QString *errorString) const;
};

// Although not actually used for constructing a specialized loader, this is required for
// casting a QObject * into QGeoServiceProviderFactoryV2 *
Q_DECLARE_INTERFACE(QGeoServiceProviderFactoryV2,
                    "org.qt-project.qt.geoservice.serviceproviderfactoryV2/5.0")

class Q_LOCATION_EXPORT QGeoServiceProviderFactoryV3 : public QGeoServiceProviderFactoryV2
{
public:
    virtual ~QGeoServiceProviderFactoryV3() {}

    virtual void setQmlEngine(QQmlEngine * engine);
};

// Although not actually used for constructing a specialized loader, this is required for
// casting a QObject * into QGeoServiceProviderFactoryV3 *
Q_DECLARE_INTERFACE(QGeoServiceProviderFactoryV3,
                    "org.qt-project.qt.geoservice.serviceproviderfactoryV3/5.0")

QT_END_NAMESPACE

#endif
