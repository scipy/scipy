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

#ifndef QGEOROUTINGMANAGER_H
#define QGEOROUTINGMANAGER_H

#include <QtCore/QObject>
#include <QtCore/QLocale>
#include <QtLocation/QGeoRouteRequest>
#include <QtLocation/QGeoRouteReply>

QT_BEGIN_NAMESPACE

class QGeoRoutingManagerEngine;
class QGeoRoutingManagerPrivate;

class Q_LOCATION_EXPORT QGeoRoutingManager : public QObject
{
    Q_OBJECT

public:
    ~QGeoRoutingManager();

    QString managerName() const;
    int managerVersion() const;

    QGeoRouteReply *calculateRoute(const QGeoRouteRequest &request);
    QGeoRouteReply *updateRoute(const QGeoRoute &route, const QGeoCoordinate &position);

    QGeoRouteRequest::TravelModes supportedTravelModes() const;
    QGeoRouteRequest::FeatureTypes supportedFeatureTypes() const;
    QGeoRouteRequest::FeatureWeights supportedFeatureWeights() const;
    QGeoRouteRequest::RouteOptimizations supportedRouteOptimizations() const;
    QGeoRouteRequest::SegmentDetails supportedSegmentDetails() const;
    QGeoRouteRequest::ManeuverDetails supportedManeuverDetails() const;

    void setLocale(const QLocale &locale);
    QLocale locale() const;
    void setMeasurementSystem(QLocale::MeasurementSystem system);
    QLocale::MeasurementSystem measurementSystem() const;

Q_SIGNALS:
    void finished(QGeoRouteReply *reply);
    void error(QGeoRouteReply *reply, QGeoRouteReply::Error error, QString errorString = QString());

private:
    explicit QGeoRoutingManager(QGeoRoutingManagerEngine *engine, QObject *parent = nullptr);

    QGeoRoutingManagerPrivate *d_ptr;
    Q_DISABLE_COPY(QGeoRoutingManager)

    friend class QGeoServiceProvider;
    friend class QGeoServiceProviderPrivate;
};

QT_END_NAMESPACE

#endif
