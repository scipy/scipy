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

#ifndef QGEOROUTE_H
#define QGEOROUTE_H

#include <QtPositioning/QGeoCoordinate>
#include <QtLocation/QGeoRouteRequest>

#include <QtCore/QExplicitlySharedDataPointer>
#include <QtCore/QList>
#include <QtCore/QHash>

QT_BEGIN_NAMESPACE

class QGeoRectangle;
class QGeoRouteSegment;

class QGeoRoutePrivate;
class QGeoRouteLeg;
class Q_LOCATION_EXPORT QGeoRoute
{
public:
    QGeoRoute();
    QGeoRoute(const QGeoRoute &other);
    ~QGeoRoute(); // ### Qt6: make this virtual

    QGeoRoute &operator = (const QGeoRoute &other);

    bool operator == (const QGeoRoute &other) const;
    bool operator != (const QGeoRoute &other) const;

    void setRouteId(const QString &id);
    QString routeId() const;

    void setRequest(const QGeoRouteRequest &request);
    QGeoRouteRequest request() const;

    void setBounds(const QGeoRectangle &bounds);
    QGeoRectangle bounds() const;

    void setFirstRouteSegment(const QGeoRouteSegment &routeSegment);
    QGeoRouteSegment firstRouteSegment() const;

    void setTravelTime(int secs);
    int travelTime() const;

    void setDistance(qreal distance);
    qreal distance() const;

    void setTravelMode(QGeoRouteRequest::TravelMode mode);
    QGeoRouteRequest::TravelMode travelMode() const;

    void setPath(const QList<QGeoCoordinate> &path);
    QList<QGeoCoordinate> path() const;

    void setRouteLegs(const QList<QGeoRouteLeg> &legs);
    QList<QGeoRouteLeg> routeLegs() const;

    void setExtendedAttributes(const QVariantMap &extendedAttributes);
    QVariantMap extendedAttributes() const;

protected:
    QGeoRoute(const QExplicitlySharedDataPointer<QGeoRoutePrivate> &dd);
    QExplicitlySharedDataPointer<QGeoRoutePrivate> &d();
    const QExplicitlySharedDataPointer<QGeoRoutePrivate> &const_d() const;

private:
    QExplicitlySharedDataPointer<QGeoRoutePrivate> d_ptr;
    friend class QDeclarativeGeoRoute;
    friend class QGeoRoutePrivate;
};

class Q_LOCATION_EXPORT QGeoRouteLeg: public QGeoRoute
{
public:
    QGeoRouteLeg();
    QGeoRouteLeg(const QGeoRouteLeg &other);
    QGeoRouteLeg &operator=(const QGeoRouteLeg &other) = default;
    ~QGeoRouteLeg();

    void setLegIndex(int idx);
    int legIndex() const;

    void setOverallRoute(const QGeoRoute &route);
    QGeoRoute overallRoute() const;

protected:
    QGeoRouteLeg(const QExplicitlySharedDataPointer<QGeoRoutePrivate> &dd);

    friend class QDeclarativeGeoRoute;
    friend class QGeoRoutePrivate;
};

QT_END_NAMESPACE

#endif
