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

#ifndef QGEOROUTE_P_H
#define QGEOROUTE_P_H

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
#include "qgeoroute.h"
#include "qgeorouterequest.h"
#include "qgeorectangle.h"
#include "qgeoroutesegment.h"

#include <QSharedData>
#include <QVariantMap>
#include <QScopedPointer>

QT_BEGIN_NAMESPACE

class QGeoCoordinate;

class Q_LOCATION_PRIVATE_EXPORT QGeoRoutePrivate : public QSharedData
{
public:
    QGeoRoutePrivate();
    QGeoRoutePrivate(const QGeoRoutePrivate &other);
    virtual ~QGeoRoutePrivate();
    virtual QGeoRoutePrivate *clone() = 0;

    bool operator == (const QGeoRoutePrivate &other) const;

    virtual void setId(const QString &id);
    virtual QString id() const;

    virtual void setRequest(const QGeoRouteRequest &request);
    virtual QGeoRouteRequest request() const;

    virtual void setBounds(const QGeoRectangle &bounds);
    virtual QGeoRectangle bounds() const;

    virtual void setTravelTime(int travelTime);
    virtual int travelTime() const;

    virtual void setDistance(qreal distance);
    virtual qreal distance() const;

    virtual void setTravelMode(QGeoRouteRequest::TravelMode mode);
    virtual QGeoRouteRequest::TravelMode travelMode() const;

    virtual void setPath(const QList<QGeoCoordinate> &path);
    virtual QList<QGeoCoordinate> path() const;

    virtual void setFirstSegment(const QGeoRouteSegment &firstSegment);
    virtual QGeoRouteSegment firstSegment() const;

    virtual QVariantMap metadata() const;

    virtual void setRouteLegs(const QList<QGeoRouteLeg> &legs);
    virtual QList<QGeoRouteLeg> routeLegs() const;

    virtual void setExtendedAttributes(const QVariantMap &extendedAttributes);
    virtual QVariantMap extendedAttributes() const;

    virtual QString engineName() const = 0;
    virtual int segmentsCount() const = 0;

    // QGeoRouteLeg API
    virtual void setLegIndex(int idx);
    virtual int legIndex() const;
    virtual void setContainingRoute(const QGeoRoute &route);
    virtual QGeoRoute containingRoute() const;

    static const QGeoRoutePrivate *routePrivateData(const QGeoRoute &route);

protected:
    virtual bool equals(const QGeoRoutePrivate &other) const;
};

class Q_LOCATION_PRIVATE_EXPORT  QGeoRoutePrivateDefault : public QGeoRoutePrivate
{
public:
    QGeoRoutePrivateDefault();
    QGeoRoutePrivateDefault(const QGeoRoutePrivateDefault &other);
    ~QGeoRoutePrivateDefault() override;
    virtual QGeoRoutePrivate *clone() override;

    virtual void setId(const QString &id) override;
    virtual QString id() const override;

    virtual void setRequest(const QGeoRouteRequest &request) override;
    virtual QGeoRouteRequest request() const override;

    virtual void setBounds(const QGeoRectangle &bounds) override;
    virtual QGeoRectangle bounds() const override;

    virtual void setTravelTime(int travelTime) override;
    virtual int travelTime() const override;

    virtual void setDistance(qreal distance) override;
    virtual qreal distance() const override;

    virtual void setTravelMode(QGeoRouteRequest::TravelMode mode) override;
    virtual QGeoRouteRequest::TravelMode travelMode() const override;

    virtual void setPath(const QList<QGeoCoordinate> &path) override;
    virtual QList<QGeoCoordinate> path() const override;

    virtual void setFirstSegment(const QGeoRouteSegment &firstSegment) override;
    virtual QGeoRouteSegment firstSegment() const override;

    virtual QString engineName() const override;
    virtual int segmentsCount() const override;

    virtual void setRouteLegs(const QList<QGeoRouteLeg> &legs) override;
    virtual QList<QGeoRouteLeg> routeLegs() const override;

    void setExtendedAttributes(const QVariantMap &extendedAttributes) override;
    QVariantMap extendedAttributes() const override;

    // QGeoRouteLeg API
    virtual void setLegIndex(int idx) override;
    virtual int legIndex() const override;
    virtual void setContainingRoute(const QGeoRoute &route) override;
    virtual QGeoRoute containingRoute() const override;

    QString m_id;
    QGeoRouteRequest m_request;

    QGeoRectangle m_bounds;
    mutable QList<QGeoRouteSegment> m_routeSegments;

    int m_travelTime;
    qreal m_distance;

    QGeoRouteRequest::TravelMode m_travelMode;

    QList<QGeoCoordinate> m_path;
    QList<QGeoRouteLeg> m_legs;
    QGeoRouteSegment m_firstSegment;
    mutable int m_numSegments;
    QScopedPointer<QGeoRoute> m_containingRoute;
    QVariantMap m_extendedAttributes;
    int m_legIndex = 0;
};

QT_END_NAMESPACE

#endif
