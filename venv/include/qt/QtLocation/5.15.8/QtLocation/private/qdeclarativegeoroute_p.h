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

#ifndef QDECLARATIVEGEOROUTE_H
#define QDECLARATIVEGEOROUTE_H

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
#include <QtLocation/private/qdeclarativegeoroutesegment_p.h>

#include <QtCore/QObject>
#include <QtQml/QQmlListProperty>
#include <QtLocation/QGeoRoute>

QT_BEGIN_NAMESPACE
class QDeclarativeGeoRouteQuery;

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeGeoRoute : public QObject
{
    Q_OBJECT

    Q_PROPERTY(QGeoRectangle bounds READ bounds CONSTANT)
    Q_PROPERTY(int travelTime READ travelTime CONSTANT)
    Q_PROPERTY(qreal distance READ distance CONSTANT)
    Q_PROPERTY(QJSValue path READ path WRITE setPath NOTIFY pathChanged)
    Q_PROPERTY(QQmlListProperty<QDeclarativeGeoRouteSegment> segments READ segments CONSTANT)
    Q_PROPERTY(QDeclarativeGeoRouteQuery *routeQuery READ routeQuery REVISION 11)
    Q_PROPERTY(QList<QObject *> legs READ legs CONSTANT REVISION 12)
    Q_PROPERTY(QObject *extendedAttributes READ extendedAttributes CONSTANT REVISION 13)

public:
    explicit QDeclarativeGeoRoute(QObject *parent = 0);
    QDeclarativeGeoRoute(const QGeoRoute &route, QObject *parent = 0);
    ~QDeclarativeGeoRoute();

    QGeoRectangle bounds() const;
    int travelTime() const;
    qreal distance() const;

    QJSValue path() const;
    void setPath(const QJSValue &value);

    QQmlListProperty<QDeclarativeGeoRouteSegment> segments();

    void appendSegment(QDeclarativeGeoRouteSegment *segment);
    void clearSegments();

    int segmentsCount() const;
    const QGeoRoute &route() const;
    QDeclarativeGeoRouteQuery *routeQuery();
    QList<QObject *> legs();
    QQmlPropertyMap *extendedAttributes() const;

    Q_INVOKABLE bool equals(QDeclarativeGeoRoute *other) const;

Q_SIGNALS:
    void pathChanged();

private:
    static void segments_append(QQmlListProperty<QDeclarativeGeoRouteSegment> *prop, QDeclarativeGeoRouteSegment *segment);
    static int segments_count(QQmlListProperty<QDeclarativeGeoRouteSegment> *prop);
    static QDeclarativeGeoRouteSegment *segments_at(QQmlListProperty<QDeclarativeGeoRouteSegment> *prop, int index);
    static void segments_clear(QQmlListProperty<QDeclarativeGeoRouteSegment> *prop);

    void initSegments(unsigned int lastIndex = ~0U);
    QList<QGeoCoordinate> routePath();

    QGeoRoute route_;
    QDeclarativeGeoRouteQuery *routeQuery_ = nullptr;
    QList<QDeclarativeGeoRouteSegment *> segments_;
    QList<QObject *> legs_;
    bool segmentsDirty_ = true;
    QQmlPropertyMap *m_extendedAttributes = nullptr;

    friend class QDeclarativeRouteMapItem;
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeGeoRouteLeg : public QDeclarativeGeoRoute
{
    Q_OBJECT

    Q_PROPERTY(int legIndex READ legIndex CONSTANT)
    Q_PROPERTY(QObject * overallRoute READ overallRoute CONSTANT)

public:
    explicit QDeclarativeGeoRouteLeg(QObject *parent = nullptr);
    QDeclarativeGeoRouteLeg(const QGeoRouteLeg &routeLeg, QObject *parent = nullptr);
    ~QDeclarativeGeoRouteLeg() override;

    int legIndex() const;
    QObject *overallRoute() const;

private:
    QGeoRouteLeg m_routeLeg;
};

QT_END_NAMESPACE

#endif
