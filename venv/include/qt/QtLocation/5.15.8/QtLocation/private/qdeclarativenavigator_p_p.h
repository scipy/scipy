/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
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

#ifndef QDECLARATIVENAVIGATOR_P_P_H
#define QDECLARATIVENAVIGATOR_P_P_H

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

#include <QtCore/qlist.h>
#include <QtLocation/private/qlocationglobal_p.h>
#include <QtCore/qpointer.h>
#include <QtLocation/qgeoroute.h>
#include <QtLocation/private/qdeclarativenavigator_p.h>
#include <QAbstractListModel>
#include <QtLocation/private/qdeclarativegeoroute_p.h>
#include <QtLocation/private/qdeclarativegeoroutemodel_p.h>

QT_BEGIN_NAMESPACE

class QDeclarativeGeoServiceProvider;
class QDeclarativeGeoMap;
class QNavigationManager;
class QDeclarativeGeoRoute;
class QDeclarativeGeoRouteLeg;
class QDeclarativePositionSource;
class QGeoMapParameter;
class QDeclarativeGeoRouteSegment;
class QParameterizableObject;
class QAbstractNavigator;

template<typename T, int Role>
class ReadOnlyListModel : public QAbstractListModel
{
public:
    explicit ReadOnlyListModel(const QByteArray &dataRoleName, QObject *parent = nullptr)
        : QAbstractListModel(parent)
    {
        m_roleNames.insert(Role, dataRoleName);
    }

    int rowCount(const QModelIndex &) const override
    {
        return m_data.size();
    }

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override
    {
        const int row = index.row();
        if (!index.isValid() || row < 0 || row >= m_data.size() || role != Role)
            return QVariant();

        return QVariant::fromValue(m_data.at(row));
    }

    QHash<int, QByteArray> roleNames() const override
    {
        return m_roleNames;
    }

    void updateData(const QList<T*> &data)
    {
        beginResetModel();
        qDeleteAll(m_data);
        m_data = data;
        endResetModel();
    }

protected:
    QHash<int, QByteArray> m_roleNames;
    QList<T*> m_data;
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeNavigationBasicDirections : public QObject
{
    Q_OBJECT

    Q_PROPERTY(QVariant nextManeuverIcon READ nextManeuverIcon NOTIFY nextManeuverIconChanged)
    Q_PROPERTY(qreal distanceToNextManeuver READ distanceToNextManeuver NOTIFY progressInformationChanged)
    Q_PROPERTY(qreal remainingTravelDistance READ remainingTravelDistance NOTIFY progressInformationChanged)
    Q_PROPERTY(qreal remainingTravelDistanceToNextWaypoint READ remainingTravelDistanceToNextWaypoint NOTIFY progressInformationChanged)
    Q_PROPERTY(qreal traveledDistance READ traveledDistance NOTIFY progressInformationChanged)
    Q_PROPERTY(int timeToNextManeuver READ timeToNextManeuver NOTIFY progressInformationChanged)
    Q_PROPERTY(int remainingTravelTime READ remainingTravelTime NOTIFY progressInformationChanged)
    Q_PROPERTY(int remainingTravelTimeToNextWaypoint READ remainingTravelTimeToNextWaypoint NOTIFY progressInformationChanged)
    Q_PROPERTY(int traveledTime READ traveledTime NOTIFY progressInformationChanged)
    Q_PROPERTY(QDeclarativeGeoRoute *currentRoute READ currentRoute NOTIFY currentRouteChanged)
    Q_PROPERTY(QDeclarativeGeoRouteLeg *currentRouteLeg READ currentRouteLeg NOTIFY currentRouteChanged)
    Q_PROPERTY(int currentSegment READ currentSegment NOTIFY currentSegmentChanged)
    Q_PROPERTY(QAbstractItemModel *alternativeRoutes READ alternativeRoutes CONSTANT)

public:
    explicit QDeclarativeNavigationBasicDirections(QDeclarativeNavigator *parent);

    QVariant nextManeuverIcon() const;
    qreal distanceToNextManeuver() const;
    qreal remainingTravelDistance() const;
    qreal remainingTravelDistanceToNextWaypoint() const;
    qreal traveledDistance() const;
    int timeToNextManeuver() const;
    int remainingTravelTime() const;
    int remainingTravelTimeToNextWaypoint() const;
    int traveledTime() const;

    QDeclarativeGeoRoute *currentRoute() const;
    QDeclarativeGeoRouteLeg *currentRouteLeg() const;
    int currentSegment() const;
    QAbstractItemModel *alternativeRoutes();

Q_SIGNALS:
    void progressInformationChanged();
    void nextManeuverIconChanged();
    void currentRouteChanged();
    void currentRouteLegChanged();
    void currentSegmentChanged();
    void waypointReached(const QDeclarativeGeoWaypoint *pos);
    void destinationReached();

protected slots:
    void onCurrentRouteChanged();
    void onCurrentRouteLegChanged();
    void onAlternativeRoutesChanged();

protected:
    QDeclarativeNavigator *m_navigator;
    QDeclarativeNavigatorPrivate *m_navigatorPrivate;
    QPointer<QDeclarativeGeoRoute> m_currentRoute;
    QPointer<QDeclarativeGeoRouteLeg> m_currentRouteLeg;
    ReadOnlyListModel<QDeclarativeGeoRoute, QDeclarativeGeoRouteModel::RouteRole> m_routes;

    friend class QDeclarativeNavigator;
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeNavigatorParams
{
public:
    QPointer<QDeclarativeGeoMap> m_map;
    QPointer<QDeclarativeGeoRoute> m_route;
    QGeoRoute m_geoRoute;
    QPointer<QDeclarativePositionSource> m_positionSource;
    QList<QPointer<QGeoMapParameter>> m_parameters;
    bool m_trackPositionSource = true;
    bool m_autoRerouting = true;
};

class QDeclarativeNavigatorPrivate
{
public:
    QDeclarativeNavigatorPrivate(QParameterizableObject *q_);

    QParameterizableObject *q = nullptr;
    QSharedPointer<QDeclarativeNavigatorParams> m_params;
    QScopedPointer<QAbstractNavigator> m_navigator;
    QDeclarativeGeoServiceProvider *m_plugin = nullptr;
    QDeclarativeNavigationBasicDirections m_basicDirections;

    bool m_active = false;
    bool m_completed = false;
    bool m_ready = false;
    QDeclarativeNavigator::NavigationError m_error = QDeclarativeNavigator::NoError;
    QString m_errorString;
};

QT_END_NAMESPACE

#endif // QDECLARATIVENAVIGATOR_P_P_H
