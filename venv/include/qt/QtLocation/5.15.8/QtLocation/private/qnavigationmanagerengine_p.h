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

#ifndef QNAVIGATIONMANAGERENGINE_H
#define QNAVIGATIONMANAGERENGINE_H

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
#include <QObject>
#include <QLocale>
#include <QGeoCoordinate>

QT_BEGIN_NAMESPACE

class QAbstractNavigatorPrivate;
class QGeoMap;
class QGeoMapParameter;
class QMapRouteObject;
class QGeoRoute;
class QGeoRouteLeg;
class QNavigationManager;
class QNavigationManagerEnginePrivate;
class QDeclarativeNavigatorParams;
class QDeclarativeGeoWaypoint;
class QDeclarativeGeoRouteLeg;
class QDeclarativeGeoRoute;

/*
    This class is not supposed to react on QDeclarativeNavigator properties changes.
    This class is meant to react only on start, stop and setTrackPosition.
    Upon start(), it is supposed to fetch all info from the QDeclarativeNavigatorParams that the engine is supposed
    to inject.
*/
class Q_LOCATION_PRIVATE_EXPORT QAbstractNavigator: public QObject
{
    Q_OBJECT
public:
    QAbstractNavigator(QObject *parent = nullptr);
    ~QAbstractNavigator() override;
    virtual void setLocale(const QLocale &locale);
    virtual QLocale locale() const;
    virtual void setMeasurementSystem(QLocale::MeasurementSystem system);
    virtual QLocale::MeasurementSystem measurementSystem() const;
    virtual bool active() const = 0;
    virtual bool ready() const = 0;

    virtual QVariant nextManeuverIcon() const;
    virtual double distanceToNextManeuver() const;
    virtual int timeToNextManeuver() const;
    virtual int remainingTravelTime() const;
    virtual double remainingTravelDistance() const;
    virtual int remainingTravelTimeToNextWaypoint() const;
    virtual double remainingTravelDistanceToNextWaypoint() const;
    virtual double traveledDistance() const;
    virtual int traveledTime() const;
    virtual QGeoRoute currentRoute() const;
    virtual QGeoRouteLeg currentRouteLeg() const;
    virtual QList<QGeoRoute> alternativeRoutes() const = 0;
    virtual int currentSegment() const;
    virtual void setAutomaticReroutingEnabled(bool autoRerouting) = 0;
    virtual bool automaticReroutingEnabled() const = 0; // configured via navigation params at construction time
    virtual bool isOnRoute() = 0;
    virtual void recalculateRoutes() = 0;

public slots:
    virtual bool start() = 0;
    virtual bool stop() = 0;
    virtual void setTrackPosition(bool trackPosition) = 0;

signals:
    // These must be emitted by the engine
    void activeChanged(bool active);
    void waypointReached(const QDeclarativeGeoWaypoint *pos);
    void destinationReached();
    void currentRouteChanged();
    void currentRouteLegChanged();
    void currentSegmentChanged();

    void nextManeuverIconChanged();
    void progressInformationChanged();
    void isOnRouteChanged();
    void alternativeRoutesChanged();

private:
    QScopedPointer<QAbstractNavigatorPrivate> d;
};

class Q_LOCATION_PRIVATE_EXPORT QNavigationManagerEngine : public QObject
{
    Q_OBJECT
public:
    explicit QNavigationManagerEngine(const QVariantMap &parameters, QObject *parent = nullptr);
    ~QNavigationManagerEngine() override;

    void setManagerName(const QString &name);
    QString managerName() const;
    void setManagerVersion(int version);
    int managerVersion() const;
    virtual void setLocale(const QLocale &locale);
    virtual QLocale locale() const;
    virtual void setMeasurementSystem(QLocale::MeasurementSystem system);
    virtual QLocale::MeasurementSystem measurementSystem() const;

    virtual bool isInitialized() const;
    virtual QAbstractNavigator *createNavigator(const QSharedPointer<QDeclarativeNavigatorParams> &navigator) = 0;

signals:
    void initialized();

protected:
    /*!
        Marks the engine as initialized. Subclasses of QGeoMappingManagerEngine are to
        call this method after performing implementation-specific initialization within
        the constructor.
    */
    virtual void engineInitialized();

    QScopedPointer<QNavigationManagerEnginePrivate> d;
};

QT_END_NAMESPACE

#endif // QNAVIGATIONMANAGERENGINE_H
