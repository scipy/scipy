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

#ifndef QDECLARATIVEGEOMANEUVER_H
#define QDECLARATIVEGEOMANEUVER_H

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
#include <QtLocation/qgeomaneuver.h>
#include <QtQml/QQmlPropertyMap>
#include <QtPositioning/QGeoCoordinate>

#include <QObject>


QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeGeoManeuver : public QObject
{
    Q_OBJECT
    Q_ENUMS(Direction)

    Q_PROPERTY(bool valid READ valid CONSTANT)
    Q_PROPERTY(QGeoCoordinate position READ position CONSTANT)
    Q_PROPERTY(QString instructionText READ instructionText CONSTANT)
    Q_PROPERTY(Direction direction READ direction CONSTANT)
    Q_PROPERTY(int timeToNextInstruction READ timeToNextInstruction CONSTANT)
    Q_PROPERTY(qreal distanceToNextInstruction READ distanceToNextInstruction CONSTANT)
    Q_PROPERTY(QGeoCoordinate waypoint READ waypoint CONSTANT)
    Q_PROPERTY(bool waypointValid READ waypointValid CONSTANT)
    Q_PROPERTY(QObject *extendedAttributes READ extendedAttributes NOTIFY extendedAttributesChanged REVISION 11)

public:
    enum Direction {
        NoDirection = QGeoManeuver::NoDirection,
        DirectionForward = QGeoManeuver::DirectionForward,
        DirectionBearRight = QGeoManeuver::DirectionBearRight,
        DirectionLightRight = QGeoManeuver::DirectionLightRight,
        DirectionRight = QGeoManeuver::DirectionRight,
        DirectionHardRight = QGeoManeuver::DirectionHardRight,
        DirectionUTurnRight = QGeoManeuver::DirectionUTurnRight,
        DirectionUTurnLeft = QGeoManeuver::DirectionUTurnLeft,
        DirectionHardLeft = QGeoManeuver::DirectionHardLeft,
        DirectionLeft = QGeoManeuver::DirectionLeft,
        DirectionLightLeft = QGeoManeuver::DirectionLightLeft,
        DirectionBearLeft = QGeoManeuver::DirectionBearLeft
    };

    explicit QDeclarativeGeoManeuver(QObject *parent = 0);
    QDeclarativeGeoManeuver(const QGeoManeuver &maneuver, QObject *parent = 0);
    ~QDeclarativeGeoManeuver();

    bool valid() const;
    bool waypointValid() const;

    QGeoCoordinate position() const;
    QString instructionText() const;
    Direction direction() const;
    int timeToNextInstruction() const;
    qreal distanceToNextInstruction() const;
    QGeoCoordinate waypoint() const;
    QQmlPropertyMap *extendedAttributes() const;

Q_SIGNALS:
    void extendedAttributesChanged();   //in practice is never emitted since parameters cannot be re-assigned
                                        //the declaration is needed to avoid warnings about non-notifyable properties

private:
    QGeoManeuver maneuver_;
    QQmlPropertyMap *m_extendedAttributes = nullptr;
};

QT_END_NAMESPACE

#endif
