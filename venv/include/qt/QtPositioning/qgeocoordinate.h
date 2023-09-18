/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtPositioning module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPL3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl-3.0.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or (at your option) the GNU General
** Public license version 3 or any later version approved by the KDE Free
** Qt Foundation. The licenses are as published by the Free Software
** Foundation and appearing in the file LICENSE.GPL2 and LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-2.0.html and
** https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QGEOCOORDINATE_H
#define QGEOCOORDINATE_H

#include <QtCore/QMetaType>
#include <QtCore/QString>
#include <QtCore/QSharedDataPointer>
#include <QtPositioning/qpositioningglobal.h>

QT_BEGIN_NAMESPACE

class QDebug;
class QDataStream;

class QGeoCoordinatePrivate;
class Q_POSITIONING_EXPORT QGeoCoordinate
{
    Q_GADGET
    Q_ENUMS(CoordinateFormat)

    Q_PROPERTY(double latitude READ latitude WRITE setLatitude)
    Q_PROPERTY(double longitude READ longitude WRITE setLongitude)
    Q_PROPERTY(double altitude READ altitude WRITE setAltitude)
    Q_PROPERTY(bool isValid READ isValid)

public:

    enum CoordinateType {
        InvalidCoordinate,
        Coordinate2D,
        Coordinate3D
    };

    enum CoordinateFormat {
        Degrees,
        DegreesWithHemisphere,
        DegreesMinutes,
        DegreesMinutesWithHemisphere,
        DegreesMinutesSeconds,
        DegreesMinutesSecondsWithHemisphere
    };

    QGeoCoordinate();
    QGeoCoordinate(double latitude, double longitude);
    QGeoCoordinate(double latitude, double longitude, double altitude);
    QGeoCoordinate(const QGeoCoordinate &other);
    ~QGeoCoordinate();

    QGeoCoordinate &operator=(const QGeoCoordinate &other);

    bool operator==(const QGeoCoordinate &other) const;
    inline bool operator!=(const QGeoCoordinate &other) const {
        return !operator==(other);
    }

    bool isValid() const;
    CoordinateType type() const;

    void setLatitude(double latitude);
    double latitude() const;

    void setLongitude(double longitude);
    double longitude() const;

    void setAltitude(double altitude);
    double altitude() const;

    Q_INVOKABLE qreal distanceTo(const QGeoCoordinate &other) const;
    Q_INVOKABLE qreal azimuthTo(const QGeoCoordinate &other) const;

    Q_INVOKABLE QGeoCoordinate atDistanceAndAzimuth(qreal distance, qreal azimuth, qreal distanceUp = 0.0) const;

    Q_INVOKABLE QString toString(CoordinateFormat format = DegreesMinutesSecondsWithHemisphere) const;

private:
    QGeoCoordinate(QGeoCoordinatePrivate &dd);
    QSharedDataPointer<QGeoCoordinatePrivate> d;
    friend class QGeoCoordinatePrivate;
    friend class QQuickGeoCoordinateAnimation;
};

Q_DECLARE_TYPEINFO(QGeoCoordinate, Q_MOVABLE_TYPE);

#ifndef QT_NO_DEBUG_STREAM
Q_POSITIONING_EXPORT QDebug operator<<(QDebug, const QGeoCoordinate &);
#endif

Q_POSITIONING_EXPORT uint qHash(const QGeoCoordinate &coordinate, uint seed = 0);

#ifndef QT_NO_DATASTREAM
Q_POSITIONING_EXPORT QDataStream &operator<<(QDataStream &stream, const QGeoCoordinate &coordinate);
Q_POSITIONING_EXPORT QDataStream &operator>>(QDataStream &stream, QGeoCoordinate &coordinate);
#endif

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QGeoCoordinate)

#endif
