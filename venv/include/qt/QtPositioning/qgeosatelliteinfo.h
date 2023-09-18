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
#ifndef QGEOSATELLITEINFO_H
#define QGEOSATELLITEINFO_H

#include <QtPositioning/qpositioningglobal.h>

QT_BEGIN_NAMESPACE

class QDebug;
class QDataStream;

class QGeoSatelliteInfoPrivate;
class Q_POSITIONING_EXPORT QGeoSatelliteInfo
{
public:
    enum Attribute {
        Elevation,
        Azimuth
    };

    enum SatelliteSystem {
        Undefined = 0x00,
        GPS = 0x01,
        GLONASS = 0x02
    };

    QGeoSatelliteInfo();
    QGeoSatelliteInfo(const QGeoSatelliteInfo &other);
    QGeoSatelliteInfo(QGeoSatelliteInfoPrivate &dd);
    ~QGeoSatelliteInfo();

    QGeoSatelliteInfo &operator=(const QGeoSatelliteInfo &other);

    bool operator==(const QGeoSatelliteInfo &other) const;
    inline bool operator!=(const QGeoSatelliteInfo &other) const {
        return !operator==(other);
    }

    void setSatelliteSystem(SatelliteSystem system);
    SatelliteSystem satelliteSystem() const;

    void setSatelliteIdentifier(int satId);
    int satelliteIdentifier() const;

    void setSignalStrength(int signalStrength);
    int signalStrength() const;

    void setAttribute(Attribute attribute, qreal value);
    qreal attribute(Attribute attribute) const;
    void removeAttribute(Attribute attribute);

    bool hasAttribute(Attribute attribute) const;

private:
#ifndef QT_NO_DEBUG_STREAM
    friend Q_POSITIONING_EXPORT QDebug operator<<(QDebug dbg, const QGeoSatelliteInfo &info);
#endif
#ifndef QT_NO_DATASTREAM
    friend Q_POSITIONING_EXPORT QDataStream &operator<<(QDataStream &stream, const QGeoSatelliteInfo &info);
    friend Q_POSITIONING_EXPORT QDataStream &operator>>(QDataStream &stream, QGeoSatelliteInfo &info);
#endif
    QGeoSatelliteInfoPrivate *d;
    friend class QGeoSatelliteInfoPrivate;
};

#ifndef QT_NO_DEBUG_STREAM
Q_POSITIONING_EXPORT QDebug operator<<(QDebug dbg, const QGeoSatelliteInfo &info);
#endif

#ifndef QT_NO_DATASTREAM
Q_POSITIONING_EXPORT QDataStream &operator<<(QDataStream &stream, const QGeoSatelliteInfo &info);
Q_POSITIONING_EXPORT QDataStream &operator>>(QDataStream &stream, QGeoSatelliteInfo &info);
#endif

QT_END_NAMESPACE

#endif
