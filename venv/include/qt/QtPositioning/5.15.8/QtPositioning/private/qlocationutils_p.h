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
#ifndef QLOCATIONUTILS_P_H
#define QLOCATIONUTILS_P_H

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

#include <QtCore/QtGlobal>
#include <math.h> // needed for non-std:: versions of functions
#include <qmath.h>
#include <QtPositioning/QGeoCoordinate>
#include <QtPositioning/private/qpositioningglobal_p.h>

static const double offsetEpsilon = 1e-12; // = 0.000000000001
static const double leftOffset = -180.0 + offsetEpsilon;
static const double rightOffset = 180.0 - offsetEpsilon;

QT_BEGIN_NAMESPACE
class QTime;
class QByteArray;

class QGeoPositionInfo;
class QGeoSatelliteInfo;
class Q_POSITIONING_PRIVATE_EXPORT QLocationUtils
{
public:
    enum CardinalDirection {
        CardinalN,
        CardinalE,
        CardinalS,
        CardinalW,
        CardinalNE,
        CardinalSE,
        CardinalSW,
        CardinalNW,
        CardinalNNE,
        CardinalENE,
        CardinalESE,
        CardinalSSE,
        CardinalSSW,
        CardinalWSW,
        CardinalWNW,
        CardinalNNW
    };

    enum NmeaSentence {
        NmeaSentenceInvalid,
        NmeaSentenceGGA, // Fix information
        NmeaSentenceGSA, // Overall Satellite data, such as HDOP and VDOP
        NmeaSentenceGLL, // Lat/Lon data
        NmeaSentenceRMC, // Recommended minimum data for gps
        NmeaSentenceVTG, // Vector track an Speed over the Ground
        NmeaSentenceZDA, // Date and Time
        NmeaSentenceGSV  // Per-Satellite Info
    };

    inline static bool isValidLat(double lat) {
        return lat >= -90.0 && lat <= 90.0;
    }
    inline static bool isValidLong(double lng) {
        return lng >= -180.0 && lng <= 180.0;
    }

    inline static double clipLat(double lat, double clipValue = 90.0) {
        if (lat > clipValue)
            lat = clipValue;
        else if (lat < -clipValue)
            lat = -clipValue;
        return lat;
    }

    inline static double wrapLong(double lng) {
        if (lng > 180.0)
            lng -= 360.0;
        else if (lng < -180.0)
            lng += 360.0;
        return lng;
    }

    inline static CardinalDirection azimuthToCardinalDirection4(double azimuth)
    {
        azimuth = fmod(azimuth, 360.0);
        if (azimuth < 45.0 || azimuth > 315.0 )
            return CardinalN;
        else if (azimuth < 135.0)
            return CardinalE;
        else if (azimuth < 225.0)
            return CardinalS;
        else
            return CardinalW;
    }

    inline static CardinalDirection azimuthToCardinalDirection8(double azimuth)
    {
        azimuth = fmod(azimuth, 360.0);
        if (azimuth < 22.5 || azimuth > 337.5 )
            return CardinalN;
        else if (azimuth < 67.5)
            return CardinalNE;
        else if (azimuth < 112.5)
            return CardinalE;
        else if (azimuth < 157.5)
            return CardinalSE;
        else if (azimuth < 202.5)
            return CardinalS;

        else if (azimuth < 247.5)
            return CardinalSW;
        else if (azimuth < 292.5)
            return CardinalW;
        else
            return CardinalNW;
    }

    inline static CardinalDirection azimuthToCardinalDirection16(double azimuth)
    {
        azimuth = fmod(azimuth, 360.0);
        if (azimuth < 11.5 || azimuth > 348.75 )
            return CardinalN;
        else if (azimuth < 33.75)
            return CardinalNNE;
        else if (azimuth < 56.25)
            return CardinalNE;
        else if (azimuth < 78.75)
            return CardinalENE;
        else if (azimuth < 101.25)
            return CardinalE;
        else if (azimuth < 123.75)
            return CardinalESE;
        else if (azimuth < 146.25)
            return CardinalSE;
        else if (azimuth < 168.75)
            return CardinalSSE;
        else if (azimuth < 191.25)
            return CardinalS;

        else if (azimuth < 213.75)
            return CardinalSSW;
        else if (azimuth < 236.25)
            return CardinalSW;
        else if (azimuth < 258.75)
            return CardinalWSW;
        else if (azimuth < 281.25)
            return CardinalW;
        else if (azimuth < 303.75)
            return CardinalWNW;
        else if (azimuth < 326.25)
            return CardinalNW;
        else
            return CardinalNNW;
    }

    // For values exceeding +- 720.0
    inline static double wrapLongExt(double lng) {
        double remainder = fmod(lng + 180.0, 360.0);
        return fmod(remainder + 360.0, 360.0) - 180.0;
    }

    // Mirrors the azimuth against the X axis. Azimuth assumed to be in [0,360[
    inline static double mirrorAzimuthX(double azimuth) {
        if (azimuth <= 90.0)
             return 180.0 - azimuth;
         else
             return 180.0 + (360.0 - azimuth);
    }

    // Mirrors the azimuth against the Y axis. Azimuth assumed to be in [0,360[
    inline static double mirrorAzimuthY(double azimuth) {
        if (azimuth == 0.0)
            return 0.0;
        return 360.0 - azimuth;
    }

    inline static double radians(double degrees)
    {
        return qDegreesToRadians(degrees);
    }

    inline static double degrees(double radians)
    {
        return qRadiansToDegrees(radians);
    }

    inline static double earthMeanRadius()
    {
        return 6371007.2;
    }

    inline static double earthMeanCircumference()
    {
        return earthMeanRadius() * 2.0 * M_PI;
    }

    inline static double mercatorMaxLatitude()
    {
        return 85.05113;
    }

    inline static QGeoCoordinate antipodalPoint(const QGeoCoordinate &p)
    {
        return QGeoCoordinate(-p.latitude(), wrapLong(p.longitude() + 180.0));
    }

    // Leftmost longitude before wrapping kicks in
    inline static double mapLeftLongitude(double centerLongitude)
    {
        return wrapLong(centerLongitude + leftOffset);
    }

    // Rightmost longitude before wrapping kicks in
    inline static double  mapRightLongitude(double centerLongitude)
    {
        return wrapLong(centerLongitude - leftOffset);
    }

    inline static void split_double(double input, float *hipart, float *lopart)
    {
        *hipart = (float) input;
        double delta = input - ((double) *hipart);
        *lopart = (float) delta;
    }

    static qreal metersPerPixel(qreal zoomLevel, const QGeoCoordinate &coordinate)
    {
        const qreal metersPerTile = earthMeanCircumference() * std::cos(radians(coordinate.latitude())) / std::pow(2, zoomLevel);
        return metersPerTile / 256.0;
    }

    /*
        returns the NMEA sentence type.
    */
    static NmeaSentence getNmeaSentenceType(const char *data, int size);

    /*
        Creates a QGeoPositionInfo from a GGA, GLL, RMC, VTG or ZDA sentence.

        Note:
        - GGA and GLL sentences have time but not date so the update's
          QDateTime object will have an invalid date.
        - RMC reports date with a two-digit year so in this case the year
          is assumed to be after the year 2000.
    */
    static bool getPosInfoFromNmea(const char *data,
                                   int size,
                                   QGeoPositionInfo *info, double uere,
                                   bool *hasFix = nullptr);

    /*
        Retruns a list of QGeoSatelliteInfo in the view.

        Note: this function has to be called repeatedly until it returns true.
        Reason being that GSV sentences can be split into multiple samples, so getting the full data
        requires parsing multiple sentences.
     */
    enum GSVParseStatus {
        GSVNotParsed,
        GSVPartiallyParsed,
        GSVFullyParsed
    };
    static GSVParseStatus getSatInfoFromNmea(const char *data,
                                   int size,
                                   QList<QGeoSatelliteInfo> &infos);

    /*
        Parses GSA for satellites in use.
     */
    static bool getSatInUseFromNmea(const char *data,
                                    int size,
                                    QList<int> &pnrsInUse);

    /*
        Returns true if the given NMEA sentence has a valid checksum.
    */
    static bool hasValidNmeaChecksum(const char *data, int size);

    /*
        Returns time from a string in hhmmss or hhmmss.z+ format.
    */
    static bool getNmeaTime(const QByteArray &bytes, QTime *time);

    /*
        Accepts for example ("2734.7964", 'S', "15306.0124", 'E') and returns the
        lat-long values. Fails if lat or long fail isValidLat() or isValidLong().
    */
    static bool getNmeaLatLong(const QByteArray &latString,
                               char latDirection,
                               const QByteArray &lngString,
                               char lngDirection,
                               double *lat,
                               double *lon);
};

QT_END_NAMESPACE

#endif
