/****************************************************************************
**
** Qt adaptation of geosimplify.js, https://github.com/mapbox/geosimplify-js, (c) 2017, Mapbox
**
** Copyright (C) 2020 Paolo Angelelli <paolo.angelelli@gmail.com>
** Copyright (C) 2020 The Qt Company Ltd.
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

#ifndef QGEOSIMPLIFY_P_H
#define QGEOSIMPLIFY_P_H

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

#include <QtPositioning/QGeoCoordinate>
#include <QtPositioning/private/qdoublevector2d_p.h>
#include <QtPositioning/private/qwebmercator_p.h>

QT_BEGIN_NAMESPACE

class QGeoSimplify {
protected:
    // Distance between two points in metres
    static double getDist(const QGeoCoordinate &p1, const QGeoCoordinate &p2);

    // p, a and b are intended as "unwrapped" around the left bound
    static QDoubleVector2D closestPoint( const QDoubleVector2D &p,
                                  const QDoubleVector2D &a,
                                  const QDoubleVector2D &b);

    static QGeoCoordinate closestPoint( const QGeoCoordinate &pc,
                                 const QGeoCoordinate &ac,
                                 const QGeoCoordinate &bc,
                                 const double &leftBound);

    // Distance from a point to a segment (line between two points) in metres
    static double getSegDist(const QGeoCoordinate &pc,
                      const QGeoCoordinate &ac,
                      const QGeoCoordinate &bc,
                      const double &leftBound);

    // doublevectors Intended as wrapped
    static double getSegDist(const QDoubleVector2D &p,
                      const QDoubleVector2D &a,
                      const QDoubleVector2D &b,
                      const double &leftBound);

    static void simplifyDPStep(const QList<QGeoCoordinate> &points,
                        const double &leftBound,
                        int first,
                        int last,
                        double offsetTolerance,
                        QList<QGeoCoordinate> &simplified);

    static double getDist(QDoubleVector2D a,
                   QDoubleVector2D b,
                   const double &leftBound);

    static void simplifyDPStep(const QList<QDoubleVector2D> &points,
                        const double &leftBound,
                        int first,
                        int last,
                        double offsetTolerance,
                        QList<QDoubleVector2D> &simplified);

    static void simplifyDPStepZL(const QList<QDoubleVector2D> &points,
                        const double &leftBound,
                        int first,
                        int last,
                        int zoomLevel,
                        QList<QDoubleVector2D> &simplified);

    // simplification using Ramer-Douglas-Peucker algorithm
    static QList<QGeoCoordinate> simplifyDouglasPeucker(const QList<QGeoCoordinate> &points,
                                                 const double &leftBound,
                                                 double offsetTolerance);

    static QList<QDoubleVector2D> simplifyDouglasPeucker(const QList<QDoubleVector2D> &points,
                                                 const double &leftBound,
                                                 double offsetTolerance);

    static QList<QDoubleVector2D> simplifyDouglasPeuckerZL(const QList<QDoubleVector2D> &points,
                                                 const double &leftBound,
                                                 int zoomLevel);

public:
    /*
        offsetTolerance - how far outside the straight line a point
                          needs to be for it to be "kept"
    */
    static QList<QGeoCoordinate> geoSimplify(const QList<QGeoCoordinate> &points,
                                   const double &leftBound,
                                   double offsetTolerance); // in meters

    static QList<QDoubleVector2D> geoSimplify(const QList<QDoubleVector2D> &points,
                                   const double &leftBound,
                                   double offsetTolerance); // in meters

    // This overload tries to be adaptive in the offsetTolerance across latitudes,
    // and return a simplification adequate for the given zoomLevel.
    static QList<QDoubleVector2D> geoSimplifyZL(const QList<QDoubleVector2D> &points,
                                   const double &leftBound,
                                   int zoomLevel); // in meters
};

QT_END_NAMESPACE

#endif // QGEOSIMPLIFY_P_H
