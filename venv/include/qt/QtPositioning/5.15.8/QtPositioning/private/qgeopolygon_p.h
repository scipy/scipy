/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
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

#ifndef QGEOPOLYGON_P_H
#define QGEOPOLYGON_P_H

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

#include <QtPositioning/private/qgeopath_p.h>
#include <QtPositioning/qgeopolygon.h>
#include <QtPositioning/private/qclipperutils_p.h>

QT_BEGIN_NAMESPACE

class Q_POSITIONING_PRIVATE_EXPORT QGeoPolygonPrivate : public QGeoPathPrivate
{
public:
    QGeoPolygonPrivate();
    QGeoPolygonPrivate(const QList<QGeoCoordinate> &path);
    ~QGeoPolygonPrivate();

// QGeoShape API
    virtual QGeoShapePrivate *clone() const override;
    virtual bool isValid() const override;
    virtual bool contains(const QGeoCoordinate &coordinate) const override;
    virtual void translate(double degreesLatitude, double degreesLongitude) override;
    virtual bool operator==(const QGeoShapePrivate &other) const override;

// QGeoPath API
    virtual void markDirty() override;

// QGeoPolygonPrivate API
    int holesCount() const;
    bool polygonContains(const QGeoCoordinate &coordinate) const;
    const QList<QGeoCoordinate> holePath(int index) const;

    virtual void addHole(const QList<QGeoCoordinate> &holePath);
    virtual void removeHole(int index);
    virtual void updateClipperPath();

// data members
    bool m_clipperDirty = true;
    QList<QList<QGeoCoordinate>> m_holesList;
    QtClipperLib::Path m_clipperPath;
};

class Q_POSITIONING_PRIVATE_EXPORT QGeoPolygonPrivateEager : public QGeoPolygonPrivate
{
public:
    QGeoPolygonPrivateEager();
    QGeoPolygonPrivateEager(const QList<QGeoCoordinate> &path);
    ~QGeoPolygonPrivateEager();

// QGeoShape API
    virtual QGeoShapePrivate *clone() const override;
    virtual void translate(double degreesLatitude, double degreesLongitude) override;

// QGeoPath API
    virtual void markDirty() override;
    virtual void addCoordinate(const QGeoCoordinate &coordinate) override;
    virtual void computeBoundingBox() override;

// QGeoPolygonPrivate API

// *Eager API
    void updateBoundingBox();

// data members
    QVector<double> m_deltaXs;      // longitude deltas from m_path[0]
    double m_minX = 0;              // minimum value inside deltaXs
    double m_maxX = 0;              // maximum value inside deltaXs
    double m_minLati = 0;           // minimum latitude. paths do not wrap around through the poles
    double m_maxLati = 0;           // minimum latitude. paths do not wrap around through the poles
};

// This is a mean of creating a QGeoPolygonPrivateEager and injecting it into QGeoPolygons via operator=
class Q_POSITIONING_PRIVATE_EXPORT QGeoPolygonEager : public QGeoPolygon
{
    Q_GADGET
public:

    QGeoPolygonEager();
    QGeoPolygonEager(const QList<QGeoCoordinate> &path);
    QGeoPolygonEager(const QGeoPolygon &other);
    QGeoPolygonEager(const QGeoShape &other);
    ~QGeoPolygonEager();
};

QT_END_NAMESPACE

#endif // QGEOPOLYGON_P_H
