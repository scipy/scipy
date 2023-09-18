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

#ifndef QGEOCIRCLE_H
#define QGEOCIRCLE_H

#include <QtPositioning/QGeoRectangle>

QT_BEGIN_NAMESPACE

class QGeoCoordinate;
class QGeoCirclePrivate;

class Q_POSITIONING_EXPORT QGeoCircle : public QGeoShape
{
    Q_GADGET
    Q_PROPERTY(QGeoCoordinate center READ center WRITE setCenter)
    Q_PROPERTY(qreal radius READ radius WRITE setRadius)

public:
    QGeoCircle();
    QGeoCircle(const QGeoCoordinate &center, qreal radius = -1.0);
    QGeoCircle(const QGeoCircle &other);
    QGeoCircle(const QGeoShape &other);

    ~QGeoCircle();

    QGeoCircle &operator=(const QGeoCircle &other);

    using QGeoShape::operator==;
    bool operator==(const QGeoCircle &other) const;

    using QGeoShape::operator!=;
    bool operator!=(const QGeoCircle &other) const;

    void setCenter(const QGeoCoordinate &center);
    QGeoCoordinate center() const;

    void setRadius(qreal radius);
    qreal radius() const;

    Q_INVOKABLE void translate(double degreesLatitude, double degreesLongitude);
    Q_INVOKABLE QGeoCircle translated(double degreesLatitude, double degreesLongitude) const;
    Q_INVOKABLE void extendCircle(const QGeoCoordinate &coordinate);

    Q_INVOKABLE QString toString() const;

private:
    inline QGeoCirclePrivate *d_func();
    inline const QGeoCirclePrivate *d_func() const;
};

Q_DECLARE_TYPEINFO(QGeoCircle, Q_MOVABLE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QGeoCircle)

#endif

