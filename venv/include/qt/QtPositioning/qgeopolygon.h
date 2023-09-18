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

#ifndef QGEOPOLYGON_H
#define QGEOPOLYGON_H

#include <QtPositioning/QGeoShape>
#include <QtCore/QVariantList>

QT_BEGIN_NAMESPACE

class QGeoCoordinate;
class QGeoPolygonPrivate;

class Q_POSITIONING_EXPORT QGeoPolygon : public QGeoShape
{
    Q_GADGET
    Q_PROPERTY(QVariantList perimeter READ perimeter WRITE setPerimeter REVISION 12)

public:
    QGeoPolygon();
    QGeoPolygon(const QList<QGeoCoordinate> &path);
    QGeoPolygon(const QGeoPolygon &other);
    QGeoPolygon(const QGeoShape &other);

    ~QGeoPolygon();

    QGeoPolygon &operator=(const QGeoPolygon &other);

    using QGeoShape::operator==;
    bool operator==(const QGeoPolygon &other) const;

    using QGeoShape::operator!=;
    bool operator!=(const QGeoPolygon &other) const;

    void setPath(const QList<QGeoCoordinate> &path); // ### Qt6: rename into setPerimeter
    const QList<QGeoCoordinate> &path() const;

    Q_INVOKABLE void addHole(const QVariant &holePath);
                void addHole(const QList<QGeoCoordinate> &holePath);
    Q_INVOKABLE const QVariantList hole(int index) const;
                const QList<QGeoCoordinate> holePath(int index) const;
    Q_INVOKABLE void removeHole(int index);
    Q_INVOKABLE int holesCount() const;
    Q_INVOKABLE void translate(double degreesLatitude, double degreesLongitude);
    Q_INVOKABLE QGeoPolygon translated(double degreesLatitude, double degreesLongitude) const;
    Q_INVOKABLE double length(int indexFrom = 0, int indexTo = -1) const;
    Q_INVOKABLE int size() const;
    Q_INVOKABLE void addCoordinate(const QGeoCoordinate &coordinate);
    Q_INVOKABLE void insertCoordinate(int index, const QGeoCoordinate &coordinate);
    Q_INVOKABLE void replaceCoordinate(int index, const QGeoCoordinate &coordinate);
    Q_INVOKABLE QGeoCoordinate coordinateAt(int index) const;
    Q_INVOKABLE bool containsCoordinate(const QGeoCoordinate &coordinate) const;
    Q_INVOKABLE void removeCoordinate(const QGeoCoordinate &coordinate);
    Q_INVOKABLE void removeCoordinate(int index);

    Q_INVOKABLE QString toString() const;

protected:
    void setPerimeter(const QVariantList &path);
    QVariantList perimeter() const;

private:
    inline QGeoPolygonPrivate *d_func();
    inline const QGeoPolygonPrivate *d_func() const;
};

Q_DECLARE_TYPEINFO(QGeoPolygon, Q_MOVABLE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QGeoPolygon)

#endif // QGEOPOLYGON_H
