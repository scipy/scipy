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

#ifndef QDECLARATIVEGEOMAPTYPE_H
#define QDECLARATIVEGEOMAPTYPE_H

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

#include <QtCore/QObject>
#include <QtQml/qqml.h>
#include <QtLocation/private/qgeomaptype_p.h>

QT_BEGIN_NAMESPACE

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeGeoCameraCapabilities: public QObject
{
    Q_OBJECT
    Q_PROPERTY(qreal minimumZoomLevel READ minimumZoomLevel CONSTANT)
    Q_PROPERTY(qreal maximumZoomLevel READ maximumZoomLevel CONSTANT)
    Q_PROPERTY(qreal minimumTilt READ minimumTilt CONSTANT)
    Q_PROPERTY(qreal maximumTilt READ maximumTilt CONSTANT)
    Q_PROPERTY(qreal minimumFieldOfView READ minimumFieldOfView CONSTANT)
    Q_PROPERTY(qreal maximumFieldOfView READ maximumFieldOfView CONSTANT)

public:
    QDeclarativeGeoCameraCapabilities(const QGeoCameraCapabilities &cameraCaps, QObject *parent = 0);
    ~QDeclarativeGeoCameraCapabilities();

    qreal minimumZoomLevel() const;
    qreal maximumZoomLevel() const;
    qreal minimumTilt() const;
    qreal maximumTilt() const;
    qreal minimumFieldOfView() const;
    qreal maximumFieldOfView() const;

private:
    QGeoCameraCapabilities cameraCaps_;
};

class Q_LOCATION_PRIVATE_EXPORT QDeclarativeGeoMapType : public QObject
{
    Q_OBJECT
    Q_ENUMS(MapStyle)

    Q_PROPERTY(MapStyle style READ style CONSTANT)
    Q_PROPERTY(QString name READ name CONSTANT)
    Q_PROPERTY(QString description READ description CONSTANT)
    Q_PROPERTY(bool mobile READ mobile CONSTANT)
    Q_PROPERTY(bool night READ night CONSTANT REVISION 1)
    Q_PROPERTY(QDeclarativeGeoCameraCapabilities *cameraCapabilities READ cameraCapabilities CONSTANT)
    Q_PROPERTY(QVariantMap metadata READ metadata CONSTANT)

public:
    enum MapStyle {
        NoMap = QGeoMapType::NoMap,
        StreetMap = QGeoMapType::StreetMap,
        SatelliteMapDay = QGeoMapType::SatelliteMapDay,
        SatelliteMapNight = QGeoMapType::SatelliteMapNight,
        TerrainMap = QGeoMapType::TerrainMap,
        HybridMap = QGeoMapType::HybridMap,
        TransitMap = QGeoMapType::TransitMap,
        GrayStreetMap = QGeoMapType::GrayStreetMap,
        PedestrianMap = QGeoMapType::PedestrianMap,
        CarNavigationMap = QGeoMapType::CarNavigationMap,
        CycleMap = QGeoMapType::CycleMap,
        CustomMap = 100
    };

    QDeclarativeGeoMapType(const QGeoMapType mapType, QObject *parent = 0);
    ~QDeclarativeGeoMapType();

    MapStyle style() const;
    QString name() const;
    QString description() const;
    bool mobile() const;
    bool night() const;
    QDeclarativeGeoCameraCapabilities *cameraCapabilities() const;
    QVariantMap metadata() const;

    const QGeoMapType mapType() { return mapType_; }

private:
    QGeoMapType mapType_;
    QDeclarativeGeoCameraCapabilities *cameraCapabilities_;
};

QT_END_NAMESPACE

QML_DECLARE_TYPE(QDeclarativeGeoMapType)

#endif
