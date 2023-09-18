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

#ifndef QGEOPOSITIONINFOSOURCEFACTORY_H
#define QGEOPOSITIONINFOSOURCEFACTORY_H

#include <QtPositioning/QGeoPositionInfoSource>
#include <QtPositioning/QGeoSatelliteInfoSource>
#include <QtPositioning/QGeoAreaMonitorSource>
#include <QtCore/QList>

QT_BEGIN_NAMESPACE

class Q_POSITIONING_EXPORT QGeoPositionInfoSourceFactory
{
public:
    virtual ~QGeoPositionInfoSourceFactory();

    virtual QGeoPositionInfoSource *positionInfoSource(QObject *parent) = 0;
    virtual QGeoSatelliteInfoSource *satelliteInfoSource(QObject *parent) = 0;
    virtual QGeoAreaMonitorSource *areaMonitor(QObject *parent) = 0;
};

#define QT_POSITION_SOURCE_INTERFACE
Q_DECLARE_INTERFACE(QGeoPositionInfoSourceFactory,
                    "org.qt-project.qt.position.sourcefactory/5.0")

class Q_POSITIONING_EXPORT QGeoPositionInfoSourceFactoryV2 : public QGeoPositionInfoSourceFactory
{
public:
    virtual ~QGeoPositionInfoSourceFactoryV2();

    virtual QGeoPositionInfoSource *positionInfoSourceWithParameters(QObject *parent, const QVariantMap &parameters) = 0;
    virtual QGeoSatelliteInfoSource *satelliteInfoSourceWithParameters(QObject *parent, const QVariantMap &parameters) = 0;
    virtual QGeoAreaMonitorSource *areaMonitorWithParameters(QObject *parent, const QVariantMap &parameters) = 0;
};

// Although not actually used for constructing a specialized loader, this is required for
// casting a QObject * into QGeoPositionInfoSourceFactoryV2 *
Q_DECLARE_INTERFACE(QGeoPositionInfoSourceFactoryV2,
                    "org.qt-project.qt.position.sourcefactoryV2/5.0")

QT_END_NAMESPACE

#endif // QGEOPOSITIONINFOSOURCEFACTORY_H
