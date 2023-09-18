/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtSensors module of the Qt Toolkit.
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

#ifndef QMAGNETOMETER_H
#define QMAGNETOMETER_H

#include <QtSensors/qsensor.h>

QT_BEGIN_NAMESPACE

class QMagnetometerReadingPrivate;

class Q_SENSORS_EXPORT QMagnetometerReading : public QSensorReading
{
    Q_OBJECT
    Q_PROPERTY(qreal x READ x)
    Q_PROPERTY(qreal y READ y)
    Q_PROPERTY(qreal z READ z)
    Q_PROPERTY(qreal calibrationLevel READ calibrationLevel)
    DECLARE_READING(QMagnetometerReading)
public:
    qreal x() const;
    void setX(qreal x);

    qreal y() const;
    void setY(qreal y);

    qreal z() const;
    void setZ(qreal z);

    qreal calibrationLevel() const;
    void setCalibrationLevel(qreal calibrationLevel);
};

class Q_SENSORS_EXPORT QMagnetometerFilter : public QSensorFilter
{
public:
    virtual bool filter(QMagnetometerReading *reading) = 0;
private:
    bool filter(QSensorReading *reading) override;
};

class QMagnetometerPrivate;

class Q_SENSORS_EXPORT QMagnetometer : public QSensor
{
    Q_OBJECT
    Q_PROPERTY(bool returnGeoValues READ returnGeoValues WRITE setReturnGeoValues NOTIFY returnGeoValuesChanged)
public:
    explicit QMagnetometer(QObject *parent = Q_NULLPTR);
    virtual ~QMagnetometer();
    QMagnetometerReading *reading() const;
    static char const * const type;

    bool returnGeoValues() const;
    void setReturnGeoValues(bool returnGeoValues);

Q_SIGNALS:
    void returnGeoValuesChanged(bool returnGeoValues);

private:
    Q_DECLARE_PRIVATE(QMagnetometer)
    Q_DISABLE_COPY(QMagnetometer)
};

QT_END_NAMESPACE

#endif

