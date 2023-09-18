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

#ifndef QTAPSENSOR_H
#define QTAPSENSOR_H

#include <QtSensors/qsensor.h>

QT_BEGIN_NAMESPACE

class QTapReadingPrivate;

class Q_SENSORS_EXPORT QTapReading : public QSensorReading
{
    Q_OBJECT
    Q_ENUMS(TapDirection)
    Q_PROPERTY(TapDirection tapDirection READ tapDirection)
    Q_PROPERTY(bool doubleTap READ isDoubleTap)
    DECLARE_READING(QTapReading)
public:
    enum TapDirection {
        Undefined = 0,
        X      = 0x0001,
        Y      = 0x0002,
        Z      = 0x0004,
        X_Pos  = 0x0011,
        Y_Pos  = 0x0022,
        Z_Pos  = 0x0044,
        X_Neg  = 0x0101,
        Y_Neg  = 0x0202,
        Z_Neg  = 0x0404,
        X_Both = 0x0111,
        Y_Both = 0x0222,
        Z_Both = 0x0444
    };

    TapDirection tapDirection() const;
    void setTapDirection(TapDirection tapDirection);

    bool isDoubleTap() const;
    void setDoubleTap(bool doubleTap);
};

class Q_SENSORS_EXPORT QTapFilter : public QSensorFilter
{
public:
    virtual bool filter(QTapReading *reading) = 0;
private:
    bool filter(QSensorReading *reading) override;
};

class QTapSensorPrivate;

class Q_SENSORS_EXPORT QTapSensor : public QSensor
{
    Q_OBJECT
    Q_PROPERTY(bool returnDoubleTapEvents READ returnDoubleTapEvents WRITE setReturnDoubleTapEvents
               NOTIFY returnDoubleTapEventsChanged)
public:
    explicit QTapSensor(QObject *parent = Q_NULLPTR);
    virtual ~QTapSensor();
    QTapReading *reading() const;
    static char const * const type;

    bool returnDoubleTapEvents() const;
    void setReturnDoubleTapEvents(bool returnDoubleTapEvents);

Q_SIGNALS:
    void returnDoubleTapEventsChanged(bool returnDoubleTapEvents);

private:
    Q_DECLARE_PRIVATE(QTapSensor)
    Q_DISABLE_COPY(QTapSensor)
};

QT_END_NAMESPACE

#endif

