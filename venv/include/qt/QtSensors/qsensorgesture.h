/****************************************************************************
**
** Copyright (C) 2019 The Qt Company Ltd.
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

#ifndef QSENSORGESTURE_H
#define QSENSORGESTURE_H

#include <QtCore/QObject>
#include <QtCore/QStringList>
#include <QtSensors/qsensorsglobal.h>

#include <QtCore/QList>
#include <QtCore/QMap>
#include <QtCore/QVector>

#include <QtCore/qmetatype.h>

QT_BEGIN_NAMESPACE

class QSensorGesturePrivate;

class Q_SENSORS_EXPORT QSensorGesture : public QObject
{
    //Do not use Q_OBJECT here
public:
    explicit QSensorGesture(const QStringList &ids, QObject *parent = Q_NULLPTR);
    ~QSensorGesture();

    bool isActive();

    QStringList validIds() const;
    QStringList invalidIds() const;

    QStringList gestureSignals() const;

    void startDetection();
    void stopDetection();

private:
    QSensorGesturePrivate * d_ptr;

    // ### fixme: Qt 6: Make public to enable Qt for Python bindings
private:
    // Pretend to be a Q_OBJECT
    const QMetaObject *metaObject() const override;
    int qt_metacall(QMetaObject::Call, int, void **) override;

Q_SIGNALS:
    // these signals are created at runtime, along with
    // gesture recognizer specific signals.
     void detected(QString);
};

QT_END_NAMESPACE


#endif // QSENSORGESTURE_H
