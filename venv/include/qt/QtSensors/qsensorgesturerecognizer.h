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

#ifndef QSENSORGESTURERECOGNIZER_H
#define QSENSORGESTURERECOGNIZER_H

#include <QtCore/QDebug>
#include <QtCore/QTimer>
#include <QtCore/QStringList>

#include <QtSensors/qsensorgesture.h>

QT_BEGIN_NAMESPACE

class QSensorGestureRecognizerPrivate;
class Q_SENSORS_EXPORT QSensorGestureRecognizer : public QObject
{
    Q_OBJECT
public:
    explicit QSensorGestureRecognizer(QObject *parent = Q_NULLPTR);
    virtual ~QSensorGestureRecognizer();

    virtual QString id() const = 0;

    virtual bool isActive() = 0;

    void startBackend();
    void stopBackend();
    void createBackend();

    QStringList gestureSignals() const;

Q_SIGNALS:
    void detected(const QString &);

protected:
    virtual void create() = 0;
    virtual bool start() = 0;
    virtual bool stop() = 0;

private:
        QSensorGestureRecognizerPrivate * d_ptr;
};

QT_END_NAMESPACE

#endif // QSENSORGESTURERECOGNIZER_H
