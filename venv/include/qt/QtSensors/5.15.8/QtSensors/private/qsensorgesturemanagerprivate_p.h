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

#ifndef QSENSORGESTUREMANAGERPRIVATE_P_H
#define QSENSORGESTUREMANAGERPRIVATE_P_H

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

#include <QObject>
#include <QMap>
#include <QStringList>
#include <QDebug>
#include <QSharedPointer>
#include <QPluginLoader>

#include "qsensorgesture.h"
#include "qsensorgesturerecognizer.h"

QT_BEGIN_NAMESPACE

class QFactoryLoader;

class QSensorGestureManagerPrivate : public QObject
{
    Q_OBJECT
public:
    explicit QSensorGestureManagerPrivate(QObject *parent = 0);
    ~QSensorGestureManagerPrivate();

    QMap<QString, QSensorGestureRecognizer *> registeredSensorGestures;

    QList <QObject *> plugins;

    QFactoryLoader *loader;
    void loadPlugins();
    bool loadRecognizer(const QString &id);

    QSensorGestureRecognizer *sensorGestureRecognizer(const QString &id);

    bool registerSensorGestureRecognizer(QSensorGestureRecognizer *recognizer);
    QStringList gestureIds();
    QStringList knownIds;
    void initPlugin(QObject *o);
#ifdef SIMULATOR_BUILD
    void recognizerStarted(const QSensorGestureRecognizer *);
    void recognizerStopped(const QSensorGestureRecognizer *);
#endif

    static QSensorGestureManagerPrivate * instance();
Q_SIGNALS:
        void newSensorGestureAvailable();

#ifdef SIMULATOR_BUILD
Q_SIGNALS:
    void newSensorGestures(QStringList);
    void removeSensorGestures(QStringList);

private slots:
    void sensorGestureDetected();

#endif
};

QT_END_NAMESPACE

#endif // QSENSORGESTUREMANAGERPRIVATE_P_H
