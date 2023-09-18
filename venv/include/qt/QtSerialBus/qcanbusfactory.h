/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtSerialBus module of the Qt Toolkit.
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

#ifndef QCANBUSFACTORY_H
#define QCANBUSFACTORY_H

#include <QtCore/qstringlist.h>
#include <QtSerialBus/qtserialbusglobal.h>
#include <QtSerialBus/qcanbusdevice.h>
#include <QtSerialBus/qcanbusdeviceinfo.h>

QT_BEGIN_NAMESPACE

class Q_SERIALBUS_EXPORT QCanBusFactory
{
public:
    virtual QCanBusDevice *createDevice(const QString &interfaceName,
                                        QString *errorMessage) const = 0;
protected:
    virtual ~QCanBusFactory() {}
};

Q_DECLARE_INTERFACE(QCanBusFactory, "org.qt-project.Qt.QCanBusFactory")

class Q_SERIALBUS_EXPORT QCanBusFactoryV2 : public QCanBusFactory
{
public:
    virtual QCanBusDevice *createDevice(const QString &interfaceName,
                                        QString *errorMessage) const = 0;
    virtual QList<QCanBusDeviceInfo> availableDevices(QString *errorMessage) const = 0;

protected:
    virtual ~QCanBusFactoryV2();
};

Q_DECLARE_INTERFACE(QCanBusFactoryV2, "org.qt-project.Qt.QCanBusFactoryV2")

QT_END_NAMESPACE

#endif // QCANBUSFACTORY_H

