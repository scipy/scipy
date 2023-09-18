/****************************************************************************
**
** Copyright (C) 2017 Andre Hartmann <aha_1980@gmx.de>
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

#ifndef QCANBUSDEVICEINFO_H
#define QCANBUSDEVICEINFO_H

#include <QtCore/qshareddata.h>
#include <QtCore/qstring.h>
#include <QtSerialBus/qtserialbusglobal.h>

QT_BEGIN_NAMESPACE

class QCanBusDeviceInfoPrivate;

class Q_SERIALBUS_EXPORT QCanBusDeviceInfo
{
public:
    QCanBusDeviceInfo() = delete;
    QCanBusDeviceInfo(const QCanBusDeviceInfo &other);
    ~QCanBusDeviceInfo();

    void swap(QCanBusDeviceInfo &other) Q_DECL_NOTHROW
    {
         qSwap(d_ptr, other.d_ptr);
    }

    QCanBusDeviceInfo &operator=(const QCanBusDeviceInfo &other);
    QCanBusDeviceInfo &operator=(QCanBusDeviceInfo &&other) Q_DECL_NOTHROW
    {
        swap(other);
        return *this;
    }

    QString name() const;
    QString description() const;
    QString serialNumber() const;
    int channel() const;

    bool hasFlexibleDataRate() const;
    bool isVirtual() const;

private:
    friend class QCanBusDevice;

    explicit QCanBusDeviceInfo(QCanBusDeviceInfoPrivate &dd);

    QSharedDataPointer<QCanBusDeviceInfoPrivate> d_ptr;
};

Q_DECLARE_SHARED(QCanBusDeviceInfo)

QT_END_NAMESPACE

#endif // QCANBUSDEVICEINFO_H
