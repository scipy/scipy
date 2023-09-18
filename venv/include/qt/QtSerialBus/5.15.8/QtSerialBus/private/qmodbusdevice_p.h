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

#ifndef QMODBUSDEVICE_P_H
#define QMODBUSDEVICE_P_H

#include <QtCore/qvariant.h>
#include <QtSerialBus/qmodbusdevice.h>
#if QT_CONFIG(modbus_serialport)
#include <QtSerialPort/qserialport.h>
#endif

#include <private/qobject_p.h>

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API. It exists purely as an
// implementation detail. This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

QT_BEGIN_NAMESPACE

class QModbusDevicePrivate : public QObjectPrivate
{
    Q_DECLARE_PUBLIC(QModbusDevice)

public:
    QModbusDevice::State state = QModbusDevice::UnconnectedState;
    QModbusDevice::Error error = QModbusDevice::NoError;
    QString errorString;

#if QT_CONFIG(modbus_serialport)
    QString m_comPort;
    QSerialPort::DataBits m_dataBits = QSerialPort::Data8;
    QSerialPort::Parity m_parity = QSerialPort::EvenParity;
    QSerialPort::StopBits m_stopBits = QSerialPort::OneStop;
    QSerialPort::BaudRate m_baudRate = QSerialPort::Baud19200;
#endif

    int m_networkPort = 502;
    QString m_networkAddress = QStringLiteral("127.0.0.1");

    QHash<int, QVariant> m_userConnectionParams; // ### Qt6: remove

    virtual QIODevice *device() const { return nullptr; }
};

QT_END_NAMESPACE

#endif // QMODBUSDEVICE_P_H
