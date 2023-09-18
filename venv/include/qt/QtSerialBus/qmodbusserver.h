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

#ifndef QMODBUSERVER_H
#define QMODBUSERVER_H

#include <QtCore/qobject.h>
#include <QtCore/qvariant.h>
#include <QtSerialBus/qmodbusdataunit.h>
#include <QtSerialBus/qmodbusdevice.h>
#include <QtSerialBus/qmodbuspdu.h>

QT_BEGIN_NAMESPACE

class QModbusServerPrivate;

class Q_SERIALBUS_EXPORT QModbusServer : public QModbusDevice
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QModbusServer)

public:
    enum Option {
        DiagnosticRegister,
        ExceptionStatusOffset,
        DeviceBusy,
        AsciiInputDelimiter,
        ListenOnlyMode,
        ServerIdentifier,
        RunIndicatorStatus,
        AdditionalData,
        DeviceIdentification,
        // Reserved
        UserOption = 0x100
    };
    Q_ENUM(Option)

    explicit QModbusServer(QObject *parent = nullptr);
    ~QModbusServer();

    int serverAddress() const;
    void setServerAddress(int serverAddress);

    virtual bool setMap(const QModbusDataUnitMap &map);
    virtual bool processesBroadcast() const { return false; }

    virtual QVariant value(int option) const;
    virtual bool setValue(int option, const QVariant &value);

    bool data(QModbusDataUnit *newData) const;
    bool setData(const QModbusDataUnit &unit);

    bool setData(QModbusDataUnit::RegisterType table, quint16 address, quint16 data);
    bool data(QModbusDataUnit::RegisterType table, quint16 address, quint16 *data) const;

Q_SIGNALS:
    void dataWritten(QModbusDataUnit::RegisterType table, int address, int size);

protected:
    QModbusServer(QModbusServerPrivate &dd, QObject *parent = nullptr);

    virtual bool writeData(const QModbusDataUnit &unit);
    virtual bool readData(QModbusDataUnit *newData) const;

    virtual QModbusResponse processRequest(const QModbusPdu &request);
    virtual QModbusResponse processPrivateRequest(const QModbusPdu &request);
};

Q_DECLARE_TYPEINFO(QModbusServer::Option, Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

#endif // QMODBUSERVER_H
