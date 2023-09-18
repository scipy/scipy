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

#ifndef QMODBUSREPLY_H
#define QMODBUSREPLY_H

#include <QtSerialBus/qmodbusdataunit.h>
#include <QtSerialBus/qmodbusdevice.h>
#include <QtSerialBus/qmodbuspdu.h>

QT_BEGIN_NAMESPACE

class QModbusReplyPrivate;

class Q_SERIALBUS_EXPORT QModbusReply : public QObject
{
    Q_OBJECT
    Q_DECLARE_PRIVATE(QModbusReply)

public:
    enum ReplyType {
        Raw,
        Common,
        Broadcast
    };
    Q_ENUM(ReplyType)

    QModbusReply(ReplyType type, int serverAddress, QObject *parent = nullptr);

    ReplyType type() const;
    int serverAddress() const;

    bool isFinished() const;

    QModbusDataUnit result() const;
    QModbusResponse rawResult() const;

    QString errorString() const;
    QModbusDevice::Error error() const;

    void setResult(const QModbusDataUnit &unit);
    void setRawResult(const QModbusResponse &unit);

    void setFinished(bool isFinished);
    void setError(QModbusDevice::Error error, const QString &errorText);

Q_SIGNALS:
    void finished();
    void errorOccurred(QModbusDevice::Error error);
};
Q_DECLARE_TYPEINFO(QModbusReply::ReplyType, Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

#endif // QMODBUSREPLY_H
