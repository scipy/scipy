/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtBluetooth module of the Qt Toolkit.
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

#ifndef QBLUETOOTHTRANSFERREQUEST_H
#define QBLUETOOTHTRANSFERREQUEST_H

#include <QtBluetooth/qtbluetoothglobal.h>
#include <QtBluetooth/QBluetoothAddress>

#include <QtCore/QtGlobal>
#include <QtCore/QVariant>

QT_BEGIN_NAMESPACE

class QBluetoothAddress;
class QBluetoothTransferRequestPrivate;

class Q_BLUETOOTH_EXPORT QBluetoothTransferRequest
{
public:
    enum Attribute {
        DescriptionAttribute,
        TimeAttribute,
        TypeAttribute,
        LengthAttribute,
        NameAttribute
    };

    explicit QBluetoothTransferRequest(const QBluetoothAddress &address = QBluetoothAddress());
    QBluetoothTransferRequest(const QBluetoothTransferRequest &other);
    ~QBluetoothTransferRequest();

    QVariant attribute(Attribute code, const QVariant &defaultValue = QVariant()) const;
    void setAttribute(Attribute code, const QVariant &value);

    QBluetoothAddress address() const;

    bool operator!=(const QBluetoothTransferRequest &other) const;
    QBluetoothTransferRequest &operator=(const QBluetoothTransferRequest &other);
    bool operator==(const QBluetoothTransferRequest &other) const;

protected:
    QBluetoothTransferRequestPrivate *d_ptr;

private:
    Q_DECLARE_PRIVATE(QBluetoothTransferRequest)

};

QT_END_NAMESPACE

#endif // QBLUETOOTHTRANSFERREQUEST_H
