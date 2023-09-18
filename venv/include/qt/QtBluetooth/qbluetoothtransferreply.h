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

#ifndef QBLUETOOTHTRANSFERREPLY_H
#define QBLUETOOTHTRANSFERREPLY_H

#include <QtCore/QIODevice>

#include <QtBluetooth/QBluetoothTransferRequest>
#include <QtBluetooth/QBluetoothTransferManager>

QT_BEGIN_NAMESPACE

class QBluetoothTransferReplyPrivate;

class Q_BLUETOOTH_EXPORT QBluetoothTransferReply : public QObject
{
    Q_OBJECT

public:
    enum TransferError {
        NoError = 0,
        UnknownError,
        FileNotFoundError,
        HostNotFoundError,
        UserCanceledTransferError,
        IODeviceNotReadableError,
        ResourceBusyError,
        SessionError
    };
    Q_ENUM(TransferError)

    ~QBluetoothTransferReply();

    virtual bool isFinished() const = 0;
    virtual bool isRunning() const = 0;

    QBluetoothTransferManager *manager() const;

    virtual TransferError error() const = 0;
    virtual QString errorString() const = 0;

    QBluetoothTransferRequest request() const;

public Q_SLOTS:
    void abort();

Q_SIGNALS:
    //TODO Remove QBluetoothTransferReply* parameter in Qt 6
    void finished(QBluetoothTransferReply *);
    void transferProgress(qint64 bytesTransferred, qint64 bytesTotal);
    void error(QBluetoothTransferReply::TransferError lastError);

protected:
    explicit QBluetoothTransferReply(QObject *parent = nullptr);
    void setManager(QBluetoothTransferManager *manager);
    void setRequest(const QBluetoothTransferRequest &request);

protected:
    QBluetoothTransferReplyPrivate *d_ptr;

private:
    Q_DECLARE_PRIVATE(QBluetoothTransferReply)

};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QBluetoothTransferReply::TransferError)

#endif // QBLUETOOTHTRANSFERREPLY_H
