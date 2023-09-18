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

#ifndef QMODBUSTCPCLIENT_P_H
#define QMODBUSTCPCLIENT_P_H

#include <QtCore/qloggingcategory.h>
#include <QtNetwork/qhostaddress.h>
#include <QtNetwork/qtcpsocket.h>
#include "QtSerialBus/qmodbustcpclient.h"

#include "private/qmodbusclient_p.h"

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

Q_DECLARE_LOGGING_CATEGORY(QT_MODBUS)
Q_DECLARE_LOGGING_CATEGORY(QT_MODBUS_LOW)

class QModbusTcpClientPrivate : public QModbusClientPrivate
{
    Q_DECLARE_PUBLIC(QModbusTcpClient)

public:
    void setupTcpSocket()
    {
        Q_Q(QModbusTcpClient);

        m_socket = new QTcpSocket(q);

        QObject::connect(m_socket, &QAbstractSocket::connected, q, [this]() {
            qCDebug(QT_MODBUS) << "(TCP client) Connected to" << m_socket->peerAddress()
                               << "on port" << m_socket->peerPort();
            Q_Q(QModbusTcpClient);
            responseBuffer.clear();
            q->setState(QModbusDevice::ConnectedState);
        });

        QObject::connect(m_socket, &QAbstractSocket::disconnected, q, [this]() {
           qCDebug(QT_MODBUS)  << "(TCP client) Connection closed.";
           Q_Q(QModbusTcpClient);
           q->setState(QModbusDevice::UnconnectedState);
           cleanupTransactionStore();
        });

        QObject::connect(m_socket,
                         &QAbstractSocket::errorOccurred, q,
                         [this](QAbstractSocket::SocketError /*error*/)
        {
            Q_Q(QModbusTcpClient);

            if (m_socket->state() == QAbstractSocket::UnconnectedState) {
                cleanupTransactionStore();
                q->setState(QModbusDevice::UnconnectedState);
            }
            q->setError(QModbusClient::tr("TCP socket error (%1).").arg(m_socket->errorString()),
                        QModbusDevice::ConnectionError);
        });

        QObject::connect(m_socket, &QIODevice::readyRead, q, [this](){
            responseBuffer += m_socket->read(m_socket->bytesAvailable());
            qCDebug(QT_MODBUS_LOW) << "(TCP client) Response buffer:" << responseBuffer.toHex();

            while (!responseBuffer.isEmpty()) {
                // can we read enough for Modbus ADU header?
                if (responseBuffer.size() < mbpaHeaderSize) {
                    qCDebug(QT_MODBUS_LOW) << "(TCP client) Modbus ADU not complete";
                    return;
                }

                quint8 serverAddress;
                quint16 transactionId, bytesPdu, protocolId;
                QDataStream input(responseBuffer);
                input >> transactionId >> protocolId >> bytesPdu >> serverAddress;

                // stop the timer as soon as we know enough about the transaction
                const bool knownTransaction = m_transactionStore.contains(transactionId);
                if (knownTransaction && m_transactionStore[transactionId].timer)
                    m_transactionStore[transactionId].timer->stop();

                qCDebug(QT_MODBUS) << "(TCP client) tid:" << Qt::hex << transactionId << "size:"
                    << bytesPdu << "server address:" << serverAddress;

                // The length field is the byte count of the following fields, including the Unit
                // Identifier and the PDU, so we remove on byte.
                bytesPdu--;

                int tcpAduSize = mbpaHeaderSize + bytesPdu;
                if (responseBuffer.size() < tcpAduSize) {
                    qCDebug(QT_MODBUS) << "(TCP client) PDU too short. Waiting for more data";
                    return;
                }

                QModbusResponse responsePdu;
                input >> responsePdu;
                qCDebug(QT_MODBUS) << "(TCP client) Received PDU:" << responsePdu.functionCode()
                                   << responsePdu.data().toHex();

                responseBuffer.remove(0, tcpAduSize);

                if (!knownTransaction) {
                    qCDebug(QT_MODBUS) << "(TCP client) No pending request for response with "
                        "given transaction ID, ignoring response message.";
                } else {
                    processQueueElement(responsePdu, m_transactionStore[transactionId]);
                }
            }
        });
    }

    QModbusReply *enqueueRequest(const QModbusRequest &request, int serverAddress,
                                 const QModbusDataUnit &unit,
                                 QModbusReply::ReplyType type) override
    {
        auto writeToSocket = [this](quint16 tId, const QModbusRequest &request, int address) {
            QByteArray buffer;
            QDataStream output(&buffer, QIODevice::WriteOnly);
            output << tId << quint16(0) << quint16(request.size() + 1) << quint8(address) << request;

            int writtenBytes = m_socket->write(buffer);
            if (writtenBytes == -1 || writtenBytes < buffer.size()) {
                Q_Q(QModbusTcpClient);
                qCDebug(QT_MODBUS) << "(TCP client) Cannot write request to socket.";
                q->setError(QModbusTcpClient::tr("Could not write request to socket."),
                            QModbusDevice::WriteError);
                return false;
            }
            qCDebug(QT_MODBUS_LOW) << "(TCP client) Sent TCP ADU:" << buffer.toHex();
            qCDebug(QT_MODBUS) << "(TCP client) Sent TCP PDU:" << request << "with tId:" <<Qt:: hex
                << tId;
            return true;
        };

        const int tId = transactionId();
        if (!writeToSocket(tId, request, serverAddress))
            return nullptr;

        Q_Q(QModbusTcpClient);
        auto reply = new QModbusReply(type, serverAddress, q);
        const auto element = QueueElement{ reply, request, unit, m_numberOfRetries,
            m_responseTimeoutDuration };
        m_transactionStore.insert(tId, element);

        q->connect(reply, &QObject::destroyed, q, [this, tId](QObject *) {
            if (!m_transactionStore.contains(tId))
                return;
            const QueueElement element = m_transactionStore.take(tId);
            if (element.timer)
                element.timer->stop();
        });

        if (element.timer) {
            q->connect(q, &QModbusClient::timeoutChanged,
                       element.timer.data(), QOverload<int>::of(&QTimer::setInterval));
            QObject::connect(element.timer.data(), &QTimer::timeout, q, [this, writeToSocket, tId]() {
                if (!m_transactionStore.contains(tId))
                    return;

                QueueElement elem = m_transactionStore.take(tId);
                if (elem.reply.isNull())
                    return;

                if (elem.numberOfRetries > 0) {
                    elem.numberOfRetries--;
                    if (!writeToSocket(tId, elem.requestPdu, elem.reply->serverAddress()))
                        return;
                    m_transactionStore.insert(tId, elem);
                    elem.timer->start();
                    qCDebug(QT_MODBUS) << "(TCP client) Resend request with tId:" << Qt::hex << tId;
                } else {
                    qCDebug(QT_MODBUS) << "(TCP client) Timeout of request with tId:" <<Qt::hex << tId;
                    elem.reply->setError(QModbusDevice::TimeoutError,
                        QModbusClient::tr("Request timeout."));
                }
            });
            element.timer->start();
        } else {
            qCWarning(QT_MODBUS) << "(TCP client) No response timeout timer for request with tId:"
                << Qt::hex << tId << ". Expected timeout:" << m_responseTimeoutDuration;
        }
        incrementTransactionId();

        return reply;
    }

    // TODO: Review once we have a transport layer in place.
    bool isOpen() const override
    {
        if (m_socket)
            return m_socket->isOpen();
        return false;
    }

    void cleanupTransactionStore()
    {
        if (m_transactionStore.isEmpty())
            return;

        qCDebug(QT_MODBUS) << "(TCP client) Cleanup of pending requests";

        for (const auto &elem : qAsConst(m_transactionStore)) {
            if (elem.reply.isNull())
                continue;
            elem.reply->setError(QModbusDevice::ReplyAbortedError,
                                 QModbusClient::tr("Reply aborted due to connection closure."));
        }
        m_transactionStore.clear();
    }

    // This doesn't overflow, it rather "wraps around". Expected.
    inline void incrementTransactionId() { m_transactionId++; }
    inline int transactionId() const { return m_transactionId; }

    QIODevice *device() const override { return m_socket; }

    QTcpSocket *m_socket = nullptr;
    QByteArray responseBuffer;
    QHash<quint16, QueueElement> m_transactionStore;
    int mbpaHeaderSize = 7;

private:   // Private to avoid using the wrong id inside the timer lambda,
    quint16 m_transactionId = 0; // capturing 'this' will not copy the id.
};

QT_END_NAMESPACE

#endif // QMODBUSTCPCLIENT_P_H
