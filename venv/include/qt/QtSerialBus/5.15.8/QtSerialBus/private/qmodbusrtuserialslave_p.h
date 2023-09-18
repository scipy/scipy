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

#ifndef QMODBUSRTUSERIALSLAVE_P_H
#define QMODBUSRTUSERIALSLAVE_P_H

#include <QtCore/qbytearray.h>
#include <QtCore/qdebug.h>
#include <QtCore/qelapsedtimer.h>
#include <QtCore/qloggingcategory.h>
#include <QtCore/qmath.h>
#include <QtSerialBus/qmodbusrtuserialslave.h>
#include <QtSerialPort/qserialport.h>

#include <private/qmodbusadu_p.h>
#include <private/qmodbusserver_p.h>

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

class QModbusRtuSerialSlavePrivate : public QModbusServerPrivate
{
    Q_DECLARE_PUBLIC(QModbusRtuSerialSlave)

public:
    void setupSerialPort()
    {
        Q_Q(QModbusRtuSerialSlave);

        m_serialPort = new QSerialPort(q);
        QObject::connect(m_serialPort, &QSerialPort::readyRead, q, [this]() {

            if (m_interFrameTimer.isValid()
                    && m_interFrameTimer.elapsed() > m_interFrameDelayMilliseconds
                    && !m_requestBuffer.isEmpty()) {
                // This permits response buffer clearing if it contains garbage
                // but still permits cases where very slow baud rates can cause
                // chunked and delayed packets
                qCDebug(QT_MODBUS_LOW) << "(RTU server) Dropping older ADU fragments due to larger than 3.5 char delay (expected:"
                                       << m_interFrameDelayMilliseconds << ", max:"
                                       << m_interFrameTimer.elapsed() << ")";
                m_requestBuffer.clear();
            }

            m_interFrameTimer.start();

            const qint64 size = m_serialPort->size();
            m_requestBuffer += m_serialPort->read(size);

            const QModbusSerialAdu adu(QModbusSerialAdu::Rtu, m_requestBuffer);
            qCDebug(QT_MODBUS_LOW) << "(RTU server) Received ADU:" << adu.rawData().toHex();

            // Index                         -> description
            // Server address                -> 1 byte
            // FunctionCode                  -> 1 byte
            // FunctionCode specific content -> 0-252 bytes
            // CRC                           -> 2 bytes
            Q_Q(QModbusRtuSerialSlave);
            QModbusCommEvent event = QModbusCommEvent::ReceiveEvent;
            if (q->value(QModbusServer::ListenOnlyMode).toBool())
                event |= QModbusCommEvent::ReceiveFlag::CurrentlyInListenOnlyMode;

            // We expect at least the server address, function code and CRC.
            if (adu.rawSize() < 4) { // TODO: LRC should be 3 bytes.
                qCWarning(QT_MODBUS) << "(RTU server) Incomplete ADU received, ignoring";

                // The quantity of CRC errors encountered by the remote device since its last
                // restart, clear counters operation, or power-up. In case of a message
                // length < 4 bytes, the receiving device is not able to calculate the CRC.
                incrementCounter(QModbusServerPrivate::Counter::BusCommunicationError);
                storeModbusCommEvent(event | QModbusCommEvent::ReceiveFlag::CommunicationError);
                return;
            }

            // Server address is set to 0, this is a broadcast.
            m_processesBroadcast = (adu.serverAddress() == 0);
            if (q->processesBroadcast())
                event |= QModbusCommEvent::ReceiveFlag::BroadcastReceived;

            const int pduSizeWithoutFcode = QModbusRequest::calculateDataSize(adu.pdu());

            // server address byte + function code byte + PDU size + 2 bytes CRC
            if ((pduSizeWithoutFcode < 0) || ((2 + pduSizeWithoutFcode + 2) != adu.rawSize())) {
                qCWarning(QT_MODBUS) << "(RTU server) ADU does not match expected size, ignoring";
                // The quantity of messages addressed to the remote device that it could not
                // handle due to a character overrun condition, since its last restart, clear
                // counters operation, or power-up. A character overrun is caused by data
                // characters arriving at the port faster than they can be stored, or by the loss
                // of a character due to a hardware malfunction.
                incrementCounter(QModbusServerPrivate::Counter::BusCharacterOverrun);
                storeModbusCommEvent(event | QModbusCommEvent::ReceiveFlag::CharacterOverrun);
                return;
            }

            // We received the full message, including checksum. We do not expect more bytes to
            // arrive, so clear the buffer. All new bytes are considered part of the next message.
            m_requestBuffer.resize(0);

            if (!adu.matchingChecksum()) {
                qCWarning(QT_MODBUS) << "(RTU server) Discarding request with wrong CRC, received:"
                                     << adu.checksum<quint16>() << ", calculated CRC:"
                                     << QModbusSerialAdu::calculateCRC(adu.data(), adu.size());
                // The quantity of CRC errors encountered by the remote device since its last
                // restart, clear counters operation, or power-up.
                incrementCounter(QModbusServerPrivate::Counter::BusCommunicationError);
                storeModbusCommEvent(event | QModbusCommEvent::ReceiveFlag::CommunicationError);
                return;
            }

            // The quantity of messages that the remote device has detected on the communications
            // system since its last restart, clear counters operation, or power-up.
            incrementCounter(QModbusServerPrivate::Counter::BusMessage);

            // If we do not process a Broadcast ...
            if (!q->processesBroadcast()) {
                // check if the server address matches ...
                if (q->serverAddress() != adu.serverAddress()) {
                    // no, not our address! Ignore!
                    qCDebug(QT_MODBUS) << "(RTU server) Wrong server address, expected"
                        << q->serverAddress() << "got" << adu.serverAddress();
                    return;
                }
            } // else { Broadcast -> Server address will never match, deliberately ignore }

            storeModbusCommEvent(event); // store the final event before processing

            const QModbusRequest req = adu.pdu();
            qCDebug(QT_MODBUS) << "(RTU server) Request PDU:" << req;
            QModbusResponse response; // If the device ...
            if (q->value(QModbusServer::DeviceBusy).value<quint16>() == 0xffff) {
                // is busy, update the quantity of messages addressed to the remote device for
                // which it returned a Server Device Busy exception response, since its last
                // restart, clear counters operation, or power-up.
                incrementCounter(QModbusServerPrivate::Counter::ServerBusy);
                response = QModbusExceptionResponse(req.functionCode(),
                    QModbusExceptionResponse::ServerDeviceBusy);
            } else {
                // is not busy, update the quantity of messages addressed to the remote device,
                // or broadcast, that the remote device has processed since its last restart,
                // clear counters operation, or power-up.
                incrementCounter(QModbusServerPrivate::Counter::ServerMessage);
                response = q->processRequest(req);
            }
            qCDebug(QT_MODBUS) << "(RTU server) Response PDU:" << response;

            event = QModbusCommEvent::SentEvent; // reset event after processing
            if (q->value(QModbusServer::ListenOnlyMode).toBool())
                event |= QModbusCommEvent::SendFlag::CurrentlyInListenOnlyMode;

            if ((!response.isValid())
                || q->processesBroadcast()
                || q->value(QModbusServer::ListenOnlyMode).toBool()) {
                // The quantity of messages addressed to the remote device for which it has
                // returned no response (neither a normal response nor an exception response),
                // since its last restart, clear counters operation, or power-up.
                incrementCounter(QModbusServerPrivate::Counter::ServerNoResponse);
                storeModbusCommEvent(event);
                return;
            }

            const QByteArray result = QModbusSerialAdu::create(QModbusSerialAdu::Rtu,
                                                               q->serverAddress(), response);

            qCDebug(QT_MODBUS_LOW) << "(RTU server) Response ADU:" << result.toHex();

            if (!m_serialPort->isOpen()) {
                qCDebug(QT_MODBUS) << "(RTU server) Requesting serial port has closed.";
                q->setError(QModbusRtuSerialSlave::tr("Requesting serial port is closed"),
                            QModbusDevice::WriteError);
                incrementCounter(QModbusServerPrivate::Counter::ServerNoResponse);
                storeModbusCommEvent(event);
                return;
            }

            qint64 writtenBytes = m_serialPort->write(result);
            if ((writtenBytes == -1) || (writtenBytes < result.size())) {
                qCDebug(QT_MODBUS) << "(RTU server) Cannot write requested response to serial port.";
                q->setError(QModbusRtuSerialSlave::tr("Could not write response to client"),
                            QModbusDevice::WriteError);
                incrementCounter(QModbusServerPrivate::Counter::ServerNoResponse);
                storeModbusCommEvent(event);
                m_serialPort->clear(QSerialPort::Output);
                return;
            }

            if (response.isException()) {
                switch (response.exceptionCode()) {
                case QModbusExceptionResponse::IllegalFunction:
                case QModbusExceptionResponse::IllegalDataAddress:
                case QModbusExceptionResponse::IllegalDataValue:
                    event |= QModbusCommEvent::SendFlag::ReadExceptionSent;
                    break;

                case QModbusExceptionResponse::ServerDeviceFailure:
                    event |= QModbusCommEvent::SendFlag::ServerAbortExceptionSent;
                    break;

                case QModbusExceptionResponse::ServerDeviceBusy:
                    // The quantity of messages addressed to the remote device for which it
                    // returned a server device busy exception response, since its last restart,
                    // clear counters operation, or power-up.
                    incrementCounter(QModbusServerPrivate::Counter::ServerBusy);
                    event |= QModbusCommEvent::SendFlag::ServerBusyExceptionSent;
                    break;

                case  QModbusExceptionResponse::NegativeAcknowledge:
                    // The quantity of messages addressed to the remote device for which it
                    // returned a negative acknowledge (NAK) exception response, since its last
                    // restart, clear counters operation, or power-up.
                    incrementCounter(QModbusServerPrivate::Counter::ServerNAK);
                    event |= QModbusCommEvent::SendFlag::ServerProgramNAKExceptionSent;
                    break;

                default:
                    break;
                }
                // The quantity of Modbus exception responses returned by the remote device since
                // its last restart, clear counters operation, or power-up.
                incrementCounter(QModbusServerPrivate::Counter::BusExceptionError);
            } else {
                switch (quint16(req.functionCode())) {
                case 0x0a: // Poll 484 (not in the official Modbus specification) *1
                case 0x0e: // Poll Controller (not in the official Modbus specification) *1
                case QModbusRequest::GetCommEventCounter: // fall through and bail out
                    break;
                default:
                    // The device's event counter is incremented once for each successful message
                    // completion. Do not increment for exception responses, poll commands, or fetch
                    // event counter commands.            *1 but mentioned here ^^^
                    incrementCounter(QModbusServerPrivate::Counter::CommEvent);
                    break;
                }
            }
            storeModbusCommEvent(event); // store the final event after processing
        });

        QObject::connect(m_serialPort, &QSerialPort::errorOccurred, q,
                         [this](QSerialPort::SerialPortError error) {
            if (error == QSerialPort::NoError)
                return;

            qCDebug(QT_MODBUS) << "(RTU server) QSerialPort error:" << error
                               << (m_serialPort ? m_serialPort->errorString() : QString());

            Q_Q(QModbusRtuSerialSlave);

            switch (error) {
            case QSerialPort::DeviceNotFoundError:
                q->setError(QModbusDevice::tr("Referenced serial device does not exist."),
                            QModbusDevice::ConnectionError);
                break;
            case QSerialPort::PermissionError:
                q->setError(QModbusDevice::tr("Cannot open serial device due to permissions."),
                            QModbusDevice::ConnectionError);
                break;
            case QSerialPort::OpenError:
            case QSerialPort::NotOpenError:
                q->setError(QModbusDevice::tr("Cannot open serial device."),
                            QModbusDevice::ConnectionError);
                break;
            case QSerialPort::WriteError:
                q->setError(QModbusDevice::tr("Write error."), QModbusDevice::WriteError);
                break;
            case QSerialPort::ReadError:
                q->setError(QModbusDevice::tr("Read error."), QModbusDevice::ReadError);
                break;
            case QSerialPort::ResourceError:
                q->setError(QModbusDevice::tr("Resource error."), QModbusDevice::ConnectionError);
                break;
            case QSerialPort::UnsupportedOperationError:
                q->setError(QModbusDevice::tr("Device operation is not supported error."),
                            QModbusDevice::ConfigurationError);
                break;
            case QSerialPort::TimeoutError:
                q->setError(QModbusDevice::tr("Timeout error."), QModbusDevice::TimeoutError);
                break;
            case QSerialPort::UnknownError:
                q->setError(QModbusDevice::tr("Unknown error."), QModbusDevice::UnknownError);
                break;
            default:
                qCDebug(QT_MODBUS) << "(RTU server) Unhandled QSerialPort error" << error;
                break;
            }
        });

        QObject::connect(m_serialPort, &QSerialPort::aboutToClose, q, [this]() {
            Q_Q(QModbusRtuSerialSlave);
            // update state if socket closure was caused by remote side
            if (q->state() != QModbusDevice::ClosingState)
                q->setState(QModbusDevice::UnconnectedState);
        });
    }

    void setupEnvironment() {
        if (m_serialPort) {
            m_serialPort->setPortName(m_comPort);
            m_serialPort->setParity(m_parity);
            m_serialPort->setBaudRate(m_baudRate);
            m_serialPort->setDataBits(m_dataBits);
            m_serialPort->setStopBits(m_stopBits);
        }

        // for calculation details see
        // QModbusRtuSerialMasterPrivate::calculateInterFrameDelay()
        m_interFrameDelayMilliseconds = qMax(m_interFrameDelayMilliseconds,
                                             qCeil(3500. / (qreal(m_baudRate) / 11.)));

        m_requestBuffer.clear();
    }

    QIODevice *device() const override { return m_serialPort; }

    QByteArray m_requestBuffer;
    bool m_processesBroadcast = false;
    QSerialPort *m_serialPort = nullptr;
    QElapsedTimer m_interFrameTimer;
    int m_interFrameDelayMilliseconds = 2;
};

QT_END_NAMESPACE

#endif // QMODBUSRTUSERIALSLAVE_P_H
