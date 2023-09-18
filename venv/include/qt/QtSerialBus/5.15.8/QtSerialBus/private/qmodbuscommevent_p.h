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

#ifndef QMODBUSCOMMEVENT_P_H
#define QMODBUSCOMMEVENT_P_H

#include <QtSerialBus/qtserialbusglobal.h>

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

QT_BEGIN_NAMESPACE

class QModbusCommEvent
{
public:
    enum struct SendFlag : quint8 {
        ReadExceptionSent = 0x01,
        ServerAbortExceptionSent = 0x02,
        ServerBusyExceptionSent = 0x04,
        ServerProgramNAKExceptionSent = 0x08,
        WriteTimeoutErrorOccurred = 0x10,
        CurrentlyInListenOnlyMode = 0x20,
    };

    enum struct ReceiveFlag : quint8 {
        /* Unused */
        CommunicationError = 0x02,
        /* Unused */
        /* Unused */
        CharacterOverrun = 0x10,
        CurrentlyInListenOnlyMode = 0x20,
        BroadcastReceived = 0x40
    };

    enum EventByte {
        SentEvent = 0x40,
        ReceiveEvent = 0x80,
        EnteredListenOnlyMode = 0x04,
        InitiatedCommunicationRestart = 0x00
    };

    inline QModbusCommEvent(QModbusCommEvent::EventByte byte)
        : m_eventByte(byte) {}

    operator quint8() const { return m_eventByte; }
    operator QModbusCommEvent::EventByte() const {
        return static_cast<QModbusCommEvent::EventByte> (m_eventByte);
    }

    inline QModbusCommEvent operator=(QModbusCommEvent::EventByte byte) {
        m_eventByte = byte;
        return *this;
    }
    inline QModbusCommEvent operator|=(QModbusCommEvent::SendFlag sf) {
        m_eventByte |= quint8(sf);
        return *this;
    }
    inline QModbusCommEvent operator|=(QModbusCommEvent::ReceiveFlag rf) {
        m_eventByte |= quint8(rf);
        return *this;
    }

private:
    quint8 m_eventByte;
};

inline QModbusCommEvent::EventByte operator|(QModbusCommEvent::EventByte b,
    QModbusCommEvent::SendFlag sf) { return QModbusCommEvent::EventByte(quint8(b) | quint8(sf)); }
inline QModbusCommEvent::EventByte operator|(QModbusCommEvent::SendFlag sf,
    QModbusCommEvent::EventByte b) { return operator|(b, sf); }
inline QModbusCommEvent::EventByte operator|(QModbusCommEvent::EventByte b,
    QModbusCommEvent::ReceiveFlag rf) { return QModbusCommEvent::EventByte(quint8(b) | quint8(rf)); }
inline QModbusCommEvent::EventByte operator|(QModbusCommEvent::ReceiveFlag rf,
    QModbusCommEvent::EventByte b) { return operator|(b, rf); }

QT_END_NAMESPACE

#endif // QMODBUSCOMMEVENT_P_H
