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

#ifndef QCANBUSFRAME_H
#define QCANBUSFRAME_H

#include <QtCore/qmetatype.h>
#include <QtCore/qobject.h>
#include <QtSerialBus/qtserialbusglobal.h>

QT_BEGIN_NAMESPACE

class QDataStream;

class Q_SERIALBUS_EXPORT QCanBusFrame
{
public:
    class TimeStamp {
    public:
        Q_DECL_CONSTEXPR TimeStamp(qint64 s = 0, qint64 usec = 0) Q_DECL_NOTHROW
            : secs(s), usecs(usec) {}

        Q_DECL_CONSTEXPR static TimeStamp fromMicroSeconds(qint64 usec) Q_DECL_NOTHROW
        { return TimeStamp(usec / 1000000, usec % 1000000); }

        Q_DECL_CONSTEXPR qint64 seconds() const Q_DECL_NOTHROW { return secs; }
        Q_DECL_CONSTEXPR qint64 microSeconds() const Q_DECL_NOTHROW { return usecs; }

    private:
        qint64 secs;
        qint64 usecs;
    };

    enum FrameType {
        UnknownFrame        = 0x0,
        DataFrame           = 0x1,
        ErrorFrame          = 0x2,
        RemoteRequestFrame  = 0x3,
        InvalidFrame        = 0x4
    };

    explicit QCanBusFrame(FrameType type = DataFrame) Q_DECL_NOTHROW :
        isExtendedFrame(0x0),
        version(Qt_5_10),
        isFlexibleDataRate(0x0),
        isBitrateSwitch(0x0),
        isErrorStateIndicator(0x0),
        isLocalEcho(0x0),
        reserved0(0x0)
    {
        Q_UNUSED(reserved0);
        ::memset(reserved, 0, sizeof(reserved));
        setFrameId(0x0);
        setFrameType(type);
    }

    enum FrameError {
        NoError                     = 0,
        TransmissionTimeoutError    = (1 << 0),
        LostArbitrationError        = (1 << 1),
        ControllerError             = (1 << 2),
        ProtocolViolationError      = (1 << 3),
        TransceiverError            = (1 << 4),
        MissingAcknowledgmentError  = (1 << 5),
        BusOffError                 = (1 << 6),
        BusError                    = (1 << 7),
        ControllerRestartError      = (1 << 8),
        UnknownError                = (1 << 9),
        AnyError                    = 0x1FFFFFFFU
        //only 29 bits usable
    };
    Q_DECLARE_FLAGS(FrameErrors, FrameError)
    Q_FLAGS(FrameErrors)

    explicit QCanBusFrame(quint32 identifier, const QByteArray &data) :
        format(DataFrame),
        isExtendedFrame(0x0),
        version(Qt_5_10),
        isFlexibleDataRate(data.length() > 8 ? 0x1 : 0x0),
        isBitrateSwitch(0x0),
        isErrorStateIndicator(0x0),
        isLocalEcho(0x0),
        reserved0(0x0),
        load(data)
    {
        ::memset(reserved, 0, sizeof(reserved));
        setFrameId(identifier);
    }

    bool isValid() const Q_DECL_NOTHROW
    {
        if (format == InvalidFrame)
            return false;

        // long id used, but extended flag not set
        if (!isExtendedFrame && (canId & 0x1FFFF800U))
            return false;

        if (!isValidFrameId)
            return false;

        // maximum permitted payload size in CAN or CAN FD
        const int length = load.length();
        if (isFlexibleDataRate) {
            if (format == RemoteRequestFrame)
                return false;

            return length <= 8 || length == 12 || length == 16 || length == 20
                    || length == 24 || length == 32 || length == 48 || length == 64;
        }

        return length <= 8;
    }

    FrameType frameType() const Q_DECL_NOTHROW
    {
        switch (format) {
        case 0x1: return DataFrame;
        case 0x2: return ErrorFrame;
        case 0x3: return RemoteRequestFrame;
        case 0x4: return InvalidFrame;
        // no default to trigger warning
        }

        return UnknownFrame;
    }

    void setFrameType(FrameType newFormat) Q_DECL_NOTHROW
    {
        switch (newFormat) {
        case DataFrame:
            format = 0x1; return;
        case ErrorFrame:
            format = 0x2; return;
        case RemoteRequestFrame:
            format = 0x3; return;
        case UnknownFrame:
            format = 0x0; return;
        case InvalidFrame:
            format = 0x4; return;
        }
    }

    bool hasExtendedFrameFormat() const Q_DECL_NOTHROW { return (isExtendedFrame & 0x1); }
    void setExtendedFrameFormat(bool isExtended) Q_DECL_NOTHROW
    {
        isExtendedFrame = (isExtended & 0x1);
    }

    quint32 frameId() const Q_DECL_NOTHROW
    {
        if (Q_UNLIKELY(format == ErrorFrame))
            return 0;
        return (canId & 0x1FFFFFFFU);
    }
    void setFrameId(quint32 newFrameId)
    {
        if (Q_LIKELY(newFrameId < 0x20000000U)) {
            isValidFrameId = true;
            canId = newFrameId;
            setExtendedFrameFormat(isExtendedFrame || (newFrameId & 0x1FFFF800U));
        } else {
            isValidFrameId = false;
            canId = 0;
        }
    }

    void setPayload(const QByteArray &data)
    {
        load = data;
        if (data.length() > 8)
            isFlexibleDataRate = 0x1;
    }
    void setTimeStamp(TimeStamp ts) Q_DECL_NOTHROW { stamp = ts; }

    QByteArray payload() const { return load; }
    TimeStamp timeStamp() const Q_DECL_NOTHROW { return stamp; }

    FrameErrors error() const Q_DECL_NOTHROW
    {
        if (format != ErrorFrame)
            return NoError;

        return FrameErrors(canId & 0x1FFFFFFFU);
    }
    void setError(FrameErrors e)
    {
        if (format != ErrorFrame)
            return;
        canId = (e & AnyError);
    }

    QString toString() const;

    bool hasFlexibleDataRateFormat() const Q_DECL_NOTHROW { return (isFlexibleDataRate & 0x1); }
    void setFlexibleDataRateFormat(bool isFlexibleData) Q_DECL_NOTHROW
    {
        isFlexibleDataRate = (isFlexibleData & 0x1);
        if (!isFlexibleData) {
            isBitrateSwitch = 0x0;
            isErrorStateIndicator = 0x0;
        }
    }

    bool hasBitrateSwitch() const Q_DECL_NOTHROW { return (isBitrateSwitch & 0x1); }
    void setBitrateSwitch(bool bitrateSwitch) Q_DECL_NOTHROW
    {
        isBitrateSwitch = (bitrateSwitch & 0x1);
        if (bitrateSwitch)
            isFlexibleDataRate = 0x1;
    }

    bool hasErrorStateIndicator() const Q_DECL_NOTHROW { return (isErrorStateIndicator & 0x1); }
    void setErrorStateIndicator(bool errorStateIndicator) Q_DECL_NOTHROW
    {
        isErrorStateIndicator = (errorStateIndicator & 0x1);
        if (errorStateIndicator)
            isFlexibleDataRate = 0x1;
    }
    bool hasLocalEcho() const Q_DECL_NOTHROW { return (isLocalEcho & 0x1); }
    void setLocalEcho(bool localEcho) Q_DECL_NOTHROW
    {
        isLocalEcho = (localEcho & 0x1);
    }

#ifndef QT_NO_DATASTREAM
    friend Q_SERIALBUS_EXPORT QDataStream &operator<<(QDataStream &, const QCanBusFrame &);
    friend Q_SERIALBUS_EXPORT QDataStream &operator>>(QDataStream &, QCanBusFrame &);
#endif

private:
    enum Version {
        Qt_5_8 = 0x0,
        Qt_5_9 = 0x1,
        Qt_5_10 = 0x2
    };

    quint32 canId:29; // acts as container for error codes too
    quint8 format:3; // max of 8 frame types

    quint8 isExtendedFrame:1;
    quint8 version:5;
    quint8 isValidFrameId:1;
    quint8 isFlexibleDataRate:1;

    quint8 isBitrateSwitch:1;
    quint8 isErrorStateIndicator:1;
    quint8 isLocalEcho:1;
    quint8 reserved0:5;

    // reserved for future use
    quint8 reserved[2];

    QByteArray load;
    TimeStamp stamp;
};

Q_DECLARE_TYPEINFO(QCanBusFrame, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(QCanBusFrame::FrameError, Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QCanBusFrame::FrameType, Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QCanBusFrame::TimeStamp, Q_PRIMITIVE_TYPE);

Q_DECLARE_OPERATORS_FOR_FLAGS(QCanBusFrame::FrameErrors)

#ifndef QT_NO_DATASTREAM
Q_SERIALBUS_EXPORT QDataStream &operator<<(QDataStream &, const QCanBusFrame &);
Q_SERIALBUS_EXPORT QDataStream &operator>>(QDataStream &, QCanBusFrame &);
#endif

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QCanBusFrame::FrameType)
Q_DECLARE_METATYPE(QCanBusFrame::FrameErrors)

#endif // QCANBUSFRAME_H
