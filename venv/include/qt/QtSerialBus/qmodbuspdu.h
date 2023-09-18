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
#ifndef QMODBUSPDU_H
#define QMODBUSPDU_H

#include <QtCore/qdatastream.h>
#include <QtCore/qmetatype.h>
#include <QtCore/qvector.h>
#include <QtSerialBus/qtserialbusglobal.h>

QT_BEGIN_NAMESPACE

class QModbusPdu
{
public:
    enum ExceptionCode {
        IllegalFunction = 0x01,
        IllegalDataAddress = 0x02,
        IllegalDataValue = 0x03,
        ServerDeviceFailure = 0x04,
        Acknowledge = 0x05,
        ServerDeviceBusy = 0x06,
        NegativeAcknowledge = 0x07,
        MemoryParityError = 0x08,
        GatewayPathUnavailable = 0x0A,
        GatewayTargetDeviceFailedToRespond = 0x0B,
        ExtendedException = 0xFF,
    };

    enum FunctionCode {
        Invalid = 0x00,
        ReadCoils = 0x01,
        ReadDiscreteInputs = 0x02,
        ReadHoldingRegisters = 0x03,
        ReadInputRegisters = 0x04,
        WriteSingleCoil = 0x05,
        WriteSingleRegister = 0x06,
        ReadExceptionStatus = 0x07,
        Diagnostics = 0x08,
        GetCommEventCounter = 0x0B,
        GetCommEventLog = 0x0C,
        WriteMultipleCoils = 0x0F,
        WriteMultipleRegisters = 0x10,
        ReportServerId = 0x11,
        ReadFileRecord = 0x14,
        WriteFileRecord = 0x15,
        MaskWriteRegister = 0x16,
        ReadWriteMultipleRegisters = 0x17,
        ReadFifoQueue = 0x18,
        EncapsulatedInterfaceTransport = 0x2B,
        UndefinedFunctionCode = 0x100
    };

    QModbusPdu() = default;
    virtual ~QModbusPdu() = default;

    bool isValid() const {
        return (m_code >= ReadCoils && m_code < UndefinedFunctionCode)
                && (m_data.size() < 253);
    }

    static const quint8 ExceptionByte = 0x80;
    ExceptionCode exceptionCode() const {
        if (!m_data.size() || !isException())
            return ExtendedException;
        return static_cast<ExceptionCode>(m_data.at(0));
    }
    bool isException() const { return m_code & ExceptionByte; }

    qint16 size() const { return dataSize() + 1; }
    qint16 dataSize() const { return qint16(m_data.size()); }

    FunctionCode functionCode() const {
        return FunctionCode(quint8(m_code) &~ ExceptionByte);
    }
    virtual void setFunctionCode(FunctionCode code) { m_code = code; }

    QByteArray data() const { return m_data; }
    void setData(const QByteArray &newData) { m_data = newData; }

    template <typename ... Args> void encodeData(Args ... newData) {
        encode(std::forward<Args>(newData)...);
    }

    template <typename ... Args> void decodeData(Args && ... newData) const {
        decode(std::forward<Args>(newData)...);
    }

protected:
    QModbusPdu(FunctionCode code, const QByteArray &newData)
        : m_code(code)
        , m_data(newData)
    {}

    QModbusPdu(const QModbusPdu &) = default;
    QModbusPdu &operator=(const QModbusPdu &) = default;

    template <typename ... Args>
    QModbusPdu(FunctionCode code, Args ... newData)
        : m_code(code)
    {
        encode(std::forward<Args>(newData)...);
    }

private:
    template <typename T, typename ... Ts> struct IsType { enum { value = false }; };
    template <typename T, typename T1, typename ... Ts> struct IsType<T, T1, Ts...> {
        enum { value = std::is_same<T, T1>::value || IsType<T, Ts...>::value };
    };

    template <typename T>
    using is_pod = std::integral_constant<bool, std::is_trivial<T>::value && std::is_standard_layout<T>::value>;

    template <typename T> void encode(QDataStream *stream, const T &t) {
        static_assert(is_pod<T>::value, "Only POD types supported.");
        static_assert(IsType<T, quint8, quint16>::value, "Only quint8 and quint16 supported.");
        (*stream) << t;
    }
    template <typename T> void decode(QDataStream *stream, T &t) const {
        static_assert(is_pod<T>::value, "Only POD types supported.");
        static_assert(IsType<T, quint8 *, quint16 *>::value, "Only quint8* and quint16* supported.");
        (*stream) >> *t;
    }
    template <typename T> void encode(QDataStream *stream, const QVector<T> &vector) {
        static_assert(is_pod<T>::value, "Only POD types supported.");
        static_assert(IsType<T, quint8, quint16>::value, "Only quint8 and quint16 supported.");
        for (int i = 0; i < vector.count(); ++i)
            (*stream) << vector[i];
    }

    template<typename ... Args> void encode(Args ... newData) {
        m_data.clear();
        Q_CONSTEXPR quint32 argCount = sizeof...(Args);
        if (argCount > 0) {
            QDataStream stream(&m_data, QIODevice::WriteOnly);
            char tmp[argCount] = { (encode(&stream, newData), void(), '0')... };
            Q_UNUSED(tmp)
        }
    }
    template<typename ... Args> void decode(Args ... newData) const {
        Q_CONSTEXPR quint32 argCount = sizeof...(Args);
        if (argCount > 0 && !m_data.isEmpty()) {
            QDataStream stream(m_data);
            char tmp[argCount] = { (decode(&stream, newData), void(), '0')... };
            Q_UNUSED(tmp)
        }
    }

private:
    FunctionCode m_code = Invalid;
    QByteArray m_data;
    friend class QModbusSerialAdu;
    friend struct QModbusPduPrivate;
};
Q_SERIALBUS_EXPORT QDebug operator<<(QDebug debug, const QModbusPdu &pdu);
Q_SERIALBUS_EXPORT QDataStream &operator<<(QDataStream &stream, const QModbusPdu &pdu);

class QModbusRequest : public QModbusPdu
{
public:
    QModbusRequest() = default;
    QModbusRequest(const QModbusPdu &pdu)
        : QModbusPdu(pdu)
    {}

    explicit QModbusRequest(FunctionCode code, const QByteArray &newData = QByteArray())
        : QModbusPdu(code, newData)
    {}

    Q_SERIALBUS_EXPORT static int minimumDataSize(const QModbusRequest &pdu);
    Q_SERIALBUS_EXPORT static int calculateDataSize(const QModbusRequest &pdu);

    using CalcFuncPtr = decltype(&calculateDataSize);
    Q_SERIALBUS_EXPORT static void registerDataSizeCalculator(FunctionCode fc, CalcFuncPtr func);

    template <typename ... Args>
    QModbusRequest(FunctionCode code, Args ... newData)
        : QModbusPdu(code, newData...)
    {}
};
Q_SERIALBUS_EXPORT QDataStream &operator>>(QDataStream &stream, QModbusRequest &pdu);
inline QDataStream &operator<<(QDataStream &stream, const QModbusRequest &pdu)
{ return stream << static_cast<const QModbusPdu &>(pdu); }

class QModbusResponse : public QModbusPdu
{
public:
    QModbusResponse() = default;
    QModbusResponse(const QModbusPdu &pdu)
        : QModbusPdu(pdu)
    {}

    explicit QModbusResponse(FunctionCode code, const QByteArray &newData = QByteArray())
        : QModbusPdu(code, newData)
    {}

    Q_SERIALBUS_EXPORT static int minimumDataSize(const QModbusResponse &pdu);
    Q_SERIALBUS_EXPORT static int calculateDataSize(const QModbusResponse &pdu);

    using CalcFuncPtr = decltype(&calculateDataSize);
    Q_SERIALBUS_EXPORT static void registerDataSizeCalculator(FunctionCode fc, CalcFuncPtr func);

    template <typename ... Args>
    QModbusResponse(FunctionCode code, Args ... newData)
        : QModbusPdu(code, newData...)
    {}
};

class QModbusExceptionResponse : public QModbusResponse
{
public:
    QModbusExceptionResponse() = default;
    QModbusExceptionResponse(const QModbusPdu &pdu)
        : QModbusResponse(pdu)
    {}
    QModbusExceptionResponse(FunctionCode fc, ExceptionCode ec)
        : QModbusResponse(FunctionCode(quint8(fc) | ExceptionByte), static_cast<quint8> (ec))
    {}

    void setFunctionCode(FunctionCode c) {
        QModbusPdu::setFunctionCode(FunctionCode(quint8(c) | ExceptionByte));
    }
    void setExceptionCode(ExceptionCode ec) { QModbusPdu::encodeData(quint8(ec)); }
};
Q_SERIALBUS_EXPORT QDataStream &operator>>(QDataStream &stream, QModbusResponse &pdu);
inline QDataStream &operator<<(QDataStream &stream, const QModbusResponse &pdu)
{ return stream << static_cast<const QModbusPdu &>(pdu); }

Q_DECLARE_TYPEINFO(QModbusPdu, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(QModbusPdu::ExceptionCode, Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QModbusPdu::FunctionCode, Q_PRIMITIVE_TYPE);

Q_DECLARE_TYPEINFO(QModbusRequest, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(QModbusResponse, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(QModbusExceptionResponse, Q_MOVABLE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QModbusPdu::ExceptionCode)
Q_DECLARE_METATYPE(QModbusPdu::FunctionCode)

#endif // QMODBUSPDU_H
