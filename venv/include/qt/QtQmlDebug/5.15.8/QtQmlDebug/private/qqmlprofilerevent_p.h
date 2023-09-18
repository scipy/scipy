/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQml module of the Qt Toolkit.
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

#ifndef QQMLPROFILEREVENT_P_H
#define QQMLPROFILEREVENT_P_H

#include "qqmlprofilerclientdefinitions_p.h"

#include <QtCore/qstring.h>
#include <QtCore/qbytearray.h>
#include <QtCore/qvarlengtharray.h>
#include <QtCore/qmetatype.h>

#include <initializer_list>
#include <limits>
#include <type_traits>

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

struct QQmlProfilerEvent {
    QQmlProfilerEvent() :
        m_timestamp(-1), m_typeIndex(-1), m_dataType(Inline8Bit), m_dataLength(0)
    {}

    template<typename Number>
    QQmlProfilerEvent(qint64 timestamp, int typeIndex, std::initializer_list<Number> list)
        : m_timestamp(timestamp), m_typeIndex(typeIndex)
    {
        assignNumbers<std::initializer_list<Number>, Number>(list);
    }

    QQmlProfilerEvent(qint64 timestamp, int typeIndex, const QString &data)
        : m_timestamp(timestamp), m_typeIndex(typeIndex)
    {
        assignNumbers<QByteArray, qint8>(data.toUtf8());
    }

    template<typename Number>
    QQmlProfilerEvent(qint64 timestamp, int typeIndex, const QVector<Number> &data)
        : m_timestamp(timestamp), m_typeIndex(typeIndex)
    {
        assignNumbers<QVector<Number>, Number>(data);
    }

    QQmlProfilerEvent(const QQmlProfilerEvent &other)
        : m_timestamp(other.m_timestamp), m_typeIndex(other.m_typeIndex),
          m_dataType(other.m_dataType), m_dataLength(other.m_dataLength)
    {
        assignData(other);
    }

    QQmlProfilerEvent(QQmlProfilerEvent &&other)
    {
        memcpy(static_cast<void *>(this), static_cast<const void *>(&other), sizeof(QQmlProfilerEvent));
        other.m_dataType = Inline8Bit; // prevent dtor from deleting the pointer
    }

    QQmlProfilerEvent &operator=(const QQmlProfilerEvent &other)
    {
        if (this != &other) {
            clearPointer();
            m_timestamp = other.m_timestamp;
            m_typeIndex = other.m_typeIndex;
            m_dataType = other.m_dataType;
            m_dataLength = other.m_dataLength;
            assignData(other);
        }
        return *this;
    }

    QQmlProfilerEvent &operator=(QQmlProfilerEvent &&other)
    {
        if (this != &other) {
            memcpy(static_cast<void *>(this), static_cast<const void *>(&other), sizeof(QQmlProfilerEvent));
            other.m_dataType = Inline8Bit;
        }
        return *this;
    }

    ~QQmlProfilerEvent()
    {
        clearPointer();
    }

    qint64 timestamp() const { return m_timestamp; }
    void setTimestamp(qint64 timestamp) { m_timestamp = timestamp; }

    int typeIndex() const { return m_typeIndex; }
    void setTypeIndex(int typeIndex) { m_typeIndex = typeIndex; }

    template<typename Number>
    Number number(int i) const
    {
        // Trailing zeroes can be omitted, for example for SceneGraph events
        if (i >= m_dataLength)
            return 0;
        switch (m_dataType) {
        case Inline8Bit:
            return m_data.internal8bit[i];
QT_WARNING_PUSH
QT_WARNING_DISABLE_GCC("-Warray-bounds") // Mingw 5.3 gcc doesn't get the type/length logic.
        case Inline16Bit:
            return m_data.internal16bit[i];
        case Inline32Bit:
            return m_data.internal32bit[i];
        case Inline64Bit:
            return m_data.internal64bit[i];
QT_WARNING_POP
        case External8Bit:
            return static_cast<const qint8 *>(m_data.external)[i];
        case External16Bit:
            return static_cast<const qint16 *>(m_data.external)[i];
        case External32Bit:
            return static_cast<const qint32 *>(m_data.external)[i];
        case External64Bit:
            return static_cast<const qint64 *>(m_data.external)[i];
        default:
            return 0;
        }
    }

    template<typename Number>
    void setNumber(int i, Number number)
    {
        QVarLengthArray<Number> nums = numbers<QVarLengthArray<Number>, Number>();
        int prevSize = nums.size();
        if (i >= prevSize) {
            nums.resize(i + 1);
            // Fill with zeroes. We don't want to accidentally prevent squeezing.
            while (prevSize < i)
                nums[prevSize++] = 0;
        }
        nums[i] = number;
        setNumbers<QVarLengthArray<Number>, Number>(nums);
    }

    template<typename Container, typename Number>
    void setNumbers(const Container &numbers)
    {
        clearPointer();
        assignNumbers<Container, Number>(numbers);
    }

    template<typename Number>
    void setNumbers(std::initializer_list<Number> numbers)
    {
        setNumbers<std::initializer_list<Number>, Number>(numbers);
    }

    template<typename Container, typename Number = qint64>
    Container numbers() const
    {
        Container container;
        for (int i = 0; i < m_dataLength; ++i)
            container.append(number<Number>(i));
        return container;
    }

    QString string() const
    {
        switch (m_dataType) {
        case External8Bit:
            return QString::fromUtf8(static_cast<const char *>(m_data.external), m_dataLength);
        case Inline8Bit:
            return QString::fromUtf8(m_data.internalChar, m_dataLength);
        default:
            Q_UNREACHABLE();
            return QString();
        }
    }

    void setString(const QString &data)
    {
        clearPointer();
        assignNumbers<QByteArray, char>(data.toUtf8());
    }

    Message rangeStage() const
    {
        Q_ASSERT(m_dataType == Inline8Bit);
        return static_cast<Message>(m_data.internal8bit[0]);
    }

    void setRangeStage(Message stage)
    {
        clearPointer();
        m_dataType = Inline8Bit;
        m_dataLength = 1;
        m_data.internal8bit[0] = stage;
    }

    bool isValid() const
    {
        return m_timestamp != -1;
    }

private:
    enum Type: quint16 {
        External = 1,
        Inline8Bit = 8,
        External8Bit = Inline8Bit | External,
        Inline16Bit = 16,
        External16Bit = Inline16Bit | External,
        Inline32Bit = 32,
        External32Bit = Inline32Bit | External,
        Inline64Bit = 64,
        External64Bit = Inline64Bit | External
    };

    qint64 m_timestamp;

    static const int s_internalDataLength = 8;
    union {
        void  *external;
        char   internalChar [s_internalDataLength];
        qint8  internal8bit [s_internalDataLength];
        qint16 internal16bit[s_internalDataLength / 2];
        qint32 internal32bit[s_internalDataLength / 4];
        qint64 internal64bit[s_internalDataLength / 8];
    } m_data;

    qint32 m_typeIndex;
    Type m_dataType;
    quint16 m_dataLength;

    void assignData(const QQmlProfilerEvent &other)
    {
        if (m_dataType & External) {
            uint length = m_dataLength * (other.m_dataType / 8);
            m_data.external = malloc(length);
            memcpy(m_data.external, other.m_data.external, length);
        } else {
            memcpy(&m_data, &other.m_data, sizeof(m_data));
        }
    }

    template<typename Big, typename Small>
    bool squeezable(Big source)
    {
        return static_cast<Small>(source) == source;
    }

    template<typename Container, typename Number>
    typename std::enable_if<(sizeof(Number) > 1), bool>::type
    squeeze(const Container &numbers)
    {
        typedef typename QIntegerForSize<sizeof(Number) / 2>::Signed Small;
        for (Number item : numbers) {
            if (!squeezable<Number, Small>(item))
                return false;
        }
        assignNumbers<Container, Small>(numbers);
        return true;
    }

    template<typename Container, typename Number>
    typename std::enable_if<(sizeof(Number) <= 1), bool>::type
    squeeze(const Container &)
    {
        return false;
    }

    template<typename Container, typename Number>
    void assignNumbers(const Container &numbers)
    {
        Number *data;
        m_dataLength = squeezable<size_t, quint16>(static_cast<size_t>(numbers.size())) ?
                    static_cast<quint16>(numbers.size()) : std::numeric_limits<quint16>::max();
        if (m_dataLength > sizeof(m_data) / sizeof(Number)) {
            if (squeeze<Container, Number>(numbers))
                return;
            m_dataType = static_cast<Type>((sizeof(Number) * 8) | External);
            m_data.external = malloc(m_dataLength * sizeof(Number));
            data = static_cast<Number *>(m_data.external);
        } else {
            m_dataType = static_cast<Type>(sizeof(Number) * 8);
            data = static_cast<Number *>(m_dataType & External ? m_data.external : &m_data);
        }
        quint16 i = 0;
        for (Number item : numbers) {
            if (i >= m_dataLength)
                break;
            data[i++] = item;
        }
    }

    void clearPointer()
    {
        if (m_dataType & External)
            free(m_data.external);
    }

    friend QDataStream &operator>>(QDataStream &stream, QQmlProfilerEvent &event);
    friend QDataStream &operator<<(QDataStream &stream, const QQmlProfilerEvent &event);
};

bool operator==(const QQmlProfilerEvent &event1, const QQmlProfilerEvent &event2);
bool operator!=(const QQmlProfilerEvent &event1, const QQmlProfilerEvent &event2);

QDataStream &operator>>(QDataStream &stream, QQmlProfilerEvent &event);
QDataStream &operator<<(QDataStream &stream, const QQmlProfilerEvent &event);

Q_DECLARE_TYPEINFO(QQmlProfilerEvent, Q_MOVABLE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QQmlProfilerEvent)

#endif // QQMLPROFILEREVENT_P_H
