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

#ifndef QMODBUSDATAUNIT_H
#define QMODBUSDATAUNIT_H

#include <QtCore/qmap.h>
#include <QtCore/qmetatype.h>
#include <QtCore/qvector.h>

QT_BEGIN_NAMESPACE

class QModbusDataUnit
{
public:
    enum RegisterType {
        Invalid,
        DiscreteInputs,
        Coils,
        InputRegisters,
        HoldingRegisters
    };

    QModbusDataUnit() = default;

    explicit QModbusDataUnit(RegisterType type)
        : QModbusDataUnit(type, 0, 0)
    {}

    QModbusDataUnit(RegisterType type, int newStartAddress, quint16 newValueCount)
        : QModbusDataUnit(type, newStartAddress, QVector<quint16>(newValueCount))
    {}

    QModbusDataUnit(RegisterType type, int newStartAddress, const QVector<quint16> &newValues)
        : m_type(type)
        , m_startAddress(newStartAddress)
        , m_values(newValues)
        , m_valueCount(newValues.size())
    {}

    RegisterType registerType() const { return m_type; }
    void setRegisterType(RegisterType type) { m_type = type; }

    inline int startAddress() const { return m_startAddress; }
    inline void setStartAddress(int newAddress) { m_startAddress = newAddress; }

    inline QVector<quint16> values() const { return m_values; }
    inline void setValues(const QVector<quint16> &newValues)
    {
        m_values = newValues;
        m_valueCount = newValues.size();
    }

    inline uint valueCount() const { return m_valueCount; }
    inline void setValueCount(uint newCount) { m_valueCount = newCount; }

    inline void setValue(int index, quint16 newValue)
    {
        if (m_values.isEmpty() || index >= m_values.size())
            return;
        m_values[index] = newValue;
    }
    inline quint16 value(int index) const { return m_values.value(index); }

    bool isValid() const { return m_type != Invalid && m_startAddress != -1; }

private:
    RegisterType m_type = Invalid;
    int m_startAddress = -1;
    QVector<quint16> m_values;
    uint m_valueCount = 0;
};
typedef QMap<QModbusDataUnit::RegisterType, QModbusDataUnit> QModbusDataUnitMap;

Q_DECLARE_TYPEINFO(QModbusDataUnit, Q_MOVABLE_TYPE);
Q_DECLARE_TYPEINFO(QModbusDataUnit::RegisterType, Q_PRIMITIVE_TYPE);

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QModbusDataUnit::RegisterType)

#endif // QMODBUSDATAUNIT_H
