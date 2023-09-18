/***************************************************************************
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

#ifndef QLOWENERGYCONNECTIONPARAMETERS_H
#define QLOWENERGYCONNECTIONPARAMETERS_H

#include <QtBluetooth/qtbluetoothglobal.h>
#include <QtCore/qmetatype.h>
#include <QtCore/qshareddata.h>

QT_BEGIN_NAMESPACE

class QLowEnergyConnectionParametersPrivate;

class Q_BLUETOOTH_EXPORT QLowEnergyConnectionParameters
{
    friend Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyConnectionParameters &p1,
                                              const QLowEnergyConnectionParameters &p2);
public:
    QLowEnergyConnectionParameters();
    QLowEnergyConnectionParameters(const QLowEnergyConnectionParameters &other);
    ~QLowEnergyConnectionParameters();

    QLowEnergyConnectionParameters &operator=(const QLowEnergyConnectionParameters &other);

    void setIntervalRange(double minimum, double maximum);
    double minimumInterval() const;
    double maximumInterval() const;

    void setLatency(int latency);
    int latency() const;

    void setSupervisionTimeout(int timeout);
    int supervisionTimeout() const;

    void swap(QLowEnergyConnectionParameters &other) Q_DECL_NOTHROW { qSwap(d, other.d); }

private:
    QSharedDataPointer<QLowEnergyConnectionParametersPrivate> d;
};

Q_BLUETOOTH_EXPORT bool operator==(const QLowEnergyConnectionParameters &p1,
                                   const QLowEnergyConnectionParameters &p2);
inline bool operator!=(const QLowEnergyConnectionParameters &p1,
                       const QLowEnergyConnectionParameters &p2)
{
    return !(p1 == p2);
}

Q_DECLARE_SHARED(QLowEnergyConnectionParameters)

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QLowEnergyConnectionParameters)

#endif // Include guard
