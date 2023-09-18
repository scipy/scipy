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

#ifndef QLEADVERTISER_P_H
#define QLEADVERTISER_P_H

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

#include "qlowenergyadvertisingdata.h"
#include "qlowenergyadvertisingparameters.h"

#if QT_CONFIG(bluez)
#include "bluez/bluez_data_p.h"
#endif

#include <QtCore/qobject.h>
#include <QtCore/qvector.h>

QT_BEGIN_NAMESPACE

class QLeAdvertiser : public QObject
{
    Q_OBJECT
public:
    void startAdvertising() { doStartAdvertising(); }
    void stopAdvertising() { doStopAdvertising(); }

signals:
    void errorOccurred();

public:
    QLeAdvertiser(const QLowEnergyAdvertisingParameters &params,
                  const QLowEnergyAdvertisingData &advData,
                  const QLowEnergyAdvertisingData &responseData, QObject *parent)
        : QObject(parent), m_params(params), m_advData(advData), m_responseData(responseData) {}
    virtual ~QLeAdvertiser() { }

protected:
    const QLowEnergyAdvertisingParameters &parameters() const { return m_params; }
    const QLowEnergyAdvertisingData &advertisingData() const { return m_advData; }
    const QLowEnergyAdvertisingData &scanResponseData() const { return m_responseData; }

private:
    virtual void doStartAdvertising() = 0;
    virtual void doStopAdvertising() = 0;

    const QLowEnergyAdvertisingParameters m_params;
    const QLowEnergyAdvertisingData m_advData;
    const QLowEnergyAdvertisingData m_responseData;
};


#if QT_CONFIG(bluez)
struct AdvData;
struct AdvParams;
class HciManager;

class QLeAdvertiserBluez : public QLeAdvertiser
{
public:
    QLeAdvertiserBluez(const QLowEnergyAdvertisingParameters &params,
                       const QLowEnergyAdvertisingData &advertisingData,
                       const QLowEnergyAdvertisingData &scanResponseData, HciManager &hciManager,
                       QObject *parent = nullptr);
    ~QLeAdvertiserBluez() override;

private:
    void doStartAdvertising() override;
    void doStopAdvertising() override;

    void setPowerLevel(AdvData &advData);
    void setFlags(AdvData &advData);
    void setServicesData(const QLowEnergyAdvertisingData &src, AdvData &dest);
    void setManufacturerData(const QLowEnergyAdvertisingData &src, AdvData &dest);
    void setLocalNameData(const QLowEnergyAdvertisingData &src, AdvData &dest);

    void queueCommand(OpCodeCommandField ocf, const QByteArray &advertisingData);
    void sendNextCommand();
    void queueAdvertisingCommands();
    void queueReadTxPowerLevelCommand();
    void toggleAdvertising(bool enable);
    void setAdvertisingParams();
    void setAdvertisingInterval(AdvParams &params);
    void setData(bool isScanResponseData);
    void setAdvertisingData();
    void setScanResponseData();
    void setWhiteList();

    void handleCommandCompleted(quint16 opCode, quint8 status, const QByteArray &advertisingData);
    void handleError();

    HciManager &m_hciManager;

    struct Command {
        Command() {}
        Command(OpCodeCommandField ocf, const QByteArray &data) : ocf(ocf), data(data) { }
        OpCodeCommandField ocf;
        QByteArray data;
    };
    QVector<Command> m_pendingCommands;

    quint8 m_powerLevel;
    bool m_sendPowerLevel;
};
#endif // QT_CONFIG(bluez)

QT_END_NAMESPACE

#endif // Include guard.
