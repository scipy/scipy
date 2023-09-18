/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
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

#ifndef QLOWENERGYCONTROLLERBLUEZ_P_H
#define QLOWENERGYCONTROLLERBLUEZ_P_H

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

#include <qglobal.h>
#include <QtCore/QQueue>
#include <QtCore/QVector>
#include <QtBluetooth/qbluetooth.h>
#include <QtBluetooth/qlowenergycharacteristic.h>
#include "qlowenergycontroller.h"
#include "qlowenergycontrollerbase_p.h"

#include <QtBluetooth/QBluetoothSocket>
#include <functional>

QT_BEGIN_NAMESPACE

class QLowEnergyServiceData;
class QTimer;

class HciManager;
class LeCmacCalculator;
class QSocketNotifier;
class RemoteDeviceManager;

extern void registerQLowEnergyControllerMetaType();

class QLeAdvertiser;

class QLowEnergyControllerPrivateBluez final: public QLowEnergyControllerPrivate
{
    Q_OBJECT
public:
    QLowEnergyControllerPrivateBluez();
    ~QLowEnergyControllerPrivateBluez() override;

    void init() override;

    void connectToDevice() override;
    void disconnectFromDevice() override;

    void discoverServices() override;
    void discoverServiceDetails(const QBluetoothUuid &service) override;

    void startAdvertising(const QLowEnergyAdvertisingParameters &params,
                          const QLowEnergyAdvertisingData &advertisingData,
                          const QLowEnergyAdvertisingData &scanResponseData) override;
    void stopAdvertising() override;

    void requestConnectionUpdate(const QLowEnergyConnectionParameters &params) override;

    // read data
    void readCharacteristic(const QSharedPointer<QLowEnergyServicePrivate> service,
                            const QLowEnergyHandle charHandle) override;
    void readDescriptor(const QSharedPointer<QLowEnergyServicePrivate> service,
                        const QLowEnergyHandle charHandle,
                        const QLowEnergyHandle descriptorHandle) override;

    // write data
    void writeCharacteristic(const QSharedPointer<QLowEnergyServicePrivate> service,
                             const QLowEnergyHandle charHandle,
                             const QByteArray &newValue, QLowEnergyService::WriteMode mode) override;
    void writeDescriptor(const QSharedPointer<QLowEnergyServicePrivate> service,
                         const QLowEnergyHandle charHandle,
                         const QLowEnergyHandle descriptorHandle,
                         const QByteArray &newValue) override;

    void addToGenericAttributeList(const QLowEnergyServiceData &service,
                                   QLowEnergyHandle startHandle) override;

    struct Attribute {
        Attribute() : handle(0) {}

        QLowEnergyHandle handle;
        QLowEnergyHandle groupEndHandle;
        QLowEnergyCharacteristic::PropertyTypes properties;
        QBluetooth::AttAccessConstraints readConstraints;
        QBluetooth::AttAccessConstraints writeConstraints;
        QBluetoothUuid type;
        QByteArray value;
        int minLength;
        int maxLength;
    };
    QVector<Attribute> localAttributes;

private:
    quint16 connectionHandle = 0;
    QBluetoothSocket *l2cpSocket = nullptr;
    struct Request {
        quint8 command;
        QByteArray payload;
        // TODO reference below is ugly but until we know all commands and their
        // requirements this is WIP
        QVariant reference;
        QVariant reference2;
    };
    QQueue<Request> openRequests;

    struct WriteRequest {
        WriteRequest() {}
        WriteRequest(quint16 h, quint16 o, const QByteArray &v)
            : handle(h), valueOffset(o), value(v) {}
        quint16 handle;
        quint16 valueOffset;
        QByteArray value;
    };
    QVector<WriteRequest> openPrepareWriteRequests;

    // Invariant: !scheduledIndications.isEmpty => indicationInFlight == true
    QVector<QLowEnergyHandle> scheduledIndications;
    bool indicationInFlight = false;

    struct TempClientConfigurationData {
        TempClientConfigurationData(QLowEnergyServicePrivate::DescData *dd = nullptr,
                                    QLowEnergyHandle chHndl = 0, QLowEnergyHandle coHndl = 0)
            : descData(dd), charValueHandle(chHndl), configHandle(coHndl) {}

        QLowEnergyServicePrivate::DescData *descData;
        QLowEnergyHandle charValueHandle;
        QLowEnergyHandle configHandle;
    };

    struct ClientConfigurationData {
        ClientConfigurationData(QLowEnergyHandle chHndl = 0, QLowEnergyHandle coHndl = 0,
                                quint16 val = 0)
            : charValueHandle(chHndl), configHandle(coHndl), configValue(val) {}

        QLowEnergyHandle charValueHandle;
        QLowEnergyHandle configHandle;
        quint16 configValue;
        bool charValueWasUpdated = false;
    };
    QHash<quint64, QVector<ClientConfigurationData>> clientConfigData;

    struct SigningData {
        SigningData() = default;
        SigningData(const quint128 &csrk, quint32 signCounter = quint32(-1))
            : key(csrk), counter(signCounter) {}

        quint128 key;
        quint32 counter = quint32(-1);
    };
    QHash<quint64, SigningData> signingData;
    LeCmacCalculator *cmacCalculator = nullptr;

    bool requestPending;
    quint16 mtuSize;
    int securityLevelValue;
    bool encryptionChangePending;
    bool receivedMtuExchangeRequest = false;

    HciManager *hciManager = nullptr;
    QLeAdvertiser *advertiser = nullptr;
    QSocketNotifier *serverSocketNotifier = nullptr;
    QTimer *requestTimer = nullptr;
    RemoteDeviceManager* device1Manager = nullptr;

    /*
      Defines the maximum number of milliseconds the implementation will
      wait for requests that require a response.

      This addresses the problem that some non-conformant BTLE devices
      do not implement the request/response system properly. In such cases
      the queue system would hang forever.

      Once timeout has been triggered we gracefully continue with the next request.
      Depending on the type of the timed out ATT command we either ignore it
      or artifically trigger an error response to ensure the API gives the
      appropriate response. Potentially this can cause problems when the
      response for the dropped requests arrives very late. That's why a big warning
      is printed about the compromised state when a timeout is triggered.
     */
    int gattRequestTimeout = 20000;

    void handleConnectionRequest();
    void closeServerSocket();

    bool isBonded() const;
    QVector<TempClientConfigurationData> gatherClientConfigData();
    void storeClientConfigurations();
    void restoreClientConfigurations();

    enum SigningKeyType { LocalSigningKey, RemoteSigningKey };
    void loadSigningDataIfNecessary(SigningKeyType keyType);
    void storeSignCounter(SigningKeyType keyType) const;
    QString signingKeySettingsGroup(SigningKeyType keyType) const;
    QString keySettingsFilePath() const;

    void sendPacket(const QByteArray &packet);
    void sendNextPendingRequest();
    void processReply(const Request &request, const QByteArray &reply);

    void sendReadByGroupRequest(QLowEnergyHandle start, QLowEnergyHandle end,
                                quint16 type);
    void sendReadByTypeRequest(QSharedPointer<QLowEnergyServicePrivate> serviceData,
                               QLowEnergyHandle nextHandle, quint16 attributeType);
    void sendReadValueRequest(QLowEnergyHandle attributeHandle, bool isDescriptor);
    void readServiceValues(const QBluetoothUuid &service,
                           bool readCharacteristics);
    void readServiceValuesByOffset(uint handleData, quint16 offset,
                                   bool isLastValue);

    void discoverServiceDescriptors(const QBluetoothUuid &serviceUuid);
    void discoverNextDescriptor(QSharedPointer<QLowEnergyServicePrivate> serviceData,
                                const QList<QLowEnergyHandle> pendingCharHandles,
                                QLowEnergyHandle startingHandle);
    void processUnsolicitedReply(const QByteArray &msg);
    void exchangeMTU();
    bool setSecurityLevel(int level);
    int securityLevel() const;
    void sendExecuteWriteRequest(const QLowEnergyHandle attrHandle,
                                 const QByteArray &newValue,
                                 bool isCancelation);
    void sendNextPrepareWriteRequest(const QLowEnergyHandle handle,
                                     const QByteArray &newValue, quint16 offset);
    bool increaseEncryptLevelfRequired(quint8 errorCode);

    void resetController();

    void handleAdvertisingError();

    bool checkPacketSize(const QByteArray &packet, int minSize, int maxSize = -1);
    bool checkHandle(const QByteArray &packet, QLowEnergyHandle handle);
    bool checkHandlePair(quint8 request, QLowEnergyHandle startingHandle,
                         QLowEnergyHandle endingHandle);

    void handleExchangeMtuRequest(const QByteArray &packet);
    void handleFindInformationRequest(const QByteArray &packet);
    void handleFindByTypeValueRequest(const QByteArray &packet);
    void handleReadByTypeRequest(const QByteArray &packet);
    void handleReadRequest(const QByteArray &packet);
    void handleReadBlobRequest(const QByteArray &packet);
    void handleReadMultipleRequest(const QByteArray &packet);
    void handleReadByGroupTypeRequest(const QByteArray &packet);
    void handleWriteRequestOrCommand(const QByteArray &packet);
    void handlePrepareWriteRequest(const QByteArray &packet);
    void handleExecuteWriteRequest(const QByteArray &packet);

    void sendErrorResponse(quint8 request, quint16 handle, quint8 code);

    using ElemWriter = std::function<void(const Attribute &, char *&)>;
    void sendListResponse(const QByteArray &packetStart, int elemSize,
                          const QVector<Attribute> &attributes, const ElemWriter &elemWriter);

    void sendNotification(QLowEnergyHandle handle);
    void sendIndication(QLowEnergyHandle handle);
    void sendNotificationOrIndication(quint8 opCode, QLowEnergyHandle handle);
    void sendNextIndication();

    void ensureUniformAttributes(QVector<Attribute> &attributes, const std::function<int(const Attribute &)> &getSize);
    void ensureUniformUuidSizes(QVector<Attribute> &attributes);
    void ensureUniformValueSizes(QVector<Attribute> &attributes);

    using AttributePredicate = std::function<bool(const Attribute &)>;
    QVector<Attribute> getAttributes(QLowEnergyHandle startHandle, QLowEnergyHandle endHandle,
            const AttributePredicate &attributePredicate = [](const Attribute &) { return true; });

    int checkPermissions(const Attribute &attr, QLowEnergyCharacteristic::PropertyType type);
    int checkReadPermissions(const Attribute &attr);
    int checkReadPermissions(QVector<Attribute> &attributes);

    bool verifyMac(const QByteArray &message, const quint128 &csrk, quint32 signCounter,
                   quint64 expectedMac);

    void updateLocalAttributeValue(
            QLowEnergyHandle handle,
            const QByteArray &value,
            QLowEnergyCharacteristic &characteristic,
            QLowEnergyDescriptor &descriptor);

    void writeCharacteristicForPeripheral(
            QLowEnergyServicePrivate::CharData &charData,
            const QByteArray &newValue);
    void writeCharacteristicForCentral(const QSharedPointer<QLowEnergyServicePrivate> &service,
            QLowEnergyHandle charHandle,
            QLowEnergyHandle valueHandle,
            const QByteArray &newValue,
            QLowEnergyService::WriteMode mode);

    void writeDescriptorForPeripheral(
            const QSharedPointer<QLowEnergyServicePrivate> &service,
            const QLowEnergyHandle charHandle,
            const QLowEnergyHandle descriptorHandle,
            const QByteArray &newValue);
    void writeDescriptorForCentral(
            const QLowEnergyHandle charHandle,
            const QLowEnergyHandle descriptorHandle,
            const QByteArray &newValue);

    void restartRequestTimer();
    void establishL2cpClientSocket();
    void createServicesForCentralIfRequired();

private slots:
    void l2cpConnected();
    void l2cpDisconnected();
    void l2cpErrorChanged(QBluetoothSocket::SocketError);
    void l2cpReadyRead();
    void encryptionChangedEvent(const QBluetoothAddress&, bool);
    void handleGattRequestTimeout();
    void activeConnectionTerminationDone();
};

Q_DECLARE_TYPEINFO(QLowEnergyControllerPrivateBluez::Attribute, Q_MOVABLE_TYPE);

QT_END_NAMESPACE

#endif //QLOWENERGYCONTROLLERBLUEZ_P_H
