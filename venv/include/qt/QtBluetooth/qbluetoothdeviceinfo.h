/****************************************************************************
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

#ifndef QBLUETOOTHDEVICEINFO_H
#define QBLUETOOTHDEVICEINFO_H

#include <QtBluetooth/qtbluetoothglobal.h>

#include <QtCore/qstring.h>
#include <QtCore/qmetatype.h>
#include <QtCore/qbytearray.h>
#include <QtCore/qvector.h>

QT_BEGIN_NAMESPACE

class QBluetoothDeviceInfoPrivate;
class QBluetoothAddress;
class QBluetoothUuid;

class Q_BLUETOOTH_EXPORT QBluetoothDeviceInfo
{
public:
    enum MajorDeviceClass {
        MiscellaneousDevice = 0,
        ComputerDevice = 1,
        PhoneDevice = 2,
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
        LANAccessDevice = 3,
#endif
        NetworkDevice = 3,
        AudioVideoDevice = 4,
        PeripheralDevice = 5,
        ImagingDevice = 6,
        WearableDevice = 7,
        ToyDevice = 8,
        HealthDevice = 9,
        UncategorizedDevice = 31
    };

    enum MinorMiscellaneousClass {
        UncategorizedMiscellaneous = 0
    };

    enum MinorComputerClass {
        UncategorizedComputer = 0,
        DesktopComputer = 1,
        ServerComputer = 2,
        LaptopComputer = 3,
        HandheldClamShellComputer = 4,
        HandheldComputer = 5,
        WearableComputer = 6
    };

    enum MinorPhoneClass {
        UncategorizedPhone = 0,
        CellularPhone = 1,
        CordlessPhone = 2,
        SmartPhone = 3,
        WiredModemOrVoiceGatewayPhone = 4,
        CommonIsdnAccessPhone = 5
    };

    enum MinorNetworkClass {
        NetworkFullService = 0x00,
        NetworkLoadFactorOne = 0x08,
        NetworkLoadFactorTwo = 0x10,
        NetworkLoadFactorThree = 0x18,
        NetworkLoadFactorFour = 0x20,
        NetworkLoadFactorFive = 0x28,
        NetworkLoadFactorSix = 0x30,
        NetworkNoService = 0x38
    };

    enum MinorAudioVideoClass {
        UncategorizedAudioVideoDevice = 0,
        WearableHeadsetDevice = 1,
        HandsFreeDevice = 2,
        // reserved = 3,
        Microphone = 4,
        Loudspeaker = 5,
        Headphones = 6,
        PortableAudioDevice = 7,
        CarAudio = 8,
        SetTopBox = 9,
        HiFiAudioDevice = 10,
        Vcr = 11,
        VideoCamera = 12,
        Camcorder = 13,
        VideoMonitor = 14,
        VideoDisplayAndLoudspeaker = 15,
        VideoConferencing = 16,
        // reserved = 17,
        GamingDevice = 18
    };

    enum MinorPeripheralClass {
        UncategorizedPeripheral = 0,
        KeyboardPeripheral = 0x10,
        PointingDevicePeripheral = 0x20,
        KeyboardWithPointingDevicePeripheral = 0x30,

        JoystickPeripheral = 0x01,
        GamepadPeripheral = 0x02,
        RemoteControlPeripheral = 0x03,
        SensingDevicePeripheral = 0x04,
        DigitizerTabletPeripheral = 0x05,
        CardReaderPeripheral = 0x06
    };

    enum MinorImagingClass {
        UncategorizedImagingDevice = 0,
        ImageDisplay = 0x04,
        ImageCamera = 0x08,
        ImageScanner = 0x10,
        ImagePrinter = 0x20
    };

    enum MinorWearableClass {
        UncategorizedWearableDevice = 0,
        WearableWristWatch = 1,
        WearablePager = 2,
        WearableJacket = 3,
        WearableHelmet = 4,
        WearableGlasses = 5
    };

    enum MinorToyClass {
        UncategorizedToy = 0,
        ToyRobot = 1,
        ToyVehicle = 2,
        ToyDoll = 3,
        ToyController = 4,
        ToyGame = 5
    };

    enum MinorHealthClass {
        UncategorizedHealthDevice = 0,
        HealthBloodPressureMonitor = 0x1,
        HealthThermometer = 0x2,
        HealthWeightScale = 0x3,
        HealthGlucoseMeter = 0x4,
        HealthPulseOximeter = 0x5,
        HealthDataDisplay = 0x7,
        HealthStepCounter = 0x8
    };

    enum ServiceClass {
        NoService = 0x0000,
        PositioningService = 0x0001,
        NetworkingService = 0x0002,
        RenderingService = 0x0004,
        CapturingService = 0x0008,
        ObjectTransferService = 0x0010,
        AudioService = 0x0020,
        TelephonyService = 0x0040,
        InformationService = 0x0080,
        AllServices = 0x07ff
    };
    Q_DECLARE_FLAGS(ServiceClasses, ServiceClass)

#if QT_DEPRECATED_SINCE(5, 13)
    // adding QT_DEPRECATED causes compile failure with gcc 7
    enum DataCompleteness {
        DataComplete,
        DataIncomplete,
        DataUnavailable
    };
#endif

    enum class Field {
        None = 0x0000,
        RSSI = 0x0001,
        ManufacturerData = 0x0002,
        All = 0x7fff
    };
    Q_DECLARE_FLAGS(Fields, Field)

    enum CoreConfiguration {
        UnknownCoreConfiguration = 0x0,
        LowEnergyCoreConfiguration = 0x01,
        BaseRateCoreConfiguration = 0x02,
        BaseRateAndLowEnergyCoreConfiguration = 0x03
    };
    Q_DECLARE_FLAGS(CoreConfigurations, CoreConfiguration)

    QBluetoothDeviceInfo();
    QBluetoothDeviceInfo(const QBluetoothAddress &address, const QString &name,
                         quint32 classOfDevice);
    QBluetoothDeviceInfo(const QBluetoothUuid &uuid, const QString &name,
                         quint32 classOfDevice);
    QBluetoothDeviceInfo(const QBluetoothDeviceInfo &other);
    ~QBluetoothDeviceInfo();

    bool isValid() const;
    bool isCached() const;

    void setCached(bool cached);

    QBluetoothDeviceInfo &operator=(const QBluetoothDeviceInfo &other);
    bool operator==(const QBluetoothDeviceInfo &other) const;
    bool operator!=(const QBluetoothDeviceInfo &other) const;

    QBluetoothAddress address() const;
    QString name() const;

    ServiceClasses serviceClasses() const;
    MajorDeviceClass majorDeviceClass() const;
    quint8 minorDeviceClass() const;

    qint16 rssi() const;
    void setRssi(qint16 signal);

#if QT_DEPRECATED_SINCE(5, 13)
    QT_DEPRECATED void setServiceUuids(const QList<QBluetoothUuid> &uuids, DataCompleteness completeness);
    QT_DEPRECATED DataCompleteness serviceUuidsCompleteness() const;
#endif

#if QT_VERSION >= QT_VERSION_CHECK(6, 0, 0)
#ifndef Q_QDOC //suppress qdoc warnings
    QVector<QBluetoothUuid> serviceUuids() const;
#endif // Q_QDOC
#elif QT_DEPRECATED_SINCE(5, 13)
    QList<QBluetoothUuid> serviceUuids(DataCompleteness *completeness = nullptr) const;
#else
    QList<QBluetoothUuid> serviceUuids() const;
#endif
    void setServiceUuids(const QVector<QBluetoothUuid> &uuids);

    // TODO Qt6 manufacturerData() need to be changed to return
    // QMultiHash<quint16, QByteArray>
    QVector<quint16> manufacturerIds() const;
    QByteArray manufacturerData(quint16 manufacturerId) const;
    bool setManufacturerData(quint16 manufacturerId, const QByteArray &data);
    QHash<quint16, QByteArray> manufacturerData() const;

    void setCoreConfigurations(QBluetoothDeviceInfo::CoreConfigurations coreConfigs);
    QBluetoothDeviceInfo::CoreConfigurations coreConfigurations() const;

    void setDeviceUuid(const QBluetoothUuid &uuid);
    QBluetoothUuid deviceUuid() const;

protected:
    QBluetoothDeviceInfoPrivate *d_ptr;

private:
    Q_DECLARE_PRIVATE(QBluetoothDeviceInfo)
};

Q_DECLARE_OPERATORS_FOR_FLAGS(QBluetoothDeviceInfo::CoreConfigurations)
Q_DECLARE_OPERATORS_FOR_FLAGS(QBluetoothDeviceInfo::ServiceClasses)

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QBluetoothDeviceInfo)
#ifdef QT_WINRT_BLUETOOTH
Q_DECLARE_METATYPE(QBluetoothDeviceInfo::Fields)
#endif

#endif
