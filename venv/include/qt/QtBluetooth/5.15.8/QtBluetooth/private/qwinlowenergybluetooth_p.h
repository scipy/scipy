/****************************************************************************
**
** Copyright (C) 2018 The Qt Company Ltd.
** Copyright (C) 2014 Denis Shienkov <denis.shienkov@gmail.com>
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

#ifndef QWINLOWENERGYBLUETOOTH_P_H
#define QWINLOWENERGYBLUETOOTH_P_H

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

#include <QtCore/qlibrary.h>

#include <qt_windows.h>

#define WIN32_FROM_HRESULT(hr)  \
    (SUCCEEDED(hr) ? ERROR_SUCCESS : \
    (HRESULT_FACILITY(hr) == FACILITY_WIN32 ? HRESULT_CODE(hr) : (hr)))

#define BLUETOOTH_GATT_FLAG_NONE                        0x00000000
#define BLUETOOTH_GATT_FLAG_CONNECTION_ENCRYPTED        0x00000001
#define BLUETOOTH_GATT_FLAG_CONNECTION_AUTHENTICATED    0x00000002
#define BLUETOOTH_GATT_FLAG_FORCE_READ_FROM_DEVICE      0x00000004
#define BLUETOOTH_GATT_FLAG_FORCE_READ_FROM_CACHE       0x00000008
#define BLUETOOTH_GATT_FLAG_SIGNED_WRITE                0x00000010
#define BLUETOOTH_GATT_FLAG_WRITE_WITHOUT_RESPONSE      0x00000020
#define BLUETOOTH_GATT_FLAG_RETURN_ALL                  0x00000040

typedef enum _BTH_LE_GATT_DESCRIPTOR_TYPE {
    CharacteristicExtendedProperties,
    CharacteristicUserDescription,
    ClientCharacteristicConfiguration,
    ServerCharacteristicConfiguration,
    CharacteristicFormat,
    CharacteristicAggregateFormat,
    CustomDescriptor
} BTH_LE_GATT_DESCRIPTOR_TYPE, *PBTH_LE_GATT_DESCRIPTOR_TYPE;

typedef enum _BTH_LE_GATT_EVENT_TYPE {
    CharacteristicValueChangedEvent
} BTH_LE_GATT_EVENT_TYPE;

typedef struct _BTH_LE_UUID {
    BOOLEAN IsShortUuid;
    union {
        USHORT ShortUuid;
        GUID LongUuid;
    } Value;
} BTH_LE_UUID, *PBTH_LE_UUID;

typedef struct _BTH_LE_GATT_SERVICE {
    BTH_LE_UUID ServiceUuid;
    USHORT AttributeHandle;
} BTH_LE_GATT_SERVICE, *PBTH_LE_GATT_SERVICE;

typedef struct _BTH_LE_GATT_CHARACTERISTIC {
    USHORT ServiceHandle;
    BTH_LE_UUID CharacteristicUuid;
    USHORT AttributeHandle;
    USHORT CharacteristicValueHandle;
    BOOLEAN IsBroadcastable;
    BOOLEAN IsReadable;
    BOOLEAN IsWritable;
    BOOLEAN IsWritableWithoutResponse;
    BOOLEAN IsSignedWritable;
    BOOLEAN IsNotifiable;
    BOOLEAN IsIndicatable;
    BOOLEAN HasExtendedProperties;
} BTH_LE_GATT_CHARACTERISTIC, *PBTH_LE_GATT_CHARACTERISTIC;

typedef struct _BTH_LE_GATT_CHARACTERISTIC_VALUE {
    ULONG DataSize;
    UCHAR Data[1];
} BTH_LE_GATT_CHARACTERISTIC_VALUE, *PBTH_LE_GATT_CHARACTERISTIC_VALUE;

typedef struct _BTH_LE_GATT_DESCRIPTOR {
    USHORT ServiceHandle;
    USHORT CharacteristicHandle;
    BTH_LE_GATT_DESCRIPTOR_TYPE DescriptorType;
    BTH_LE_UUID DescriptorUuid;
    USHORT AttributeHandle;
} BTH_LE_GATT_DESCRIPTOR, *PBTH_LE_GATT_DESCRIPTOR;

typedef struct _BTH_LE_GATT_DESCRIPTOR_VALUE {
    BTH_LE_GATT_DESCRIPTOR_TYPE DescriptorType;
    BTH_LE_UUID DescriptorUuid;
    union {
        struct {
            BOOLEAN IsReliableWriteEnabled;
            BOOLEAN IsAuxiliariesWritable;
        } CharacteristicExtendedProperties;
        struct {
            BOOLEAN IsSubscribeToNotification;
            BOOLEAN IsSubscribeToIndication;
        } ClientCharacteristicConfiguration;
        struct {
            BOOLEAN IsBroadcast;
        } ServerCharacteristicConfiguration;
        struct {
            UCHAR Format;
            UCHAR Exponent;
            BTH_LE_UUID Unit;
            UCHAR NameSpace;
            BTH_LE_UUID Description;
        } CharacteristicFormat;
    };
    ULONG DataSize;
    UCHAR Data[1];
} BTH_LE_GATT_DESCRIPTOR_VALUE, *PBTH_LE_GATT_DESCRIPTOR_VALUE;

typedef struct _BLUETOOTH_GATT_VALUE_CHANGED_EVENT {
    USHORT ChangedAttributeHandle;
    size_t CharacteristicValueDataSize;
    PBTH_LE_GATT_CHARACTERISTIC_VALUE CharacteristicValue;
} BLUETOOTH_GATT_VALUE_CHANGED_EVENT, *PBLUETOOTH_GATT_VALUE_CHANGED_EVENT;

typedef struct _BLUETOOTH_GATT_VALUE_CHANGED_EVENT_REGISTRATION {
    USHORT NumCharacteristics;
    BTH_LE_GATT_CHARACTERISTIC Characteristics[1];
} BLUETOOTH_GATT_VALUE_CHANGED_EVENT_REGISTRATION, *PBLUETOOTH_GATT_VALUE_CHANGED_EVENT_REGISTRATION;

typedef VOID (CALLBACK *PFNBLUETOOTH_GATT_EVENT_CALLBACK)(
        BTH_LE_GATT_EVENT_TYPE EventType,
        PVOID EventOutParameter,
        PVOID Context
        );

typedef ULONG64 BTH_LE_GATT_RELIABLE_WRITE_CONTEXT, *PBTH_LE_GATT_RELIABLE_WRITE_CONTEXT;

#define DEFINEFUNC(ret, func, ...) \
    typedef ret (WINAPI *fp_##func)(__VA_ARGS__); \
    static fp_##func func;

#define RESOLVEFUNC(func) \
    func = (fp_##func)resolveFunction(library, #func); \
    if (!func) \
        return false;

DEFINEFUNC(HRESULT, BluetoothGATTGetServices, HANDLE, USHORT, PBTH_LE_GATT_SERVICE, PUSHORT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTGetIncludedServices, HANDLE, PBTH_LE_GATT_SERVICE, USHORT, PBTH_LE_GATT_SERVICE, PUSHORT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTGetCharacteristics, HANDLE, PBTH_LE_GATT_SERVICE, USHORT, PBTH_LE_GATT_CHARACTERISTIC, PUSHORT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTGetDescriptors, HANDLE, PBTH_LE_GATT_CHARACTERISTIC, USHORT, PBTH_LE_GATT_DESCRIPTOR, PUSHORT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTGetCharacteristicValue, HANDLE, PBTH_LE_GATT_CHARACTERISTIC, ULONG, PBTH_LE_GATT_CHARACTERISTIC_VALUE, PUSHORT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTGetDescriptorValue, HANDLE, PBTH_LE_GATT_DESCRIPTOR, ULONG, PBTH_LE_GATT_DESCRIPTOR_VALUE, PUSHORT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTBeginReliableWrite, HANDLE, PBTH_LE_GATT_RELIABLE_WRITE_CONTEXT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTEndReliableWrite, HANDLE, BTH_LE_GATT_RELIABLE_WRITE_CONTEXT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTAbortReliableWrite, HANDLE, BTH_LE_GATT_RELIABLE_WRITE_CONTEXT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTSetCharacteristicValue, HANDLE, PBTH_LE_GATT_CHARACTERISTIC, PBTH_LE_GATT_CHARACTERISTIC_VALUE, BTH_LE_GATT_RELIABLE_WRITE_CONTEXT, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTSetDescriptorValue, HANDLE, PBTH_LE_GATT_DESCRIPTOR, PBTH_LE_GATT_DESCRIPTOR_VALUE, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTRegisterEvent, HANDLE, BTH_LE_GATT_EVENT_TYPE, PVOID, PFNBLUETOOTH_GATT_EVENT_CALLBACK, PVOID, PHANDLE, ULONG)
DEFINEFUNC(HRESULT, BluetoothGATTUnregisterEvent, HANDLE, ULONG)

static inline QFunctionPointer resolveFunction(QLibrary *library, const char *func)
{
    QFunctionPointer symbolFunctionPointer = library->resolve(func);
    if (!symbolFunctionPointer)
        qWarning("Cannot resolve '%s' in '%s'.", func, qPrintable(library->fileName()));
    return symbolFunctionPointer;
}

static inline bool resolveFunctions(QLibrary *library)
{
    if (!library->isLoaded()) {
        library->setFileName(QStringLiteral("bluetoothapis"));
        if (!library->load()) {
            qWarning("Unable to load '%s' library.", qPrintable(library->fileName()));
            return false;
        }
    }

    RESOLVEFUNC(BluetoothGATTGetServices)
    RESOLVEFUNC(BluetoothGATTGetIncludedServices)
    RESOLVEFUNC(BluetoothGATTGetCharacteristics)
    RESOLVEFUNC(BluetoothGATTGetDescriptors)
    RESOLVEFUNC(BluetoothGATTGetCharacteristicValue)
    RESOLVEFUNC(BluetoothGATTGetDescriptorValue)
    RESOLVEFUNC(BluetoothGATTBeginReliableWrite)
    RESOLVEFUNC(BluetoothGATTEndReliableWrite)
    RESOLVEFUNC(BluetoothGATTAbortReliableWrite)
    RESOLVEFUNC(BluetoothGATTSetCharacteristicValue)
    RESOLVEFUNC(BluetoothGATTSetDescriptorValue)
    RESOLVEFUNC(BluetoothGATTRegisterEvent)
    RESOLVEFUNC(BluetoothGATTUnregisterEvent)

    return true;
}

#endif // QWINLOWENERGYBLUETOOTH_P_H
