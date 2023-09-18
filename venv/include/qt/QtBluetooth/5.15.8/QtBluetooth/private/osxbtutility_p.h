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

#ifndef OSXBTUTILITY_P_H
#define OSXBTUTILITY_P_H

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

#include "osxbluetooth_p.h"

#include <QtCore/qloggingcategory.h>
#include <QtCore/qscopedpointer.h>
#include <QtCore/qglobal.h>

#include <Foundation/Foundation.h>

QT_BEGIN_NAMESPACE

class QLowEnergyCharacteristicData;
class QBluetoothAddress;
class QBluetoothUuid;

namespace OSXBluetooth {

struct NSObjectDeleter {
    static void cleanup(NSObject *obj)
    {
        [obj release];
    }
};

template<class T>
class ObjCScopedPointer : public QScopedPointer<NSObject, NSObjectDeleter>
{
public:
    // TODO: remove default argument, add 'retain' parameter,
    // add a default ctor??? This will make the semantics more
    // transparent + will simplify the future transition to ARC
    // (if it will ever happen).
    explicit ObjCScopedPointer(T *ptr = nullptr) : QScopedPointer(ptr){}
    operator T*() const
    {
        return data();
    }

    T *data()const
    {
        return static_cast<T *>(QScopedPointer::data());
    }

    T *take()
    {
        return static_cast<T *>(QScopedPointer::take());
    }
};

#define QT_BT_MAC_AUTORELEASEPOOL const QMacAutoReleasePool pool;

template<class T>
class ObjCStrongReference {
public:
    ObjCStrongReference()
        : m_ptr(nil)
    {
    }
    ObjCStrongReference(T *obj, bool retain)
    {
        if (retain)
            m_ptr = [obj retain];
        else
            m_ptr = obj; // For example, created with initWithXXXX.
    }
    ObjCStrongReference(const ObjCStrongReference &rhs)
    {
        m_ptr = [rhs.m_ptr retain];
    }
    ObjCStrongReference &operator = (const ObjCStrongReference &rhs)
    {
        // "Old-style" implementation:
        if (this != &rhs && m_ptr != rhs.m_ptr) {
            [m_ptr release];
            m_ptr = [rhs.m_ptr retain];
        }

        return *this;
    }

#ifdef Q_COMPILER_RVALUE_REFS
    ObjCStrongReference(ObjCStrongReference &&xval)
    {
        m_ptr = xval.m_ptr;
        xval.m_ptr = nil;
    }

    ObjCStrongReference &operator = (ObjCStrongReference &&xval)
    {
        m_ptr = xval.m_ptr;
        xval.m_ptr = nil;
        return *this;
    }
#endif

    ~ObjCStrongReference()
    {
        [m_ptr release];
    }

    void reset(T *newVal)
    {
        if (m_ptr != newVal) {
            [m_ptr release];
            m_ptr = [newVal retain];
        }
    }

    void resetWithoutRetain(T *newVal)
    {
        if (m_ptr != newVal) {
            [m_ptr release];
            m_ptr = newVal;
        }
    }

    operator T *() const
    {
        return m_ptr;
    }

    T *data() const
    {
        return m_ptr;
    }

    T *take()
    {
        T * p = m_ptr;
        m_ptr = nil;
        return p;
    }
private:
    T *m_ptr;
};

// The type 'T' is some XXXRef from CoreFoundation and co.
// In principle, we can do a trick removing a pointer from a type
// when template is instantiated, but it's quite a lot of ugly pp-tokens
// like OSXBluetooth::CFStrongReference<OSXBluetooth::remove_pointer<CFUUIDRref> > strongReference;
// so instead we use 'T' everywhere, not 'T *' as can expected
// from a smart pointer.
template<class T>
class CFStrongReference {
public:
    CFStrongReference()
        : m_ptr(nullptr)
    {
    }

    CFStrongReference(T obj, bool retain)
        : m_ptr(obj)
    {
        if (m_ptr && retain)
            CFRetain(m_ptr);
    }

    CFStrongReference(const CFStrongReference &rhs)
    {
        if ((m_ptr = rhs.m_ptr))
            CFRetain(m_ptr);
    }

    CFStrongReference &operator = (const CFStrongReference &rhs)
    {
        // "Old-style" implementation:
        if (this != &rhs && m_ptr != rhs.m_ptr) {
            if (m_ptr)
                CFRelease(m_ptr);
            if ((m_ptr = rhs.m_ptr))
                CFRetain(m_ptr);
        }

        return *this;
    }

#ifdef Q_COMPILER_RVALUE_REFS
    CFStrongReference(CFStrongReference &&xval)
    {
        m_ptr = xval.m_ptr;
        xval.m_ptr = nullptr;
    }

    CFStrongReference &operator = (CFStrongReference &&xval)
    {
        m_ptr = xval.m_ptr;
        xval.m_ptr = nullptr;
        return *this;
    }
#endif

    ~CFStrongReference()
    {
        if (m_ptr)
            CFRelease(m_ptr);
    }

    void reset(T newVal)
    {
        if (m_ptr != newVal) {
            if (m_ptr)
                CFRelease(m_ptr);
            if ((m_ptr = newVal))
                CFRetain(m_ptr);
        }
    }

    operator T() const
    {
        return m_ptr;
    }

    T data() const
    {
        return m_ptr;
    }

    T take()
    {
        T p = m_ptr;
        m_ptr = nullptr;
        return p;
    }
private:
    T m_ptr;
};

QString qt_address(NSString *address);

#ifndef QT_IOS_BLUETOOTH

QBluetoothAddress qt_address(const BluetoothDeviceAddress *address);
BluetoothDeviceAddress iobluetooth_address(const QBluetoothAddress &address);

ObjCStrongReference<IOBluetoothSDPUUID> iobluetooth_uuid(const QBluetoothUuid &uuid);
QBluetoothUuid qt_uuid(IOBluetoothSDPUUID *uuid);
QString qt_error_string(IOReturn errorCode);
void qt_test_iobluetooth_runloop();

#endif // !QT_IOS_BLUETOOTH

QBluetoothUuid qt_uuid(CBUUID *uuid);
CFStrongReference<CFUUIDRef> cf_uuid(const QBluetoothUuid &qtUuid);
ObjCStrongReference<CBUUID> cb_uuid(const QBluetoothUuid &qtUuid);
bool equal_uuids(const QBluetoothUuid &qtUuid, CBUUID *cbUuid);
bool equal_uuids(CBUUID *cbUuid, const QBluetoothUuid &qtUuid);
QByteArray qt_bytearray(NSData *data);
QByteArray qt_bytearray(NSObject *data);

ObjCStrongReference<NSData> data_from_bytearray(const QByteArray &qtData);
ObjCStrongReference<NSMutableData> mutable_data_from_bytearray(const QByteArray &qtData);

dispatch_queue_t qt_LE_queue();

extern const int defaultLEScanTimeoutMS;
extern const int maxValueLength;

// Add more keys if needed, for now this one is enough:
extern NSString *const bluetoothUsageKey;

bool qt_appNeedsBluetoothUsageDescription();
bool qt_appPlistContainsDescription(NSString *key);

} // namespace OSXBluetooth

// Logging category for both OS X and iOS.
Q_DECLARE_LOGGING_CATEGORY(QT_BT_OSX)

QT_END_NAMESPACE

#if QT_MACOS_PLATFORM_SDK_EQUAL_OR_ABOVE(101300) && QT_MACOS_DEPLOYMENT_TARGET_BELOW(101300)

 // In the macOS 10.13 SDK, the identifier property was moved from the CBPeripheral
 // and CBCentral classes to a new CBPeer base class. Because CBPeer is only available
 // on macOS 10.13 and above, the same is true for -[CBPeer identifier]. However,
 // since we know that the derived classes have always had this property,
 // we'll explicitly mark its availability here. This will not adversely affect
 // using the identifier through the CBPeer base class, which will still require macOS 10.13.

@interface CBPeripheral (UnguardedWorkaround)
@property (readonly, nonatomic) NSUUID *identifier NS_AVAILABLE(10_7, 5_0);
@end

@interface CBCentral (UnguardedWorkaround)
@property (readonly, nonatomic) NSUUID *identifier NS_AVAILABLE(10_7, 5_0);
@end

#endif

#endif
