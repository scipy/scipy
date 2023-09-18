/****************************************************************************
**
** Copyright (C) 2019 Denis Shienkov <denis.shienkov@gmail.com>
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtSerialPort module of the Qt Toolkit.
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

#ifndef QTNTDLL_P_H
#define QTNTDLL_P_H

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
#include <QtCore/qstring.h>
#include <QtCore/qdebug.h>

#include <qt_windows.h>

// Internal control codes.

#ifndef CTL_CODE
#  define CTL_CODE(DeviceType, Function, Method, Access) ( \
    ((DeviceType) << 16) | ((Access) << 14) | ((Function) << 2) | (Method) \
    )
#endif

#ifndef FILE_DEVICE_SERIAL_PORT
#  define FILE_DEVICE_SERIAL_PORT  27
#endif

#ifndef METHOD_BUFFERED
#  define METHOD_BUFFERED  0
#endif

#ifndef FILE_ANY_ACCESS
#  define FILE_ANY_ACCESS  0x00000000
#endif

#ifndef IOCTL_SERIAL_GET_DTRRTS
#  define IOCTL_SERIAL_GET_DTRRTS \
    CTL_CODE(FILE_DEVICE_SERIAL_PORT, 30, METHOD_BUFFERED, FILE_ANY_ACCESS)
#endif

#ifndef SERIAL_DTR_STATE
#  define SERIAL_DTR_STATE  0x00000001
#endif

#ifndef SERIAL_RTS_STATE
#  define SERIAL_RTS_STATE  0x00000002
#endif

#ifndef IOCTL_SERIAL_WAIT_ON_MASK
#  define IOCTL_SERIAL_WAIT_ON_MASK \
    CTL_CODE(FILE_DEVICE_SERIAL_PORT, 18, METHOD_BUFFERED, FILE_ANY_ACCESS)
#endif

// Internal NT-based data types.

#ifndef NT_SUCCESS
#define NT_SUCCESS(status) (((NTSTATUS)(status)) >= 0)
#endif

typedef struct _IO_STATUS_BLOCK {
    union {
        NTSTATUS Status;
        PVOID Pointer;
    } DUMMYUNIONNAME;

    ULONG_PTR Information;
} IO_STATUS_BLOCK, *PIO_STATUS_BLOCK;

typedef VOID (WINAPI *PIO_APC_ROUTINE) (
        PVOID ApcContext,
        PIO_STATUS_BLOCK IoStatusBlock,
        ULONG Reserved
        );

// Resolving macros.

#define GENERATE_SYMBOL_VARIABLE(returnType, symbolName, ...) \
    typedef returnType (WINAPI *fp_##symbolName)(__VA_ARGS__); \
    static fp_##symbolName symbolName;

#define RESOLVE_SYMBOL(symbolName) \
    symbolName = reinterpret_cast<fp_##symbolName>(resolveNtdllSymbol(ntLibrary, #symbolName)); \
    if (!symbolName) \
        return false;

GENERATE_SYMBOL_VARIABLE(ULONG, RtlNtStatusToDosError, NTSTATUS)
GENERATE_SYMBOL_VARIABLE(NTSTATUS, NtDeviceIoControlFile, HANDLE, HANDLE, PIO_APC_ROUTINE, PVOID, PIO_STATUS_BLOCK, ULONG, PVOID, ULONG, PVOID, ULONG)

inline QFunctionPointer resolveNtdllSymbol(QLibrary *ntLibrary, const char *symbolName)
{
    QFunctionPointer symbolFunctionPointer = ntLibrary->resolve(symbolName);
    if (!symbolFunctionPointer)
        qWarning("Failed to resolve the symbol: %s", symbolName);

    return symbolFunctionPointer;
}

inline bool resolveNtdllSymbols(QLibrary *ntLibrary)
{
    if (!ntLibrary->isLoaded()) {
        ntLibrary->setFileName(QStringLiteral("ntdll"));
        if (!ntLibrary->load()) {
            qWarning("Failed to load the library: %s", qPrintable(ntLibrary->fileName()));
            return false;
        }
    }

    RESOLVE_SYMBOL(RtlNtStatusToDosError)
    RESOLVE_SYMBOL(NtDeviceIoControlFile)

    return true;
}

#endif // QTNTDLL_P_H
