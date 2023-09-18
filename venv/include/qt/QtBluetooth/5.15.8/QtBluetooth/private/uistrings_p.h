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

#ifndef TRANSLATIONS_H
#define TRANSLATIONS_H

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

#include <QtCore/QCoreApplication>
#include <QtCore/QString>

QT_BEGIN_NAMESPACE

// QBluetoothDeviceDiscoveryAgent related strings
extern const char DEV_DISCOVERY[];
extern const char DD_POWERED_OFF[];
extern const char DD_INVALID_ADAPTER[];
extern const char DD_IO[];
extern const char DD_NOTSUPPORTED[];
extern const char DD_UNKNOWN_ERROR[];
extern const char DD_NOT_STARTED[];
extern const char DD_NOT_STARTED_LE[];
extern const char DD_NOT_STOPPED[];

// QBluetoothServiceDiscoveryAgent related strings
extern const char SERVICE_DISCOVERY[];
extern const char SD_LOCAL_DEV_OFF[];
extern const char SD_MINIMAL_FAILED[];
extern const char SD_INVALID_ADDRESS[];

// QBluetoothSocket related strings
extern const char SOCKET[];
extern const char SOC_NETWORK_ERROR[];
extern const char SOC_NOWRITE[];
extern const char SOC_CONNECT_IN_PROGRESS[];
extern const char SOC_SERVICE_NOT_FOUND[];
extern const char SOC_INVAL_DATASIZE[];
extern const char SOC_NOREAD[];

// QBluetoothTransferReply related strings
extern const char TRANSFER_REPLY[];
extern const char TR_INVAL_TARGET[];
extern const char TR_SESSION_NO_START[];
extern const char TR_CONNECT_FAILED[];
extern const char TR_FILE_NOT_EXIST[];
extern const char TR_NOT_READ_IODEVICE[];
extern const char TR_SESSION_FAILED[];
extern const char TR_INVALID_DEVICE[];
extern const char TR_OP_CANCEL[];
extern const char TR_IN_PROGRESS[];
extern const char TR_SERVICE_NO_FOUND[];

// QLowEnergyController related strings
extern const char LE_CONTROLLER[];
extern const char LEC_RDEV_NO_FOUND[];
extern const char LEC_NO_LOCAL_DEV[];
extern const char LEC_IO_ERROR[];
extern const char LEC_UNKNOWN_ERROR[];

QT_END_NAMESPACE

#endif // TRANSLATIONS_H

