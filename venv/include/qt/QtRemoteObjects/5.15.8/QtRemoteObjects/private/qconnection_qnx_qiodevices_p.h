/****************************************************************************
**
** Copyright (C) 2017-2016 Ford Motor Company
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtRemoteObjects module of the Qt Toolkit.
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

#ifndef QQNXNATIVEIO_P_H
#define QQNXNATIVEIO_P_H

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

#include "qconnection_qnx_qiodevices.h"
#include "qconnection_qnx_global_p.h"

#include <QtCore/qreadwritelock.h>
#include <QtCore/qscopedpointer.h>

#include "private/qiodevice_p.h"
#include "private/qringbuffer_p.h"

QT_BEGIN_NAMESPACE

class QQnxNativeIoPrivate : public QIODevicePrivate
{
    Q_DECLARE_PUBLIC(QQnxNativeIo)

    mutable QReadWriteLock ibLock;
    mutable QReadWriteLock obLock;

public:
    QQnxNativeIoPrivate();
    ~QQnxNativeIoPrivate();
    void thread_func();
    bool establishConnection();
    void teardownConnection();
    void stopThread();
    QString serverName;
    int serverId, channelId, connectionId;
    sigevent tx_pulse;
    QAbstractSocket::SocketState state;
    QScopedPointer<QRingBuffer> obuffer;
    MsgType msgType;
    iov_t tx_iov[3], rx_iov[2];
    Thread<QQnxNativeIoPrivate> thread;
};

class QIOQnxSourcePrivate : public QIODevicePrivate
{
    Q_DECLARE_PUBLIC(QIOQnxSource)

    friend class QQnxNativeServerPrivate;

    mutable QReadWriteLock ibLock;
    mutable QReadWriteLock obLock;

public:
    QIOQnxSourcePrivate(int _rcvid);
    int rcvid;
    QAbstractSocket::SocketState state;
    MsgType msgType;
    sigevent m_event;
    QRingBuffer obuffer;
    QAtomicInt m_serverClosing;
};

QT_END_NAMESPACE

#endif // QQNXNATIVEIO_P_H

