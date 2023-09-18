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

#ifndef QQNXNATIVESERVER_P_H
#define QQNXNATIVESERVER_P_H

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

#include "private/qobject_p.h"
#include "qconnection_qnx_server.h"
#include "qconnection_qnx_global_p.h"

#include <QtCore/qatomic.h>
#include <QtCore/qmutex.h>
#include <QtCore/qsharedpointer.h>

QT_BEGIN_NAMESPACE

class QQnxNativeServerPrivate : public QObjectPrivate
{
    Q_DECLARE_PUBLIC(QQnxNativeServer)

public:
    QQnxNativeServerPrivate();

    ~QQnxNativeServerPrivate();

    void thread_func();

    void cleanupIOSource(QIOQnxSource *conn);
    void teardownServer();
    void createSource(int rcvid, uint64_t uid, pid_t toPid);
#ifdef USE_HAM
    bool initializeHam();
    void configureHamDeath(int sentChannelId, pid_t toPid, uint64_t uid);
    void closeHamResources();
#endif

    bool listen(const QString &name);
    QString errorString;
    QAbstractSocket::SocketError error;
    QString serverName;
    name_attach_t *attachStruct;
    QHash<int, QSet<int> > connections;
    QHash<uint64_t, QSharedPointer<QIOQnxSource>> sources;
    QList<QSharedPointer<QIOQnxSource>> pending;
    QAtomicInt running;
    Thread<QQnxNativeServerPrivate> thread;
    mutable QMutex mutex;
    int terminateCoid;
#ifdef USE_HAM
    ham_entity_t *hamEntityHandle;
    ham_condition_t *hamConditionHandle;
    QHash<uint64_t, ham_action_t*> hamActions;
    bool hamAvailable = false;
    bool hamInitialized = false;
#endif
};

QT_END_NAMESPACE

#endif // QQNXNATIVESERVER_P_H

