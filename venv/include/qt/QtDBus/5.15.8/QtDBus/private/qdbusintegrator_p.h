/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Copyright (C) 2016 Intel Corporation.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtDBus module of the Qt Toolkit.
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

//
//  W A R N I N G
//  -------------
//
// This file is not part of the public API.  This header file may
// change from version to version without notice, or even be
// removed.
//
// We mean it.
//
//

#ifndef QDBUSINTEGRATOR_P_H
#define QDBUSINTEGRATOR_P_H

#include <QtDBus/private/qtdbusglobal_p.h>
#include "qdbus_symbols_p.h"

#include "qcoreevent.h"
#include "qeventloop.h"
#include "qobject.h"
#include "private/qobject_p.h"
#include "qlist.h"
#include "qpointer.h"
#include "qsemaphore.h"

#include "qdbusconnection.h"
#include "qdbusmessage.h"
#include "qdbusconnection_p.h"

#ifndef QT_NO_DBUS

QT_BEGIN_NAMESPACE

class QDBusConnectionPrivate;
class QDBusMessage;

// Really private structs used by qdbusintegrator.cpp
// Things that aren't used by any other file

struct QDBusSlotCache
{
    struct Data
    {
        int flags;
        int slotIdx;
        QVector<int> metaTypes;

        void swap(Data &other) noexcept
        {
            qSwap(flags,     other.flags);
            qSwap(slotIdx,   other.slotIdx);
            qSwap(metaTypes, other.metaTypes);
        }
    };
    typedef QMultiHash<QString, Data> Hash;
    Hash hash;

    void swap(QDBusSlotCache &other) noexcept { qSwap(hash, other.hash); }
};
Q_DECLARE_SHARED(QDBusSlotCache::Data)
Q_DECLARE_SHARED(QDBusSlotCache)

class QDBusCallDeliveryEvent: public QAbstractMetaCallEvent
{
public:
    QDBusCallDeliveryEvent(const QDBusConnection &c, int id, QObject *sender,
                           const QDBusMessage &msg, const QVector<int> &types, int f = 0)
        : QAbstractMetaCallEvent(sender, -1), connection(c), message(msg), metaTypes(types), id(id), flags(f)
        { }

    void placeMetaCall(QObject *object) override
    {
        QDBusConnectionPrivate::d(connection)->deliverCall(object, flags, message, metaTypes, id);
    }

private:
    QDBusConnection connection; // just for refcounting
    QDBusMessage message;
    QVector<int> metaTypes;
    int id;
    int flags;
};

class QDBusActivateObjectEvent: public QAbstractMetaCallEvent
{
public:
    QDBusActivateObjectEvent(const QDBusConnection &c, QObject *sender,
                             const QDBusConnectionPrivate::ObjectTreeNode &n,
                             int p, const QDBusMessage &m, QSemaphore *s = nullptr)
        : QAbstractMetaCallEvent(sender, -1, s), connection(c), node(n),
          pathStartPos(p), message(m), handled(false)
        { }
    ~QDBusActivateObjectEvent() override;

    void placeMetaCall(QObject *) override;

private:
    QDBusConnection connection; // just for refcounting
    QDBusConnectionPrivate::ObjectTreeNode node;
    int pathStartPos;
    QDBusMessage message;
    bool handled;
};

class QDBusSpyCallEvent : public QAbstractMetaCallEvent
{
public:
    typedef void (*Hook)(const QDBusMessage&);
    QDBusSpyCallEvent(QDBusConnectionPrivate *cp, const QDBusConnection &c, const QDBusMessage &msg,
                      const Hook *hooks, int count)
        : QAbstractMetaCallEvent(cp, 0), conn(c), msg(msg), hooks(hooks), hookCount(count)
    {}
    ~QDBusSpyCallEvent() override;
    void placeMetaCall(QObject *) override;
    static inline void invokeSpyHooks(const QDBusMessage &msg, const Hook *hooks, int hookCount);

    QDBusConnection conn;   // keeps the refcount in QDBusConnectionPrivate up
    QDBusMessage msg;
    const Hook *hooks;
    int hookCount;
};

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QDBusSlotCache)

#endif // QT_NO_DBUS
#endif
