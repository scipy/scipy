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

#ifndef QDBUSCONNECTION_P_H
#define QDBUSCONNECTION_P_H

#include <QtDBus/private/qtdbusglobal_p.h>
#include <qdbuserror.h>
#include <qdbusconnection.h>

#include <QtCore/qatomic.h>
#include <QtCore/qhash.h>
#include <QtCore/qobject.h>
#include <QtCore/qpointer.h>
#include <QtCore/qreadwritelock.h>
#include <QtCore/qstringlist.h>
#include <QtCore/qvarlengtharray.h>
#include <QtCore/qvector.h>

#include "qdbus_symbols_p.h"

#include <qdbusmessage.h>
#include <qdbusservicewatcher.h>    // for the WatchMode enum

#ifndef QT_NO_DBUS

QT_BEGIN_NAMESPACE

class QDBusMessage;
class QSocketNotifier;
class QTimerEvent;
class QDBusObjectPrivate;
class QDBusCallDeliveryEvent;
class QDBusActivateObjectEvent;
class QMetaMethod;
class QDBusInterfacePrivate;
struct QDBusMetaObject;
class QDBusAbstractInterface;
class QDBusConnectionInterface;
class QDBusPendingCallPrivate;
class QDBusServer;

#ifndef QT_BOOTSTRAPPED

class QDBusErrorInternal
{
    mutable DBusError error;
    Q_DISABLE_COPY_MOVE(QDBusErrorInternal)
public:
    inline QDBusErrorInternal() { q_dbus_error_init(&error); }
    inline ~QDBusErrorInternal() { q_dbus_error_free(&error); }
    inline bool operator !() const { return !q_dbus_error_is_set(&error); }
    inline operator DBusError *() { q_dbus_error_free(&error); return &error; }
    inline operator QDBusError() const { QDBusError err(&error); q_dbus_error_free(&error); return err; }
};

// QDBusConnectionPrivate holds the DBusConnection and
// can have many QDBusConnection objects referring to it

class Q_AUTOTEST_EXPORT QDBusConnectionPrivate: public QObject
{
    Q_OBJECT
public:
    // structs and enums
    enum ConnectionMode { InvalidMode, ServerMode, ClientMode, PeerMode }; // LocalMode

    struct Watcher
    {
        Watcher(): watch(nullptr), read(nullptr), write(nullptr) {}
        DBusWatch *watch;
        QSocketNotifier *read;
        QSocketNotifier *write;
    };

    struct ArgMatchRules {
        QStringList args;
        QString arg0namespace;
        bool operator==(const ArgMatchRules &other) const {
            return args == other.args &&
                   arg0namespace == other.arg0namespace;
        }
    };

    struct SignalHook
    {
        inline SignalHook() : obj(nullptr), midx(-1) { }
        QString service, path, signature;
        QObject* obj;
        int midx;
        QVector<int> params;
        ArgMatchRules argumentMatch;
        QByteArray matchRule;
    };

    enum TreeNodeType {
        Object = 0x0,
        VirtualObject = 0x01000000
    };

    struct ObjectTreeNode
    {
        typedef QVector<ObjectTreeNode> DataList;

        inline ObjectTreeNode() : obj(nullptr), flags(0) { }
        inline ObjectTreeNode(const QString &n) // intentionally implicit
            : name(n), obj(nullptr), flags(0) { }
        inline bool operator<(const QString &other) const
            { return name < other; }
        inline bool operator<(const QStringRef &other) const
            { return QStringRef(&name) < other; }
        inline bool isActive() const
        { return obj || !children.isEmpty(); }

        QString name;
        QString interfaceName;
        union {
            QObject *obj;
            QDBusVirtualObject *treeNode;
        };
        int flags;

        DataList children;
    };

public:
    // typedefs
    typedef QMultiHash<qintptr, Watcher> WatcherHash;
    typedef QHash<int, DBusTimeout *> TimeoutHash;
    typedef QVector<QDBusMessage> PendingMessageList;

    typedef QMultiHash<QString, SignalHook> SignalHookHash;
    typedef QHash<QString, QDBusMetaObject* > MetaObjectHash;
    typedef QHash<QByteArray, int> MatchRefCountHash;
    typedef QVector<QDBusPendingCallPrivate*> PendingCallList;

    struct WatchedServiceData {
        WatchedServiceData() : refcount(0) {}
        WatchedServiceData(const QString &owner, int refcount = 0)
            : owner(owner), refcount(refcount)
        {}
        QString owner;
        int refcount;
    };
    typedef QHash<QString, WatchedServiceData> WatchedServicesHash;

public:
    // public methods are entry points from other objects
    explicit QDBusConnectionPrivate(QObject *parent = nullptr);
    ~QDBusConnectionPrivate();

    void createBusService();
    void setPeer(DBusConnection *connection, const QDBusErrorInternal &error);
    void setConnection(DBusConnection *connection, const QDBusErrorInternal &error);
    void setServer(QDBusServer *object, DBusServer *server, const QDBusErrorInternal &error);
    void closeConnection();

    QString getNameOwner(const QString &service);

    bool shouldWatchService(const QString &service);
    void watchService(const QString &service, QDBusServiceWatcher::WatchMode mode,
                      QObject *obj, const char *member);
    void unwatchService(const QString &service, QDBusServiceWatcher::WatchMode mode,
                        QObject *obj, const char *member);

    bool send(const QDBusMessage &message);
    QDBusMessage sendWithReply(const QDBusMessage &message, int mode, int timeout = -1);
    QDBusMessage sendWithReplyLocal(const QDBusMessage &message);
    QDBusPendingCallPrivate *sendWithReplyAsync(const QDBusMessage &message, QObject *receiver,
                                                const char *returnMethod, const char *errorMethod,int timeout = -1);

    bool connectSignal(const QString &service, const QString &path, const QString& interface,
                       const QString &name, const QStringList &argumentMatch, const QString &signature,
                       QObject *receiver, const char *slot);
    bool disconnectSignal(const QString &service, const QString &path, const QString& interface,
                          const QString &name, const QStringList &argumentMatch, const QString &signature,
                          QObject *receiver, const char *slot);
    bool connectSignal(const QString &service, const QString &path, const QString& interface,
                       const QString &name, const ArgMatchRules &argumentMatch, const QString &signature,
                       QObject *receiver, const char *slot);
    bool disconnectSignal(const QString &service, const QString &path, const QString& interface,
                          const QString &name, const ArgMatchRules &argumentMatch, const QString &signature,
                          QObject *receiver, const char *slot);
    void registerObject(const ObjectTreeNode *node);
    void unregisterObject(const QString &path, QDBusConnection::UnregisterMode mode);
    void connectRelay(const QString &service,
                      const QString &path, const QString &interface,
                      QDBusAbstractInterface *receiver, const QMetaMethod &signal);
    void disconnectRelay(const QString &service,
                         const QString &path, const QString &interface,
                         QDBusAbstractInterface *receiver, const QMetaMethod &signal);
    void registerService(const QString &serviceName);
    void unregisterService(const QString &serviceName);

    bool handleMessage(const QDBusMessage &msg);

    QDBusMetaObject *findMetaObject(const QString &service, const QString &path,
                                    const QString &interface, QDBusError &error);

    void postEventToThread(int action, QObject *target, QEvent *event);

private:
    void checkThread();
    bool handleError(const QDBusErrorInternal &error);

    void handleSignal(const QString &key, const QDBusMessage &msg);
    void handleSignal(const QDBusMessage &msg);
    void handleObjectCall(const QDBusMessage &message);

    void activateSignal(const SignalHook& hook, const QDBusMessage &msg);
    void activateObject(ObjectTreeNode &node, const QDBusMessage &msg, int pathStartPos);
    bool activateInternalFilters(const ObjectTreeNode &node, const QDBusMessage &msg);
    bool activateCall(QObject *object, int flags, const QDBusMessage &msg);

    void sendInternal(QDBusPendingCallPrivate *pcall, void *msg, int timeout);
    void sendError(const QDBusMessage &msg, QDBusError::ErrorType code);
    void deliverCall(QObject *object, int flags, const QDBusMessage &msg,
                     const QVector<int> &metaTypes, int slotIdx);

    SignalHookHash::Iterator removeSignalHookNoLock(SignalHookHash::Iterator it);
    void collectAllObjects(ObjectTreeNode &node, QSet<QObject *> &set);

    bool isServiceRegisteredByThread(const QString &serviceName);

    QString getNameOwnerNoCache(const QString &service);

    void watchForDBusDisconnection();

    void _q_newConnection(QDBusConnectionPrivate *newConnection);

    void handleAuthentication();

protected:
    void timerEvent(QTimerEvent *e) override;

public slots:
    // public slots
    void setDispatchEnabled(bool enable);
    void doDispatch();
    void socketRead(qintptr);
    void socketWrite(qintptr);
    void objectDestroyed(QObject *o);
    void relaySignal(QObject *obj, const QMetaObject *, int signalId, const QVariantList &args);
    bool addSignalHook(const QString &key, const SignalHook &hook);
    bool removeSignalHook(const QString &key, const SignalHook &hook);

private slots:
    void serviceOwnerChangedNoLock(const QString &name, const QString &oldOwner, const QString &newOwner);
    void registerServiceNoLock(const QString &serviceName);
    void unregisterServiceNoLock(const QString &serviceName);
    void handleDBusDisconnection();

signals:
    void dispatchStatusChanged();
    void spyHooksFinished(const QDBusMessage &msg);
    void messageNeedsSending(QDBusPendingCallPrivate *pcall, void *msg, int timeout = -1);
    bool signalNeedsConnecting(const QString &key, const QDBusConnectionPrivate::SignalHook &hook);
    bool signalNeedsDisconnecting(const QString &key, const QDBusConnectionPrivate::SignalHook &hook);
    void serviceOwnerChanged(const QString &name, const QString &oldOwner, const QString &newOwner);
    void callWithCallbackFailed(const QDBusError &error, const QDBusMessage &message);
    void newServerConnection(QDBusConnectionPrivate *newConnection);

public:
    QAtomicInt ref;
    QAtomicInt capabilities;
    QDBusConnection::ConnectionCapabilities connectionCapabilities() const
    {
        return (QDBusConnection::ConnectionCapabilities)capabilities.loadRelaxed();
    }
    QString name;               // this connection's name
    QString baseService;        // this connection's base service
    QStringList serverConnectionNames;

    ConnectionMode mode;
    union {
        QDBusConnectionInterface *busService;
        QDBusServer *serverObject;
    };

    union {
        DBusConnection *connection;
        DBusServer *server;
    };
    WatcherHash watchers;
    TimeoutHash timeouts;
    PendingMessageList pendingMessages;

    // the master lock protects our own internal state
    QReadWriteLock lock;
    QDBusError lastError;

    QStringList serviceNames;
    WatchedServicesHash watchedServices;
    SignalHookHash signalHooks;
    MatchRefCountHash matchRefCounts;
    ObjectTreeNode rootNode;
    MetaObjectHash cachedMetaObjects;
    PendingCallList pendingCalls;

    bool anonymousAuthenticationAllowed;
    bool dispatchEnabled;               // protected by the dispatch lock, not the main lock
    bool isAuthenticated;

public:
    // static methods
    static int findSlot(QObject *obj, const QByteArray &normalizedName, QVector<int> &params);
    static bool prepareHook(QDBusConnectionPrivate::SignalHook &hook, QString &key,
                            const QString &service,
                            const QString &path, const QString &interface, const QString &name,
                            const ArgMatchRules &argMatch,
                            QObject *receiver, const char *signal, int minMIdx,
                            bool buildSignature);
    static DBusHandlerResult messageFilter(DBusConnection *, DBusMessage *, void *);
    static QDBusCallDeliveryEvent *prepareReply(QDBusConnectionPrivate *target, QObject *object,
                                                int idx, const QVector<int> &metaTypes,
                                                const QDBusMessage &msg);
    static void processFinishedCall(QDBusPendingCallPrivate *call);

    static QDBusConnectionPrivate *d(const QDBusConnection& q) { return q.d; }
    static QDBusConnection q(QDBusConnectionPrivate *connection) { return QDBusConnection(connection); }

    friend class QDBusActivateObjectEvent;
    friend class QDBusCallDeliveryEvent;
    friend class QDBusServer;
};

// in qdbusmisc.cpp
extern int qDBusParametersForMethod(const QMetaMethod &mm, QVector<int> &metaTypes, QString &errorMsg);
#endif // QT_BOOTSTRAPPED
extern Q_DBUS_EXPORT int qDBusParametersForMethod(const QList<QByteArray> &parameters, QVector<int>& metaTypes, QString &errorMsg);
extern Q_DBUS_EXPORT bool qDBusCheckAsyncTag(const char *tag);
#ifndef QT_BOOTSTRAPPED
extern bool qDBusInterfaceInObject(QObject *obj, const QString &interface_name);
extern QString qDBusInterfaceFromMetaObject(const QMetaObject *mo);

// in qdbusinternalfilters.cpp
extern QString qDBusIntrospectObject(const QDBusConnectionPrivate::ObjectTreeNode &node, const QString &path);
extern QDBusMessage qDBusPropertyGet(const QDBusConnectionPrivate::ObjectTreeNode &node,
                                     const QDBusMessage &msg);
extern QDBusMessage qDBusPropertySet(const QDBusConnectionPrivate::ObjectTreeNode &node,
                                     const QDBusMessage &msg);
extern QDBusMessage qDBusPropertyGetAll(const QDBusConnectionPrivate::ObjectTreeNode &node,
                                        const QDBusMessage &msg);

// can be replaced with a lambda in Qt 5.7
class QDBusConnectionDispatchEnabler : public QObject
{
    Q_OBJECT
    QDBusConnectionPrivate *con;
public:
    QDBusConnectionDispatchEnabler(QDBusConnectionPrivate *con) : con(con) {}

public slots:
    void execute()
    {
        // This call cannot race with something disabling dispatch only because dispatch is
        // never re-disabled from Qt code on an in-use connection once it has been enabled.
        QMetaObject::invokeMethod(con, "setDispatchEnabled", Qt::QueuedConnection, Q_ARG(bool, true));
        if (!con->ref.deref())
            con->deleteLater();
        deleteLater();
    }
};

#endif // QT_BOOTSTRAPPED

QT_END_NAMESPACE

#endif // QT_NO_DBUS
#endif
