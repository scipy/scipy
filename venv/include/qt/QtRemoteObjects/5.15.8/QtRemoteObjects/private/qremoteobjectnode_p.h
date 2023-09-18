/****************************************************************************
**
** Copyright (C) 2017 Ford Motor Company
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

#ifndef QREMOTEOBJECTNODE_P_H
#define QREMOTEOBJECTNODE_P_H

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

#include <QtCore/private/qobject_p.h>
#include "qremoteobjectsourceio_p.h"
#include "qremoteobjectreplica.h"
#include "qremoteobjectnode.h"

#include <QtCore/qbasictimer.h>
#include <QtCore/qmutex.h>

QT_BEGIN_NAMESPACE

#define qRODebug(x) qCDebug(QT_REMOTEOBJECT) << qPrintable(QtPrivate::deref_for_methodcall(x).objectName())
#define qROWarning(x) qCWarning(QT_REMOTEOBJECT) << qPrintable(QtPrivate::deref_for_methodcall(x).objectName())
#define qROCritical(x) qCCritical(QT_REMOTEOBJECT) << qPrintable(QtPrivate::deref_for_methodcall(x).objectName())
#define qROFatal(x) qCFatal(QT_REMOTEOBJECT) << qPrintable(QtPrivate::deref_for_methodcall(x).objectName())
#define qROPrivDebug() qCDebug(QT_REMOTEOBJECT) << qPrintable(q_ptr->objectName())
#define qROPrivWarning() qCWarning(QT_REMOTEOBJECT) << qPrintable(q_ptr->objectName())
#define qROPrivCritical() qCCritical(QT_REMOTEOBJECT) << qPrintable(q_ptr->objectName())
#define qROPrivFatal() qCFatal(QT_REMOTEOBJECT) << qPrintable(q_ptr->objectName())

class QRemoteObjectRegistry;
class QRegistrySource;
class QConnectedReplicaImplementation;

class QRemoteObjectAbstractPersistedStorePrivate : public QObjectPrivate
{
public:
    QRemoteObjectAbstractPersistedStorePrivate();
    virtual ~QRemoteObjectAbstractPersistedStorePrivate();

    Q_DECLARE_PUBLIC(QRemoteObjectAbstractPersistedStore)
};

class QRemoteObjectMetaObjectManager
{
public:
    QRemoteObjectMetaObjectManager() {}
    ~QRemoteObjectMetaObjectManager();

    const QMetaObject *metaObjectForType(const QString &type);
    QMetaObject *addDynamicType(IoDeviceBase* connection, QDataStream &in);
    void addFromMetaObject(const QMetaObject *);

private:
    QHash<QString, QMetaObject*> dynamicTypes;
    QHash<QString, const QMetaObject*> staticTypes;
};

struct ProxyReplicaInfo;
class ProxyInfo : public QObject
{
    Q_OBJECT
public:
    ProxyInfo(QRemoteObjectNode *node, QRemoteObjectHostBase *parent, QRemoteObjectHostBase::RemoteObjectNameFilter filter);
    ~ProxyInfo() override;
    enum class ProxyDirection { Forward, Reverse };

    bool setReverseProxy(QRemoteObjectHostBase::RemoteObjectNameFilter filter);
    void proxyObject(const QRemoteObjectSourceLocation &entry, ProxyDirection direction = ProxyDirection::Forward);
    void unproxyObject(const QRemoteObjectSourceLocation &entry);

    QRemoteObjectNode *proxyNode;
    QRemoteObjectHostBase *parentNode;
    QRemoteObjectHostBase::RemoteObjectNameFilter proxyFilter;
    QRemoteObjectHostBase::RemoteObjectNameFilter reverseFilter;
    QHash<QString, ProxyReplicaInfo*> proxiedReplicas;

private:
    void disableAndDeleteObject(ProxyReplicaInfo* info);
};

struct ProxyReplicaInfo
{
    // We need QObject, so we can hold Dynamic Replicas and QAIM Adapters
    QObject* replica;
    ProxyInfo::ProxyDirection direction;
    ~ProxyReplicaInfo() { delete replica; }
};

class QRemoteObjectNodePrivate : public QObjectPrivate
{
public:
    QRemoteObjectNodePrivate();
    ~QRemoteObjectNodePrivate() override;

    virtual QRemoteObjectSourceLocations remoteObjectAddresses() const;

    void setReplicaImplementation(const QMetaObject *, QRemoteObjectReplica *, const QString &);

    void setLastError(QRemoteObjectNode::ErrorCode errorCode);

    void connectReplica(QObject *object, QRemoteObjectReplica *instance);
    void openConnectionIfNeeded(const QString &name);

    bool initConnection(const QUrl &address);
    bool hasInstance(const QString &name);
    void setRegistry(QRemoteObjectRegistry *);
    QVariant handlePointerToQObjectProperty(QConnectedReplicaImplementation *rep, int index, const QVariant &property);
    void handlePointerToQObjectProperties(QConnectedReplicaImplementation *rep, QVariantList &properties);

    void onClientRead(QObject *obj);
    void onRemoteObjectSourceAdded(const QRemoteObjectSourceLocation &entry);
    void onRemoteObjectSourceRemoved(const QRemoteObjectSourceLocation &entry);
    void onRegistryInitialized();
    void onShouldReconnect(ClientIoDevice *ioDevice);

    virtual QReplicaImplementationInterface *handleNewAcquire(const QMetaObject *meta, QRemoteObjectReplica *instance, const QString &name);
    void handleReplicaConnection(const QString &name);
    void handleReplicaConnection(const QByteArray &sourceSignature, QConnectedReplicaImplementation *rep, IoDeviceBase *connection);
    void initialize();
private:
    bool checkSignatures(const QByteArray &a, const QByteArray &b);

public:
    struct SourceInfo
    {
        IoDeviceBase* device;
        QString typeName;
        QByteArray objectSignature;
    };

    QMutex mutex;
    QUrl registryAddress;
    QHash<QString, QWeakPointer<QReplicaImplementationInterface> > replicas;
    QMap<QString, SourceInfo> connectedSources;
    QMap<QString, QRemoteObjectNode::RemoteObjectSchemaHandler> schemaHandlers;
    QSet<ClientIoDevice*> pendingReconnect;
    QSet<QUrl> requestedUrls;
    QRemoteObjectRegistry *registry;
    int retryInterval;
    QBasicTimer reconnectTimer;
    QRemoteObjectNode::ErrorCode lastError;
    QString rxName;
    QRemoteObjectPackets::ObjectInfoList rxObjects;
    QVariantList rxArgs;
    QVariant rxValue;
    QRemoteObjectAbstractPersistedStore *persistedStore;
    bool m_handshakeReceived = false;
    int m_heartbeatInterval = 0;
    QRemoteObjectMetaObjectManager dynamicTypeManager;
    Q_DECLARE_PUBLIC(QRemoteObjectNode)
};

class QRemoteObjectHostBasePrivate : public QRemoteObjectNodePrivate
{
public:
    QRemoteObjectHostBasePrivate();
    ~QRemoteObjectHostBasePrivate() override;
    QReplicaImplementationInterface *handleNewAcquire(const QMetaObject *meta, QRemoteObjectReplica *instance, const QString &name) override;

public:
    QRemoteObjectSourceIo *remoteObjectIo;
    ProxyInfo *proxyInfo = nullptr;
    Q_DECLARE_PUBLIC(QRemoteObjectHostBase);
};

class QRemoteObjectHostPrivate : public QRemoteObjectHostBasePrivate
{
public:
    QRemoteObjectHostPrivate();
    ~QRemoteObjectHostPrivate() override;
    Q_DECLARE_PUBLIC(QRemoteObjectHost);
};

class QRemoteObjectRegistryHostPrivate : public QRemoteObjectHostBasePrivate
{
public:
    QRemoteObjectRegistryHostPrivate();
    ~QRemoteObjectRegistryHostPrivate() override;
    QRemoteObjectSourceLocations remoteObjectAddresses() const override;
    QRegistrySource *registrySource;
    Q_DECLARE_PUBLIC(QRemoteObjectRegistryHost);
};

QT_END_NAMESPACE

#endif
