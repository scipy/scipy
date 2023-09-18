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

#ifndef QREMOTEOBJECTNODE_H
#define QREMOTEOBJECTNODE_H

#include <QtCore/qsharedpointer.h>
#include <QtCore/qmetaobject.h>
#include <QtRemoteObjects/qtremoteobjectglobal.h>
#include <QtRemoteObjects/qremoteobjectregistry.h>
#include <QtRemoteObjects/qremoteobjectdynamicreplica.h>

#include <functional>

QT_BEGIN_NAMESPACE

class QRemoteObjectReplica;
class SourceApiMap;
class QAbstractItemModel;
class QAbstractItemModelReplica;
class QItemSelectionModel;
class QRemoteObjectAbstractPersistedStorePrivate;
class QRemoteObjectNodePrivate;
class QRemoteObjectHostBasePrivate;
class QRemoteObjectHostPrivate;
class QRemoteObjectRegistryHostPrivate;

class Q_REMOTEOBJECTS_EXPORT QRemoteObjectAbstractPersistedStore : public QObject
{
    Q_OBJECT

public:
    QRemoteObjectAbstractPersistedStore (QObject *parent = nullptr);
    virtual ~QRemoteObjectAbstractPersistedStore();

    virtual void saveProperties(const QString &repName, const QByteArray &repSig, const QVariantList &values) = 0;
    virtual QVariantList restoreProperties(const QString &repName, const QByteArray &repSig) = 0;

protected:
    QRemoteObjectAbstractPersistedStore(QRemoteObjectAbstractPersistedStorePrivate &, QObject *parent);
    Q_DECLARE_PRIVATE(QRemoteObjectAbstractPersistedStore)
};

class Q_REMOTEOBJECTS_EXPORT QRemoteObjectNode : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QUrl registryUrl READ registryUrl WRITE setRegistryUrl)
    Q_PROPERTY(QRemoteObjectAbstractPersistedStore* persistedStore READ persistedStore WRITE setPersistedStore)
    Q_PROPERTY(int heartbeatInterval READ heartbeatInterval WRITE setHeartbeatInterval NOTIFY heartbeatIntervalChanged)

public:
    enum ErrorCode{
        NoError,
        RegistryNotAcquired,
        RegistryAlreadyHosted,
        NodeIsNoServer,
        ServerAlreadyCreated,
        UnintendedRegistryHosting,
        OperationNotValidOnClientNode,
        SourceNotRegistered,
        MissingObjectName,
        HostUrlInvalid,
        ProtocolMismatch,
        ListenFailed
    };
    Q_ENUM(ErrorCode)

    QRemoteObjectNode(QObject *parent = nullptr);
    QRemoteObjectNode(const QUrl &registryAddress, QObject *parent = nullptr);
    ~QRemoteObjectNode() override;

    Q_INVOKABLE bool connectToNode(const QUrl &address);
    void addClientSideConnection(QIODevice *ioDevice);
    virtual void setName(const QString &name);
    template < class ObjectType >
    ObjectType *acquire(const QString &name = QString())
    {
        return new ObjectType(this, name);
    }

    template<typename T>
    QStringList instances() const
    {
        const QMetaObject *mobj = &T::staticMetaObject;
        const int index = mobj->indexOfClassInfo(QCLASSINFO_REMOTEOBJECT_TYPE);
        if (index == -1)
            return QStringList();

        const QString typeName = QString::fromLatin1(mobj->classInfo(index).value());
        return instances(typeName);
    }
    QStringList instances(const QString &typeName) const;

    QRemoteObjectDynamicReplica *acquireDynamic(const QString &name);
    QAbstractItemModelReplica *acquireModel(const QString &name, QtRemoteObjects::InitialAction action = QtRemoteObjects::FetchRootSize, const QVector<int> &rolesHint = {});
    QUrl registryUrl() const;
    virtual bool setRegistryUrl(const QUrl &registryAddress);
    bool waitForRegistry(int timeout = 30000);
    const QRemoteObjectRegistry *registry() const;

    QRemoteObjectAbstractPersistedStore *persistedStore() const;
    void setPersistedStore(QRemoteObjectAbstractPersistedStore *persistedStore);

    ErrorCode lastError() const;

    int heartbeatInterval() const;
    void setHeartbeatInterval(int interval);

    typedef std::function<void (QUrl)> RemoteObjectSchemaHandler;
    void registerExternalSchema(const QString &schema, RemoteObjectSchemaHandler handler);

Q_SIGNALS:
    void remoteObjectAdded(const QRemoteObjectSourceLocation &);
    void remoteObjectRemoved(const QRemoteObjectSourceLocation &);

    void error(QRemoteObjectNode::ErrorCode errorCode);
    void heartbeatIntervalChanged(int heartbeatInterval);

protected:
    QRemoteObjectNode(QRemoteObjectNodePrivate &, QObject *parent);

    void timerEvent(QTimerEvent*) override;

private:
    void initializeReplica(QRemoteObjectReplica *instance, const QString &name = QString());
    void persistProperties(const QString &repName, const QByteArray &repSig, const QVariantList &props);
    QVariantList retrieveProperties(const QString &repName, const QByteArray &repSig);

    Q_DECLARE_PRIVATE(QRemoteObjectNode)
    friend class QRemoteObjectReplica;
};

class Q_REMOTEOBJECTS_EXPORT QRemoteObjectHostBase : public QRemoteObjectNode
{
    Q_OBJECT
public:
    enum AllowedSchemas { BuiltInSchemasOnly, AllowExternalRegistration };
    Q_ENUM(AllowedSchemas)
    ~QRemoteObjectHostBase() override;
    void setName(const QString &name) override;

    template <template <typename> class ApiDefinition, typename ObjectType>
    bool enableRemoting(ObjectType *object)
    {
        ApiDefinition<ObjectType> *api = new ApiDefinition<ObjectType>(object);
        return enableRemoting(object, api);
    }
    Q_INVOKABLE bool enableRemoting(QObject *object, const QString &name = QString());
    bool enableRemoting(QAbstractItemModel *model, const QString &name, const QVector<int> roles, QItemSelectionModel *selectionModel = nullptr);
    Q_INVOKABLE bool disableRemoting(QObject *remoteObject);
    void addHostSideConnection(QIODevice *ioDevice);

    typedef std::function<bool(const QString &, const QString &)> RemoteObjectNameFilter;
    bool proxy(const QUrl &registryUrl, const QUrl &hostUrl={},
               RemoteObjectNameFilter filter=[](const QString &, const QString &) {return true; });
    // ### Qt 6: Fix -> This should only be part of the QRemoteObjectRegistryHost type, since the
    // reverse aspect requires the registry.
    bool reverseProxy(RemoteObjectNameFilter filter=[](const QString &, const QString &) {return true; });

protected:
    virtual QUrl hostUrl() const;
    virtual bool setHostUrl(const QUrl &hostAddress, AllowedSchemas allowedSchemas=BuiltInSchemasOnly);
    QRemoteObjectHostBase(QRemoteObjectHostBasePrivate &, QObject *);

private:
    bool enableRemoting(QObject *object, const SourceApiMap *, QObject *adapter = nullptr);
    Q_DECLARE_PRIVATE(QRemoteObjectHostBase)
};

class Q_REMOTEOBJECTS_EXPORT QRemoteObjectHost : public QRemoteObjectHostBase
{
    Q_OBJECT
    Q_PROPERTY(QUrl hostUrl READ hostUrl WRITE setHostUrl NOTIFY hostUrlChanged)

public:
    QRemoteObjectHost(QObject *parent = nullptr);
    QRemoteObjectHost(const QUrl &address, const QUrl &registryAddress = QUrl(),
                      AllowedSchemas allowedSchemas=BuiltInSchemasOnly, QObject *parent = nullptr);
    QRemoteObjectHost(const QUrl &address, QObject *parent);
    ~QRemoteObjectHost() override;
    QUrl hostUrl() const override;
    bool setHostUrl(const QUrl &hostAddress, AllowedSchemas allowedSchemas=BuiltInSchemasOnly) override;

Q_SIGNALS:
    void hostUrlChanged();

protected:
    QRemoteObjectHost(QRemoteObjectHostPrivate &, QObject *);

private:
    Q_DECLARE_PRIVATE(QRemoteObjectHost)
};

class Q_REMOTEOBJECTS_EXPORT QRemoteObjectRegistryHost : public QRemoteObjectHostBase
{
    Q_OBJECT
public:
    QRemoteObjectRegistryHost(const QUrl &registryAddress = QUrl(), QObject *parent = nullptr);
    ~QRemoteObjectRegistryHost() override;
    bool setRegistryUrl(const QUrl &registryUrl) override;

protected:
    QRemoteObjectRegistryHost(QRemoteObjectRegistryHostPrivate &, QObject *);

private:
    Q_DECLARE_PRIVATE(QRemoteObjectRegistryHost)
};

QT_END_NAMESPACE

#endif
