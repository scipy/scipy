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

#ifndef QREMOTEOBJECTREPLICA_P_H
#define QREMOTEOBJECTREPLICA_P_H

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

#include "qremoteobjectreplica.h"

#include "qremoteobjectpendingcall.h"

#include "qremoteobjectpacket_p.h"

#include <QtCore/qpointer.h>
#include <QtCore/qvector.h>
#include <QtCore/qdatastream.h>
#include <QtCore/qcompilerdetection.h>
#include <QtCore/qtimer.h>

QT_BEGIN_NAMESPACE

class QRemoteObjectReplica;
class QRemoteObjectSource;
class IoDeviceBase;

class QReplicaImplementationInterface
{
public:
    virtual ~QReplicaImplementationInterface() {}
    virtual const QVariant getProperty(int i) const = 0;
    virtual void setProperties(const QVariantList &) = 0;
    virtual void setProperty(int i, const QVariant &) = 0;
    virtual bool isInitialized() const = 0;
    virtual QRemoteObjectReplica::State state() const = 0;
    virtual bool waitForSource(int) = 0;
    virtual QRemoteObjectNode *node() const = 0;

    virtual void _q_send(QMetaObject::Call call, int index, const QVariantList &args) = 0;
    virtual QRemoteObjectPendingCall _q_sendWithReply(QMetaObject::Call call, int index, const QVariantList &args) = 0;
};

class QStubReplicaImplementation final : public QReplicaImplementationInterface
{
public:
    explicit QStubReplicaImplementation();
    ~QStubReplicaImplementation() override;

    const QVariant getProperty(int i) const override;
    void setProperties(const QVariantList &) override;
    void setProperty(int i, const QVariant &) override;
    bool isInitialized() const override { return false; }
    QRemoteObjectReplica::State state() const override { return QRemoteObjectReplica::State::Uninitialized;}
    bool waitForSource(int) override { return false; }
    QRemoteObjectNode *node() const override { return nullptr; }

    void _q_send(QMetaObject::Call call, int index, const QVariantList &args) override;
    QRemoteObjectPendingCall _q_sendWithReply(QMetaObject::Call call, int index, const QVariantList &args) override;
    QVariantList m_propertyStorage;
};

class QRemoteObjectReplicaImplementation : public QObject, public QReplicaImplementationInterface
{
public:
    explicit QRemoteObjectReplicaImplementation(const QString &name, const QMetaObject *, QRemoteObjectNode *);
    ~QRemoteObjectReplicaImplementation() override;

    bool needsDynamicInitialization() const;

    const QVariant getProperty(int i) const override = 0;
    void setProperties(const QVariantList &) override = 0;
    void setProperty(int i, const QVariant &) override = 0;
    virtual bool isShortCircuit() const = 0;
    bool isInitialized() const override { return true; }
    QRemoteObjectReplica::State state() const override { return QRemoteObjectReplica::State(m_state.loadRelaxed()); }
    void setState(QRemoteObjectReplica::State state);
    bool waitForSource(int) override { return true; }
    virtual bool waitForFinished(const QRemoteObjectPendingCall &, int) { return true; }
    virtual void notifyAboutReply(int, const QVariant &) {}
    virtual void configurePrivate(QRemoteObjectReplica *);
    void emitInitialized();
    void emitNotified();
    QRemoteObjectNode *node() const override { return m_node; }

    void _q_send(QMetaObject::Call call, int index, const QVariantList &args) override = 0;
    QRemoteObjectPendingCall _q_sendWithReply(QMetaObject::Call call, int index, const QVariantList &args) override = 0;

    //Dynamic replica functions
    virtual void setDynamicMetaObject(const QMetaObject *meta);
    virtual void setDynamicProperties(const QVariantList &values);

    QString m_objectName;
    const QMetaObject *m_metaObject;

    //Dynamic Replica data
    int m_numSignals;//TODO maybe here too
    int m_methodOffset;
    int m_signalOffset;
    int m_propertyOffset;
    QRemoteObjectNode *m_node;
    QByteArray m_objectSignature;
    QAtomicInt m_state;
};

class QConnectedReplicaImplementation final : public QRemoteObjectReplicaImplementation
{
public:
    explicit QConnectedReplicaImplementation(const QString &name, const QMetaObject *, QRemoteObjectNode *);
    ~QConnectedReplicaImplementation() override;
    const QVariant getProperty(int i) const override;
    void setProperties(const QVariantList &) override;
    void setProperty(int i, const QVariant &) override;
    bool isShortCircuit() const final { return false; }
    bool isInitialized() const override;
    bool waitForSource(int timeout) override;
    QVector<int> childIndices() const;
    void initialize(QVariantList &values);
    void configurePrivate(QRemoteObjectReplica *) override;
    void requestRemoteObjectSource();
    bool sendCommand();
    QRemoteObjectPendingCall sendCommandWithReply(int serialId);
    bool waitForFinished(const QRemoteObjectPendingCall &call, int timeout) override;
    void notifyAboutReply(int ackedSerialId, const QVariant &value) override;
    void setConnection(IoDeviceBase *conn);
    void setDisconnected();

    void _q_send(QMetaObject::Call call, int index, const QVariantList &args) override;
    QRemoteObjectPendingCall _q_sendWithReply(QMetaObject::Call call, int index, const QVariantList& args) override;

    void setDynamicMetaObject(const QMetaObject *meta) override;
    void setDynamicProperties(const QVariantList&) override;
    QVector<QRemoteObjectReplica *> m_parentsNeedingConnect;
    QVariantList m_propertyStorage;
    QVector<int> m_childIndices;
    QPointer<IoDeviceBase> connectionToSource;

    // pending call data
    int m_curSerialId = 1; // 0 is reserved for heartbeat signals
    QHash<int, QRemoteObjectPendingCall> m_pendingCalls;
    QRemoteObjectPackets::DataStreamPacket m_packet;
    QTimer m_heartbeatTimer;
};

class QInProcessReplicaImplementation final : public QRemoteObjectReplicaImplementation
{
public:
    explicit QInProcessReplicaImplementation(const QString &name, const QMetaObject *, QRemoteObjectNode *);
    ~QInProcessReplicaImplementation() override;

    const QVariant getProperty(int i) const override;
    void setProperties(const QVariantList &) override;
    void setProperty(int i, const QVariant &) override;
    bool isShortCircuit() const final { return true; }

    void _q_send(QMetaObject::Call call, int index, const QVariantList &args) override;
    QRemoteObjectPendingCall _q_sendWithReply(QMetaObject::Call call, int index, const QVariantList& args) override;

    QPointer<QRemoteObjectSourceBase> connectionToSource;
};

QT_END_NAMESPACE

#endif
