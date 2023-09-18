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

#ifndef QQREMOTEOBJECTREPLICA_H
#define QQREMOTEOBJECTREPLICA_H

#include <QtRemoteObjects/qtremoteobjectglobal.h>

#include <QtCore/qsharedpointer.h>

QT_BEGIN_NAMESPACE

class QObjectPrivate;
class QRemoteObjectPendingCall;
class QRemoteObjectReplicaImplementation;
class QReplicaImplementationInterface;
class QRemoteObjectNode;

class Q_REMOTEOBJECTS_EXPORT QRemoteObjectReplica : public QObject
{
    Q_OBJECT
    Q_PROPERTY(QRemoteObjectNode *node READ node WRITE setNode)
    Q_PROPERTY(State state READ state NOTIFY stateChanged)
public:
    enum State {
        Uninitialized,
        Default,
        Valid,
        Suspect,
        SignatureMismatch
    };
    Q_ENUM(State)

public:
    ~QRemoteObjectReplica() override;

    bool isReplicaValid() const;
    bool waitForSource(int timeout = 30000);
    bool isInitialized() const;
    State state() const;
    QRemoteObjectNode *node() const;
    virtual void setNode(QRemoteObjectNode *node);

Q_SIGNALS:
    void initialized();
    void notified();
    void stateChanged(State state, State oldState);

protected:
    enum ConstructorType {DefaultConstructor, ConstructWithNode};
    explicit QRemoteObjectReplica(ConstructorType t = DefaultConstructor);
    QRemoteObjectReplica(QObjectPrivate &dptr, QObject *parent);

    virtual void initialize();
    void send(QMetaObject::Call call, int index, const QVariantList &args);
    QRemoteObjectPendingCall sendWithReply(QMetaObject::Call call, int index, const QVariantList &args);

protected:
    void setProperties(const QVariantList &);
    void setChild(int i, const QVariant &);
    const QVariant propAsVariant(int i) const;
    void persistProperties(const QString &repName, const QByteArray &repSig, const QVariantList &props) const;
    QVariantList retrieveProperties(const QString &repName, const QByteArray &repSig) const;
    void initializeNode(QRemoteObjectNode *node, const QString &name = QString());
    QSharedPointer<QReplicaImplementationInterface> d_impl;
private:
    friend class QRemoteObjectNodePrivate;
    friend class QConnectedReplicaImplementation;
};

QT_END_NAMESPACE

#endif
