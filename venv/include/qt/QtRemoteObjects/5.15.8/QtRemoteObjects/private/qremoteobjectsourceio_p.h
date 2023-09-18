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

#ifndef QREMOTEOBJECTSOURCEIO_P_H
#define QREMOTEOBJECTSOURCEIO_P_H

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

#include "qconnectionfactories_p.h"
#include "qtremoteobjectglobal.h"
#include "qremoteobjectpacket_p.h"

#include <QtCore/qiodevice.h>
#include <QtCore/qscopedpointer.h>

QT_BEGIN_NAMESPACE

class QRemoteObjectSourceBase;
class QRemoteObjectRootSource;
class SourceApiMap;
class QRemoteObjectHostBase;

class QRemoteObjectSourceIo : public QObject
{
    Q_OBJECT
public:
    explicit QRemoteObjectSourceIo(const QUrl &address, QObject *parent = nullptr);
    explicit QRemoteObjectSourceIo(QObject *parent = nullptr);
    ~QRemoteObjectSourceIo() override;

    bool startListening();
    bool enableRemoting(QObject *object, const QMetaObject *meta, const QString &name,
                        const QString &typeName);
    bool enableRemoting(QObject *object, const SourceApiMap *api, QObject *adapter = nullptr);
    bool disableRemoting(QObject *object);
    void newConnection(IoDeviceBase *conn);

    QUrl serverAddress() const;

public Q_SLOTS:
    void handleConnection();
    void onServerDisconnect(QObject *obj = nullptr);
    void onServerRead(QObject *obj);

Q_SIGNALS:
    void remoteObjectAdded(const QRemoteObjectSourceLocation &);
    void remoteObjectRemoved(const QRemoteObjectSourceLocation &);
    void serverRemoved(const QUrl& url);

public:
    void registerSource(QRemoteObjectSourceBase *source);
    void unregisterSource(QRemoteObjectSourceBase *source);

    QHash<QIODevice*, quint32> m_readSize;
    QSet<IoDeviceBase*> m_connections;
    QHash<QObject *, QRemoteObjectRootSource*> m_objectToSourceMap;
    QMap<QString, QRemoteObjectSourceBase*> m_sourceObjects;
    QMap<QString, QRemoteObjectRootSource*> m_sourceRoots;
    QHash<IoDeviceBase*, QUrl> m_registryMapping;
    QScopedPointer<QConnectionAbstractServer> m_server;
    QRemoteObjectPackets::DataStreamPacket m_packet;
    QString m_rxName;
    QVariantList m_rxArgs;
    QUrl m_address;
};

QT_END_NAMESPACE

#endif
