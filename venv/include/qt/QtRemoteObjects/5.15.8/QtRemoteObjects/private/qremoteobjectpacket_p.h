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

#ifndef QTREMOTEOBJECTPACKET_P_H
#define QTREMOTEOBJECTPACKET_P_H

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

#include "qtremoteobjectglobal.h"
#include "qremoteobjectsource.h"
#include "qconnectionfactories_p.h"

#include <QtCore/qhash.h>
#include <QtCore/qmap.h>
#include <QtCore/qpair.h>
#include <QtCore/qurl.h>
#include <QtCore/qvariant.h>
#include <QtCore/qloggingcategory.h>

#include <cstdlib>

QT_BEGIN_NAMESPACE

class QMetaObjectBuilder;
class QRemoteObjectSourceBase;
class QRemoteObjectRootSource;

namespace QRemoteObjectPackets {

Q_NAMESPACE

class DataStreamPacket;

struct ObjectInfo
{
    QString name;
    QString typeName;
    QByteArray signature;
};

inline QDebug operator<<(QDebug dbg, const ObjectInfo &info)
{
    dbg.nospace() << "ObjectInfo(" << info.name << ", " << info.typeName << ", " << info.signature <<")";
    return dbg.space();
}

inline QDataStream& operator<<(QDataStream &stream, const ObjectInfo &info)
{
    return stream << info.name << info.typeName << info.signature;
}

inline QDataStream& operator>>(QDataStream &stream, ObjectInfo &info)
{
    return stream >> info.name >> info.typeName >> info.signature;
}

typedef QVector<ObjectInfo> ObjectInfoList;

enum class ObjectType : quint8 { CLASS, MODEL, GADGET };
Q_ENUM_NS(ObjectType)

// Use a short name, as QVariant::save writes the name every time a qvariant of
// this type is serialized
class QRO_
{
public:
    QRO_() : type(ObjectType::CLASS), isNull(true) {}
    explicit QRO_(QRemoteObjectSourceBase *source);
    explicit QRO_(const QVariant &value);
    QString name, typeName;
    ObjectType type;
    bool isNull;
    QByteArray classDefinition;
    QByteArray parameters;
};

inline QDebug operator<<(QDebug dbg, const QRO_ &info)
{
    dbg.nospace() << "QRO_(name: " << info.name << ", typeName: " << info.typeName << ", type: " << info.type
                  << ", valid: " << (info.isNull ? "true" : "false") << ", paremeters: {" << info.parameters <<")"
                  << (info.classDefinition.isEmpty() ? " no definitions)" : " with definitions)");
    return dbg.space();
}

QDataStream& operator<<(QDataStream &stream, const QRO_ &info);

QDataStream& operator>>(QDataStream &stream, QRO_ &info);

void serializeObjectListPacket(DataStreamPacket&, const ObjectInfoList&);
void deserializeObjectListPacket(QDataStream&, ObjectInfoList&);

//Helper class for creating a QByteArray from a QRemoteObjectPacket
class DataStreamPacket : public QDataStream
{
public:
    DataStreamPacket(quint16 id = QtRemoteObjects::InvokePacket)
        : QDataStream(&array, QIODevice::WriteOnly)
        , baseAddress(0)
        , size(0)
    {
        this->setVersion(QtRemoteObjects::dataStreamVersion);
        *this << quint32(0);
        *this << id;
    }
    void setId(quint16 id)
    {
        device()->seek(baseAddress);
        *this << quint32(0);
        *this << id;
    }

    void finishPacket()
    {
        size = device()->pos();
        device()->seek(baseAddress);
        *this << quint32(size - baseAddress - sizeof(quint32));
    }
    QByteArray array;
    int baseAddress;
    int size;

private:
    Q_DISABLE_COPY(DataStreamPacket)
};

const QVariant encodeVariant(const QVariant &value);
QVariant &decodeVariant(QVariant &value, int type);

void serializeProperty(QDataStream &, const QRemoteObjectSourceBase *source, int internalIndex);

void serializeHandshakePacket(DataStreamPacket &);
void serializeInitPacket(DataStreamPacket &, const QRemoteObjectRootSource*);
void serializeProperties(DataStreamPacket &, const QRemoteObjectSourceBase*);
void deserializeInitPacket(QDataStream &, QVariantList&);

void serializeInitDynamicPacket(DataStreamPacket &, const QRemoteObjectRootSource*);
void serializeDefinition(QDataStream &, const QRemoteObjectSourceBase*);

void serializeAddObjectPacket(DataStreamPacket &, const QString &name, bool isDynamic);
void deserializeAddObjectPacket(QDataStream &, bool &isDynamic);

void serializeRemoveObjectPacket(DataStreamPacket&, const QString &name);
//There is no deserializeRemoveObjectPacket - no parameters other than id and name

void serializeInvokePacket(DataStreamPacket&, const QString &name, int call, int index, const QVariantList &args, int serialId = -1, int propertyIndex = -1);
void deserializeInvokePacket(QDataStream& in, int &call, int &index, QVariantList &args, int &serialId, int &propertyIndex);

void serializeInvokeReplyPacket(DataStreamPacket&, const QString &name, int ackedSerialId, const QVariant &value);
void deserializeInvokeReplyPacket(QDataStream& in, int &ackedSerialId, QVariant &value);

//TODO do we need the object name or could we go with an id in backend code, this could be a costly allocation
void serializePropertyChangePacket(QRemoteObjectSourceBase *source, int signalIndex);
void deserializePropertyChangePacket(QDataStream& in, int &index, QVariant &value);

// Heartbeat packets
void serializePingPacket(DataStreamPacket &ds, const QString &name);
void serializePongPacket(DataStreamPacket &ds, const QString &name);


} // namespace QRemoteObjectPackets

QT_END_NAMESPACE

Q_DECLARE_METATYPE(QRemoteObjectPackets::QRO_)

#endif
