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

#ifndef QTREMOTEOBJECTGLOBAL_H
#define QTREMOTEOBJECTGLOBAL_H

#include <QtCore/qglobal.h>
#include <QtCore/qhash.h>
#include <QtCore/qurl.h>
#include <QtCore/qloggingcategory.h>

QT_BEGIN_NAMESPACE

struct QRemoteObjectSourceLocationInfo
{
    QRemoteObjectSourceLocationInfo() = default;
    QRemoteObjectSourceLocationInfo(const QString &typeName_, const QUrl &hostUrl_)
        : typeName(typeName_), hostUrl(hostUrl_) {}

    inline bool operator==(const QRemoteObjectSourceLocationInfo &other) const Q_DECL_NOTHROW
    {
        return other.typeName == typeName && other.hostUrl == hostUrl;
    }
    inline bool operator!=(const QRemoteObjectSourceLocationInfo &other) const Q_DECL_NOTHROW
    {
        return !(*this == other);
    }

    QString typeName;
    QUrl hostUrl;
};

inline QDebug operator<<(QDebug dbg, const QRemoteObjectSourceLocationInfo &info)
{
    dbg.nospace() << "SourceLocationInfo(" << info.typeName << ", " << info.hostUrl << ")";
    return dbg.space();
}

inline QDataStream& operator<<(QDataStream &stream, const QRemoteObjectSourceLocationInfo &info)
{
    return stream << info.typeName << info.hostUrl;
}

inline QDataStream& operator>>(QDataStream &stream, QRemoteObjectSourceLocationInfo &info)
{
    return stream >> info.typeName >> info.hostUrl;
}

typedef QPair<QString, QRemoteObjectSourceLocationInfo> QRemoteObjectSourceLocation;
typedef QHash<QString, QRemoteObjectSourceLocationInfo> QRemoteObjectSourceLocations;
typedef QHash<int, QByteArray> QIntHash;

QT_END_NAMESPACE
Q_DECLARE_METATYPE(QRemoteObjectSourceLocation)
Q_DECLARE_METATYPE(QRemoteObjectSourceLocations)
Q_DECLARE_METATYPE(QIntHash)
QT_BEGIN_NAMESPACE

#ifndef QT_STATIC
#  if defined(QT_BUILD_REMOTEOBJECTS_LIB)
#    define Q_REMOTEOBJECTS_EXPORT Q_DECL_EXPORT
#  else
#    define Q_REMOTEOBJECTS_EXPORT Q_DECL_IMPORT
#  endif
#else
#  define Q_REMOTEOBJECTS_EXPORT
#endif

#define QCLASSINFO_REMOTEOBJECT_TYPE "RemoteObject Type"
#define QCLASSINFO_REMOTEOBJECT_SIGNATURE "RemoteObject Signature"

class QDataStream;

namespace QRemoteObjectStringLiterals {

// when QStringLiteral is used with the same string in different functions,
// it creates duplicate static data. Wrapping it in inline functions prevents it.

inline QString local() { return QStringLiteral("local"); }
inline QString tcp() { return QStringLiteral("tcp"); }
inline QString CLASS() { return QStringLiteral("Class::%1"); }
inline QString MODEL() { return QStringLiteral("Model::%1"); }
inline QString QAIMADAPTER() { return QStringLiteral("QAbstractItemModelAdapter"); }

}

Q_DECLARE_LOGGING_CATEGORY(QT_REMOTEOBJECT)
Q_DECLARE_LOGGING_CATEGORY(QT_REMOTEOBJECT_MODELS)
Q_DECLARE_LOGGING_CATEGORY(QT_REMOTEOBJECT_IO)

namespace QtRemoteObjects {

Q_NAMESPACE

Q_REMOTEOBJECTS_EXPORT void copyStoredProperties(const QMetaObject *mo, const void *src, void *dst);
Q_REMOTEOBJECTS_EXPORT void copyStoredProperties(const QMetaObject *mo, const void *src, QDataStream &dst);
Q_REMOTEOBJECTS_EXPORT void copyStoredProperties(const QMetaObject *mo, QDataStream &src, void *dst);

QString getTypeNameAndMetaobjectFromClassInfo(const QMetaObject *& meta);

template <typename T>
void copyStoredProperties(const T *src, T *dst)
{
    copyStoredProperties(&T::staticMetaObject, src, dst);
}

template <typename T>
void copyStoredProperties(const T *src, QDataStream &dst)
{
    copyStoredProperties(&T::staticMetaObject, src, dst);
}

template <typename T>
void copyStoredProperties(QDataStream &src, T *dst)
{
    copyStoredProperties(&T::staticMetaObject, src, dst);
}

enum QRemoteObjectPacketTypeEnum
{
    Invalid = 0,
    Handshake,
    InitPacket,
    InitDynamicPacket,
    AddObject,
    RemoveObject,
    InvokePacket,
    InvokeReplyPacket,
    PropertyChangePacket,
    ObjectList,
    Ping,
    Pong
};
Q_ENUM_NS(QRemoteObjectPacketTypeEnum)

enum InitialAction {
    FetchRootSize,
    PrefetchData
};
Q_ENUM_NS(InitialAction)

}

QT_END_NAMESPACE

#endif // QTREMOTEOBJECTSGLOBAL_H
