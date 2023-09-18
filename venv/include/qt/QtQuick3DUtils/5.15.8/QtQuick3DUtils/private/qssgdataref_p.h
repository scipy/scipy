/****************************************************************************
**
** Copyright (C) 2008-2012 NVIDIA Corporation.
** Copyright (C) 2019 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of Qt Quick 3D.
**
** $QT_BEGIN_LICENSE:GPL$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see https://www.qt.io/terms-conditions. For further
** information use the contact form at https://www.qt.io/contact-us.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 3 or (at your option) any later version
** approved by the KDE Free Qt Foundation. The licenses are as published by
** the Free Software Foundation and appearing in the file LICENSE.GPL3
** included in the packaging of this file. Please review the following
** information to ensure the GNU General Public License requirements will
** be met: https://www.gnu.org/licenses/gpl-3.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QSSGDATAREF_H
#define QSSGDATAREF_H

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

#include <QtQuick3DUtils/private/qtquick3dutilsglobal_p.h>

#include <QtCore/qvector.h>
#include <QtCore/qbytearray.h>

QT_BEGIN_NAMESPACE

template<typename T>
struct QSSGDataView
{
    const T *mData;
    int mSize;

    explicit QSSGDataView(const QVector<T> &data) : mData(data.constBegin()), mSize(data.size()) { Q_ASSERT(mSize >= 0); }
    QSSGDataView(const T *inData, qint32 inSize) : mData(inData), mSize(inSize) { Q_ASSERT(mSize >= 0); }
    constexpr QSSGDataView() : mData(nullptr), mSize(0) {}

    qint32 size() const { return mSize; }

    const T *begin() const { return mData; }
    const T *end() const { return mData + mSize; }

    const T &operator[](int index) const
    {
        Q_ASSERT(index > -1);
        Q_ASSERT(index < mSize);
        return mData[index];
    }

    void clear()
    {
        mData = nullptr;
        mSize = 0;
    }

    operator const void *() { return reinterpret_cast<const void *>(mData); }
};

template<>
struct QSSGDataView<quint8>
{
    const quint8 *mData;
    int mSize;

    explicit QSSGDataView(const QByteArray &data)
        : mData(reinterpret_cast<const quint8 *>(data.constBegin())), mSize(data.size())
    { Q_ASSERT(mSize >= 0); }
    template<typename T>
    explicit QSSGDataView(const QVector<T> &data)
        : mData(reinterpret_cast<const quint8 *>(data.constBegin())), mSize(data.size()*sizeof(T))
    { Q_ASSERT(mSize >= 0); }
    QSSGDataView(const quint8 *inData, qint32 inSize) : mData(inData), mSize(inSize) { Q_ASSERT(mSize >= 0); }
    template<typename T>
    QSSGDataView(const T *inData, qint32 inSize)
        : mData(reinterpret_cast<const quint8 *>(inData)), mSize(inSize*sizeof(T))
    { Q_ASSERT(mSize >= 0); }
    constexpr QSSGDataView() : mData(nullptr), mSize(0) {}

    qint32 size() const { return mSize; }

    const quint8 *begin() const { return mData; }
    const quint8 *end() const { return mData + mSize; }

    const quint8 &operator[](int index) const
    {
        Q_ASSERT(index > -1);
        Q_ASSERT(index < mSize);
        return mData[index];
    }

    void clear()
    {
        mData = nullptr;
        mSize = 0;
    }

    operator const void *() { return reinterpret_cast<const void *>(mData); }
};

using QSSGByteView = QSSGDataView<quint8>;

template<typename T>
inline QSSGDataView<T> toDataView(const T &type)
{
    return QSSGDataView<T>(&type, 1);
}

template<typename T>
inline QSSGDataView<T> toDataView(const QVector<T> &type)
{
    return QSSGDataView<T>(type);
}

template<typename T>
inline QSSGByteView toByteView(const T &type)
{
    return QSSGByteView(&type, 1);
}

template<typename T>
inline QSSGByteView toByteView(const QVector<T> &type)
{
    return QSSGByteView(type);
}

template<>
inline QSSGByteView toByteView(const QByteArray &type)
{
    return QSSGByteView(type);
}

inline QSSGByteView toByteView(const char *str)
{
    return QSSGByteView(str, qstrlen(str));
}

template<typename T>
inline QSSGDataView<T> toDataView(const T *type, quint32 count)
{
    return QSSGDataView<T>(type, count);
}

template<typename T>
inline QSSGByteView toByteView(const T *type, quint32 count)
{
    return QSSGByteView(type, count);
}

template<typename T>
struct QSSGDataRef
{
    T *mData;
    qint32 mSize;

    QSSGDataRef(T *inData, qint32 inSize) : mData(inData), mSize(inSize) { Q_ASSERT(inSize >= 0); }
    QSSGDataRef() : mData(nullptr), mSize(0) {}
    qint32 size() const { return mSize; }

    T *begin() { return mData; }
    T *end() { return mData + mSize; }

    T *begin() const { return mData; }
    T *end() const { return mData + mSize; }

    T &operator[](qint32 index)
    {
        Q_ASSERT(index >= 0);
        Q_ASSERT(index < mSize);
        return mData[index];
    }

    const T &operator[](qint32 index) const
    {
        Q_ASSERT(index >= 0);
        Q_ASSERT(index < mSize);
        return mData[index];
    }

    void clear()
    {
        mData = nullptr;
        mSize = 0;
    }

    operator QSSGDataView<T>() const { return QSSGDataView<T>(mData, mSize); }
    operator void *() { return reinterpret_cast<void *>(mData); }
};

using QSSGByteRef = QSSGDataRef<quint8>;

template<typename T>
inline QSSGDataRef<T> toDataRef(T &type)
{
    return QSSGDataRef<T>(&type, 1);
}

template<typename T>
inline QSSGByteRef toByteRef(T &type)
{
    return QSSGByteRef(reinterpret_cast<quint8 *>(&type), sizeof(T));
}

template<typename T>
inline QSSGDataRef<T> toDataRef(T *type, quint32 count)
{
    return QSSGDataRef<T>(type, count);
}

template<typename T>
inline QSSGByteRef toByteRef(T *type, quint32 count)
{
    return QSSGByteRef(reinterpret_cast<quint8 *>(type), sizeof(T) * count);
}

QT_END_NAMESPACE

#endif // QSSGDATAREF_H
