/****************************************************************************
**
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

#ifndef QTSSGGLOBAL_P_H
#define QTSSGGLOBAL_P_H

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

#include <QtCore/qglobal.h>

QT_BEGIN_NAMESPACE

#ifndef QT_STATIC
#if defined(QT_BUILD_QUICK3DUTILS_LIB)
#define Q_QUICK3DUTILS_EXPORT Q_DECL_EXPORT
#else
#define Q_QUICK3DUTILS_EXPORT Q_DECL_IMPORT
#endif
#else
#define Q_QUICK3DUTILS_EXPORT
#endif

#define Q_QUICK3DUTILS_PRIVATE_EXPORT Q_QUICK3DUTILS_EXPORT

template<typename T>
class QSSGRef
{
    T *d;

public:
    T *data() const { return d; }
    T *get() const { return d; }
    T *take()
    {
        T *t = d;
        d = nullptr;
        return t;
    }
    bool isNull() const { return !d; }
    operator bool() const { return d; }
    bool operator!() const { return !d; }
    T &operator*() const { return *d; }
    T *operator->() const { return d; }

    // constructors
    constexpr QSSGRef() : d(nullptr) {}
    constexpr QSSGRef(std::nullptr_t) : d(nullptr) {}

    QSSGRef(T *ptr) : d(ptr)
    {
        if (d)
            d->ref.ref();
    }
    template<typename X>
    QSSGRef(X *ptr) : d(ptr)
    {
        if (d)
            d->ref.ref();
    }
    QSSGRef(const QSSGRef<T> &other) : d(other.d)
    {
        if (d)
            d->ref.ref();
    }
    template<typename X>
    QSSGRef(const QSSGRef<X> &other) : d(other.get())
    {
        if (d)
            d->ref.ref();
    }

    QSSGRef(QSSGRef<T> &&other) : d(other.take()) {}
    template<typename X>
    QSSGRef(QSSGRef<X> &&other) : d(other.take())
    {
    }

    ~QSSGRef()
    {
        if (d && !d->ref.deref())
            delete d;
    }

    template<typename X>
    QSSGRef<T> &operator=(const QSSGRef<X> &other)
    {
        if (d != other.get()) {
            if (d && !d->ref.deref())
                delete d;
            d = other.get();
            if (d)
                d->ref.ref();
        }
        return *this;
    }
    template<typename X>
    QSSGRef<T> &operator=(QSSGRef<X> &&other)
    {
        clear();
        d = other.take();
        return *this;
    }
    QSSGRef<T> &operator=(const QSSGRef<T> &other)
    {
        if (d != other.get()) {
            if (d && !d->ref.deref())
                delete d;
            d = other.get();
            if (d)
                d->ref.ref();
        }
        return *this;
    }
    QSSGRef<T> &operator=(QSSGRef<T> &&other)
    {
        qSwap(d, other.d);
        return *this;
    }

    void swap(QSSGRef<T> &other) { qSwap(d, other.d); }

    void clear()
    {
        if (d && !d->ref.deref())
            delete d;
        d = nullptr;
    }
};

template<class T, class X>
bool operator==(const QSSGRef<T> &ptr1, const QSSGRef<X> &ptr2)
{
    return ptr1.get() == ptr2.get();
}
template<class T, class X>
bool operator!=(const QSSGRef<T> &ptr1, const QSSGRef<X> &ptr2)
{
    return ptr1.get() != ptr2.get();
}
template<class T, class X>
bool operator==(const QSSGRef<T> &ptr1, const X *ptr2)
{
    return ptr1.get() == ptr2;
}
template<class T, class X>
bool operator!=(const QSSGRef<T> &ptr1, const X *ptr2)
{
    return ptr1.get() != ptr2;
}
template<class T, class X>
bool operator==(const T *ptr1, const QSSGRef<X> &ptr2)
{
    return ptr1 == ptr2.get();
}
template<class T, class X>
bool operator!=(const T *ptr1, const QSSGRef<X> &ptr2)
{
    return ptr1 != ptr2.get();
}
template<class T>
bool operator==(const QSSGRef<T> &lhs, std::nullptr_t)
{
    return !lhs.get();
}
template<class T>
bool operator!=(const QSSGRef<T> &lhs, std::nullptr_t)
{
    return lhs.get();
}
template<class T>
bool operator==(std::nullptr_t, const QSSGRef<T> &rhs)
{
    return !rhs.get();
}
template<class T>
bool operator!=(std::nullptr_t, const QSSGRef<T> &rhs)
{
    return rhs.get();
}

QT_END_NAMESPACE

#endif // QTSSGGLOBAL_P_H
