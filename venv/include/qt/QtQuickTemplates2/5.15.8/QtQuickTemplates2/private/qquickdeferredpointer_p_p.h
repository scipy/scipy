/****************************************************************************
**
** Copyright (C) 2017 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the Qt Quick Templates 2 module of the Qt Toolkit.
**
** $QT_BEGIN_LICENSE:LGPL3$
** Commercial License Usage
** Licensees holding valid commercial Qt licenses may use this file in
** accordance with the commercial license agreement provided with the
** Software or, alternatively, in accordance with the terms contained in
** a written agreement between you and The Qt Company. For licensing terms
** and conditions see http://www.qt.io/terms-conditions. For further
** information use the contact form at http://www.qt.io/contact-us.
**
** GNU Lesser General Public License Usage
** Alternatively, this file may be used under the terms of the GNU Lesser
** General Public License version 3 as published by the Free Software
** Foundation and appearing in the file LICENSE.LGPLv3 included in the
** packaging of this file. Please review the following information to
** ensure the GNU Lesser General Public License version 3 requirements
** will be met: https://www.gnu.org/licenses/lgpl.html.
**
** GNU General Public License Usage
** Alternatively, this file may be used under the terms of the GNU
** General Public License version 2.0 or later as published by the Free
** Software Foundation and appearing in the file LICENSE.GPL included in
** the packaging of this file. Please review the following information to
** ensure the GNU General Public License version 2.0 requirements will be
** met: http://www.gnu.org/licenses/gpl-2.0.html.
**
** $QT_END_LICENSE$
**
****************************************************************************/

#ifndef QQUICKDEFERREDPOINTER_P_P_H
#define QQUICKDEFERREDPOINTER_P_P_H

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

template<typename T>
class QQuickDeferredPointer
{
public:
    inline QQuickDeferredPointer();
    inline QQuickDeferredPointer(T *);
    inline QQuickDeferredPointer(const QQuickDeferredPointer<T> &o);

    inline bool isNull() const;

    inline bool wasExecuted() const;
    inline void setExecuted();

    inline bool isExecuting() const;
    inline void setExecuting(bool);

    inline operator T*() const;
    inline operator bool() const;

    inline T *data() const;
    inline T *operator*() const;
    inline T *operator->() const;

    inline QQuickDeferredPointer<T> &operator=(T *);
    inline QQuickDeferredPointer<T> &operator=(const QQuickDeferredPointer &o);

private:
    quintptr ptr_value = 0;

    static const quintptr WasExecutedBit = 0x1;
    static const quintptr IsExecutingBit = 0x2;
    static const quintptr FlagsMask = WasExecutedBit | IsExecutingBit;
};

template<typename T>
QQuickDeferredPointer<T>::QQuickDeferredPointer()
{
}

template<typename T>
QQuickDeferredPointer<T>::QQuickDeferredPointer(T *v)
: ptr_value(quintptr(v))
{
    Q_ASSERT((ptr_value & FlagsMask) == 0);
}

template<typename T>
QQuickDeferredPointer<T>::QQuickDeferredPointer(const QQuickDeferredPointer<T> &o)
: ptr_value(o.ptr_value)
{
}

template<typename T>
bool QQuickDeferredPointer<T>::isNull() const
{
    return 0 == (ptr_value & (~FlagsMask));
}

template<typename T>
bool QQuickDeferredPointer<T>::wasExecuted() const
{
    return ptr_value & WasExecutedBit;
}

template<typename T>
void QQuickDeferredPointer<T>::setExecuted()
{
    ptr_value |= WasExecutedBit;
}

template<typename T>
bool QQuickDeferredPointer<T>::isExecuting() const
{
    return ptr_value & IsExecutingBit;
}

template<typename T>
void QQuickDeferredPointer<T>::setExecuting(bool b)
{
    if (b)
        ptr_value |= IsExecutingBit;
    else
        ptr_value &= ~IsExecutingBit;
}

template<typename T>
QQuickDeferredPointer<T>::operator T*() const
{
    return data();
}

template<typename T>
QQuickDeferredPointer<T>::operator bool() const
{
    return !isNull();
}

template<typename T>
T *QQuickDeferredPointer<T>::data() const
{
    return (T *)(ptr_value & ~FlagsMask);
}

template<typename T>
T *QQuickDeferredPointer<T>::operator*() const
{
    return (T *)(ptr_value & ~FlagsMask);
}

template<typename T>
T *QQuickDeferredPointer<T>::operator->() const
{
    return (T *)(ptr_value & ~FlagsMask);
}

template<typename T>
QQuickDeferredPointer<T> &QQuickDeferredPointer<T>::operator=(T *o)
{
    Q_ASSERT((quintptr(o) & FlagsMask) == 0);

    ptr_value = quintptr(o) | (ptr_value & FlagsMask);
    return *this;
}

template<typename T>
QQuickDeferredPointer<T> &QQuickDeferredPointer<T>::operator=(const QQuickDeferredPointer &o)
{
    ptr_value = o.ptr_value;
    return *this;
}

QT_END_NAMESPACE

#endif // QQUICKDEFERREDPOINTER_P_P_H
