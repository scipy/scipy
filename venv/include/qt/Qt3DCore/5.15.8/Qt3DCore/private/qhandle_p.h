/****************************************************************************
**
** Copyright (C) 2014 Klaralvdalens Datakonsult AB (KDAB).
** Contact: https://www.qt.io/licensing/
**
** This file is part of the Qt3D module of the Qt Toolkit.
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

#ifndef QT3DCORE_QHANDLE_P_H
#define QT3DCORE_QHANDLE_P_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists for the convenience
// of other Qt classes.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.
//

#include <Qt3DCore/qt3dcore_global.h>
#include <QtCore/QDebug>
#include <QtCore/qhashfunctions.h>

QT_BEGIN_NAMESPACE

namespace Qt3DCore {

template <typename T>
class QHandle
{
public:
    struct Data {
        union {
            quintptr counter;
            Data *nextFree;
        };
    };
    QHandle()
        : d(nullptr),
          counter(0)
    {}
    QHandle(Data *d)
        : d(d),
          counter(d->counter)
    {
    }
    QHandle(const QHandle &other)
        : d(other.d),
          counter(other.counter)
    {
    }
    QHandle &operator=(const QHandle &other)
    {
        d = other.d;
        counter = other.counter;
        return *this;
    }

    inline T *operator->() const;
    T *data() const;

    quintptr handle() const { return reinterpret_cast<quintptr>(d); }
    bool isNull() const { return !d; }

    Data *data_ptr() const { return d; }

    bool operator==(const QHandle &other) const { return d == other.d && counter == other.counter; }
    bool operator!=(const QHandle &other) const { return !operator==(other); }
private:
    Data *d;
    quintptr counter;
};


template <typename T>
QDebug operator<<(QDebug dbg, const QHandle<T> &h)
{
    QDebugStateSaver saver(dbg);
    QString binNumber = QString::number(h.handle(), 2).rightJustified(32, QChar::fromLatin1('0'));
    dbg.nospace() << " m_handle = " << h.handle()
                  << " = " << binNumber;
    return dbg;
}

template <typename T>
uint qHash(const QHandle<T> &h, uint seed)
{
    using QT_PREPEND_NAMESPACE(qHash);
    return qHash(h.handle(), seed);
}

} // Qt3DCore

// simpler than fighting the Q_DECLARE_TYPEINFO macro, use QString as a dummy to get movable semantics
template <typename T>
class QTypeInfo<Qt3DCore::QHandle<T> >
    : public QTypeInfoMerger<Qt3DCore::QHandle<T>, QString> {};

QT_END_NAMESPACE

#endif // QT3DCORE_QRHANDLE_H
