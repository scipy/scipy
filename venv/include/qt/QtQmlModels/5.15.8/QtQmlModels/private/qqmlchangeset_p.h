/****************************************************************************
**
** Copyright (C) 2016 The Qt Company Ltd.
** Contact: https://www.qt.io/licensing/
**
** This file is part of the QtQml module of the Qt Toolkit.
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

#ifndef QQMLCHANGESET_P_H
#define QQMLCHANGESET_P_H

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

#include <QtCore/qdebug.h>
#include <QtCore/qvector.h>
#include <QtQmlModels/private/qtqmlmodelsglobal_p.h>

QT_BEGIN_NAMESPACE

class Q_QMLMODELS_PRIVATE_EXPORT QQmlChangeSet
{
public:
    struct MoveKey
    {
        MoveKey() {}
        MoveKey(int moveId, int offset) : moveId(moveId), offset(offset) {}
        int moveId = -1;
        int offset = 0;
    };

    // The storrage for Change (below). This struct is trivial, which it has to be in order to store
    // it in a QV4::Heap::Base object. The Change struct doesn't add any storage fields, so it is
    // safe to cast ChangeData to/from Change.
    struct ChangeData
    {
        int index;
        int count;
        int moveId;
        int offset;
    };

    struct Change: ChangeData
    {
        Change() {
            index = 0;
            count = 0;
            moveId = -1;
            offset = 0;
        }
        Change(int index, int count, int moveId = -1, int offset = 0) {
            this->index = index;
            this->count = count;
            this->moveId = moveId;
            this->offset = offset;
        }

        bool isMove() const { return moveId >= 0; }

        MoveKey moveKey(int index) const {
            return MoveKey(moveId, index - Change::index + offset); }

        int start() const { return index; }
        int end() const { return index + count; }
    };

    QQmlChangeSet();
    QQmlChangeSet(const QQmlChangeSet &changeSet);
    ~QQmlChangeSet();

    QQmlChangeSet &operator =(const QQmlChangeSet &changeSet);

    const QVector<Change> &removes() const { return m_removes; }
    const QVector<Change> &inserts() const { return m_inserts; }
    const QVector<Change> &changes() const { return m_changes; }

    void insert(int index, int count);
    void remove(int index, int count);
    void move(int from, int to, int count, int moveId);
    void change(int index, int count);

    void insert(const QVector<Change> &inserts);
    void remove(const QVector<Change> &removes, QVector<Change> *inserts = nullptr);
    void move(const QVector<Change> &removes, const QVector<Change> &inserts);
    void change(const QVector<Change> &changes);
    void apply(const QQmlChangeSet &changeSet);

    bool isEmpty() const { return m_removes.empty() && m_inserts.empty() && m_changes.isEmpty(); }

    void clear()
    {
        m_removes.clear();
        m_inserts.clear();
        m_changes.clear();
        m_difference = 0;
    }

    int difference() const { return m_difference; }

private:
    void remove(QVector<Change> *removes, QVector<Change> *inserts);
    void change(QVector<Change> *changes);

    QVector<Change> m_removes;
    QVector<Change> m_inserts;
    QVector<Change> m_changes;
    int m_difference;
};

Q_DECLARE_TYPEINFO(QQmlChangeSet::Change, Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QQmlChangeSet::MoveKey, Q_PRIMITIVE_TYPE);

inline uint qHash(const QQmlChangeSet::MoveKey &key) { return qHash(qMakePair(key.moveId, key.offset)); }
inline bool operator ==(const QQmlChangeSet::MoveKey &l, const QQmlChangeSet::MoveKey &r) {
    return l.moveId == r.moveId && l.offset == r.offset; }

Q_QMLMODELS_PRIVATE_EXPORT QDebug operator <<(QDebug debug, const QQmlChangeSet::Change &change);
Q_QMLMODELS_PRIVATE_EXPORT QDebug operator <<(QDebug debug, const QQmlChangeSet &change);

QT_END_NAMESPACE

#endif
