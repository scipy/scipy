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

#ifndef QQMLLISTCOMPOSITOR_P_H
#define QQMLLISTCOMPOSITOR_P_H

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
#include <QtCore/qvector.h>

#include <private/qqmlchangeset_p.h>

#include <QtCore/qdebug.h>

QT_BEGIN_NAMESPACE

class Q_AUTOTEST_EXPORT QQmlListCompositor
{
public:
    enum { MinimumGroupCount = 3, MaximumGroupCount = 11 };

    enum Group
    {
        Cache   = 0,
        Default = 1,
        Persisted = 2
    };

    enum Flag
    {
        CacheFlag       = 1 << Cache,
        DefaultFlag     = 1 << Default,
        PersistedFlag   = 1 << Persisted,
        PrependFlag     = 0x10000000,
        AppendFlag      = 0x20000000,
        UnresolvedFlag  = 0x40000000,
        MovedFlag       = 0x80000000,
        GroupMask       = ~(PrependFlag | AppendFlag | UnresolvedFlag | MovedFlag | CacheFlag)
    };

    class Range
    {
    public:
        Range() : next(this), previous(this) {}
        Range(Range *next, void *list, int index, int count, uint flags)
            : next(next), previous(next->previous), list(list), index(index), count(count), flags(flags) {
            next->previous = this; previous->next = this; }

        Range *next;
        Range *previous;
        void *list = nullptr;
        int index = 0;
        int count = 0;
        uint flags = 0;

        inline int start() const { return index; }
        inline int end() const { return index + count; }

        inline int groups() const { return flags & GroupMask; }

        inline bool inGroup() const { return flags & GroupMask; }
        inline bool inCache() const { return flags & CacheFlag; }
        inline bool inGroup(int group) const { return flags & (1 << group); }
        inline bool isUnresolved() const { return flags & UnresolvedFlag; }

        inline bool prepend() const { return flags & PrependFlag; }
        inline bool append() const { return flags & AppendFlag; }
    };

    class Q_AUTOTEST_EXPORT iterator
    {
    public:
        inline iterator();
        inline iterator(const iterator &it) = default;
        inline iterator(Range *range, int offset, Group group, int groupCount);
        inline ~iterator() {}

        bool operator ==(const iterator &it) const { return range == it.range && offset == it.offset; }
        bool operator !=(const iterator &it) const { return range != it.range || offset != it.offset; }

        bool operator ==(Group group) const { return range->flags & (1 << group); }
        bool operator !=(Group group) const { return !(range->flags & (1 << group)); }

        Range *&operator *() { return range; }
        Range * const &operator *() const { return range; }
        Range *operator ->() { return range; }
        const Range *operator ->() const { return range; }

        iterator &operator=(const iterator &) = default;
        iterator &operator +=(int difference);

        template<typename T> T *list() const { return static_cast<T *>(range->list); }
        int modelIndex() const { return range->index + offset; }

        void incrementIndexes(int difference) { incrementIndexes(difference, range->flags); }
        void decrementIndexes(int difference) { decrementIndexes(difference, range->flags); }

        inline void incrementIndexes(int difference, uint flags);
        inline void decrementIndexes(int difference, uint flags);

        void setGroup(Group g) { group = g; groupFlag = 1 << g; }

        Range *range = nullptr;
        int offset = 0;
        Group group = Default;
        int groupFlag;
        int groupCount = 0;
        union {
            struct {
                int cacheIndex;
            };
            int index[MaximumGroupCount];
        };
    };

    class Q_AUTOTEST_EXPORT insert_iterator : public iterator
    {
    public:
        inline insert_iterator() {}
        inline insert_iterator(const iterator &it) : iterator(it) {}
        inline insert_iterator(Range *, int, Group, int);
        inline ~insert_iterator() {}

        insert_iterator &operator +=(int difference);
    };

    struct Change
    {
        inline Change() {}
        inline Change(const iterator &it, int count, uint flags, int moveId = -1);
        int count;
        uint flags;
        int moveId;
        union {
            struct {
                int cacheIndex;
            };
            int index[MaximumGroupCount];
        };

        inline bool isMove() const { return moveId >= 0; }
        inline bool inCache() const { return flags & CacheFlag; }
        inline bool inGroup() const { return flags & GroupMask; }
        inline bool inGroup(int group) const { return flags & (CacheFlag << group); }

        inline int groups() const { return flags & GroupMask; }
    };

    struct Insert : public Change
    {
        Insert() {}
        Insert(const iterator &it, int count, uint flags, int moveId = -1)
            : Change(it, count, flags, moveId) {}
    };

    struct Remove : public Change
    {
        Remove() {}
        Remove(const iterator &it, int count, uint flags, int moveId = -1)
            : Change(it, count, flags, moveId) {}
    };

    QQmlListCompositor();
    ~QQmlListCompositor();

    int defaultGroups() const { return m_defaultFlags & ~PrependFlag; }
    void setDefaultGroups(int groups) { m_defaultFlags = groups | PrependFlag; }
    void setDefaultGroup(Group group) { m_defaultFlags |= (1 << group); }
    void clearDefaultGroup(Group group) { m_defaultFlags &= ~(1 << group); }
    void setRemoveGroups(int groups) { m_removeFlags = PrependFlag | AppendFlag | groups; }
    void setGroupCount(int count);

    int count(Group group) const;
    iterator find(Group group, int index);
    iterator find(Group group, int index) const;
    insert_iterator findInsertPosition(Group group, int index);

    const iterator &end() { return m_end; }

    void append(void *list, int index, int count, uint flags, QVector<Insert> *inserts = nullptr);
    void insert(Group group, int before, void *list, int index, int count, uint flags, QVector<Insert> *inserts = nullptr);
    iterator insert(iterator before, void *list, int index, int count, uint flags, QVector<Insert> *inserts = nullptr);

    void setFlags(Group fromGroup, int from, int count, Group group, int flags, QVector<Insert> *inserts = nullptr);
    void setFlags(iterator from, int count, Group group, uint flags, QVector<Insert> *inserts = nullptr);
    void setFlags(Group fromGroup, int from, int count, uint flags, QVector<Insert> *inserts = nullptr) {
        setFlags(fromGroup, from, count, fromGroup, flags, inserts); }
    void setFlags(const iterator from, int count, uint flags, QVector<Insert> *inserts = nullptr) {
        setFlags(from, count, from.group, flags, inserts); }

    void clearFlags(Group fromGroup, int from, int count, Group group, uint flags, QVector<Remove> *removals = nullptr);
    void clearFlags(iterator from, int count, Group group, uint flags, QVector<Remove> *removals = nullptr);
    void clearFlags(Group fromGroup, int from, int count, uint flags, QVector<Remove> *removals = nullptr) {
        clearFlags(fromGroup, from, count, fromGroup, flags, removals); }
    void clearFlags(const iterator &from, int count, uint flags, QVector<Remove> *removals = nullptr) {
        clearFlags(from, count, from.group, flags, removals); }

    bool verifyMoveTo(Group fromGroup, int from, Group toGroup, int to, int count, Group group) const;

    void move(
            Group fromGroup,
            int from,
            Group toGroup,
            int to,
            int count,
            Group group,
            QVector<Remove> *removals = nullptr,
            QVector<Insert> *inserts = nullptr);
    void clear();

    void listItemsInserted(void *list, int index, int count, QVector<Insert> *inserts);
    void listItemsRemoved(void *list, int index, int count, QVector<Remove> *removals);
    void listItemsMoved(void *list, int from, int to, int count, QVector<Remove> *removals, QVector<Insert> *inserts);
    void listItemsChanged(void *list, int index, int count, QVector<Change> *changes);

    void transition(
            Group from,
            Group to,
            QVector<QQmlChangeSet::Change> *removes,
            QVector<QQmlChangeSet::Change> *inserts);

private:
    Range m_ranges;
    iterator m_end;
    iterator m_cacheIt;
    int m_groupCount;
    int m_defaultFlags;
    int m_removeFlags;
    int m_moveId;

    inline Range *insert(Range *before, void *list, int index, int count, uint flags);
    inline Range *erase(Range *range);

    struct MovedFlags
    {
        MovedFlags() {}
        MovedFlags(int moveId, uint flags) : moveId(moveId), flags(flags) {}

        int moveId;
        uint flags;
    };

    void listItemsRemoved(
            QVector<Remove> *translatedRemovals,
            void *list,
            QVector<QQmlChangeSet::Change> *removals,
            QVector<QQmlChangeSet::Change> *insertions = nullptr,
            QVector<MovedFlags> *movedFlags = nullptr);
    void listItemsInserted(
            QVector<Insert> *translatedInsertions,
            void *list,
            const QVector<QQmlChangeSet::Change> &insertions,
            const QVector<MovedFlags> *movedFlags = nullptr);
    void listItemsChanged(
            QVector<Change> *translatedChanges,
            void *list,
            const QVector<QQmlChangeSet::Change> &changes);

    friend Q_AUTOTEST_EXPORT QDebug operator <<(QDebug debug, const QQmlListCompositor &list);
};

Q_DECLARE_TYPEINFO(QQmlListCompositor::Change, Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QQmlListCompositor::Remove, Q_PRIMITIVE_TYPE);
Q_DECLARE_TYPEINFO(QQmlListCompositor::Insert, Q_PRIMITIVE_TYPE);

inline QQmlListCompositor::iterator::iterator() {}

inline QQmlListCompositor::iterator::iterator(
        Range *range, int offset, Group group, int groupCount)
    : range(range)
    , offset(offset)
    , group(group)
    , groupFlag(1 << group)
    , groupCount(groupCount)
{
    for (int i = 0; i < groupCount; ++i)
        index[i] = 0;
}

inline void QQmlListCompositor::iterator::incrementIndexes(int difference, uint flags)
{
    for (int i = 0; i < groupCount; ++i) {
        if (flags & (1 << i))
            index[i] += difference;
    }
}

inline void QQmlListCompositor::iterator::decrementIndexes(int difference, uint flags)
{
    for (int i = 0; i < groupCount; ++i) {
        if (flags & (1 << i))
            index[i] -= difference;
    }
}

inline QQmlListCompositor::insert_iterator::insert_iterator(
        Range *range, int offset, Group group, int groupCount)
    : iterator(range, offset, group, groupCount) {}

inline QQmlListCompositor::Change::Change(const iterator &it, int count, uint flags, int moveId)
    : count(count), flags(flags), moveId(moveId)
{
    for (int i = 0; i < MaximumGroupCount; ++i)
        index[i] = it.index[i];
}

Q_AUTOTEST_EXPORT QDebug operator <<(QDebug debug, const QQmlListCompositor::Group &group);
Q_AUTOTEST_EXPORT QDebug operator <<(QDebug debug, const QQmlListCompositor::Range &range);
Q_AUTOTEST_EXPORT QDebug operator <<(QDebug debug, const QQmlListCompositor::iterator &it);
Q_AUTOTEST_EXPORT QDebug operator <<(QDebug debug, const QQmlListCompositor::Change &change);
Q_AUTOTEST_EXPORT QDebug operator <<(QDebug debug, const QQmlListCompositor::Remove &remove);
Q_AUTOTEST_EXPORT QDebug operator <<(QDebug debug, const QQmlListCompositor::Insert &insert);
Q_AUTOTEST_EXPORT QDebug operator <<(QDebug debug, const QQmlListCompositor &list);

QT_END_NAMESPACE

#endif
