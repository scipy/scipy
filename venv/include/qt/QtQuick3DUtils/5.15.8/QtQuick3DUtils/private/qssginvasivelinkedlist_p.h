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

#ifndef QSSGINVASIVELINKEDLIST_H
#define QSSGINVASIVELINKEDLIST_H

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

QT_BEGIN_NAMESPACE

// Used for singly linked list where
// items have either no head or tail ptr.
template<typename T>
struct QSSGNullOp
{
    static void set(T &, T *) {}
    static T *get(const T &) { return nullptr; }
};

template <typename T, T *T::*n>
struct QSSGListAccessorNext
{
    static inline T *get(T &o) { return o.*n; }
    static inline const T *get(const T &o) { return o.*n; }
    static inline void set(T &o, T *next) { o.*n = next; }
};

template <typename T, T *T::*p>
struct QSSGListAccessorPrevious
{
    static inline T *get(T &o) { return o.*p; }
    static inline const T *get(const T &o) { return o.*p; }
    static inline void set(T &o, T *next) { o.*p = next; }
};

// Base linked list without an included head or tail member.
template<typename T, typename HeadOp, typename TailOp>
struct QSSGInvasiveLinkListBase
{
    inline T *tail(T *inObj) { return inObj ? TailOp::get(inObj) : nullptr; }
    inline T *head(T *inObj) { return inObj ? HeadOp::get(inObj) : nullptr; }
    inline const T *tail(const T *inObj) { return inObj ? TailOp::get(inObj) : nullptr; }
    inline const T *head(const T *inObj) { return inObj ? HeadOp::get(inObj) : nullptr; }

    void remove(T &inObj)
    {
        T *theHead = HeadOp::get(inObj);
        T *theTail = TailOp::get(inObj);
        if (theHead)
            TailOp::set(*theHead, theTail);
        if (theTail)
            HeadOp::set(*theTail, theHead);
        HeadOp::set(inObj, nullptr);
        TailOp::set(inObj, nullptr);
    }

    void insert_after(T &inPosition, T &inObj)
    {
        T *theHead = &inPosition;
        T *theTail = TailOp::get(inPosition);
        insert_unsafe(theHead, theTail, inObj);
    }

    void insert_before(T &inPosition, T &inObj)
    {
        T *theHead = HeadOp::get(inPosition);
        T *theTail = &inPosition;
        insert_unsafe(theHead, theTail, inObj);
    }

    // The name here is intentionally named "unsafe" to discourage use
    // with out knowing what it implies to call this function.
    // In most cases this will be used to insert before and after a neighboring pair
    // (see: insert_before()/_after), but it can also be convenient if you want to split
    // up a list and retain the chain for the removed section.
    void insert_unsafe(T *inHead, T *inTail, T &inObj)
    {
        if (inHead)
            TailOp::set(*inHead, &inObj);
        if (inTail)
            HeadOp::set(*inTail, &inObj);
        HeadOp::set(inObj, inHead);
        TailOp::set(inObj, inTail);
    }
};

template<typename T, typename TailOp>
struct QSSGLinkedListIterator
{
    using Iterator = QSSGLinkedListIterator<T, TailOp>;
    T *m_obj;
    explicit QSSGLinkedListIterator(T *inObj = nullptr) : m_obj(inObj) {}

    inline bool operator!=(const Iterator &inIter) const { return m_obj != inIter.m_obj; }
    inline bool operator==(const Iterator &inIter) const { return m_obj == inIter.m_obj; }

    Iterator &operator++()
    {
        if (m_obj)
            m_obj = TailOp::get(*m_obj);
        return *this;
    }

    Iterator &operator++(int)
    {
        Iterator retval(*this);
        ++(*this);
        return retval;
    }

    T &operator*() { return *m_obj; }
    T *operator->() { return m_obj; }
};

template<typename T, T *T::*Next>
struct QSSGInvasiveSingleLinkedList : public QSSGInvasiveLinkListBase<T, QSSGNullOp<T>, QSSGListAccessorNext<T, Next>>
{
    using TailOp = QSSGListAccessorNext<T, Next>;
    using List = QSSGInvasiveSingleLinkedList<T, Next>;
    using BaseList = QSSGInvasiveLinkListBase<T, QSSGNullOp<T>, TailOp>;
    using iterator = QSSGLinkedListIterator<T, TailOp>;
    using const_iterator = iterator;
    T *m_head = nullptr;
    QSSGInvasiveSingleLinkedList() = default;
    QSSGInvasiveSingleLinkedList(const List &inOther) : m_head(inOther.m_head) {}
    List &operator=(const List &inOther)
    {
        m_head = inOther.m_head;
        return *this;
    }

    inline  T &front() const { return *m_head; }

    void push_front(T &inObj)
    {
        if (m_head != nullptr)
            BaseList::insert_before(*m_head, inObj);
        m_head = &inObj;
    }

    void push_back(T &inObj)
    {
        if (m_head == nullptr)
            m_head = &inObj;
        else {
            T *lastObj = nullptr;
            for (iterator iter = begin(), endIter = end(); iter != endIter; ++iter)
                lastObj = &(*iter);

            Q_ASSERT(lastObj);
            if (lastObj)
                TailOp::set(*lastObj, &inObj);
        }
    }

    void remove(T &inObj)
    {
        if (m_head == &inObj)
            m_head = TailOp::get(inObj);
        BaseList::remove(inObj);
    }

    inline bool empty() const { return m_head == nullptr; }

    inline iterator begin() { return iterator(m_head); }
    inline iterator end() { return iterator(nullptr); }
    inline const_iterator begin() const { return iterator(m_head); }
    inline const_iterator end() const { return iterator(nullptr); }
};

template<typename T, T *T::*Previous, T *T::*Next>
struct QSSGInvasiveLinkedList : public QSSGInvasiveLinkListBase<T, QSSGListAccessorPrevious<T, Previous>, QSSGListAccessorNext<T, Next>>
{
    using HeadOp = QSSGListAccessorPrevious<T, Previous>;
    using TailOp = QSSGListAccessorNext<T, Next>;
    using List = QSSGInvasiveLinkedList<T, Previous, Next>;
    using BaseList = QSSGInvasiveLinkListBase<T, HeadOp, TailOp>;
    using iterator = QSSGLinkedListIterator<T, TailOp>;
    using const_iterator = iterator;
    using reverse_iterator = QSSGLinkedListIterator<T, HeadOp>;
    using const_reverse_iterator = reverse_iterator;

    T *m_head = nullptr;
    T *m_tail = nullptr;

    QSSGInvasiveLinkedList() = default;
    QSSGInvasiveLinkedList(const List &inOther) : m_head(inOther.m_head), m_tail(inOther.m_tail) {}
    List &operator=(const List &inOther)
    {
        m_head = inOther.m_head;
        m_tail = inOther.m_tail;
        return *this;
    }

    inline T &front() const
    {
        Q_ASSERT(m_head);
        return *m_head;
    }
    inline T &back() const
    {
        Q_ASSERT(m_tail);
        return *m_tail;
    }

    inline T *front_ptr() const { return m_head; }
    inline T *back_ptr() const { return m_tail; }

    void push_front(T &inObj)
    {
        if (m_head != nullptr)
            BaseList::insert_before(*m_head, inObj);
        m_head = &inObj;

        if (m_tail == nullptr)
            m_tail = &inObj;
    }

    void push_back(T &inObj)
    {
        if (m_tail != nullptr)
            BaseList::insert_after(*m_tail, inObj);
        m_tail = &inObj;

        if (m_head == nullptr)
            m_head = &inObj;
    }

    void remove(T &inObj)
    {
        if (m_head == &inObj)
            m_head = TailOp::get(inObj);
        if (m_tail == &inObj)
            m_tail = HeadOp::get(inObj);

        BaseList::remove(inObj);
    }

    inline bool empty() const { return m_head == nullptr; }

    inline iterator begin() { return iterator(m_head); }
    inline iterator end() { return iterator(nullptr); }

    inline const_iterator begin() const { return iterator(m_head); }
    inline const_iterator end() const { return iterator(nullptr); }

    inline reverse_iterator rbegin() { return reverse_iterator(m_tail); }
    inline reverse_iterator rend() { return reverse_iterator(nullptr); }

    inline const_reverse_iterator rbegin() const { return reverse_iterator(m_tail); }
    inline const_reverse_iterator rend() const { return reverse_iterator(nullptr); }
};

QT_END_NAMESPACE

#endif // QSSGINVASIVELINKEDLIST_H
