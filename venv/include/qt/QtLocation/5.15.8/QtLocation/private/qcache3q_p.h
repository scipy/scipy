/****************************************************************************
**
** Copyright (C) 2015 The Qt Company Ltd.
** Contact: http://www.qt.io/licensing/
**
** This file is part of the QtCore module of the Qt Toolkit.
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

#ifndef QCACHE3Q_H
#define QCACHE3Q_H

//
//  W A R N I N G
//  -------------
//
// This file is not part of the Qt API.  It exists purely as an
// implementation detail.  This header file may change from version to
// version without notice, or even be removed.
//
// We mean it.

#include <QtCore/qlinkedlist.h>
#include <QtCore/qhash.h>
#include <QtCore/qcache.h>
#include <QtCore/qsharedpointer.h>
#include <QDebug>

QT_BEGIN_NAMESPACE

template <class Key, class T>
class QCache3QDefaultEvictionPolicy
{
protected:
    /* called just before a key/value pair is about to be _evicted_ */
    void aboutToBeEvicted(const Key &key, QSharedPointer<T> obj);
    /* called just before a key/value pair is about to be removed, by
     * clear(), remove() or by the destructor (which calls clear) */
    void aboutToBeRemoved(const Key &key, QSharedPointer<T> obj);
};

template <class Key, class T>
void QCache3QDefaultEvictionPolicy<Key,T>::aboutToBeEvicted(const Key &key, QSharedPointer<T> obj)
{
    Q_UNUSED(key);
    Q_UNUSED(obj);
}

template <class Key, class T>
void QCache3QDefaultEvictionPolicy<Key,T>::aboutToBeRemoved(const Key &key, QSharedPointer<T> obj)
{
    Q_UNUSED(key);
    Q_UNUSED(obj);
}

/*
 * QCache3Q
 *
 * A cache template class for managing QSharedPointers to objects with an
 * associated key. It's a lot like QCache, but uses an alternative algorithm
 * '3Q' -- which is the '2Q' algorithm plus an extra queue for previously popular
 * but evicted nodes, and a 'ghost' list of recent evictions to make a better
 * placement choice if they are requested again.
 *
 * New nodes enter the cache on the "newbies" queue, which is evicted LRA
 * (least-recently-added). If a newbie is popular enough (it has been requested
 * more than promoteAt times), it will be promoted to a "regular". Regulars
 * are evicted LRU (least-recently-used). If a regular is under consideration
 * for eviction, its popularity is compared to the mean popularity of the whole
 * regulars queue. If it is greater, it is instead moved to the "hobos" queue.
 * The "hobos" queue is also evicted LRU, but has a maximum size constraint
 * so eviction from it is less likely than from the regulars.
 *
 * Tweakables:
 *  * maxCost = maximum total cost for the whole cache
 *  * minRecent = minimum size that q1 ("newbies") has to be before eviction
 *                from it takes place
 *  * maxOldPopular = maximum size that q3 ("hobos") can reach before eviction
 *                    from it takes place
 *  * promoteAt = minimum popularity necessary to promote a node from
 *                "newbie" to "regular"
 */
template <class Key, class T, class EvPolicy = QCache3QDefaultEvictionPolicy<Key,T> >
class QCache3Q : public EvPolicy
{
private:
    class Queue;
    class Node
    {
    public:
        inline explicit Node() : q(0), n(0), p(0), pop(0), cost(0) {}

        Queue *q;
        Node *n;
        Node *p;
        Key k;
        QSharedPointer<T> v;
        quint64 pop;                // popularity, incremented each ping
        int cost;
    };

    class Queue
    {
    public:
        inline explicit Queue() : f(0), l(0), cost(0), pop(0), size(0) {}

        Node *f;
        Node *l;
        int cost;               // total cost of nodes on the queue
        quint64 pop;            // sum of popularity values on the queue
        int size;               // size of the queue
    };

    Queue *q1_;          // "newbies": seen only once, evicted LRA (least-recently-added)
    Queue *q2_;          // regular nodes, promoted from newbies, evicted LRU
    Queue *q3_;          // "hobos": evicted from q2 but were very popular (above mean)
    Queue *q1_evicted_;  // ghosts of recently evicted newbies and regulars
    QHash<Key, Node *> lookup_;

public:
    explicit QCache3Q(int maxCost = 0, int minRecent = -1, int maxOldPopular = -1);
    inline ~QCache3Q() { clear(); delete q1_; delete q2_; delete q3_; delete q1_evicted_; }

    inline int maxCost() const { return maxCost_; }
    void setMaxCost(int maxCost, int minRecent = -1, int maxOldPopular = -1);

    inline int promoteAt() const { return promote_; }
    inline void setPromoteAt(int p) { promote_ = p; }

    inline int totalCost() const { return q1_->cost + q2_->cost + q3_->cost; }

    void clear();
    bool insert(const Key &key, QSharedPointer<T> object, int cost = 1);
    QSharedPointer<T> object(const Key &key) const;
    QSharedPointer<T> operator[](const Key &key) const;

    void remove(const Key &key, bool force = false);
    QList<Key> keys() const;
    void printStats();

    // Copy data directly into a queue. Designed for single use after construction
    void deserializeQueue(int queueNumber, const QList<Key> &keys,
                          const QList<QSharedPointer<T> > &values, const QList<int> &costs);
    // Copy data from specific queue into list
    void serializeQueue(int queueNumber, QList<QSharedPointer<T> > &buffer);

private:
    int maxCost_, minRecent_, maxOldPopular_;
    int hitCount_, missCount_, promote_;

    void rebalance();
    void unlink(Node *n);
    void link_front(Node *n, Queue *q);

private:
    // make these private so they can't be used
    inline QCache3Q(const QCache3Q<Key,T,EvPolicy> &) {}
    inline QCache3Q<Key,T,EvPolicy> &operator=(const QCache3Q<Key,T,EvPolicy> &) {}
};

template <class Key, class T, class EvPolicy>
void QCache3Q<Key,T,EvPolicy>::printStats()
{
    qDebug("\n=== cache %p ===", this);
    qDebug("hits: %d (%.2f%%)\tmisses: %d\tfill: %.2f%%", hitCount_,
           100.0 * float(hitCount_) / (float(hitCount_ + missCount_)),
           missCount_,
           100.0 * float(totalCost()) / float(maxCost()));
    qDebug("q1g: size=%d, pop=%llu", q1_evicted_->size, q1_evicted_->pop);
    qDebug("q1:  cost=%d, size=%d, pop=%llu", q1_->cost, q1_->size, q1_->pop);
    qDebug("q2:  cost=%d, size=%d, pop=%llu", q2_->cost, q2_->size, q2_->pop);
    qDebug("q3:  cost=%d, size=%d, pop=%llu", q3_->cost, q3_->size, q3_->pop);
}

template <class Key, class T, class EvPolicy>
QCache3Q<Key,T,EvPolicy>::QCache3Q(int maxCost, int minRecent, int maxOldPopular)
    : q1_(new Queue), q2_(new Queue), q3_(new Queue), q1_evicted_(new Queue),
      maxCost_(maxCost), minRecent_(minRecent), maxOldPopular_(maxOldPopular),
      hitCount_(0), missCount_(0), promote_(0)
{
    if (minRecent_ < 0)
        minRecent_ = maxCost_ / 3;
    if (maxOldPopular_ < 0)
        maxOldPopular_ = maxCost_ / 5;
}

template <class Key, class T, class EvPolicy>
void QCache3Q<Key,T,EvPolicy>::serializeQueue(int queueNumber, QList<QSharedPointer<T> > &buffer)
{
    Q_ASSERT(queueNumber >= 1 && queueNumber <= 4);
    Queue *queue = queueNumber == 1 ? q1_ :
                   queueNumber == 2 ? q2_ :
                   queueNumber == 3 ? q3_ :
                                      q1_evicted_;
    for (Node *node = queue->f; node; node = node->n)
        buffer.append(node->v);
}

template <class Key, class T, class EvPolicy>
void QCache3Q<Key,T,EvPolicy>::deserializeQueue(int queueNumber, const QList<Key> &keys,
                       const QList<QSharedPointer<T> > &values, const QList<int> &costs)
{
    Q_ASSERT(queueNumber >= 1 && queueNumber <= 4);
    int bufferSize = keys.size();
    if (bufferSize == 0)
        return;
    clear();
    Queue *queue = queueNumber == 1 ? q1_ :
                   queueNumber == 2 ? q2_ :
                   queueNumber == 3 ? q3_ :
                                      q1_evicted_;
    for (int i = 0; i<bufferSize; ++i) {
        Node *node = new Node;
        node->v = values[i];
        node->k = keys[i];
        node->cost = costs[i];
        link_front(node, queue);
        lookup_[keys[i]] = node;
    }
}


template <class Key, class T, class EvPolicy>
inline void QCache3Q<Key,T,EvPolicy>::setMaxCost(int maxCost, int minRecent, int maxOldPopular)
{
    maxCost_ = maxCost;
    minRecent_ = minRecent;
    maxOldPopular_ = maxOldPopular;
    if (minRecent_ < 0)
        minRecent_ = maxCost_ / 3;
    if (maxOldPopular_ < 0)
        maxOldPopular_ = maxCost_ / 5;
    rebalance();
}

template <class Key, class T, class EvPolicy>
bool QCache3Q<Key,T,EvPolicy>::insert(const Key &key, QSharedPointer<T> object, int cost)
{
    if (cost > maxCost_) {
        return false;
    }

    if (lookup_.contains(key)) {
        Node *n = lookup_[key];
        n->v = object;
        n->q->cost -= n->cost;
        n->cost = cost;
        n->q->cost += cost;

        if (n->q == q1_evicted_) {
            if (n->pop > (uint)promote_) {
                unlink(n);
                link_front(n, q2_);
                rebalance();
            }
        } else if (n->q != q1_) {
            Queue *q = n->q;
            unlink(n);
            link_front(n, q);
            rebalance();
        }

        return true;
    }

    Node *n = new Node;
    n->v = object;
    n->k = key;
    n->cost = cost;
    link_front(n, q1_);
    lookup_[key] = n;

    rebalance();

    return true;
}

template <class Key, class T, class EvPolicy>
void QCache3Q<Key,T,EvPolicy>::clear()
{
    while (q1_evicted_->f) {
        Node *n = q1_evicted_->f;
        unlink(n);
        delete n;
    }

    while (q1_->f) {
        Node *n = q1_->f;
        unlink(n);
        EvPolicy::aboutToBeRemoved(n->k, n->v);
        delete n;
    }

    while (q2_->f) {
        Node *n = q2_->f;
        unlink(n);
        EvPolicy::aboutToBeRemoved(n->k, n->v);
        delete n;
    }

    while (q3_->f) {
        Node *n = q3_->f;
        unlink(n);
        EvPolicy::aboutToBeRemoved(n->k, n->v);
        delete n;
    }

    lookup_.clear();
}

template <class Key, class T, class EvPolicy>
void QCache3Q<Key,T,EvPolicy>::unlink(Node *n)
{
    if (n->n)
        n->n->p = n->p;
    if (n->p)
        n->p->n = n->n;
    if (n->q->f == n)
        n->q->f = n->n;
    if (n->q->l == n)
        n->q->l = n->p;
    n->n = 0;
    n->p = 0;
    n->q->pop -= n->pop;
    n->q->cost -= n->cost;
    n->q->size--;
    n->q = 0;
}

template <class Key, class T, class EvPolicy>
void QCache3Q<Key,T,EvPolicy>::link_front(Node *n, Queue *q)
{
    n->n = q->f;
    n->p = 0;
    n->q = q;
    if (q->f)
        q->f->p = n;
    q->f = n;
    if (!q->l)
        q->l = n;

    q->pop += n->pop;
    q->cost += n->cost;
    q->size++;
}

template <class Key, class T, class EvPolicy>
void QCache3Q<Key,T,EvPolicy>::rebalance()
{
    while (q1_evicted_->size > (q1_->size + q2_->size + q3_->size) * 4) {
        Node *n = q1_evicted_->l;
        unlink(n);
        lookup_.remove(n->k);
        delete n;
    }

    while ((q1_->cost + q2_->cost + q3_->cost) > maxCost_) {
        if (q3_->cost > maxOldPopular_) {
            Node *n = q3_->l;
            unlink(n);
            EvPolicy::aboutToBeEvicted(n->k, n->v);
            lookup_.remove(n->k);
            delete n;
        } else if (q1_->cost > minRecent_) {
            Node *n = q1_->l;
            unlink(n);
            EvPolicy::aboutToBeEvicted(n->k, n->v);
            n->v.clear();
            n->cost = 0;
            link_front(n, q1_evicted_);
        } else {
            Node *n = q2_->l;
            unlink(n);
            if (q2_->size && n->pop > (q2_->pop / q2_->size)) {
                link_front(n, q3_);
            } else {
                EvPolicy::aboutToBeEvicted(n->k, n->v);
                n->v.clear();
                n->cost = 0;
                link_front(n, q1_evicted_);
            }
        }
    }
}

template <class Key, class T, class EvPolicy>
void QCache3Q<Key,T,EvPolicy>::remove(const Key &key, bool force)
{
    if (!lookup_.contains(key)) {
        return;
    }
    Node *n = lookup_[key];
    unlink(n);
    if (n->q != q1_evicted_ && !force)
        EvPolicy::aboutToBeRemoved(n->k, n->v);
    lookup_.remove(key);
    delete n;
}

template <class Key, class T, class EvPolicy>
QList<Key> QCache3Q<Key,T,EvPolicy>::keys() const
{
    return lookup_.keys();
}

template <class Key, class T, class EvPolicy>
QSharedPointer<T> QCache3Q<Key,T,EvPolicy>::object(const Key &key) const
{
    if (!lookup_.contains(key)) {
        const_cast<QCache3Q<Key,T,EvPolicy> *>(this)->missCount_++;
        return QSharedPointer<T>(0);
    }

    QCache3Q<Key,T,EvPolicy> *me = const_cast<QCache3Q<Key,T,EvPolicy> *>(this);

    Node *n = me->lookup_[key];
    n->pop++;
    n->q->pop++;

    if (n->q == q1_) {
        me->hitCount_++;

        if (n->pop > (quint64)promote_) {
            me->unlink(n);
            me->link_front(n, q2_);
            me->rebalance();
        }
    } else if (n->q != q1_evicted_) {
        me->hitCount_++;

        Queue *q = n->q;
        me->unlink(n);
        me->link_front(n, q);
        me->rebalance();
    } else {
        me->missCount_++;
    }

    return n->v;
}

template <class Key, class T, class EvPolicy>
inline QSharedPointer<T> QCache3Q<Key,T,EvPolicy>::operator[](const Key &key) const
{
    return object(key);
}

QT_END_NAMESPACE

#endif // QCACHE3Q_H
